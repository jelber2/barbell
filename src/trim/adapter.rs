//! Adapter trimming for Oxford Nanopore ligation-based (kit14) libraries.
//!
//! Reads from kit14 libraries (e.g. `SQK-LSK114`, `SQK-ULK114`, and the other
//! "114" kits) carry the same ligation adapter at both ends of the insert:
//!
//! ```text
//!   5' [adapter] ------------------ insert ------------------ [adapter] 3'
//!                \                                            /
//!                 both are the same double-stranded AMX114 adapter
//! ```
//!
//! Depending on which strand was captured, the first bases of the read are either
//! the adapter itself or its reverse complement:
//!
//! * forward adapter  : `CCTGTACTTCGTTCAGTTACGTATTGC`
//! * reverse adapter  : `GCAATACGTAACTGAACGAAGTACAGG`
//!
//! (These are the sequences Dorado uses to trim `SQK-LSK114`, `SQK-ULK114`,
//! `SQK-PCS114`, `SQK-RAD114`, `SQK-16S114-24`, `SQK-NBD114-*`, `SQK-PCB114-24`,
//! `SQK-RBK114-*`, `SQK-RPB114-24`, `SQK-HTB114-96`, ... see
//! `dorado/demux/adapter_primer_kits.cpp`.)
//!
//! This module finds (partial) adapter matches at the very start/end of a read and
//! reports the interval of the read that should be kept, in the same spirit as
//! [Porechop](https://github.com/rrwick/Porechop). The read is *not* reverse
//! complemented: trimming is always relative to the left (5') end of the read as
//! it is stored in the FASTQ.
//!
//! Adapters in the middle of a read are deliberately ignored, because those cannot
//! be resolved by plain end trimming (that requires splitting/filtering, which is
//! out of scope here).

use anyhow::anyhow;
use paraseq::Record;
use sassy::profiles::Iupac;
use sassy::{EncodedPatterns, Match, Searcher};
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};

use crate::config::AdapterTrimConfig;
use crate::io::io::validate_fastq_paths;

/// The forward (5') AMX114 ligation adapter used by all kit14 sequencing kits.
pub const ADAPTER_SEQUENCE: &[u8] = b"CCTGTACTTCGTTCAGTTACGTATTGC";

/// An adapter is only considered to sit at a read end when its outer edge is
/// within this many bases of the very start/end of the read. This tolerates the
/// few noisy bases that basecalling can leave in front of the adapter (Porechop
/// uses a similar slack via its `extra_end_trim`).
pub const ADAPTER_END_SLACK: usize = 10;

/// Overhang cost used for the adapter search. With `alpha = 0.5` a partial
/// alignment that only covers half of the adapter costs as much as a single
/// mismatch, which is what we want: a confidently matched adapter *sub*-sequence
/// at the read end is still evidence of an adapter.
const ADAPTER_SEARCH_ALPHA: f32 = 0.5;

/// Maximum number of edits allowed while searching for an adapter hit.
/// Together with the overhang cost this still allows partial hits (e.g. 12 bases
/// with 3 mismatches, or 8 perfectly matching bases), while keeping the search
/// selective.
const ADAPTER_MAX_EDITS: usize = 12;

/// A hit needs at least this many matching bases ...
///
/// This was tuned on simulated kit14 reads: with 14 matching bases and 75%
/// identity ~98.5% of reads carrying an adapter are trimmed while random reads
/// are only touched in <0.1% of the cases. Lowering it picks up a few more
/// (shorter) partial adapters at the cost of more false positives on random
/// sequence.
pub const ADAPTER_MIN_MATCH_LEN: usize = 14;

/// ... and at least this fraction of the alignment needs to match (gaps count as
/// mismatches).
pub const ADAPTER_MIN_IDENTITY: f64 = 0.75;

/// A trim that would leave less than this many bases is rejected (the read is
/// returned unchanged), to avoid emptying reads with many spurious matches.
const ADAPTER_MIN_INSERT_LEN: usize = 1;

/// A chained adapter has to be anchored within this many bases of the previous
/// cut. It allows a few bases in between while still rejecting an independent
/// adapter that merely happens to be close by.
const ADAPTER_CHAIN_OVERLAP: usize = 10;

/// Safety bound on the number of adapters trimmed in tandem at one read end.
const MAX_ADAPTERS_PER_END: usize = 8;

/// Which read end an adapter match was found on.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AdapterSide {
    /// The adapter (or its reverse complement) at the 5' end of the read.
    Start,
    /// The adapter (or its reverse complement) at the 3' end of the read.
    End,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct AdapterTrim {
    /// First base that is kept.
    pub keep_start: usize,
    /// One past the last base that is kept.
    pub keep_end: usize,
    /// Adapters that were found (left to right). Empty means nothing to trim.
    pub sides: [bool; 2],
}

impl AdapterTrim {
    pub fn nothing(read_len: usize) -> Self {
        Self {
            keep_start: 0,
            keep_end: read_len,
            sides: [false, false],
        }
    }

    pub fn is_empty(&self) -> bool {
        !self.sides[0] && !self.sides[1]
    }
}

/// Searches (partial) ligation adapters at the ends of Nanopore reads.
///
/// The internal search buffers make this cheap to reuse across reads, but it is
/// not `Sync`; instantiate one per thread.
pub struct AdapterTrimmer {
    searcher: Searcher<Iupac>,
    encoded: EncodedPatterns<Iupac>,
    min_identity: f64,
    min_match_len: usize,
    end_slack: usize,
}

impl Default for AdapterTrimmer {
    fn default() -> Self {
        Self::new(AdapterTrimConfig::default())
    }
}

impl AdapterTrimmer {
    pub fn new(config: AdapterTrimConfig) -> Self {
        let mut searcher = Searcher::<Iupac>::new_fwd_with_overhang(ADAPTER_SEARCH_ALPHA);
        // Pattern 0 is the adapter as found at the 5' end of a read, pattern 1 is
        // its reverse complement (which is what the 3' end looks like when the
        // other strand was captured).
        let fwd = ADAPTER_SEQUENCE.to_vec();
        let rc = reverse_complement(ADAPTER_SEQUENCE);
        let encoded = searcher.encode_patterns(&[fwd, rc]);
        Self {
            searcher,
            encoded,
            min_identity: config.min_identity,
            min_match_len: config.min_match_len,
            end_slack: config.end_slack,
        }
    }

    /// Find the adapters at the ends of `read` and return the interval that should
    /// be kept. `None` when there is nothing to trim.
    ///
    /// Only *terminal* adapters are considered, i.e. adapters that reach the very
    /// start/end of the read (allowing a few noisy bases). Adapters found in the
    /// middle of a read are ignored: those cannot be resolved by end trimming.
    ///
    /// Ligation artefacts can leave several adapters in tandem at a read end, so
    /// trimming continues as long as another adapter is anchored right at the
    /// current cut position.
    pub fn find_trim(&mut self, read: &[u8]) -> Option<AdapterTrim> {
        if read.is_empty() {
            return None;
        }

        let read_len = read.len();

        // Collect the hits first so that the (mutable) searcher borrow ends before
        // we inspect them.
        let hits: Vec<Match> = self
            .searcher
            .search_all_encoded_patterns(&self.encoded, read, ADAPTER_MAX_EDITS)
            .to_vec();

        let mut trim = AdapterTrim::nothing(read_len);
        if let Some((_, keep_start)) = self.resolve_side(&hits, read_len, AdapterSide::Start) {
            trim.keep_start = keep_start;
            trim.sides[0] = true;
        }
        if let Some((_, keep_end)) = self.resolve_side(&hits, read_len, AdapterSide::End) {
            trim.keep_end = keep_end;
            trim.sides[1] = true;
        }

        if trim.is_empty() {
            return None;
        }

        // Never return an empty interval; leaving the read as-is is better than
        // dropping it.
        if trim.keep_end <= trim.keep_start
            || trim.keep_end - trim.keep_start < ADAPTER_MIN_INSERT_LEN
        {
            return None;
        }
        Some(trim)
    }

    /// Walk inwards from one end of the read over consecutive adapters and return
    /// `(score, position_to_cut_at)`.
    ///
    /// The first adapter has to reach the very end of the read; each following one
    /// has to be anchored at the position of the previous cut (tandem adapters).
    fn resolve_side(
        &self,
        hits: &[Match],
        read_len: usize,
        side: AdapterSide,
    ) -> Option<(f64, usize)> {
        // Start at the very beginning/end of the read and move inwards.
        let mut pos = match side {
            AdapterSide::Start => 0,
            AdapterSide::End => read_len,
        };
        let mut result: Option<(f64, usize)> = None;

        // Each round has to advance `pos`, so this terminates; the bound is only
        // there to keep a pathological read from looping for a long time.
        for _ in 0..MAX_ADAPTERS_PER_END {
            let mut best: Option<(f64, usize)> = None;

            for m in hits {
                let Some((matches, identity)) =
                    score_match(m, self.min_match_len, self.min_identity)
                else {
                    continue;
                };
                // Weight matching bases by how clean the alignment is, so that an
                // alignment consisting mostly of insertions does not outrank a
                // clean (possibly shorter) one.
                let score = matches as f64 * identity;

                let (anchored, new_pos) = match side {
                    AdapterSide::Start => (
                        // Starts at/after the current cut (a few bases of junk in
                        // front of the adapter are allowed) ...
                        m.text_start <= pos + self.end_slack
                            // ... but anchored at it, not a re-alignment of the
                            // adapter we already walked past.
                            && m.text_start + ADAPTER_CHAIN_OVERLAP >= pos
                            // and it has to make progress.
                            && m.text_end > pos,
                        m.text_end,
                    ),
                    AdapterSide::End => (
                        m.text_end + self.end_slack >= pos
                            && m.text_end <= pos + ADAPTER_CHAIN_OVERLAP
                            && m.text_start < pos,
                        m.text_start,
                    ),
                };
                if !anchored {
                    continue;
                }

                let better = match best {
                    // On equal score walk over the outermost adapter first.
                    Some((best_score, best_pos)) => match side {
                        AdapterSide::Start => {
                            score > best_score || (score == best_score && new_pos > best_pos)
                        }
                        AdapterSide::End => {
                            score > best_score || (score == best_score && new_pos < best_pos)
                        }
                    },
                    None => true,
                };
                if better {
                    best = Some((score, new_pos));
                }
            }

            match best {
                Some((score, new_pos)) => {
                    pos = new_pos;
                    result = Some((score, new_pos));
                }
                None => break,
            }
        }

        result
    }

    /// Trim the adapters off a sequence (and quality string, if given).
    ///
    /// Returns `None` when there was nothing to trim.
    pub fn trim<'a>(
        &mut self,
        seq: &'a [u8],
        qual: Option<&'a [u8]>,
    ) -> Option<(&'a [u8], Option<&'a [u8]>)> {
        let trim = self.find_trim(seq)?;
        let seq = &seq[trim.keep_start..trim.keep_end];
        let qual = qual.map(|q| &q[trim.keep_start..trim.keep_end.min(q.len())]);
        Some((seq, qual))
    }
}

/// Number of matching bases and the identity of the alignment.
///
/// The identity is `matches / alignment_columns`, i.e. gaps count as mismatches.
/// This mirrors Porechop's "aligned region percent identity" and prevents
/// alignments full of insertions/deletions from being mistaken for an adapter.
fn score_match(m: &Match, min_match_len: usize, min_identity: f64) -> Option<(usize, f64)> {
    let mut matches = 0usize;
    let mut columns = 0usize;
    for op in &m.cigar.ops {
        let cnt = op.cnt as usize;
        columns += cnt;
        if op.op == pa_types::CigarOp::Match {
            matches += cnt;
        }
    }

    if columns == 0 || matches < min_match_len {
        return None;
    }
    let identity = matches as f64 / columns as f64;
    if identity < min_identity {
        return None;
    }
    Some((matches, identity))
}

/// Statistics of an adapter trimming run.
#[derive(Debug, Default, Clone, Copy)]
pub struct AdapterTrimStats {
    pub total: usize,
    pub trimmed_start: usize,
    pub trimmed_end: usize,
    pub trimmed_both: usize,
}

/// Standalone adapter trimming of FASTQ file(s).
///
/// Reads are written to `output` (a file, or a folder when `output` has no
/// extension) and keep their original order and orientation.
pub fn trim_adapters_in_fastq(
    input: &[PathBuf],
    output: &str,
    config: &AdapterTrimConfig,
) -> anyhow::Result<AdapterTrimStats> {
    validate_fastq_paths(input)?;
    let mut trimmer = AdapterTrimmer::new(config.clone());
    let mut stats = AdapterTrimStats::default();

    // A single input file can be written straight to `output`, multiple files need
    // one output per input, so they go into a folder.
    let single_file_output = input.len() == 1 && Path::new(output).extension().is_some();
    let writers: HashMap<PathBuf, Box<dyn Write>> = HashMap::new();
    let mut writers = writers;

    for read_path in input {
        let out_path = if single_file_output {
            PathBuf::from(output)
        } else {
            let stem = read_path
                .file_name()
                .map(|name| name.to_string_lossy().to_string())
                .unwrap_or_else(|| "reads.fastq".to_string());
            Path::new(output).join(format!("{stem}.adapter_trimmed.fastq"))
        };
        if let Some(parent) = out_path.parent()
            && !parent.as_os_str().is_empty()
        {
            std::fs::create_dir_all(parent)?;
        }

        let mut reader = paraseq::fastq::Reader::from_path(read_path)
            .map_err(|err| anyhow!("Failed to open FASTQ file '{}': {err}", read_path.display()))?;
        let mut record_set = reader.new_record_set();

        let file = File::create(&out_path).map_err(|err| {
            anyhow!(
                "Failed to create output file '{}': {err}",
                out_path.display()
            )
        })?;
        let gzip = config.gzip || out_path.extension().map(|ext| ext == "gz").unwrap_or(false);
        let writer: Box<dyn Write> = if gzip {
            Box::new(flate2::write::GzEncoder::new(
                file,
                flate2::Compression::default(),
            ))
        } else {
            Box::new(BufWriter::new(file))
        };
        writers.insert(out_path.clone(), writer);

        while record_set
            .fill(&mut reader)
            .map_err(|err| anyhow!("Error reading FASTQ file '{}': {err}", read_path.display()))?
        {
            for record in record_set.iter() {
                let record = record.map_err(|err| {
                    anyhow!(
                        "Error reading FASTQ record in '{}': {err}",
                        read_path.display()
                    )
                })?;
                stats.total += 1;

                let seq = record.seq();
                let qual = record.qual();

                let trim = trimmer
                    .find_trim(seq.as_ref())
                    .unwrap_or_else(|| AdapterTrim::nothing(seq.len()));
                let keep = trim.keep_start..trim.keep_end;
                stats.count(&trim);

                let writer = writers
                    .get_mut(&out_path)
                    .expect("writer was just inserted");
                writer.write_all(b"@")?;
                writer.write_all(record.id())?;
                writer.write_all(b"\n")?;
                writer.write_all(&seq.as_ref()[keep.clone()])?;
                writer.write_all(b"\n+\n")?;
                match qual {
                    Some(qual) => writer.write_all(&qual[keep.clone()])?,
                    None => writer.write_all(&vec![b'!'; keep.len()])?,
                }
                writer.write_all(b"\n")?;
            }
        }
    }

    for writer in writers.values_mut() {
        writer.flush()?;
    }

    Ok(stats)
}

impl AdapterTrimStats {
    fn count(&mut self, trim: &AdapterTrim) {
        match (trim.sides[0], trim.sides[1]) {
            (true, true) => self.trimmed_both += 1,
            (true, false) => self.trimmed_start += 1,
            (false, true) => self.trimmed_end += 1,
            (false, false) => {}
        }
    }

    pub fn trimmed(&self) -> usize {
        self.trimmed_start + self.trimmed_end + self.trimmed_both
    }
}

fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    seq.iter().rev().map(|&c| RC[c as usize]).collect()
}

const RC: [u8; 256] = {
    let mut rc = [0; 256];
    let mut i = 0;
    while i < 256 {
        rc[i] = i as u8;
        i += 1;
    }
    rc[b'A' as usize] = b'T';
    rc[b'C' as usize] = b'G';
    rc[b'T' as usize] = b'A';
    rc[b'G' as usize] = b'C';
    rc[b'a' as usize] = b't';
    rc[b'c' as usize] = b'g';
    rc[b't' as usize] = b'a';
    rc[b'g' as usize] = b'c';
    rc[b'N' as usize] = b'N';
    rc[b'n' as usize] = b'n';
    rc
};

#[cfg(test)]
mod tests {
    use super::*;

    fn trimmer() -> AdapterTrimmer {
        AdapterTrimmer::new(AdapterTrimConfig::default())
    }

    #[test]
    fn test_full_forward_adapter_at_start() {
        let mut t = trimmer();
        let insert = b"ACGTACGTACGTACGTACGT";
        let read = [ADAPTER_SEQUENCE, insert.as_slice()].concat();
        let trim = t.find_trim(&read).unwrap();
        assert_eq!(trim.keep_start, ADAPTER_SEQUENCE.len());
        assert_eq!(trim.keep_end, read.len());
        assert!(trim.sides[0] && !trim.sides[1]);
    }

    #[test]
    fn test_full_reverse_complement_adapter_at_start() {
        let mut t = trimmer();
        let insert = b"ACGTACGTACGTACGTACGT";
        let rc_adapter = reverse_complement(ADAPTER_SEQUENCE);
        let read = [rc_adapter.as_slice(), insert.as_slice()].concat();
        let trim = t.find_trim(&read).unwrap();
        assert_eq!(trim.keep_start, ADAPTER_SEQUENCE.len());
        assert!(trim.sides[0] && !trim.sides[1]);
    }

    #[test]
    fn test_full_forward_adapter_at_end() {
        let mut t = trimmer();
        let insert = b"ACGTACGTACGTACGTACGT";
        let read = [insert.as_slice(), ADAPTER_SEQUENCE].concat();
        let trim = t.find_trim(&read).unwrap();
        assert_eq!(trim.keep_start, 0);
        assert_eq!(trim.keep_end, insert.len());
        assert!(!trim.sides[0] && trim.sides[1]);
    }

    #[test]
    fn test_full_reverse_complement_adapter_at_end() {
        let mut t = trimmer();
        let insert = b"ACGTACGTACGTACGTACGT";
        let rc_adapter = reverse_complement(ADAPTER_SEQUENCE);
        let read = [insert.as_slice(), rc_adapter.as_slice()].concat();
        let trim = t.find_trim(&read).unwrap();
        assert_eq!(trim.keep_end, insert.len());
        assert!(!trim.sides[0] && trim.sides[1]);
    }

    #[test]
    fn test_adapters_at_both_ends() {
        let mut t = trimmer();
        let insert = b"ACGTACGTACGTACGTACGT";
        let read = [
            ADAPTER_SEQUENCE,
            insert.as_slice(),
            &reverse_complement(ADAPTER_SEQUENCE),
        ]
        .concat();
        let trim = t.find_trim(&read).unwrap();
        assert_eq!(trim.keep_start, ADAPTER_SEQUENCE.len());
        assert_eq!(trim.keep_end, ADAPTER_SEQUENCE.len() + insert.len());
        assert!(trim.sides[0] && trim.sides[1]);
    }

    #[test]
    fn test_partial_adapter_at_start() {
        let mut t = trimmer();
        // Only the last 20 bases of the adapter are present (basecalling can drop
        // the first bases of a read).
        let partial = &ADAPTER_SEQUENCE[7..];
        let insert = b"ACGTACGTACGTACGTACGT";
        let read = [partial, insert.as_slice()].concat();
        let trim = t.find_trim(&read).unwrap();
        assert_eq!(trim.keep_start, partial.len());
        assert!(trim.sides[0]);
    }

    #[test]
    fn test_adapter_with_mismatches_is_trimmed() {
        let mut t = trimmer();
        let mut adapter = ADAPTER_SEQUENCE.to_vec();
        adapter[3] = b'A';
        adapter[15] = b'C';
        let insert = b"ACGTACGTACGTACGTACGT";
        let read = [adapter.as_slice(), insert.as_slice()].concat();
        let trim = t.find_trim(&read).unwrap();
        assert_eq!(trim.keep_start, adapter.len());
    }

    #[test]
    fn test_no_adapter_is_untouched() {
        let mut t = trimmer();
        let read = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT".to_vec();
        assert!(t.find_trim(&read).is_none());
        assert!(t.trim(&read, None).is_none());
    }

    #[test]
    fn test_middle_adapter_is_not_trimmed() {
        let mut t = trimmer();
        let flank = vec![b'A'; 300];
        let read = [flank.as_slice(), ADAPTER_SEQUENCE, flank.as_slice()].concat();
        assert!(t.find_trim(&read).is_none());
    }

    #[test]
    fn test_trim_returns_slices() {
        let mut t = trimmer();
        let insert = b"ACGTACGTACGTACGTACGT";
        let read = [ADAPTER_SEQUENCE, insert.as_slice()].concat();
        let qual = vec![b'I'; read.len()];
        let (seq, qual) = t.trim(&read, Some(&qual)).unwrap();
        assert_eq!(seq, insert);
        assert_eq!(qual.unwrap().len(), insert.len());
    }

    #[test]
    fn test_trim_does_not_empty_read() {
        let mut t = trimmer();
        // A read that is (almost) nothing but adapter should not be emptied
        // entirely by an over-eager hit.
        let read = ADAPTER_SEQUENCE.to_vec();
        if let Some(trim) = t.find_trim(&read) {
            assert!(trim.keep_end > trim.keep_start);
        }
    }

    #[test]
    fn test_adapters_survive_internal_gaps() {
        let mut t = trimmer();
        // Adapter with a deletion and a mismatch in the middle.
        let mut adapter = ADAPTER_SEQUENCE.to_vec();
        adapter.remove(12);
        adapter[5] = b'A';
        let insert = b"ACGTACGTACGTACGTACGTACGTACGT";
        let read = [adapter.as_slice(), insert.as_slice()].concat();
        let trim = t.find_trim(&read).expect("should still trim");
        assert!(trim.keep_start >= ADAPTER_SEQUENCE.len() - 3);
    }

    #[test]
    fn test_terminal_junk_before_adapter_is_removed() {
        let mut t = trimmer();
        let insert = b"ACGTACGTACGTACGTACGT";
        // A few spurious bases in front of the adapter (common in ONT reads).
        let read = [b"GGGG".as_slice(), ADAPTER_SEQUENCE, insert.as_slice()].concat();
        let trim = t.find_trim(&read).unwrap();
        assert_eq!(trim.keep_start, 4 + ADAPTER_SEQUENCE.len());
    }

    #[test]
    fn test_two_adapters_in_tandem_at_start() {
        let mut t = trimmer();
        let insert = b"ACGTACGTACGTACGTACGT";
        // Ligation artefacts can leave two adapters in tandem.
        let read = [ADAPTER_SEQUENCE, ADAPTER_SEQUENCE, insert.as_slice()].concat();
        let trim = t.find_trim(&read).unwrap();
        assert_eq!(trim.keep_start, 2 * ADAPTER_SEQUENCE.len());
    }

    #[test]
    fn test_single_read_with_adapter_at_end_is_not_start_trimmed() {
        let mut t = trimmer();
        let insert = b"ACGTACGTACGTACGTACGT";
        let read = [insert.as_slice(), ADAPTER_SEQUENCE].concat();
        let trim = t.find_trim(&read).unwrap();
        assert!(!trim.sides[0], "must not trim the start of this read");
    }

    #[test]
    fn test_real_lsk114_read() {
        // First 493bp of the read used in dorado's adapter trimming tests
        // (tests/data/adapter_trim/lsk110_single_read.fastq). The 5' end is a
        // partial AMX114 adapter.
        let read = b"TACTTCGTTCCAGTTACGTATTGCTTGCAGGTCCAGACGGTGATGGTGATCGGCGGCGTGCTCATGATCGCAGGCGTGCTCCTTGGTCCGGTTCGGTTGATCGATGAAAAGCAGGAATCGCAAAAAAACAAAAACCCCGCCCGAAGGCGGGGTTTTCCCGAACAGCACAGTCGCAGGAATTAACAACATCCCCCCTATCAGGAGAGACGGCCCGCCGGCGAGCTGATCGAGCGGCCTGTGACCTGATCTG";
        let mut t = trimmer();
        let trim = t.find_trim(read).expect("adapter should be found");
        assert_eq!(trim.keep_start, 24);
        assert_eq!(trim.keep_end, read.len());
        assert!(trim.sides[0] && !trim.sides[1]);
    }
}
