// Copyright (C) 2025 [Thomas Steimlé]
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU Affero General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU Affero General Public License for more details.
//
// You should have received a copy of the GNU Affero General Public License
// along with this program. If not, see <https://www.gnu.org/licenses/>.

use duct::cmd;

pub fn bwa_aln(ref_path: &str, query_sequence: &str) -> Vec<Sam> {
    let stdout = cmd!(
        "echo",
        "-e",
        format!("\">sequence_to_query\n{}\"", query_sequence)
    )
    .pipe(cmd!("bwa", "mem", ref_path, "-"))
    .stdout_capture()
    .stderr_capture()
    .read()
    .unwrap();

    stdout
        .split("\n")
        .filter(|e| e.contains("sequence_to_query"))
        .map(Sam::new_from_bwa_output)
        .map(|mut e| {
            // if the sequence is returned by bwa with an N at first position remove it
            if e.seq.len() != query_sequence.len() {
                if e.seq.starts_with('N') {
                    let (_, seq) = e.seq.split_at(1);
                    e.seq = seq.to_string();
                    let cig = e.cigar.0.get_mut(0).unwrap();
                    cig.1 -= 1;
                }
                if e.seq.ends_with('N') {
                    let (seq, _) = e.seq.split_at(e.seq.len());
                    e.seq = seq.to_string();
                    let cl = e.cigar.0.len();
                    let cig = e.cigar.0.get_mut(cl - 1).unwrap();
                    cig.1 -= 1;
                }
            }
            e
        })
        .collect()
}

#[derive(Debug, Clone)]
pub struct Sam {
    pub qname: String,
    pub flag: u16,
    pub ref_name: String,
    pub pos: i32,
    pub mapq: i32,
    pub cigar: Cigar,
    pub rnext: String,
    pub pnext: String,
    pub tlen: String,
    pub seq: String,
}

impl Sam {
    pub fn new_from_bwa_output(sam_string: &str) -> Sam {
        let mut split = sam_string.split("\t");
        Sam {
            qname: split.next().unwrap().to_string(),
            flag: split.next().unwrap().parse().unwrap(),
            ref_name: split.next().unwrap().to_string(),
            pos: split.next().unwrap().parse().unwrap(),
            mapq: split.next().unwrap().parse().unwrap(),
            cigar: Cigar::new_from_str(split.next().unwrap()),
            rnext: split.next().unwrap().to_string(),
            pnext: split.next().unwrap().to_string(),
            tlen: split.next().unwrap().to_string(),
            seq: split.next().unwrap().to_string(),
        }
    }

    pub fn is_rc(&self) -> bool {
        self.flag & 0x10 == 0x10
    }
}

#[derive(Debug, Clone)]
pub struct Cigar(pub Vec<(String, u32)>);

impl Cigar {
    pub fn new_from_str(cigar_str: &str) -> Cigar {
        let mut res: Vec<(String, u32)> = Vec::new();
        let mut num_acc = String::new();
        for c in cigar_str.split("") {
            if !c.is_empty() {
                match c.parse::<usize>() {
                    Ok(_) => num_acc.push_str(c), // if a number added to the accumulator
                    Err(_) => {
                        let current_add = num_acc.parse::<u32>().unwrap();
                        num_acc = "".to_string();
                        res.push((c.to_string(), current_add));
                    }
                }
            }
        }
        Cigar(res)
    }

    pub fn len(&self) -> u32 {
        self.0.iter().map(|e| e.1).reduce(|acc, e| acc + e).unwrap()
    }

    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    pub fn expand(&self) -> String {
        let mut res = String::new();
        for (op, n) in self.0.iter() {
            let c = op.chars().next().unwrap();
            for _ in 0..(*n) {
                res.push(c);
            }
        }
        res
    }
}

type FullRange = ((String, i32, i32), (i32, i32));
pub struct Alignments {
    pub query_sequence: String,
    pub sam: Vec<Sam>,
    pub ranges: Vec<Vec<FullRange>>,
}

impl Alignments {
    pub fn to_reference(query_sequence: &str, ref_path: &str) -> Alignments {
        let res_sam = bwa_aln(ref_path, query_sequence);

        let mut res_ranges = Vec::new();
        for sam in &res_sam {
            res_ranges.push(
                ranges_from_cigar(sam.pos, sam.cigar.clone(), sam.is_rc())
                    .into_iter()
                    .map(|(reference, query)| {
                        ((sam.ref_name.clone(), reference.0, reference.1), query)
                    })
                    .collect(),
            );
        }

        Alignments {
            query_sequence: query_sequence.to_string(),
            sam: res_sam,
            ranges: res_ranges,
        }
    }
}

pub fn revcomp(dna: &str) -> String {
    let mut rdna: String = String::with_capacity(dna.len());
    for c in dna.chars().rev() {
        rdna.push(switch_base(c));
    }
    rdna
}

pub fn switch_base(c: char) -> char {
    match c {
        'a' => 't',
        'c' => 'g',
        't' => 'a',
        'g' => 'c',
        'u' => 'a',
        'A' => 'T',
        'C' => 'G',
        'T' => 'A',
        'G' => 'C',
        'U' => 'A',
        _ => 'N',
    }
}

pub struct BwaAlign {
    aligner: BwaAligner,
    header: rust_htslib::bam::header::Header,
}

impl BwaAlign {
    pub fn init(ref_path: &str) -> BwaAlign {
        let bwa_ref = BwaReference::open(ref_path).unwrap();
        let header = bwa_ref.create_bam_header();

        BwaAlign {
            aligner: BwaAligner::new(bwa_ref, BwaSettings::new(), PairedEndStats::default()),
            header,
        }
    }

    pub fn get_ref_positions(
        &self,
        sequence: &str,
    ) -> Vec<FullRange> {
        let mut all_ranges = Vec::new();
        let (a, _) = self.aligner.align_read_pair(
            b"read_name",
            sequence.as_bytes(),
            &vec![b'2'; sequence.len()],
            String::new().as_bytes(),
            String::new().as_bytes(),
        );

        for record in a.iter() {
            let contig = self.bwa_tid_contig(record.tid() as usize);
            if let Some(contig) = contig {
                let res = ranges_from_cigar(
                    record.pos() as i32 + 1,
                    Cigar::new_from_str(&format!("{}", record.cigar())),
                    record.is_reverse(),
                );
                res.into_iter().for_each(|(reference, query)| {
                    all_ranges.push(((contig.clone(), reference.0, reference.1), query))
                });
            }
        }
        all_ranges
    }
    pub fn bwa_tid_contig(&self, tid: usize) -> Option<String> {
        if let Some(h) = self.header.to_hashmap().get("SQ") {
            if let Some(h) = h.get(tid) {
                h.get("SN").cloned()
            } else {
                None
            }
        } else {
            None
        }
    }
}

/// return [start, stop] 1-based positions relative to reference and query
pub fn ranges_from_cigar(
    ref_start: i32,
    cigar: Cigar,
    is_rc: bool,
) -> Vec<((i32, i32), (i32, i32))> {
    let cigar = if is_rc {
        cigar.expand().chars().rev().collect::<Vec<char>>()
    } else {
        cigar.expand().chars().collect::<Vec<char>>()
    };
    let mut ref_counter = 0;
    let mut query_counter = 0;

    let mut query_matches: Vec<i32> = Vec::new();
    let mut ref_matches: Vec<i32> = Vec::new();
    for op in cigar.iter() {
        match op {
            'M' | '=' => {
                query_matches.push(query_counter);
                ref_matches.push(ref_counter);

                ref_counter += 1;
                query_counter += 1;
            }
            'S' | 'H' | 'I' => {
                query_counter += 1;
            }
            'D' | 'N' => {
                ref_counter += 1;
            }
            _ => panic!("Unknown cigar operand"),
        }
    }

    let mut query_ranges = Vec::new();
    let mut acc_range = None;
    for qp in query_matches.into_iter() {
        if let Some((start, stop)) = acc_range {
            if stop == qp - 1 {
                acc_range = Some((start, qp));
            } else {
                query_ranges.push(acc_range.unwrap());
                acc_range = Some((qp, qp));
            }
        } else {
            acc_range = Some((qp, qp));
        }
    }
    query_ranges.push(acc_range.unwrap());
    query_ranges = query_ranges.iter().map(|e| (e.0 + 1, e.1 + 1)).collect();

    let mut ref_ranges = Vec::new();
    let mut acc_range = None;
    for qp in ref_matches.into_iter() {
        if let Some((start, stop)) = acc_range {
            if stop == qp - 1 {
                acc_range = Some((start, qp));
            } else {
                query_ranges.push(acc_range.unwrap());
                acc_range = Some((qp, qp));
            }
        } else {
            acc_range = Some((qp, qp));
        }
    }
    ref_ranges.push(acc_range.unwrap());

    ref_ranges = ref_ranges
        .iter()
        .map(|(start, stop)| {
            if is_rc {
                (ref_start + stop, ref_start + start)
            } else {
                (ref_start + start, ref_start + stop)
            }
        })
        .collect();

    ref_ranges
        .into_iter()
        .zip(query_ranges)
        .collect()
}

pub fn format_seq(name: &str, sequence: &str, line_size: usize) -> String {
    // todo check fasta name compl
    let res = sequence
        .chars()
        .collect::<Vec<char>>()
        .chunks(line_size)
        .map(|c| c.iter().collect::<String>())
        .collect::<Vec<String>>()
        .join("\n");
    format!(">{name}\n{res}")
}

use bwa::{BwaAligner, BwaReference, BwaSettings, PairedEndStats};

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_1() {
        let ref_path = "/home/thomas/NGS/ref/hg19.fa";

        //2011N160703-HAOAA 1
        let sequence = "
        CCTCCCGCCCTGCCTTTCAGGTCGGGAAAGTCCCGGGGTTTGCAAAAGAGTGTCCGAGCG
        CCCTGGAGGCGGGGAGGGCGGCAAGGAGGGCGCCGGTGTCGCGGTTGAGTTTCTCCACTG
        CCGACCGCGGCCACGCTGCCCGGGGCTTCCCGGACAGGCTTCGCGCCGCCCACCTCGGCA
        GCCGGGGCGGA
        TGG
        CCACACAGGTTGGAGTGCATTAAGCCTTTGTCCAAAAACACCCAGC
        CGTGACCCGCTATGTATGTCTCAGCATTGGGAAGAGTCCTCTGAGTGTCATGGGAAAATA
        ATATATGAGTT"
            .replace("\n", "")
            .replace(" ", "");

        let alns = Alignments::to_reference(&sequence, ref_path);
        assert_eq!(
            alns.ranges,
            vec![
                vec![(("chr2".to_string(), 43453321, 43453131), (1, 191))],
                vec![(("chr14".to_string(), 22908007, 22908123), (195, 311))]
            ]
        );
    }
    #[test]
    fn test_2() {
        let ref_path = "/home/thomas/NGS/ref/hg19.fa";
        //2011N160703-HAOAA 2
        let sequence = "
        TGGAATTGTTCATCACTGATGTTTTTGGTGATATTTGTTTATGTTCTGAAGCTATTGCTG
        TAGACCAACATGGAGTAAACAGAAAATAATTGGGACACAGGGGCAATAAAACTCATATAT
        TATTTTCCCATGACACTCAGAGGACTCTTCCCAATGCTGAGACATACATAGCGGGTCACG
        GCTGGGTGTTTTTGGACAAAGGCTTAATGCACTCCAACCTGTGTGGCCATCCGCCCCGGC
        TGCCGAGGTGGGCGGCGCGAAGCCTGTCCGGGAAGCCCCGGGCAGCGTGGCCGCGGTCGG
        CAGTGGAGAAACTCAACCGCGACACC"
            .replace("\n", "")
            .replace(" ", "");

        //2011N160703-HAOAA 1
        let sequence_2 = "
        CCTCCCGCCCTGCCTTTCAGGTCGGGAAAGTCCCGGGGTTTGCAAAAGAGTGTCCGAGCG
        CCCTGGAGGCGGGGAGGGCGGCAAGGAGGGCGCCGGTGTCGCGGTTGAGTTTCTCCACTG
        CCGACCGCGGCCACGCTGCCCGGGGCTTCCCGGACAGGCTTCGCGCCGCCCACCTCGGCA
        GCCGGGGCGGA
        TGG
        CCACACAGGTTGGAGTGCATTAAGCCTTTGTCCAAAAACACCCAGC
        CGTGACCCGCTATGTATGTCTCAGCATTGGGAAGAGTCCTCTGAGTGTCATGGGAAAATA
        ATATATGAGTT"
            .replace("\n", "")
            .replace(" ", "");

        let aligner = BwaAlign::init(ref_path);
        let res_1 = aligner.get_ref_positions(&sequence);
        let res_2 = aligner.get_ref_positions(&sequence_2);

        assert_eq!(
            res_1,
            vec![
                (("chr14".to_string(), 22908232, 22908007), (1, 226)),
                (("chr2".to_string(), 43453131, 43453227), (230, 326))
            ]
        );

        assert_eq!(
            res_2,
            vec![
                (("chr2".to_string(), 43453321, 43453131), (1, 191)),
                (("chr14".to_string(), 22908007, 22908123), (195, 311))
            ]
        );
    }
}
