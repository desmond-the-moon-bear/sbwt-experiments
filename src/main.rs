#![allow(unused)]

mod statistics;
mod vodbg;

use std::io::Read;

use sbwt::LcsArray;
use sbwt::SbwtIndexVariant;
use sbwt::{BitPackedKmerSortingMem, ContractLeft, ExtendRight, SbwtIndex, SbwtIndexBuilder, StreamingIndex, SubsetMatrix};
use sbwt::vodbg::*;
use sbwt::vodbg::benchmark::benchmark_bms_separate_queries;
use sbwt::vodbg::count::Counts;
use sbwt::vodbg::pnsv::{
    self,
    ABS,
    LcsPnsvBp,
    LcsSimd,
    Pnsv,
    PnsvDynOwned,
    PnsvMatrix,
    PnsvSafe,
    PnsvTuned,
    Ranges,
    WWT,
};

const HELP: &str = r"possible actions:
    ms <-load_safe | -safe | -tuned> <SBWT_PATH> <QUERY_PATH> <PNSV_PATH | LCS_PATH> [SCAN_BOUND]
        * matching statistics benchmark; a PNSV_PATH is expected when the first flag has prefix -load,
          otherwise, the pnsv structure will be built from the lcs read from the LCS_PATH
    ser_safe <SBWT_PATH> <LCS_PATH> [OUTPUT_PATH]
        * create and serialize a PnsvSafe structure to the file at the given path; the default for
          OUTPUT_PATH is ./unknown.pnsv_safe
    count <SBWT_PATH> <LCS_PATH> <FILE_LIST_PATH>
";

fn main() {
    env_logger::init();
    let mut args = std::env::args();
    args.next().unwrap();
    let action = args.next().unwrap_or("".to_string());

    match action.as_str() {
        "ms" => { ms_benchmark(&mut args); }
        "ser_safe" => { serialize_pnsv_safe(&mut args); },
        "count" => { count(&mut args); },
        _ => {
            println!("{}", HELP);
        }
    };
}

fn ms_benchmark(args: &mut std::env::Args) {
    let flag = args.next().expect("expected pnsv type");
    let index = load_index(args);
    let SbwtIndexVariant::SubsetMatrix(sbwt) = index;
    log::info!("count: {}", sbwt.n_sets());
    let queries = load_query(args);

    if flag == "-load_safe" {
        let mut pnsv = load_pnsv_safe(args).unwrap();
        let scan_bound = parse_scan_bound(args);
        pnsv.scan_bound = scan_bound;
        log::info!("running PnsvSafe benchmark...");
        run_ms_benchmark(&sbwt, &queries, &pnsv);
        return;
    }

    let lcs = load_lcs(args);
    let scan_bound = parse_scan_bound(args);
    match flag.as_str() {
        "-safe" => {
            let pnsv = PnsvSafe::new(&sbwt, &lcs, sbwt.k(), scan_bound);
            drop(lcs);
            log::info!("running PnsvSafe benchmark...");
            run_ms_benchmark(&sbwt, &queries, &pnsv);
        },
        _ => {
            let pnsv = PnsvTuned::new(&sbwt, &lcs, sbwt.k(), scan_bound, PnsvTuned::DEFAULT_FALLBACK_OVERLAP);
            drop(lcs);
            log::info!("running PnsvTuned benchmark...");
            run_ms_benchmark(&sbwt, &queries, &pnsv);
        }
    }
}

fn parse_scan_bound(args: &mut std::env::Args) -> usize {
    if let Some(scan_bound_arg) = args.next() {
        match scan_bound_arg.parse::<usize>() {
            Ok(value) => value.max(PnsvSafe::DEFAULT_SCAN_BOUND),
            Err(_) => PnsvSafe::DEFAULT_SCAN_BOUND
        }
    } else {
        PnsvSafe::DEFAULT_SCAN_BOUND
    }
}

fn run_ms_benchmark(sbwt: &SbwtIndex<SubsetMatrix>, queries: &[Vec<u8>], pnsv: &impl Pnsv) {
    let index = StreamingIndex {
        extend_right: sbwt,
        contract_left: pnsv,
        n: sbwt.n_sets(),
        k: sbwt.k(),
    };

    for bound in 1..=31 {
        print!("{},", bound);
        benchmark_bms_separate_queries(&index, queries, bound);
        println!();
    }
}

fn load_query(args: &mut std::env::Args) -> Vec<Vec<u8>> {
    use sbwt::SeqStream;
    let path = args.next().expect("expected query path");
    let file = std::fs::File::open(path).unwrap();
    let reader = std::io::BufReader::new(file);
    sbwt::vodbg::benchmark::read_sequences(reader)
}

fn load_index(args: &mut std::env::Args) -> SbwtIndexVariant {
    let path = args.next().expect("expected sbwt index path");
    log::info!("loading sbwt index: {}", path);
    let mut reader = std::io::BufReader::new(std::fs::File::open(path).unwrap());
    sbwt::load_sbwt_index_variant(&mut reader).unwrap()
}

fn load_lcs(args: &mut std::env::Args) -> LcsArray {
    let path = args.next().expect("expected sbwt index path");
    log::info!("loading lcs array: {}", path);
    let mut reader = std::io::BufReader::new(std::fs::File::open(path).unwrap());
    LcsArray::load(&mut reader).unwrap()
}

fn serialize_pnsv_safe(args: &mut std::env::Args) -> std::io::Result<()> {
    let index = load_index(args);
    let SbwtIndexVariant::SubsetMatrix(sbwt) = index;
    let lcs = load_lcs(args);

    let path = match args.next() {
        Some(value) => value,
        None => "./unknown.pnsv_safe".to_string()
    };

    let file = std::fs::File::create_new(path)?;
    let mut writer = std::io::BufWriter::new(file);
    let pnsv_safe = PnsvSafe::new_with_default_values(&sbwt, &lcs, sbwt.k());
    pnsv_safe.serialize(&mut writer)?;
    Ok(())
}

fn load_pnsv_safe(args: &mut std::env::Args) -> std::io::Result<PnsvSafe> {
    log::info!("[load_pnsv_safe] creating ranges...");
    let path = args.next().expect("expected pnsv_safe path");
    let file = std::fs::File::open(path)?;
    let mut reader = std::io::BufReader::new(file);
    PnsvSafe::load(&mut reader)
}

fn count(args: &mut std::env::Args) {
    let index = load_index(args);
    let SbwtIndexVariant::SubsetMatrix(mut sbwt) = index;
    let lcs = load_lcs(args);
    let pnsv = PnsvTuned::new_with_default_values(&sbwt, &lcs, sbwt.k());
    drop(lcs);
    sbwt.build_select();
    let streaming_index = StreamingIndex {
        extend_right: &sbwt,
        contract_left: &pnsv,
        n: sbwt.n_sets(),
        k: sbwt.k(),
    };
    let vodbg = sbwt::vodbg::VoDbg::new(&sbwt, &pnsv);
    log::info!("so far so good...");
    let files = read_lines(&PathBuf::from(args.next().expect("expected file list path")));
    let sequence_stream = MyMultiFileSeqReader::new(files);
    let counts = sbwt::vodbg::count::Counts::try_new_concurrent(
        sequence_stream,
        &streaming_index,
        &vodbg,
        Counts::DEFAULT_SAMPLE_DISTANCE,
        4,
        8,
        Counts::DEFAULT_BATCH_SIZE_IN_BYTES
    );
    log::info!("result is ok? {}", counts.is_some());
    if let Some(counts) = counts {
        log::info!("number of elements whose count is > u8::MAX: {}", counts.large_counts.len());

        use std::io::Write;
        let stats_file = std::fs::File::create("../stats.txt").unwrap();
        let mut stats_writer = std::io::BufWriter::new(stats_file);
        for count in counts.iter() {
            writeln!(&mut stats_writer, "{}", count).unwrap();
        }
        drop(stats_writer);
    }
}

// ==================================================================== {
// ====================================================================
// ==================================================================== stolen from sbwt-cli
use jseqio::reader::*;
use std::path::PathBuf;

fn read_lines(path: &PathBuf) -> Vec<PathBuf> {
    use std::io::BufRead;
    let file = std::fs::File::open(path).unwrap();
    let reader = std::io::BufReader::new(file);
    let mut paths = Vec::<PathBuf>::new();
    for line in reader.lines() {
        let line = line.unwrap();
        paths.push(PathBuf::from(line));
    }
    paths
}

struct MyMultiFileSeqReader {
    paths: Vec<PathBuf>,
    next_idx: usize,
    current: Option<jseqio::reader::DynamicFastXReader>,
    local_buf: Vec<u8>,
}

impl MyMultiFileSeqReader {
    fn new(paths: Vec<PathBuf>) -> Self {
        Self {
            paths,
            next_idx: 0,
            current: None,
            local_buf: vec![]
        }
    }
}

impl sbwt::SeqStream for MyMultiFileSeqReader {
    fn stream_next(&mut self) -> Option<&[u8]> {
        loop {
            if let Some(current) = &mut self.current {
                if let Some(rec) = current.read_next().unwrap() {
                    self.local_buf.clear();
                    self.local_buf.extend_from_slice(rec.seq);
                    return Some(&self.local_buf);
                } else {
                    self.current = None;
                }
            }

            // Open next file if available
            if self.next_idx < self.paths.len() {
                let path = &self.paths[self.next_idx];
                self.next_idx += 1;
                self.current = Some(jseqio::reader::DynamicFastXReader::from_file(path).unwrap());
            } else {
                return None;
            }
        }
    }
}
// ==================================================================== stolen from sbwt-cli
// ====================================================================
// ==================================================================== }

fn correctness(n: usize, first: &impl Pnsv, second: &impl Pnsv, target_length_lower: usize, target_length_upper: usize) {
    let ten_percent = n / 10;
    for i in 0..n {
        for target_length in target_length_lower..=target_length_upper {
            let first_answer = first.previous(i, target_length);
            let second_answer = second.previous(i, target_length);
            assert_eq!(first_answer, second_answer, "p; i: {}, target_length: {}", i, target_length);
            let first_answer = first.next(i, target_length);
            let second_answer = second.next(i, target_length);
            assert_eq!(first_answer, second_answer, "n; i: {}, target_length: {}", i, target_length);
        }
        if i % ten_percent == ten_percent - 1 {
            println!("{}0%", 1 + i / ten_percent);
        }
    }
}

fn sbwt_count_investigation() {
    let max_k = 4;
    let mut seqs: Vec<Vec<u8>> = vec![
        b"ACGTACG".to_vec()
    ];

    let (sbwt, lcs) = SbwtIndexBuilder::<BitPackedKmerSortingMem>::new()
        .k(max_k).build_lcs(true)
        .add_all_dummy_paths(true)
        .build_select_support(true)
        .run_from_vecs(&seqs);
    let lcs = lcs.unwrap();
    let pnsv_tuned = PnsvTuned::new_with_default_values(&sbwt, &lcs, max_k);
    let sequence_stream = sbwt::VecSeqStream::new(&seqs);
    let streaming_index = StreamingIndex {
        extend_right: &sbwt,
        contract_left: &pnsv_tuned,
        n: sbwt.n_sets(),
        k: sbwt.k()
    };

    let vodbg = VoDbg::new(&sbwt, &pnsv_tuned);
    let counts = Counts::try_new_default(sequence_stream, streaming_index, &vodbg).unwrap();

    for i in 0..sbwt.n_sets() {
        let node = sbwt::vodbg::new_node(i, i + 1, sbwt.k());
        let kmer = vodbg.get_kmer(node);
        let kmer_str = str::from_utf8(&kmer).unwrap();
        println!("{} {} {}", i, kmer_str, counts.individual_counts[i]);
    }
}
