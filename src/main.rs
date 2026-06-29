#![allow(unused)]

mod statistics;
mod vodbg;

use sbwt::{ContractLeft, ExtendRight, SbwtIndex, StreamingIndex, SubsetMatrix};
use sbwt::SbwtIndexVariant;
use sbwt::LcsArray;

use sbwt::vodbg::{pnsv::Pnsv, benchmark::*};

use sbwt::vodbg::*;
use sbwt::vodbg::pnsv::{
    self,
    ABS,
    LcsPnsvBp,
    LcsSimd,
    PnsvDyn,
    PnsvDynOwned,
    PnsvMatrix,
    PnsvMatrixSux,
    PnsvSafe,
    PnsvTuned,
    Ranges,
    WWT,
};

fn main() {
    env_logger::init();
    vodbg::operation_benchmark();
}

fn ms_benchmark() {
    let pnsv_type = std::env::args().nth(1).expect("expected pnsv type");
    log::info!("loading data: {}", std::env::args().nth(2).expect("expected sbwt path"));
    let (index, lcs) = read_index_and_lcs(2);
    let SbwtIndexVariant::SubsetMatrix(sbwt) = index;
    log::info!("count: {}", lcs.len());
    let queries = read_query(4);

    if pnsv_type == "safe" {
        let pnsv = PnsvSafe::new_with_default_values(&sbwt, &lcs, sbwt.k());
        drop(lcs);
        log::info!("running PnsvSafe benchmark...");
        run_ms_benchmark(&sbwt, &queries, &pnsv);
    } else {
        let pnsv = PnsvTuned::new_with_default_values(&sbwt, &lcs, sbwt.k());
        drop(lcs);
        log::info!("running PnsvTuned benchmark...");
        run_ms_benchmark(&sbwt, &queries, &pnsv);
    };
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
