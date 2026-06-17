#![allow(unused)]

use sbwt::{ContractLeft, ExtendRight, StreamingIndex};
use sbwt::SbwtIndexVariant;
use sbwt::LcsArray;

use sbwt::vodbg::{pnsv::Pnsv, benchmark::*};

use sbwt::vodbg::pnsv::{
    self,
    ABS,
    LcsPnsvBp,
    LcsSimd,
    PnsvDyn,
    PnsvDynOwned,
    PnsvMatrix,
    PnsvMatrixSux,
    Ranges,
    WWT,
};

fn main() {
    env_logger::init();
    comparison();
}

fn comparison() {
    println!("loading data...");
    let (index, lcs) = read_index_and_lcs(1);
    let SbwtIndexVariant::SubsetMatrix(sbwt) = index;
    println!("lcs.len: {}", lcs.len());
    let queries = read_query(3);

    // println!("creating standard bp structure...");
    // let bp = LcsPnsvBp::new(&lcs, 2048);

    let pnsv_dyn = pnsv::pnsv_simd_fallback_matrix(&sbwt, &lcs);
    drop(lcs);

    print!("{}", 1);
    for pnsv in &pnsv_dyn.structures {
        print!(":{}", pnsv.max_target());
    }
    println!("(:30)");

    // println!("creating wavelet...");
    // let wavelet = WWT::from_iterator(iterator, 7, 4);

    // println!("creating matrix sux...");
    // let matrix = PnsvMatrixSux::from_iterator(iterator, lcs.len(), 8, 10);

    let pnsv_dyn_index = StreamingIndex {
        extend_right: &sbwt,
        contract_left: &pnsv_dyn,
        // contract_left: &bp,
        n: sbwt.n_sets(),
        k: sbwt.k(),
    };

    println!("running benchmarks...");

    let lower = 1;
    let upper = 31;

    for bound in lower..upper {
        print!("dyn,{},", bound);
        benchmark_bms_separate_queries(&pnsv_dyn_index, &queries, bound);
        println!();
    }
}

fn analyse_range_lengths(argument_start: usize) {
    let mut args = std::env::args().skip(argument_start);
    let lcs_path = args.next().expect("expected lcs index path");

    println!("reading data...");
    let mut lcs_reader = std::io::BufReader::new(std::fs::File::open(lcs_path).unwrap());
    let lcs = LcsArray::load(&mut lcs_reader).unwrap();

    println!("count: {}", lcs.len());

    let k: usize = 31;
    let mut range_counts = vec![1_usize; k];

    println!("counting...");
    for i in 1..lcs.len() {
        let item = lcs.access(i);
        #[allow(clippy::needless_range_loop)]
        for target_length in 1..k {
            if item < target_length {
                range_counts[target_length] += 1;
            }
        }
    }

    let mut previous_average_length = 1.0_f64;
    let mut average_length;
    let mut ratio;
    #[allow(clippy::needless_range_loop)]
    for target_length in 1..k {
        average_length = lcs.len() as f64 / range_counts[target_length] as f64;
        ratio = previous_average_length / average_length;
        previous_average_length = average_length;
        println!("tl: {} | avg: {:.3} | ratio: {:.3}", target_length, average_length, ratio);
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

fn simd_scan_compare() {
    let (index, lcs) = read_index_and_lcs(1);
    let SbwtIndexVariant::SubsetMatrix(sbwt) = index;
    let queries = read_query(3);

    let iterator = (0..lcs.len()).map(|index| lcs.access(index) as u8);
    let lcs_simd = LcsSimd::from_iterator(iterator, lcs.len());

    let lcs_index = StreamingIndex {
        extend_right: &sbwt,
        contract_left: &lcs,
        n: sbwt.n_sets(),
        k: sbwt.k(),
    };

    let lcs_simd_index = StreamingIndex {
        extend_right: &sbwt,
        contract_left: &lcs_simd,
        n: sbwt.n_sets(),
        k: sbwt.k(),
    };

    let bound = 10;

    print!("scan,{},", bound);
    benchmark_bms_separate_queries(&lcs_index, &queries, bound);
    println!();

    print!("simd,{},", bound);
    benchmark_bms_separate_queries(&lcs_simd_index, &queries, bound);
    println!();
}

fn statistics_impl_pnsv(pnsv: &impl Pnsv, n: usize, step: usize, target_length: usize) -> (f64, f64) {
    let start_time = std::time::Instant::now();
    for i in (0..n).step_by(step) {
        let _ = std::hint::black_box(pnsv.previous(i, target_length));
    }
    let end_time = std::time::Instant::now();
    let nanos_per_previous = (end_time - start_time).as_nanos() as f64 / n as f64;

    let start_time = std::time::Instant::now();
    for i in (0..n).step_by(step) {
        let _ = std::hint::black_box(pnsv.next(i, target_length));
    }
    let end_time = std::time::Instant::now();
    let nanos_per_next = (end_time - start_time).as_nanos() as f64 / n as f64;

    (nanos_per_previous, nanos_per_next)
}

fn statistics_lcs_simd(lcs_simd: &LcsSimd, target_length: usize, bound: usize) -> (f64, f64, f64, f64) {
    let n = lcs_simd.n;
    let target_length = target_length as u8;

    let mut successful_previous = 0;
    let start_time = std::time::Instant::now();
    for i in 0..n {
        successful_previous += if lcs_simd.scan_left_bounded(i, target_length, bound).is_ok() {
            1
        } else {
            0
        };
    }
    let end_time = std::time::Instant::now();
    let nanos_per_previous = (end_time - start_time).as_nanos() as f64 / n as f64;
    let percentage_previous = successful_previous as f64 / n as f64;

    let mut successful_next = 0;
    let start_time = std::time::Instant::now();
    for i in 0..n {
        successful_next += if lcs_simd.scan_left_bounded(i, target_length, bound).is_ok() {
            1
        } else {
            0
        };
    }
    let end_time = std::time::Instant::now();
    let nanos_per_next = (end_time - start_time).as_nanos() as f64 / n as f64;
    let percentage_next = successful_next as f64 / n as f64;

    (percentage_previous, percentage_next, nanos_per_next, nanos_per_previous)
}

fn statistics_pnsv_matrix(matrix: &PnsvMatrix, target_length: usize) -> (f64, f64) {
    let n = matrix.width;

    let start_time = std::time::Instant::now();
    for i in 0..n {
        let _ = std::hint::black_box(matrix.previous(i, target_length));
    }
    let end_time = std::time::Instant::now();
    let nanos_per_previous = (end_time - start_time).as_nanos() as f64 / n as f64;

    let start_time = std::time::Instant::now();
    for i in 0..n {
        let _ = std::hint::black_box(matrix.next(i, target_length));
    }
    let end_time = std::time::Instant::now();
    let nanos_per_next = (end_time - start_time).as_nanos() as f64 / n as f64;

    (nanos_per_previous, nanos_per_next)
}

fn simd_bounded_scan_time() {
    let (index, lcs) = read_index_and_lcs(1);
    let SbwtIndexVariant::SubsetMatrix(sbwt) = index;

    let iterator = (0..lcs.len()).map(|index| lcs.access(index) as u8);

    let lower_bound = 8;
    let upper_bound = 10;

    println!("creating lcs_simd...");
    let lcs_simd = LcsSimd::from_iterator(iterator.clone(), lcs.len());

    println!("creating matrix...");
    let matrix = PnsvMatrix::from_iterator(iterator, lcs.len(), lower_bound, upper_bound);

    println!("timing...");
    let item_bound: usize = 1000;
    let word_bound = item_bound.div_ceil(LcsSimd::LANES);

    for target_length in lower_bound..=upper_bound {
        let (
            percentage_previous,
            percentage_next,
            nanos_per_next_scan,
            nanos_per_previous_scan,
        ) = statistics_lcs_simd(&lcs_simd, target_length, word_bound);

        let (
            nanos_per_previous_matrix,
            nanos_per_next_matrix
        ) = statistics_pnsv_matrix(&matrix, target_length);

        println!("target_length: {}", target_length);
        println!(
            "%previous: {:.3} <> t_scan/t_bitvector: {:.3} ({:.3}/{:.3})",
            percentage_previous,
            nanos_per_previous_scan / nanos_per_previous_matrix,
            nanos_per_previous_scan,
            nanos_per_previous_matrix
        );

        println!(
            "%next: {:.3} <> t_scan/t_bitvector: {:.3} ({:.3}/{:.3})",
            percentage_next,
            nanos_per_next_scan / nanos_per_next_matrix,
            nanos_per_next_scan,
            nanos_per_next_matrix
        );
        println!();
    }
}

fn statistics_lcs_simd_with_matrix_fallback(lcs_simd: &LcsSimd, matrix: &PnsvMatrix, target_length: usize, bound: usize) -> (f64, f64) {
    let n = lcs_simd.n;
    let target_length_u8 = target_length as u8;

    let start_time = std::time::Instant::now();
    for i in 0..n {
        if lcs_simd.scan_left_bounded(i, target_length_u8, bound).is_err() {
            let i = i.saturating_sub(LcsSimd::LANES * bound);
            matrix.previous(i, target_length);
        }
    }
    let end_time = std::time::Instant::now();
    let nanos_per_previous = (end_time - start_time).as_nanos() as f64 / n as f64;

    let start_time = std::time::Instant::now();
    for i in 0..n {
        if lcs_simd.scan_right_bounded(i, target_length_u8, bound).is_err() {
            let i = i + LcsSimd::LANES * bound;
            matrix.next(i, target_length);
        }
    }
    let end_time = std::time::Instant::now();
    let nanos_per_next = (end_time - start_time).as_nanos() as f64 / n as f64;

    (nanos_per_next, nanos_per_previous)
}

fn statistics_augmented_bounded_scan(abs: &ABS, target_length: usize) -> (f64, f64) {
    let n = abs.lcs_simd.len();
    let target_length_u8 = target_length as u8;

    let start_time = std::time::Instant::now();
    for i in 0..n {
        let _ = abs.previous(i, target_length);
    }
    let end_time = std::time::Instant::now();
    let nanos_per_previous = (end_time - start_time).as_nanos() as f64 / n as f64;

    let start_time = std::time::Instant::now();
    for i in 0..n {
        let _ = abs.next(i, target_length);
    }
    let end_time = std::time::Instant::now();
    let nanos_per_next = (end_time - start_time).as_nanos() as f64 / n as f64;

    (nanos_per_next, nanos_per_previous)
}

fn simd_bounded_scan_with_fallback_time() {
    let (index, lcs) = read_index_and_lcs(1);
    let SbwtIndexVariant::SubsetMatrix(sbwt) = index;

    let iterator = (0..lcs.len()).map(|index| lcs.access(index) as u8);

    let target_length_lower = 8;
    let target_length_upper = 10;

    println!("creating lcs_simd...");
    let lcs_simd = LcsSimd::from_iterator(iterator.clone(), lcs.len());

    println!("creating matrix...");
    let matrix = PnsvMatrix::from_iterator(iterator.clone(), lcs.len(), target_length_lower, target_length_upper);

    let item_bound: usize = 256;
    let word_bound = item_bound.div_ceil(LcsSimd::LANES);

    println!("creating augmented bounded scan...");
    let abs = ABS::from_iterator(lcs_simd.clone(), iterator, word_bound, target_length_lower, target_length_upper);

    println!("timing...");
    for target_length in target_length_lower..=target_length_upper {
        // let (
        //     nanos_per_previous_abs,
        //     nanos_per_next_abs
        // ) = statistics_augmented_bounded_scan(&abs, target_length);

        let (
            nanos_per_previous_matrix_fallback,
            nanos_per_next_matrix_fallback,
        ) = statistics_lcs_simd_with_matrix_fallback(&lcs_simd, &matrix, target_length, word_bound);

        let (
            nanos_per_previous_matrix,
            nanos_per_next_matrix
        ) = statistics_pnsv_matrix(&matrix, target_length);

        println!("target_length: {}", target_length);
        // println!(
        //     "augmented bounded scan    | previous: {:.3} | next: {:.3}",
        //     nanos_per_previous_abs,
        //     nanos_per_next_abs
        // );
        println!(
            "scan with matrix fallback | previous: {:.3} | next: {:.3}",
            nanos_per_previous_matrix_fallback,
            nanos_per_next_matrix_fallback
        );
        println!(
            "matrix only               | previous: {:.3} | next: {:.3}",
            nanos_per_previous_matrix,
            nanos_per_next_matrix
        );
        println!();
    }
}
