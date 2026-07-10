use super::*;

pub fn analyse_range_lengths(argument_start: usize) {
    let mut args = std::env::args().skip(argument_start);
    let lcs_path = args.next().expect("expected lcs index path");

    log::info!("reading data...");
    let mut lcs_reader = std::io::BufReader::new(std::fs::File::open(lcs_path).unwrap());
    let lcs = LcsArray::load(&mut lcs_reader).unwrap();

    let count = lcs.len();
    let two_percent = count / 50;
    let mut bound = two_percent;
    let mut progress = 0;
    log::info!("count: {}", count);

    let k: usize = 31;
    let mut range_counts = vec![1_usize; k];
    let mut previous_position = vec![0_usize; k];
    let mut max_range = vec![0_usize; k];

    log::info!("counting...");
    for i in 1..lcs.len() {
        let item = lcs.access(i);
        #[allow(clippy::needless_range_loop)]
        for target_length in 1..k {
            if item < target_length {
                range_counts[target_length] += 1;
                max_range[target_length] = max_range[target_length].max(i - previous_position[target_length]);
                previous_position[target_length] = i;
            }
        }
        if i > bound {
            bound += two_percent;
            progress += 2;
            log::info!("scanning... {}%", progress);
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
        log::info!(
            "tl: {} | avg: {:.3} | ratio: {:.3} | ml: {}",
            target_length,
            average_length,
            ratio,
            max_range[target_length]
        );
    }
}

pub fn simd_scan_compare(args: &mut std::env::Args) {
    let index = load_index(args);
    let SbwtIndexVariant::SubsetMatrix(sbwt) = index;
    let lcs = load_lcs(args);
    let queries = load_query(args);

    let iterator = (0..lcs.len()).map(|index| lcs.access(index) as u8);
    let lcs_simd = LcsSimd::from_iterator(iterator, lcs.len(), sbwt.k());

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

pub fn statistics_impl_pnsv(pnsv: &impl Pnsv, n: usize, step: usize, target_length: usize) -> (f64, f64) {
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

pub fn statistics_lcs_simd(lcs_simd: &LcsSimd, target_length: usize, bound: usize) -> (f64, f64, f64, f64) {
    let n = lcs_simd.len();

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

pub fn statistics_pnsv_matrix(matrix: &PnsvMatrix, target_length: usize) -> (f64, f64) {
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

pub fn simd_bounded_scan_time(args: &mut std::env::Args) {
    let index = load_index(args);
    let SbwtIndexVariant::SubsetMatrix(sbwt) = index;
    let lcs = load_lcs(args);

    let iterator = (0..lcs.len()).map(|index| lcs.access(index) as u8);

    let lower_bound = 8;
    let upper_bound = 10;

    println!("creating lcs_simd...");
    let lcs_simd = LcsSimd::from_iterator(iterator.clone(), lcs.len(), sbwt.k());

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

pub fn statistics_lcs_simd_with_matrix_fallback(lcs_simd: &impl Scan, matrix: &PnsvMatrix, target_length: usize, bound: usize) -> (f64, f64) {
    let n = lcs_simd.len();

    let start_time = std::time::Instant::now();
    for i in 0..n {
        if lcs_simd.scan_left_bounded(i, target_length, bound).is_err() {
            let i = i.saturating_sub(LcsSimd::LANES * bound);
            matrix.previous(i, target_length);
        }
    }
    let end_time = std::time::Instant::now();
    let nanos_per_previous = (end_time - start_time).as_nanos() as f64 / n as f64;

    let start_time = std::time::Instant::now();
    for i in 0..n {
        if lcs_simd.scan_right_bounded(i, target_length, bound).is_err() {
            let i = i + LcsSimd::LANES * bound;
            matrix.next(i, target_length);
        }
    }
    let end_time = std::time::Instant::now();
    let nanos_per_next = (end_time - start_time).as_nanos() as f64 / n as f64;

    (nanos_per_next, nanos_per_previous)
}

pub fn statistics_augmented_bounded_scan(abs: &ABS, target_length: usize) -> (f64, f64) {
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

pub fn simd_bounded_scan_with_fallback_time(args: &mut std::env::Args) {
    let index = load_index(args);
    let SbwtIndexVariant::SubsetMatrix(sbwt) = index;
    let lcs = load_lcs(args);

    let iterator = (0..lcs.len()).map(|index| lcs.access(index) as u8);

    let target_length_lower = 8;
    let target_length_upper = 10;

    println!("creating lcs_simd...");
    let lcs_simd = pnsv::scan::LcsSimd8x32::from_iterator(iterator.clone(), lcs.len(), sbwt.k());

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
