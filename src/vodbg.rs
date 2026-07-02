use super::*;

pub fn operation_benchmark(args: &mut std::env::Args) {
    let index = load_index(args);
    let SbwtIndexVariant::SubsetMatrix(mut sbwt) = index;
    let lcs = load_lcs(args);
    log::info!("count: {}", lcs.len());
    let mut queries = load_query(args);
    let query = queries.swap_remove(0);

    sbwt.build_select();
    let pnsv = PnsvTuned::new(&sbwt, &lcs, sbwt.k(), 8, 0);
    let vodbg = VoDbg::new(&sbwt, &pnsv);

    // log::info!();

    let nodes: Vec<Node> = vec![];
    let mut kmer_buffer = vec![];

    let max_k = sbwt.k();
    const MIN_K: usize = 3;
    for current_k in MIN_K..=max_k {
        let mut count;
        let mut nanos;
        let q = query.windows(current_k);

        count = 0;
        let total = run_operation(q.clone(), |kmer| {
            let node = vodbg.get_node(kmer).expect("Should exist.");
            count += 1;
        });
        nanos = total / count as f64;
        println!("get_node,{},{}", current_k, nanos);

        count = 0;
        let total = run_operation(q.clone(), |kmer| {
            let node = vodbg.get_node(kmer).expect("Should exist.");
            let indegree = vodbg.indegree(node);
            count += 1;
        });
        nanos = total / count as f64;
        println!("indegree,{},{}", current_k, nanos);

        count = 0;
        let total = run_operation(q.clone(), |kmer| {
            let node = vodbg.get_node(kmer).expect("Should exist.");
            let outdegree = vodbg.outdegree(node);
            count += 1;
        });
        nanos = total / count as f64;
        println!("outdegree,{},{}", current_k, nanos);

        //     // follow_inedge
        //     let dbg_from_inedge = dbg.follow_inedge(dbg_node, i).expect("Should exist.");
        //     let vodbg_from_inedge = vodbg.follow_inedge(vodbg_node, i).expect("Should exist.");
        //     assert_eq!(dbg_in_neighbour_kmer, dbg.get_kmer(dbg_from_inedge));
        //     assert_eq!(vodbg_in_neighbour_kmer, vodbg.get_kmer(vodbg_from_inedge));

        count = 0;
        let total = run_operation(q.clone(), |kmer| {
            let node = vodbg.get_node(kmer).expect("Should exist.");
            let contracted = vodbg.contract_left(node, node.k - 1);
            count += 1;
        });
        nanos = total / count as f64;
        println!("contract_left,{},{}", current_k, nanos);

        count = 0;
        let total = run_operation(q.clone(), |kmer| {
            let node = vodbg.get_node(kmer).expect("Should exist.");
            let contracted = vodbg.contract_left(node, node.k - 1);
            let extended = vodbg.extend_left_with_character(contracted, kmer[0], &mut kmer_buffer);
            count += 1;
        });
        nanos = total / count as f64;
        println!("extend_left_character,{},{}", current_k, nanos);

        count = 0;
        let total = run_operation(q.clone(), |kmer| {
            let node = vodbg.get_node(kmer).expect("Should exist.");
            let contracted = vodbg.contract_right(node, node.k - 1);
            count += 1;
        });
        nanos = total / count as f64;
        println!("contract_right,{},{}", current_k, nanos);

        count = 0;
        let total = run_operation(q.clone(), |kmer| {
            let node = vodbg.get_node(kmer).expect("Should exist.");
            let contracted = vodbg.contract_right(node, node.k - 1);
            let extended = vodbg.extend_right(node, kmer[kmer.len() - 1]);
            count += 1;
        });
        nanos = total / count as f64;
        println!("extend_right,{},{}", current_k, nanos);
    }
}

fn run_operation<'a, I, F>(kmers: I, mut operation: F) -> f64
where 
    I: Iterator<Item = &'a [u8]>,
    F: FnMut(&[u8])
{
    let start_time = std::time::Instant::now();
    for kmer in kmers {
        operation(kmer);
    }
    let end_time = std::time::Instant::now();
    (end_time - start_time).as_nanos() as f64
}

