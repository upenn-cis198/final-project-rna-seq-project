extern crate primes;
extern crate rayon;
use rayon::prelude::*;

mod dna_hash_table;
mod dna_read_graph;
mod read_inputs;
use std::collections::HashMap;
use dna_hash_table::DNAHashTable;

use std::hash::Hash;
use std::cmp::Eq;
use std::collections::hash_map::Entry;

fn main() {
    // let segments : Vec<String> = read_inputs::read_fa_file("example_file2.txt");

    let fa_col_db = read_inputs::read_fa_file_to_cols("example_file2.txt");

    for s in &fa_col_db.seg_ids{
        println!("Seg {}", s);
    }
    let k : usize = 10;
    let thread_size : usize = 3;
    let kmer_hash_table = DNAHashTable::new(&fa_col_db.seg_strings, k);


    let mut reads : Vec<String> = Vec::new();

    let kmer : String = "AGCTGGCCCT".to_string();
    reads.push(kmer);

    let mut reads2 : Vec<String> = read_inputs::read_fq_fasta_file("example_file.txt");

    let mut partitions :Vec<Vec<String>> = Vec::new();
    let ratio : usize = (reads2.len() / thread_size) as usize;
    let mut i = 0;
    for chunk in reads2.chunks_mut(ratio) {
        for c in chunk.iter() {
            println!("CC{} {}", i, c);
        }
        i = i + 1;

        partitions.push(chunk.to_vec());
    } 

    i = 0;
    let partitions_map = partitions.par_iter().map(|ref chunk| dna_hash_table::get_segments(&kmer_hash_table, &chunk)); 
    let mut comp_result : Vec<HashMap<i32, i32>> = partitions_map.collect();

    // partitions_map.collect_into(comp_result);


    let mut global_map : HashMap<i32, i32> = HashMap::new();

    println!("Length {}", comp_result.len());
    for p in comp_result.iter() {
        for (k, v) in p.iter() {
            println!("Partition {} has {} {}", i, k, v);

            let cnt = global_map.entry(*k).or_insert(0);
            *cnt += *v;
        }
        i = i + 1;
    }




    // let segment_index_counts : HashMap<i32, i32> = dna_hash_table::get_segments(&kmer_hash_table, &reads);


    // HashMap::new();


    // match kmer_hash_table.get_kmer(&kmer) {
    // 	Some((kmers, kmer_indexes)) => {
    // 		for kmer_index in kmer_indexes {
    // 			println!("{:?}", &kmers[kmer_index]);
    //             let count_instance = &kmers[kmer_index];

    //             let count = segment_index_counts.entry(count_instance.segment_index as i32).or_insert(0);
    //             *count += 1;
    // 		}
    // 	},
    // 	None => {
    //         println!("No match");
    //     }
    // };


    read_inputs::write_output3("out.txt", global_map, fa_col_db.seg_ids);
}
