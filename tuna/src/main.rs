extern crate primes;

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
    let kmer_hash_table = DNAHashTable::new(&fa_col_db.seg_strings, k);
    let kmer : String = "AGCTGGCCCT".to_string();

    let mut segment_index_counts : HashMap<i32, i32> = HashMap::new();


    match kmer_hash_table.get_kmer(&kmer) {
    	Some((kmers, kmer_indexes)) => {
    		for kmer_index in kmer_indexes {
    			println!("{:?}", &kmers[kmer_index]);
                let count_instance = &kmers[kmer_index];

                let count = segment_index_counts.entry(count_instance.segment_index as i32).or_insert(0);
                *count += 1;
    		}
    	},
    	None => {
            println!("No match");
        }
    };


    read_inputs::write_output3("out.txt", segment_index_counts, fa_col_db.seg_ids);
}
