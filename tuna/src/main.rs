extern crate primes;

mod dna_hash_table;
mod dna_read_graph;
mod read_inputs;

use dna_hash_table::DNAHashTable;

fn main() {
    // let segments : Vec<String> = read_inputs::read_fa_file("example_file2.txt");

    let fa_col_db = read_inputs::read_fa_file_to_cols("example_file2.txt");

    for s in &fa_col_db.seg_strings{
        println!("Seg {}", s);
    }
    let k : usize = 10;
    let kmer_hash_table = DNAHashTable::new(&fa_col_db.seg_strings, k);
    let kmer : String = "AGCTGGCCCT".to_string();

    match kmer_hash_table.get_kmer(&kmer) {
    	Some((kmers, kmer_indexes)) => {
    		for kmer_index in kmer_indexes {
    			println!("{:?}", kmers[kmer_index]);
    		}
    	},
    	None => {
            println!("No match");
        }
    };
}
