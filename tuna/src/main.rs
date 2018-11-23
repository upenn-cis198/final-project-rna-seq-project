extern crate primes;

mod dna_hash_table;

use dna_hash_table::DNAHashTable;

fn main() {
    let mut segments : Vec<String> = Vec::<String>::new();
    let k : usize = 10;
    segments.push("ATGATAGATAGACATACGTACGATCG".to_string());
    let kmer_hash_table = DNAHashTable::new(&segments, k);
    let kmer : String = "ATGATAGATA".to_string();
    kmer_hash_table.get_kmer(&kmer);
}
