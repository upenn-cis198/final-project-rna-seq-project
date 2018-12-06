extern crate multiset;

use primes::PrimeSet;
use self::multiset::HashMultiSet;
use std::cmp;

pub struct DNAHashTable<'a> {
	pub hash_table : Vec<Vec<Kmer>>,
	segments : &'a Vec<String>,
	pub size : usize,
	pub k : usize,
	j : usize,
	segment_iter_index : usize,
	kmer_iter_index : usize,
}

impl<'a> DNAHashTable<'a> {
	
	//k is the size of the k-mers to be hashed, and j is the maximum index of the k-mer that can be used for hashing before overflow
	pub fn new(segments : &Vec<String>, k : usize) -> DNAHashTable {
		let j : usize = DNAHashTable::get_max_j(k);
		let size : usize = DNAHashTable::get_table_size(segments, k);
		let mut hash_table : Vec<Vec<Kmer>> = vec![Vec::<Kmer>::new(); size];

		let mut creation_time : usize = 0;
		for segment_index in 0..segments.len() {
			let segment : &String = &segments[segment_index];
			for i in 0..(segment.len() - k + 1) {
				let hash_value : usize = DNAHashTable::hash_function(&segment[i..(i+k)], j, size);
				hash_table[hash_value].push(Kmer {
					segment_index : segment_index,
					position : i,
					creation_time : creation_time,
				});
				creation_time += 1;
			}
		}

		DNAHashTable {
			hash_table : hash_table,
			segments : segments,
			size : size,
			k : k,
			j : j,
			segment_iter_index : 0,
			kmer_iter_index : 0,
		}
	}

	pub fn get_kmer(&self, kmer_string : &str) -> Option<(&Vec<Kmer>, Vec<usize>)> {
		if kmer_string.len() == self.k {
			let hash_value : usize = DNAHashTable::hash_function(kmer_string, self.j, self.size);
			let hash_element : &Vec<Kmer> = &self.hash_table[hash_value];
			match hash_element.len() {
				0 => None,
				_ => {
					let mut kmer_indexes : Vec<usize> = Vec::<usize>::new();
					kmer_indexes.reserve(hash_element.len());
					
					for i in 0..hash_element.len() {
						let kmer : &Kmer = &hash_element[i];
						if self.segments[kmer.segment_index][kmer.position..(kmer.position + self.k)] == *kmer_string {
							kmer_indexes.push(i);
						}
					}
					match kmer_indexes.len() {
						0 => None,
						_ => Some((&hash_element, kmer_indexes)),
					}
				}
			}
		} else {
			None
		}
	}

	//Get the table size
	pub fn get_table_size(segments : &Vec<String>, k : usize) -> usize {
		let mut size : usize = 0;
		let mut pset = PrimeSet::new();
		
		for segment in segments {
			size += segment.len() - k + 1;
		}
		let new_size = (size * 13) / 10; //Recommended to make size 1.3 times number of keys
		let (_index, prime_size) = pset.find(new_size as u64);
		prime_size as usize
	}

	fn get_max_j(k : usize) -> usize {
		let max_value : usize = usize::max_value();
		let j : usize = DNAHashTable::integer_log_base_4(max_value);
		cmp::min(j, k)

	}

	fn integer_log_base_4(mut value : usize) -> usize {
		let mut i = 0;
		while value > 0 {
			value /= 4;
			i += 1;
		}
		i - 1
	}

	fn dna_to_int(dna_letter : char) -> usize {
		match dna_letter {
			'A' => 0,
			'C' => 1,
			'G' => 2,
			'T' => 3,
			_ => 0,
		}
	}

	fn hash_function(kmer : &str, j : usize, size : usize) -> usize {
		let mut hash_value : usize = 0;
		let mut i = 0;
		for dna_letter in kmer[0..j].chars() {
			hash_value += DNAHashTable::dna_to_int(dna_letter) * 4_usize.pow(i);
			i += 1;
		}
		hash_value % size
	}

	pub fn get_most_likely_position(&self, segments : &Vec<String>, lmer : &str) -> Option<(usize, usize)> {
		let mut kmer_counts = HashMultiSet::<(usize, usize)>::new();
		for i in 0..(lmer.len() - self.k + 1) {
			println!("{}", i);
			let kmer_string : &str = &lmer[i..(i+self.k)];
			match DNAHashTable::get_kmer(&self, kmer_string) {
				Some((kmers, kmer_indexes)) => {
					for kmer_index in kmer_indexes {
						let kmer : &Kmer = &kmers[kmer_index];
						if !(kmer.position < i || (kmer.position - i + lmer.len()) > segments[kmer.segment_index].len()) {
							kmer_counts.insert((kmer.segment_index, kmer.position - i));
						}
					}
				},
				None => {},
			}
		}
		let mut max_count : usize = 0;
		let mut max_position : (usize, usize) = (0, 0);
		for element in kmer_counts.distinct_elements() {
			if kmer_counts.count_of(*element) > max_count {
				max_count = kmer_counts.count_of(*element);
				max_position = *element;
			}
		}
		match max_count {
			0 => None,
			_ => Some(max_position)
		}
	}
}

//impl<'a> Iterator for DNAHashTable<'a> {
//	type Item = &'a Kmer;

//	fn next(&mut self) -> Option<&'a Kmer> {
//		if self.segment_iter_index >= self.size {
//			None
//		} else {
//			let next_kmer : &'a Kmer = &self.hash_table[self.segment_iter_index][self.kmer_iter_index];
//			self.segment_iter_index += 1;
//			self.kmer_iter_index += 1;
//			Some(next_kmer)
//		}
//	}
//}

//CONSIDER USING A TRAIT HERE FOR DIFFERING KMERS ON CREATION TIME

//A locus is a location in the genome, which we represent by the segment that the k-mer mapped to and the location on the segment
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Kmer {
	pub segment_index : usize,
	pub position : usize,
	pub creation_time : usize,
}

#[cfg(test)]
mod tests {
	use super::*;

	#[test]
    fn test_hash_function_1() {
        let kmer_1 : &str = "ATG";
	    assert_eq!((0*1 + 3*4 + 2*16) % 3, DNAHashTable::hash_function(kmer_1, kmer_1.len(), 3));
    }

    #[test]
    fn test_hash_function_2() {
        let kmer_1 : &str = "TGAC";
	    assert_eq!((3*1 + 2*4 + 0*16 + 1*64) % 5, DNAHashTable::hash_function(kmer_1, kmer_1.len(), 5));
    }

    #[test]
    fn test_hash_function_3() {
        let kmer_1 : &str = "AAGTCT";
	    assert_eq!((0*1 + 0*4 + 2*16 + 3*64 + 1*256 + 3*1024) % 7, DNAHashTable::hash_function(kmer_1, kmer_1.len(), 7));
    }

    #[test]
    fn test_hash_table_create_1() {
        let mut segments : Vec<String> = Vec::<String>::new();
	    let k : usize = 1;
	    segments.push("ATGTG".to_string());
	    let kmer_hash_table = DNAHashTable::new(&segments, k);
	    let kmer_1 : Kmer = Kmer {
	    	segment_index : 0,
	    	position : 0,
	    	creation_time : 0
	    };
	    let kmer_2 : Kmer = Kmer {
	    	segment_index : 0,
	    	position : 1,
	    	creation_time : 1
	    };
	    let kmer_3 : Kmer = Kmer {
	    	segment_index : 0,
	    	position : 2,
	    	creation_time : 2
	    };
	    let kmer_4 : Kmer = Kmer {
	    	segment_index : 0,
	    	position : 3,
	    	creation_time : 3
	    };
	    let kmer_5 : Kmer = Kmer {
	    	segment_index : 0,
	    	position : 4,
	    	creation_time : 4
	    };
	    assert_eq!(7, kmer_hash_table.size);
	    assert_eq!(vec![kmer_1], kmer_hash_table.hash_table[0]);
	    assert_eq!(Vec::<Kmer>::new(), kmer_hash_table.hash_table[1]);
	    assert_eq!(vec![kmer_3, kmer_5], kmer_hash_table.hash_table[2]);
	    assert_eq!(vec![kmer_2, kmer_4], kmer_hash_table.hash_table[3]);
    }

    #[test]
    fn test_hash_table_create_2() {
        let mut segments : Vec<String> = Vec::<String>::new();
	    let k : usize = 3;
	    segments.push("ATGTG".to_string());
	    let kmer_hash_table = DNAHashTable::new(&segments, k);
	    let kmer_1 : Kmer = Kmer {
	    	segment_index : 0,
	    	position : 0,
	    	creation_time : 0
	    };
	    let kmer_2 : Kmer = Kmer {
	    	segment_index : 0,
	    	position : 1,
	    	creation_time : 1
	    };
	    let kmer_3 : Kmer = Kmer {
	    	segment_index : 0,
	    	position : 2,
	    	creation_time : 2
	    };
	    assert_eq!(3, kmer_hash_table.size);
	    assert_eq!(Vec::<Kmer>::new(), kmer_hash_table.hash_table[0]);
	    assert_eq!(vec![kmer_3], kmer_hash_table.hash_table[1]);
	    assert_eq!(vec![kmer_1, kmer_2], kmer_hash_table.hash_table[2]);
    }

    #[test]
    fn test_hash_table_create_3() {
        let mut segments : Vec<String> = Vec::<String>::new();
	    let k : usize = 3;
	    segments.push("ATGTG".to_string());
	    segments.push("GTGC".to_string());
	    let kmer_hash_table = DNAHashTable::new(&segments, k);
	    let kmer_1 : Kmer = Kmer {
	    	segment_index : 0,
	    	position : 0,
	    	creation_time : 0
	    };
	    let kmer_2 : Kmer = Kmer {
	    	segment_index : 0,
	    	position : 1,
	    	creation_time : 1
	    };
	    let kmer_3 : Kmer = Kmer {
	    	segment_index : 0,
	    	position : 2,
	    	creation_time : 2
	    };
	    let kmer_4 : Kmer = Kmer {
	    	segment_index : 1,
	    	position : 0,
	    	creation_time : 3
	    };
	    let kmer_5 : Kmer = Kmer {
	    	segment_index : 1,
	    	position : 1,
	    	creation_time : 4
	    };
	    assert_eq!(7, kmer_hash_table.size);
	    assert_eq!(Vec::<Kmer>::new(), kmer_hash_table.hash_table[0]);
	    assert_eq!(Vec::<Kmer>::new(), kmer_hash_table.hash_table[1]);
	    assert_eq!(vec![kmer_1], kmer_hash_table.hash_table[2]);
	    assert_eq!(vec![kmer_2], kmer_hash_table.hash_table[3]);
	    assert_eq!(vec![kmer_3, kmer_4], kmer_hash_table.hash_table[4]);
	    assert_eq!(Vec::<Kmer>::new(), kmer_hash_table.hash_table[5]);
	    assert_eq!(vec![kmer_5], kmer_hash_table.hash_table[6]);
    }

    #[test]
    fn test_hash_table_get_1() {
        let mut segments : Vec<String> = Vec::<String>::new();
	    let k : usize = 1;
	    segments.push("ATGTG".to_string());
	    let kmer_hash_table = DNAHashTable::new(&segments, k);
	    let kmer_1 : Kmer = Kmer {
	    	segment_index : 0,
	    	position : 0,
	    	creation_time : 0
	    };
	    let kmer_2 : Kmer = Kmer {
	    	segment_index : 0,
	    	position : 1,
	    	creation_time : 1
	    };
	    let kmer_3 : Kmer = Kmer {
	    	segment_index : 0,
	    	position : 2,
	    	creation_time : 2
	    };
	    let kmer_4 : Kmer = Kmer {
	    	segment_index : 0,
	    	position : 3,
	    	creation_time : 3
	    };
	    let kmer_5 : Kmer = Kmer {
	    	segment_index : 0,
	    	position : 4,
	    	creation_time : 4
	    };
	    assert_eq!(Some((&vec![kmer_1], vec![0])), kmer_hash_table.get_kmer("A"));
	    assert_eq!(Some((&vec![kmer_3, kmer_5], vec![0, 1])), kmer_hash_table.get_kmer("G"));
	    assert_eq!(Some((&vec![kmer_2, kmer_4], vec![0, 1])), kmer_hash_table.get_kmer("T"));
    }

    #[test]
    fn test_hash_table_get_2() {
        let mut segments : Vec<String> = Vec::<String>::new();
	    let k : usize = 3;
	    segments.push("ATGTG".to_string());
	    let kmer_hash_table = DNAHashTable::new(&segments, k);
	    let kmer_1 : Kmer = Kmer {
	    	segment_index : 0,
	    	position : 0,
	    	creation_time : 0
	    };
	    let kmer_2 : Kmer = Kmer {
	    	segment_index : 0,
	    	position : 1,
	    	creation_time : 1
	    };
	    let kmer_3 : Kmer = Kmer {
	    	segment_index : 0,
	    	position : 2,
	    	creation_time : 2
	    };
	    let vec_1 : &Vec<Kmer> = &vec![kmer_1, kmer_2];
	    assert_eq!(None, kmer_hash_table.get_kmer("TAT"));
	    assert_eq!(None, kmer_hash_table.get_kmer("TATT"));
	    assert_eq!(Some((&vec![kmer_3], vec![0])), kmer_hash_table.get_kmer("GTG"));
	    assert_eq!(Some((vec_1, vec![0])), kmer_hash_table.get_kmer("ATG"));
	    assert_eq!(Some((vec_1, vec![1])), kmer_hash_table.get_kmer("TGT"));
    }

    #[test]
    fn test_hash_table_get_3() {
        let mut segments : Vec<String> = Vec::<String>::new();
	    let k : usize = 3;
	    segments.push("ATGTG".to_string());
	    segments.push("GTGC".to_string());
	    let kmer_hash_table = DNAHashTable::new(&segments, k);
	    let kmer_1 : Kmer = Kmer {
	    	segment_index : 0,
	    	position : 0,
	    	creation_time : 0
	    };
	    let kmer_2 : Kmer = Kmer {
	    	segment_index : 0,
	    	position : 1,
	    	creation_time : 1
	    };
	    let kmer_3 : Kmer = Kmer {
	    	segment_index : 0,
	    	position : 2,
	    	creation_time : 2
	    };
	    let kmer_4 : Kmer = Kmer {
	    	segment_index : 1,
	    	position : 0,
	    	creation_time : 3
	    };
	    let kmer_5 : Kmer = Kmer {
	    	segment_index : 1,
	    	position : 1,
	    	creation_time : 4
	    };
	    assert_eq!(None, kmer_hash_table.get_kmer("TAT"));
	    assert_eq!(Some((&vec![kmer_1], vec![0])), kmer_hash_table.get_kmer("ATG"));
	    assert_eq!(Some((&vec![kmer_2], vec![0])), kmer_hash_table.get_kmer("TGT"));
	    assert_eq!(Some((&vec![kmer_3, kmer_4], vec![0, 1])), kmer_hash_table.get_kmer("GTG"));
	    assert_eq!(Some((&vec![kmer_5], vec![0])), kmer_hash_table.get_kmer("TGC"));
    }

    #[test]
    fn test_get_max_j_1() {
    	assert_eq!(31, DNAHashTable::get_max_j(100));
    	assert_eq!(3, DNAHashTable::get_max_j(3));
    }

    #[test]
    fn test_get_most_likely_position_1() {
        let mut segments : Vec<String> = Vec::<String>::new();
	    let k : usize = 3;
	    segments.push("ATGTGACGCCGATG".to_string());
	    segments.push("GTGC".to_string());
	    let kmer_hash_table = DNAHashTable::new(&segments, k);
	    assert_eq!(Some((0, 0)), kmer_hash_table.get_most_likely_position(&segments, "ATGT"));
    	assert_eq!(Some((0, 5)), kmer_hash_table.get_most_likely_position(&segments, "ACGCCGA"));
    }

    #[test]
    fn test_get_most_likely_position_2() {
        let mut segments : Vec<String> = Vec::<String>::new();
	    let k : usize = 3;
	    segments.push("ATGTGATGCCGATG".to_string());
	    segments.push("GTGCGATGATAGAG".to_string());
	    let kmer_hash_table = DNAHashTable::new(&segments, k);
	    assert_eq!(Some((1, 4)), kmer_hash_table.get_most_likely_position(&segments, "CATGATA"));
    }
}