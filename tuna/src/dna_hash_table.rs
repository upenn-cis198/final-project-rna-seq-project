use primes::PrimeSet;
use std::cmp;

pub struct DNAHashTable<'a> {
	pub hash_table : Vec<Vec<Kmer>>,
	segments : &'a Vec<String>,
	pub size : usize,
	k : usize,
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
						println!("{:?}", kmer);
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
	fn get_table_size(segments : &Vec<String>, k : usize) -> usize {
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
}

impl<'a> Iterator for DNAHashTable<'a> {
	type Item = &'a Kmer;

	fn next(&mut self) -> Option<&'a Kmer> {
		if self.segment_iter_index >= self.size {
			None
		} else {
			let next_kmer : &'a Kmer = &self.hash_table[self.segment_iter_index][self.kmer_iter_index];
			self.segment_iter_index += 1;
			self.kmer_iter_index += 1;
			Some(next_kmer)
		}
	}
}

//CONSIDER USING A TRAIT HERE FOR DIFFERING KMERS ON CREATION TIME

//A locus is a location in the genome, which we represent by the segment that the k-mer mapped to and the location on the segment
#[derive(Clone, Debug)]
pub struct Kmer {
	pub segment_index : usize,
	pub position : usize,
	pub creation_time : usize,
}