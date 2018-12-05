use dna_hash_table::DNAHashTable;
use dna_hash_table::Kmer;

pub struct DNAReadGraph<'a> {
	nodes : Vec<ReferenceRead>,
	node_hash_table : DNAHashTable<'a>,
}

impl<'a> DNAReadGraph<'a> {
	pub fn new(segments : &'a Vec<String>, kmer_hash_table : &DNAHashTable, l : usize, k : usize, d : usize) -> DNAReadGraph<'a> {
		let node_hash_table : DNAHashTable = DNAHashTable::new(segments, l);
		let default_reference_lmer : Kmer = Kmer {
			segment_index : 0,
			position : 0,
			creation_time : 0
		};
		let default_reference_read : ReferenceRead = ReferenceRead {
			lmer : default_reference_lmer,
			near_reads : Vec::<NearRead>::new(),
		};
		let mut nodes = vec![default_reference_read; node_hash_table.size];

		//CHANGE TO ITERATOR
		let mut segment_iter_index : usize = 0;
		let mut kmer_iter_index : usize = 0;

		while segment_iter_index < node_hash_table.hash_table.len() {
			while kmer_iter_index < node_hash_table.hash_table[segment_iter_index].len() {
				let lmer : &Kmer = &(node_hash_table.hash_table[segment_iter_index][kmer_iter_index]);
				//CONSIDER MAKING INTO OWN METHOD

				let mut near_reads = Vec::<NearRead>::new();

				//Hash every k-mer of the l-mer
				for i in lmer.position..(lmer.position + l) {
					let kmer_string : &str = &(segments[lmer.segment_index][i..(i + k)]);
					let (kmers, matching_indexes) = kmer_hash_table.get_kmer(kmer_string).unwrap();

					for matching_index in matching_indexes {
						let matching_kmer : &Kmer = &kmers[matching_index];
						let matching_lmer_position : usize = matching_kmer.position - i;

						//WHY IS THIS COMPARISON USELESS?
						if !(matching_lmer_position < 0 || (matching_lmer_position + l) > segments[matching_kmer.segment_index].len()) {
							let lmer_string : &str = &(segments[lmer.segment_index][lmer.position..(lmer.position + l)]);
							let matching_lmer_string : &str = &(segments[matching_kmer.segment_index][matching_lmer_position..(matching_lmer_position + l)]);
							
							if DNAReadGraph::lmers_within_distance(lmer_string, matching_lmer_string, d) {
								let (transition_positions, transition_letters) = DNAReadGraph::get_lmer_differences(lmer_string, matching_lmer_string, d);
								let matching_lmer_index : usize = DNAReadGraph::get_lmer_index(&node_hash_table, matching_kmer.segment_index, matching_lmer_position, matching_lmer_string);

								near_reads.push(NearRead {
									read_index : matching_lmer_index,
									transition_positions : transition_positions,
									transition_letters : transition_letters,
								});
							}
						}
					}
				}

				nodes[lmer.creation_time] = ReferenceRead {
					lmer : lmer.clone(),
					near_reads : near_reads,
				};
				kmer_iter_index += 1;
			}
			kmer_iter_index = 0;
			segment_iter_index += 1;
		}


		DNAReadGraph {
			nodes : nodes,
			node_hash_table : node_hash_table,
		}
	}



	fn lmers_within_distance(lmer1 : &str, lmer2 : &str, d : usize) -> bool {
		let mut differences : usize = 0;

		let mut lmer1_iterator = lmer1.chars();
		let mut lmer2_iterator = lmer2.chars();

		while let (Some(lmer1_char), Some(lmer2_char)) = (lmer1_iterator.next(), lmer2_iterator.next()) {
			if lmer1_char != lmer2_char {
				differences += 1;
			}
		}
		differences <= d
	}

	//Get it working with str
	fn get_lmer_differences(lmer1 : &str, lmer2 : &str, d : usize) -> (Vec<usize>, Vec<char>) {
		let mut transition_positions : Vec<usize> = Vec::<usize>::new();
		let mut transition_letters : Vec<char> = Vec::<char>::new();
		transition_positions.reserve(d);
		transition_letters.reserve(d);

		let mut lmer1_iterator = lmer1.chars();
		let mut lmer2_iterator = lmer2.chars();

		let mut i : usize = 0;

		while let (Some(lmer1_char), Some(lmer2_char)) = (lmer1_iterator.next(), lmer2_iterator.next()) {
			if lmer1_char != lmer2_char {
				transition_positions.push(i);
				transition_letters.push(lmer2_char);
			}
			i += 1;
		}

		(transition_positions, transition_letters)
	}

	fn get_lmer_index(node_hash_table : &DNAHashTable, lmer_segment_index : usize, lmer_position : usize, lmer_string : &str) -> usize {
		let (lmer_entries, lmer_entries_indexes) = node_hash_table.get_kmer(lmer_string).unwrap();
		let mut lmer_index : usize = 0;
		for lmer_entry_index in lmer_entries_indexes {
			let possible_lmer : &Kmer = &lmer_entries[lmer_entry_index];

			//If we found the hash table entry that corresponds to it
			if possible_lmer.segment_index == lmer_segment_index && possible_lmer.position == lmer_position {
				lmer_index = possible_lmer.creation_time;
			}
		}
		lmer_index
	}
}

//TURN INTO REFERECNE TO kmer
#[derive(Clone)]
pub struct ReferenceRead {
	pub lmer : Kmer,
	near_reads : Vec<NearRead>
}

#[derive(Clone)]
struct NearRead {
	read_index : usize,
	transition_positions : Vec<usize>,
	transition_letters : Vec<char>
}