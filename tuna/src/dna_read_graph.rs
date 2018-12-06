extern crate multiset;

use dna_hash_table::DNAHashTable;
use dna_hash_table::Kmer;
use self::multiset::HashMultiSet;
use std::collections::HashSet;
use std::collections::HashMap;
use log::*;

pub struct DNAReadGraph<'a> {
	nodes : Vec<ReferenceRead>,
	node_hash_table : DNAHashTable<'a>,
	kmer_hash_table : DNAHashTable<'a>
}

impl<'a> DNAReadGraph<'a> {
	pub fn new(segments : &'a Vec<String>, l : usize, k : usize, d : usize) -> DNAReadGraph<'a> {
		let kmer_hash_table : DNAHashTable = DNAHashTable::new(segments, k);
		let node_hash_table : DNAHashTable = DNAHashTable::new(segments, l);
		let default_reference_lmer : Kmer = Kmer {
			segment_index : 0,
			position : 0,
			creation_time : 0
		};
		let default_reference_read : ReferenceRead = ReferenceRead {
			lmer : default_reference_lmer,
			near_reads_a : Vec::<Vec<usize>>::new(),
			near_reads_c : Vec::<Vec<usize>>::new(),
			near_reads_g : Vec::<Vec<usize>>::new(),
			near_reads_t : Vec::<Vec<usize>>::new(),
			read_distances : HashMap::<usize, usize>::new(),
		};

		//Make this more efficient later - wastes space
		let mut nodes = vec![default_reference_read; DNAHashTable::get_table_size(segments, k)];

		//CHANGE TO ITERATOR
		let mut segment_iter_index : usize = 0;
		let mut kmer_iter_index : usize = 0;

		debug!("{}", node_hash_table.hash_table.len());
		while segment_iter_index < node_hash_table.hash_table.len() {
			while kmer_iter_index < node_hash_table.hash_table[segment_iter_index].len() {
				let lmer : &Kmer = &(node_hash_table.hash_table[segment_iter_index][kmer_iter_index]);
				//CONSIDER MAKING INTO OWN METHOD
				let mut traversed_matching_lmers : HashSet<(usize, usize)> = HashSet::<(usize, usize)>::new();

				let mut near_reads_a = vec![Vec::<usize>::new(); l];
				let mut near_reads_c = vec![Vec::<usize>::new(); l];
				let mut near_reads_g = vec![Vec::<usize>::new(); l];
				let mut near_reads_t = vec![Vec::<usize>::new(); l];
				let mut read_distances = HashMap::<usize, usize>::new();

				//Hash every k-mer of the l-mer
				for i in lmer.position..(lmer.position + l - k + 1) {
					let kmer_string : &str = &(segments[lmer.segment_index][i..(i + k)]);
					let (kmers, matching_indexes) = kmer_hash_table.get_kmer(kmer_string).unwrap();
					for matching_index in matching_indexes {
						let matching_kmer : &Kmer = &kmers[matching_index];
						let matching_lmer_position : usize = matching_kmer.position;
						if !(matching_lmer_position < i || (matching_lmer_position - i + l) > segments[matching_kmer.segment_index].len()) {
							let matching_lmer_position : usize = matching_lmer_position - i;
							if !(matching_lmer_position == lmer.position && matching_kmer.segment_index == lmer.segment_index) {
								
								let lmer_string : &str = &(segments[lmer.segment_index][lmer.position..(lmer.position + l)]);
								let matching_lmer_string : &str = &(segments[matching_kmer.segment_index][matching_lmer_position..(matching_lmer_position + l)]);
								let lmer_distance : usize = DNAReadGraph::lmers_distance(lmer_string, matching_lmer_string, d);
								if !traversed_matching_lmers.contains(&(matching_kmer.segment_index, matching_lmer_position)) 
								&& lmer_distance <= d {

									let (transition_positions, transition_letters) = DNAReadGraph::get_lmer_differences(lmer_string, matching_lmer_string, d);
									let matching_lmer_index : usize = DNAReadGraph::get_lmer_index(&node_hash_table, matching_kmer.segment_index, matching_lmer_position, matching_lmer_string);
									read_distances.insert(matching_lmer_index, lmer_distance);

									let mut transition_positions_iter = transition_positions.iter();
									let mut transition_letters_iter = transition_letters.iter();

									while let (Some(transition_position), Some(transition_letter)) = 
									(transition_positions_iter.next(), transition_letters_iter.next()) {
										match transition_letter {
											'A' => near_reads_a[*transition_position].push(matching_lmer_index),
											'C' => near_reads_c[*transition_position].push(matching_lmer_index),
											'G' => near_reads_g[*transition_position].push(matching_lmer_index),
											'T' => near_reads_t[*transition_position].push(matching_lmer_index),
											_ => {},
										};
									}

									traversed_matching_lmers.insert((matching_kmer.segment_index, matching_lmer_position));
								}
							}
						}
					}
				}

				nodes[lmer.creation_time] = ReferenceRead {
					lmer : lmer.clone(),
					near_reads_a : near_reads_a,
					near_reads_c : near_reads_c,
					near_reads_g : near_reads_g,
					near_reads_t : near_reads_t,
					read_distances : read_distances,
				};
				kmer_iter_index += 1;
			}
			kmer_iter_index = 0;
			segment_iter_index += 1;
			//debug!("{}", segment_iter_index);
		}


		DNAReadGraph {
			nodes : nodes,
			node_hash_table : node_hash_table,
			kmer_hash_table : kmer_hash_table,
		}
	}

	pub fn get_read_segment_indexes(&self, segments : &Vec<String>, reads : &Vec<String>) -> HashMap<i32, i32> {
		let mut segment_index_counts : HashMap<i32, i32> = HashMap::new();
		debug!("{}", segments[0]);
		for read in reads {
			debug!("{}", read);
			match self.get_read_graph_segment_index(segments, read, self.node_hash_table.k) {
		    	Some(segment_index) => {
	                let count = segment_index_counts.entry(segment_index as i32).or_insert(0);
	                *count += 1;
	                debug!("{}", *count);
	    		},
		    	None => {
		    		debug!("No match");
		    	},
		    }
	    };
		segment_index_counts
	}

	fn get_read_graph_segment_index(&self, segments: &Vec<String>, read : &str, l : usize) -> Option<usize>{
		match self.kmer_hash_table.get_most_likely_position(segments, read) {
	    	Some((segment_index, position)) => {
                let initial_node_index = DNAReadGraph::get_lmer_index(&self.node_hash_table, segment_index, position, &segments[segment_index][position..(position + l)]);
                let initial_node = &self.nodes[initial_node_index];
                let initial_lmer = &initial_node.lmer;
                let lmer_string : &str = &segments[initial_lmer.segment_index][initial_lmer.position..(initial_lmer.position + l)];
                let initial_distance = DNAReadGraph::lmers_distance(lmer_string, read, l);
                let (transition_positions, transition_letters) = DNAReadGraph::get_lmer_differences(lmer_string, read, l);
                let mut near_lmers = HashMultiSet::<usize>::new();
                for i in 0..transition_positions.len() {
                	let transition_position = transition_positions[i];
                	let transition_letter = transition_letters[i];
          			let empty_vector = Vec::<usize>::new();
                	let near_transitions = match transition_letter {
                		'A' => &initial_node.near_reads_a[transition_position],
                		'C' => &initial_node.near_reads_c[transition_position],
                		'G' => &initial_node.near_reads_g[transition_position],
                		'T' => &initial_node.near_reads_t[transition_position],
                		_ => &empty_vector,
                	};
                	for near_transition in near_transitions {
                		near_lmers.insert(*near_transition);
                	}
                }
                let mut min_read_distance = initial_distance;
                let mut min_distance_lmer = initial_node_index;
                for near_lmer in near_lmers.distinct_elements() {
                	let lmer_transition_count = near_lmers.count_of(*near_lmer);
                	let read_lmer_distance = initial_distance
                	+ initial_node.read_distances[near_lmer] - lmer_transition_count;
                	if read_lmer_distance < min_read_distance {
                		min_read_distance = read_lmer_distance;
                		min_distance_lmer = *near_lmer;
                	}
                }
                Some((&self.nodes[min_distance_lmer]).lmer.segment_index)
    		},
	    	None => None,
	    }
	}

	//Combine with get lmer differences
	fn lmers_distance(lmer1 : &str, lmer2 : &str, d : usize) -> usize {
		let mut differences : usize = 0;

		let mut lmer1_iterator = lmer1.chars();
		let mut lmer2_iterator = lmer2.chars();

		while let (Some(lmer1_char), Some(lmer2_char)) = (lmer1_iterator.next(), lmer2_iterator.next()) {
			if lmer1_char != lmer2_char {
				differences += 1;
			}
			if differences > d {
				break;
			}
		}
		differences
	}

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
	near_reads_a : Vec<Vec<usize>>,
	near_reads_c : Vec<Vec<usize>>,
	near_reads_g : Vec<Vec<usize>>,
	near_reads_t : Vec<Vec<usize>>,
	read_distances : HashMap<usize, usize>,
}

#[cfg(test)]
mod tests {
	use dna_hash_table::DNAHashTable;
	use dna_hash_table::Kmer;
	use super::*;

	#[test]
    fn test_lmers_within_distance_1() {
    	let lmer_1 : &str = "ATAGGATA";
    	let lmer_2 : &str = "ATAGGATA";
    	assert_eq!(0, DNAReadGraph::lmers_distance(lmer_1, lmer_2, 1));
    }

    #[test]
    fn test_lmers_within_distance_2() {
    	let lmer_1 : &str = "ATGGCATA";
    	let lmer_2 : &str = "ATAGGATA";
    	assert!(2, DNAReadGraph::lmers_distance(lmer_1, lmer_2, 1));
    }

    #[test]
    fn test_lmers_within_distance_3() {
    	let lmer_1 : &str = "ATGGCATA";
    	let lmer_2 : &str = "ATAGGAAA";
    	assert_eq!(3, DNAReadGraph::lmers_distance(lmer_1, lmer_2, 3));
    }

    #[test]
    fn test_get_lmer_differences_1() {
    	let lmer_1 : &str = "ATAGGATA";
    	let lmer_2 : &str = "ATAGGATA";

    	let transition_positions : Vec<usize> = Vec::<usize>::new();
    	let transition_letters : Vec<char> = Vec::<char>::new();
    	assert_eq!((transition_positions, transition_letters), 
    	DNAReadGraph::get_lmer_differences(lmer_1, lmer_2, 3));
    }

    #[test]
    fn test_get_lmer_differences_2() {
    	let lmer_1 : &str = "ATAGGATA";
    	let lmer_2 : &str = "ATAAGAAA";

    	let transition_positions : Vec<usize> = vec![3, 6];
    	let transition_letters : Vec<char> = vec!['A', 'A'];
    	assert_eq!((transition_positions, transition_letters), 
    	DNAReadGraph::get_lmer_differences(lmer_1, lmer_2, 3));
    }

    #[test]
    fn test_get_lmer_differences_3() {
    	let lmer_1 : &str = "ATAGGATA";
    	let lmer_2 : &str = "ATCATAAA";

    	let transition_positions : Vec<usize> = vec![2, 3, 4, 6];
    	let transition_letters : Vec<char> = vec!['C', 'A', 'T', 'A'];
    	assert_eq!((transition_positions, transition_letters), 
    	DNAReadGraph::get_lmer_differences(lmer_1, lmer_2, 10));
    }

    #[test]
    fn test_get_lmer_index_1() {
        let mut segments : Vec<String> = Vec::<String>::new();
	    let k : usize = 1;
	    segments.push("ATGTG".to_string());
	    let kmer_hash_table = DNAHashTable::new(&segments, k);
	    assert_eq!(0, DNAReadGraph::get_lmer_index(&kmer_hash_table, 0, 0, "A"));
	    assert_eq!(2, DNAReadGraph::get_lmer_index(&kmer_hash_table, 0, 2, "G"));
	    assert_eq!(4, DNAReadGraph::get_lmer_index(&kmer_hash_table, 0, 4, "G"));
    }

    #[test]
    fn test_get_lmer_index_2() {
        let mut segments : Vec<String> = Vec::<String>::new();
	    let k : usize = 3;
	    segments.push("ATGTG".to_string());
	    let kmer_hash_table = DNAHashTable::new(&segments, k);
	    assert_eq!(0, DNAReadGraph::get_lmer_index(&kmer_hash_table, 0, 0, "ATG"));
	    assert_eq!(1, DNAReadGraph::get_lmer_index(&kmer_hash_table, 0, 1, "TGT"));
	    assert_eq!(2, DNAReadGraph::get_lmer_index(&kmer_hash_table, 0, 2, "GTG"));
    }

    #[test]
    fn test_get_lmer_index_3() {
        let mut segments : Vec<String> = Vec::<String>::new();
	    let k : usize = 3;
	    segments.push("ATGTG".to_string());
	    segments.push("GTGC".to_string());
	    let kmer_hash_table = DNAHashTable::new(&segments, k);
	    assert_eq!(0, DNAReadGraph::get_lmer_index(&kmer_hash_table, 0, 0, "ATG"));
	    assert_eq!(1, DNAReadGraph::get_lmer_index(&kmer_hash_table, 0, 1, "TGT"));
	    assert_eq!(2, DNAReadGraph::get_lmer_index(&kmer_hash_table, 0, 2, "GTG"));
	    assert_eq!(3, DNAReadGraph::get_lmer_index(&kmer_hash_table, 1, 0, "GTG"));
	    assert_eq!(4, DNAReadGraph::get_lmer_index(&kmer_hash_table, 1, 1, "TGC"));
    }

    #[test]
    fn test_create_dna_read_graph_1() {
        let mut segments : Vec<String> = Vec::<String>::new();
	    let l : usize = 3;
	    let k : usize = 1;
	    let d : usize = 8;
	    segments.push("ATGTG".to_string());
	    let dna_read_graph = DNAReadGraph::new(&segments, l, k, d);
	    
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
	    assert_eq!(kmer_1, dna_read_graph.nodes[0].lmer);
	    assert_eq!(vec![Vec::<usize>::new(), Vec::<usize>::new(), Vec::<usize>::new()], dna_read_graph.nodes[0].near_reads_A);
	    assert_eq!(vec![vec![2], Vec::<usize>::new(), Vec::<usize>::new()], dna_read_graph.nodes[0].near_reads_G);
	    assert_eq!(kmer_2, dna_read_graph.nodes[1].lmer);
	    assert_eq!(kmer_3, dna_read_graph.nodes[2].lmer);
    }

    #[test]
    fn test_create_dna_read_graph_2() {
        let mut segments : Vec<String> = Vec::<String>::new();
	    let l : usize = 3;
	    let k : usize = 1;
	    let d : usize = 1;
	    segments.push("ATGAG".to_string());
	    let dna_read_graph = DNAReadGraph::new(&segments, l, k, d);
	    
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
	    assert_eq!(kmer_1, dna_read_graph.nodes[0].lmer);
	    assert_eq!(vec![Vec::<usize>::new(), Vec::<usize>::new(), Vec::<usize>::new()], dna_read_graph.nodes[0].near_reads_A);
	    assert_eq!(vec![Vec::<usize>::new(), Vec::<usize>::new(), Vec::<usize>::new()], dna_read_graph.nodes[0].near_reads_G);
	    
	    assert_eq!(kmer_2, dna_read_graph.nodes[1].lmer);
	    assert_eq!(kmer_3, dna_read_graph.nodes[2].lmer);
    }

    #[test]
    fn test_create_dna_read_graph_3() {
        let mut segments : Vec<String> = Vec::<String>::new();
	    let l : usize = 3;
	    let k : usize = 1;
	    let d : usize = 1;
	    segments.push("ATTTT".to_string());
	    let dna_read_graph = DNAReadGraph::new(&segments, l, k, d);
	    
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
	    assert_eq!(kmer_1, dna_read_graph.nodes[0].lmer);
	    assert_eq!(vec![Vec::<usize>::new(), Vec::<usize>::new(), Vec::<usize>::new()], dna_read_graph.nodes[0].near_reads_A);
	    assert_eq!(vec![Vec::<usize>::new(), Vec::<usize>::new(), Vec::<usize>::new()], dna_read_graph.nodes[0].near_reads_G);
	    assert_eq!(vec![vec![1, 2], Vec::<usize>::new(), Vec::<usize>::new()], dna_read_graph.nodes[0].near_reads_T);
	    assert_eq!(kmer_2, dna_read_graph.nodes[1].lmer);
	    assert_eq!(kmer_3, dna_read_graph.nodes[2].lmer);
    }
}