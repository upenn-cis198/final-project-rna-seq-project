use std::env;
use std::fs::File;
use std::io::BufReader;
use std::io::prelude::*;
use std::collections::HashMap;

pub struct FaEntry {
    pub seg_id : String,
    pub seg_string : String
}

impl FaEntry {
	pub fn new() -> Self {
		// env_logger::init();
		let fa_entry = FaEntry {
			seg_id : String::new(),
			seg_string : String::new()
		};
		fa_entry
	}

	pub fn from_read(seg_id : String, seg_string: String) -> Self {
		let fa_entry = FaEntry {
			seg_id : seg_id,
			seg_string : seg_string
		};
		fa_entry
	}
}


pub fn is_char_RNA(x : char) -> bool {
	return x == 'A' || x == 'U' || x == 'C' || x == 'G' || x == 'T';
}

/**
	read_fq_fasta_file:
	Given the name of a fq or fasta file, read it into a vector of strings consisting of
	the sequences
*/
pub fn read_fq_fasta_file(filename : &str) -> Vec<String> {
	let mut f = File::open(filename).expect("file not found");
	let mut fq_strings : Vec<String> = Vec::new();
	let mut current : String = String::new();
	for line in BufReader::new(f).lines() {
		let line_val = line.expect("line not available");
		let first_char = line_val.chars().next().unwrap();
		// println!("{}", line_val);
		// println!("{}", first_char == 'a');

		if first_char == '!' || first_char == '@' || first_char == '>' {
			// process new string
			if current != "" {
				fq_strings.push(current);
				current = String::new();
			}
			continue; // skip
		} else {
			let mapped : String = line_val.chars().filter(|x| is_char_RNA(*x)).collect();
			current.push_str(mapped.as_str());
		}
	}

	if current != "" {
		fq_strings.push(current);
		current = String::new();
	}
	return fq_strings
}

/**
	Read fa_file:
	Given the name of a fa file, read it into a vector of FaEntry consisting of a 
	sequence ID paired with the sequence. 
*/
pub fn read_fa_file(filename : &str) -> Vec<FaEntry> {
	let mut f = File::open(filename).expect("file not found");
	let mut fq_strings : Vec<FaEntry> = Vec::new();
	let mut current : String = String::new();
	let mut current_seg_id : String = String::new();
	for line in BufReader::new(f).lines() {
		let line_val = line.expect("line not available");
		let first_char = line_val.chars().next().unwrap();
		// println!("{}", line_val);
		// println!("{}", first_char == 'a');

		if first_char == '!' || first_char == '@' || first_char == '>' {
			// process new string
			if current != "" {
				let fa_entry = FaEntry::from_read(current_seg_id, current);
				fq_strings.push(fa_entry);
				current = String::new();
			}
			// split the header by whitespace, find the first word
			// assign the first word to current_seg_id
			
			let mut first_word : String = line_val.split_whitespace().next().expect("no word in header").to_string();
			first_word = first_word.chars().filter(|x| (*x).is_numeric()).collect();
			current_seg_id = first_word;
			continue; // skip
		} else {
			let mapped : String = line_val.chars().filter(|x| is_char_RNA(*x)).collect();
			current.push_str(mapped.as_str());
		}
	}

	if current != "" {
		let fa_entry = FaEntry::from_read(current_seg_id, current);
		fq_strings.push(fa_entry);
		current = String::new();
	}
	return fq_strings
}

pub fn write_output(filename : &str, counts : HashMap<String, i32>) -> bool {
	unimplemented!();
}