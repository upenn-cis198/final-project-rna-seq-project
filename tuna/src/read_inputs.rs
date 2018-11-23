use std::env;
use std::fs::File;
use std::io::BufReader;
use std::io::prelude::*;


pub fn is_char_RNA(x : char) -> bool {
	return x == 'A' || x == 'U' || x == 'C' || x == 'G' || x == 'T';
}

pub fn read_fq_file(filename : &str) -> Vec<String> {
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

