use std::env;
use std::fs::File;
use std::io::prelude::*;

pub Vec<String> read_fq_file(filename : String) {
	let mut f = File::open(filename).expect("file not found");
	for line in BufReader::new(f).lines() {
		println("{}", line?);
	}

	unimplemented!();
}