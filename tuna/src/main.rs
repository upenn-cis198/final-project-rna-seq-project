extern crate primes;
#[macro_use]
extern crate structopt;
#[macro_use]
use std::collections::HashMap;

mod read_inputs;
mod args;

fn main() {
    let input = read_inputs::read_fa_file("example_file2.txt");

    let mut input2 : HashMap<String, i32> = HashMap::new();
    input2.insert("AAAA".to_string(), 33);
    input2.insert("BBBBBB".to_string(), 1233);
    input2.insert("!#@#@!RWWEWEWERSSSS".to_string(), 132123121);

    read_inputs::write_output("aa.txt", input2);
    // for x in input {
    //     println!("{} {}", x.seg_id, x.seg_string);
    // }
}
