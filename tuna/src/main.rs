extern crate primes;
#[macro_use]
extern crate structopt;
#[macro_use]
mod read_inputs;
mod args;

fn main() {
    let input = read_inputs::read_fa_file("example_file2.txt");

    for x in input {
        println!("{} {}", x.seg_id, x.seg_string);
    }
}
