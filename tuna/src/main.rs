extern crate primes;

mod read_inputs;

fn main() {
    let input = read_inputs::read_fa_file("example_file2.txt");

    for x in input {
        println!("{} {}", x.seg_id, x.seg_string);
    }
}
