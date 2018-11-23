extern crate primes;

mod read_inputs;

fn main() {
    let input = read_inputs::read_fq_file("example_file.txt");

    for x in input {
        println!("{}", x);
    }
}
