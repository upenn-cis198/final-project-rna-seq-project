extern crate structopt;
//Command: target/debug/tuna -v -k 10 example_reads.txt example_seqs.txt out.txt

#[derive(Debug, StructOpt)]
#[structopt(name = "tuna", about = "RNA sequence aligner")]
pub struct Opt {

	/// Whether or not we want to debug
	#[structopt(short = "v", long = "verbose")]
    pub verbose: bool,

    /// The k size for the k-mers
    #[structopt(short = "k", long = "kk")]
    pub k: usize,

    /// Number of partitions, sometimes useful in determing number of threads
    #[structopt(short = "p", long = "n_partition", default_value = "4")]
    pub n_partition: usize,

    /// The name of the reads inputs
    #[structopt(parse(from_str))]
    pub read_input_filename : String,

    /// The name of the sequence inputs
    #[structopt(parse(from_str))]
    pub seq_input_filename : String,

    /// The name of the output file containing the sequence IDs mapped to their counts
    #[structopt(parse(from_str))]
    pub seqcount_output_filename : String,
}