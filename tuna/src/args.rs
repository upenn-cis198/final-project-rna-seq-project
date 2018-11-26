extern crate structopt;


#[derive(Debug, StructOpt)]
#[structopt(name = "tuna", about = "RNA sequence aligner")]
pub struct Opt {
    /// The k size for the k-mers
    #[structopt(short = "k", long = "kk")]
    pub k: i32,
   
    /// Executing program arguments
    #[structopt(parse(from_str))]
    pub exe_args: Vec<String>,
}

