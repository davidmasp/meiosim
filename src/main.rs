

mod variants;
mod recombination;
mod utils;
mod io;
mod workflows;
use workflows::wrk_generate_offspring;
use io::FamilyOut;
use noodles::vcf::header::record::value::map::info;
use recombination::RecombinationMapGenome;

use clap::Parser;
use log::{info, warn};
use simplelog;

#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
struct Args {
    #[arg(long, help = "Sets the level of verbosity")]
    verbose: bool,
    #[arg(short, long, value_name = "FOLDER", help = "Sets the folder path to recombination maps")]
    recombination: String,
    #[arg(short = 'v', long, value_name = "FOLDER", help = "Sets the folder path to VCF collection of population variants")]
    population: String,
    #[arg(short, long, value_name = "FOLDER", help = "Sets the folder path to VCF collection of denovo variants")]
    denovo: String,
    #[arg(short = 'p', long, value_name = "SAMPLE", help = "Sets the sample1 string option")]
    parent1: String,
    #[arg(short = 'P', long, value_name = "SAMPLE", help = "Sets the sample2 string option")]
    parent2: String,
    #[arg(long, value_name = "PREFIX", help = "Sets the prefix string for the output")]
    prefix: String,
    #[arg(long, value_name = "SEED", help = "Sets the seed")]
    seed: u32,
    #[arg(short = 'f', long, value_name = "SIZE", help = "Sets the family size of the generated family tree")]
    familysize: u8,
}

fn main() {
    let opts: Args = Args::parse();

    // You can then use these options in your program
    let verbose = opts.verbose;
    let recomb_maps = opts.recombination; // debug/recombmaps
    let pop_variants = opts.population; // debug/vcfcollections
    let denovo_variants = opts.denovo; // debug/dnmcollections
    let sample1 = opts.parent1;
    let sample2 = opts.parent2;
    let prefix = opts.prefix;

    let _ = simplelog::SimpleLogger::init(simplelog::LevelFilter::Info, simplelog::Config::default());
    if verbose{
        info!("Recombination maps folder: {}", recomb_maps);
        info!("Population variants folder: {}", pop_variants);
        info!("Denovo variants folder: {}", denovo_variants);
        info!("Parent1: {}", sample1);
        info!("Parent2: {}", sample2);
        info!("Prefix: {}", prefix);
    }

    // Generate the family out:
    let family = FamilyOut::new(&prefix, &sample1, &sample2, opts.familysize, verbose);
    if verbose {
        family.samples.iter().for_each(|s| {
            info!("{} -> {}", s.name, s.targetvcfout);
        });
    }

    // load recombination maps
    let genome_recomb_map = RecombinationMapGenome::from_path(&recomb_maps, "map");
    let mut popvars = variants::VCFCollection::from_path(&pop_variants, "gz", verbose);

    for i in 0..family.samples.len() {
        let sample = &family.samples[i];
        wrk_generate_offspring(&sample, &genome_recomb_map, &mut popvars, &denovo_variants, verbose);
    }
  // PARAMS
}
