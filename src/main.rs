
mod variants;
mod recombination;
mod utils;
mod io;
mod workflows;
use workflows::wrk_generate_offspring;
use workflows::wrk_format_vcf;
use io::FamilyOut;
use recombination::RecombinationMapGenome;

use clap::{Parser, Subcommand, Args};
use log::info;
// use log::warn;
use simplelog;

use rand::SeedableRng;
use rand::rngs::StdRng;
use rand::seq::SliceRandom;

#[derive(Parser)]
#[command(author, version, about, long_about = None)]
#[command(propagate_version = true)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    /// Generates files that simulate offspring from two parents
    Main(Main),
    /// Support command to convert VCF to DWGSIM format
    Vcf2dwgsim(Vcf2dwgsim)
}

#[derive(Args)]
struct Main  {
    #[arg(long, help = "Sets the level of verbosity")]
    verbose: bool,
    #[arg(long, help = "Do recombination map use header?")]
    recomheader: bool,
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
    seed: u64,
    #[arg(short = 'f', long, value_name = "SIZE", help = "Sets the family size of the generated family tree")]
    familysize: u8,
    #[arg(short = 'g', long, value_name = "GENOME", help = "Sets the genome file")]
    genome: String,
}

#[derive(Args)]
struct Vcf2dwgsim  {
    #[arg(long, help = "Sets the level of verbosity")]
    verbose: bool,
    #[arg(short, long, value_name = "FILE", help = "VCF to convert to DWGSIM format")]
    vcf: String,
    #[arg(long, value_name = "PREFIX", help = "Sets the prefix string for the output")]
    prefix: String,
}

fn main() {
    let _ = simplelog::SimpleLogger::init(simplelog::LevelFilter::Info, simplelog::Config::default());
    let cli = Cli::parse();

    match &cli.command {
        Commands::Main(opts) => {
            let verbose = opts.verbose;
            let recomb_maps = &opts.recombination; // debug/recombmaps
            let pop_variants = &opts.population; // debug/vcfcollections
            let denovo_variants = &opts.denovo; // debug/dnmcollections
            let sample1 = &opts.parent1;
            let sample2 = &opts.parent2;
            let prefix = &opts.prefix;
            let recom_header = opts.recomheader;
            let genome_file = &opts.genome;
            let seed_value: u64 = opts.seed;

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
            let genome_recomb_map = RecombinationMapGenome::from_path(
                    recomb_maps,
                    "map",
                    recom_header);
            let mut popvars = variants::VCFCollection::from_path(&pop_variants, "gz", verbose);
            let genome_hash = utils::read_genome_file(genome_file);
            
            let mut rng: StdRng = StdRng::seed_from_u64(seed_value);

            let dnm_files = utils::list_files_in_directory(
                            &denovo_variants,
                            "vcf").unwrap();

            // randomly select at least family.samples.len() files
            let selected_dnm_files: Vec<String> = dnm_files
                    .choose_multiple(&mut rng,
                            family.samples.len()).cloned().collect();

            for i in 0..family.samples.len() { 
                let sample = &family.samples[i];
                let dnm_file = &selected_dnm_files[i];
                wrk_generate_offspring(&sample,
                                    &genome_recomb_map,
                                    &mut popvars,
                                    &dnm_file,
                                    verbose,
                                    &genome_hash,
                                    &mut rng);
            }
        }
        Commands::Vcf2dwgsim(opts) => {
            let verbose = opts.verbose;
            let vcf_file = &opts.vcf;
            let prefix = &opts.prefix;

            if verbose {
                info!("Mode: VCF -> DWGSIM");
                info!("VCF file: {}", vcf_file);
                info!("Prefix: {}", prefix);
            }

            let outputfilename = format!("{}_meiosimvariants.txt", prefix);
            wrk_format_vcf(
                &outputfilename,
                vcf_file,
                verbose
            )
        }
    }
}
