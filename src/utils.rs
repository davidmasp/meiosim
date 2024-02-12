use std::fs;
use std::io;
use log::{info, warn};
use noodles::vcf::record::chromosome;
use regex::Regex;

pub fn list_files_in_directory(input_folder: &str, extension: &str) -> io::Result<Vec<String>> {
    let mut file_names = Vec::new();
    for entry in fs::read_dir(input_folder)? {
        let entry = entry?;
        // so this is only the filename, could we improve it until
        // it is the full path.
        let path = entry.path();
        if path.is_file() && path.extension() == Some(std::ffi::OsStr::new(extension)) {
            if let Ok(file_name) = entry.file_name().into_string() {
                let complete_file_name = format!("{}/{}", input_folder, file_name);
                file_names.push(complete_file_name);
            }
        }
    }
    Ok(file_names)
}

pub fn ensure_directory_exists(path: &str, verbose: bool) -> io::Result<()> {
    if !std::path::Path::new(path).exists() {
        fs::create_dir_all(path)?;
        if verbose {
            info!("Directory {} created", path);
        }
    } else {
        if verbose {
            warn!("Directory {} already exists", path);
        }
    }
    Ok(())
}

pub fn capture_chromosome_from_file_name(text: &str, verbose: bool) -> String {
    let re = Regex::new(r"chr\d+|chrX|chrY").unwrap();
    let caps = re.captures(text).unwrap();
    info!("Captures: {:?}", caps);
    if caps.len() == 0 {
        panic!("No chromosome found in file name");
    } else if caps.len() > 1 {
        panic!("More than one chromosome found in file name");
    } else {
        let chromosome = caps[0].to_string();
        if verbose {
            info!("Chromosome found: {}", chromosome);
        }
        chromosome
    }

    
}

