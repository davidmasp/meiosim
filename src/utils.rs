use std::fs;
use std::io;
use log::{info, warn};
use regex::Regex;

use std::fs::File;
use std::io::BufRead;
use std::collections::HashMap;
use std::path::Path;

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
    if verbose {
        info!("Captures: {:?}", caps);
    }
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

pub fn from_vu8_to_string(x: Vec<&[u8]>)-> Vec<String> {
    x.iter().map(|sample| {
        String::from_utf8(sample.to_vec()).unwrap_or_else(|_| String::from("Invalid UTF-8"))
    }).collect()
}


pub fn read_genome_file(file_path: String) -> HashMap<String, u64> {
    let mut contig_size = HashMap::new();
    let input_file = File::open(&Path::new(&file_path)).unwrap();
    let reader = io::BufReader::new(input_file);
    let lines = reader.lines();
    for line in lines {
        let line = line.unwrap();
        let fields: Vec<&str> = line.split_whitespace().collect();
        let contig = fields[0].to_string();
        let size = fields[1].parse::<u64>().unwrap();
        contig_size.insert(contig, size);
    }
    contig_size
}

