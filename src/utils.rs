use std::fs;
use std::io;

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