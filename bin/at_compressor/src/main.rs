use std::env;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::process;

fn main() {
    // 1. Parse command line arguments
    let args: Vec<String> = env::args().collect();
    if args.len() != 4 {
        eprintln!("Usage: {} <input.fastq> <output.fasta> <w>", args[0]);
        process::exit(1);
    }

    let input_path = &args[1];
    let output_path = &args[2];
    let w: usize = match args[3].parse() {
        Ok(val) => val,
        Err(_) => {
            eprintln!("Error: 'w' must be a valid positive integer.");
            process::exit(1);
        }
    };

    // 2. Setup fast buffered I/O
    let in_file = File::open(input_path).unwrap_or_else(|e| {
        eprintln!("Error opening input file: {}", e);
        process::exit(1);
    });
    let mut reader = BufReader::with_capacity(128 * 1024, in_file);

    let out_file = File::create(output_path).unwrap_or_else(|e| {
        eprintln!("Error creating output file: {}", e);
        process::exit(1);
    });
    let mut writer = BufWriter::with_capacity(128 * 1024, out_file);

    // 3. Pre-allocate vectors so we don't allocate memory during the loop
    let mut header = Vec::new();
    let mut seq = Vec::new();
    let mut plus = Vec::new();
    let mut qual = Vec::new();
    let mut out_seq = Vec::with_capacity(2048);

    // 4. Process the file record by record
    loop {
        header.clear();
        seq.clear();
        plus.clear();
        qual.clear();
        out_seq.clear();

        // Read 4 lines (1 FASTQ record)
        let bytes_read = reader.read_until(b'\n', &mut header).unwrap();
        if bytes_read == 0 {
            break; // End of file
        }
        reader.read_until(b'\n', &mut seq).unwrap();
        reader.read_until(b'\n', &mut plus).unwrap();
        reader.read_until(b'\n', &mut qual).unwrap();

        // Convert FASTQ header (@) to FASTA header (>)
        if !header.is_empty() && header[0] == b'@' {
            header[0] = b'>';
        }
        writer.write_all(&header).unwrap();

        // Compress the Sequence
        let mut prev_byte = b'\0';
        let mut count = 0;

        for &b in &seq {
            // Skip carriage returns and newlines from the raw byte array
            if b == b'\n' || b == b'\r' {
                continue;
            }

            // Track consecutive identical bytes
            if b == prev_byte {
                count += 1;
            } else {
                prev_byte = b;
                count = 1;
            }

            // Check if the current nucleotide is A or T (case-insensitive just in case)
            let is_at = true; //b == b'A' || b == b'T' || b == b'a' || b == b't';
            
            // Keep the base if it is NOT an A/T, OR if the streak count is within our window `w`
            if !is_at || count <= w {
                out_seq.push(b);
            }
        }
        
        // Append a newline to the compressed sequence and write to file
        out_seq.push(b'\n');
        writer.write_all(&out_seq).unwrap();
    }
}
