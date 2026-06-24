use std::env;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::process;

// Compression imports
use flate2::Compression;
use flate2::read::MultiGzDecoder;
use flate2::write::GzEncoder;

#[inline(always)]
fn process_and_write_seq(raw_seq: &[u8], out_seq: &mut Vec<u8>, writer: &mut Box<dyn Write>, w: usize) {
    // 1. Ignorer les \n et \r à la fin (au lieu de vérifier chaque byte)
    let mut len = raw_seq.len();
    while len > 0 && (raw_seq[len - 1] == b'\n' || raw_seq[len - 1] == b'\r') {
        len -= 1;
    }
    let clean_seq = &raw_seq[..len];

    // 2. Prévenir les réallocations mémoire
    out_seq.clear();
    out_seq.reserve(len + 1);

    let mut prev_byte = b'\0';
    let mut count = 0;

    // 3. Boucle sans la condition (if b == \n)
    for &b in clean_seq {
        if b == prev_byte {
            count += 1;
        } else {
            prev_byte = b;
            count = 1;
        }
        if count <= w {
            out_seq.push(b);
        }
    }

    out_seq.push(b'\n');
    writer.write_all(out_seq).unwrap();
}

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

    // --- INPUT ---
    let buf_size = 2 * 1024 * 1024;

    let in_file = File::open(input_path).expect("Error opening input");
    let mut reader: Box<dyn BufRead> = if input_path.ends_with(".gz") {
        Box::new(BufReader::with_capacity(buf_size, MultiGzDecoder::new(in_file)))
    } else {
        Box::new(BufReader::with_capacity(buf_size, in_file))
    };

    // --- OUTPUT ---
    let out_file = File::create(output_path).expect("Error creating output");
    let mut writer: Box<dyn Write> = if output_path.ends_with(".gz") {
        let encoder = GzEncoder::new(out_file, Compression::fast());
        Box::new(BufWriter::with_capacity(buf_size, encoder))
    } else {
        Box::new(BufWriter::with_capacity(buf_size, out_file))
    };

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
        process_and_write_seq(&seq, &mut out_seq, &mut writer, w);

        //Check if it not a fasta file (i.e. plus start with '>')
        if !plus.is_empty() && plus[0] == b'>' {
            writer.write_all(&plus).unwrap();
            process_and_write_seq(&qual, &mut out_seq, &mut writer, w);
        }
    }

    // Explicitly flush to ensure the Gzip stream is finished correctly
    writer.flush().unwrap();
}
