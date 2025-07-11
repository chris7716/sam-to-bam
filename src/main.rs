use std::env;
use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use flate2::Compression;
use flate2::write::GzEncoder;

const BAM_MAGIC: &[u8] = b"BAM\x01";
const BAM_EOF: [u8; 28] = [
    31, 139, 8, 4, 0, 0, 0, 0,
    0, 255, 6, 0, 66, 67, 2, 0,
    27, 0, 3, 0, 0, 0, 0, 0,
    0, 0, 0, 0,
];

fn main() -> std::io::Result<()> {
    let args: Vec<String> = env::args().collect();
    if args.len() != 3 {
        eprintln!("Usage: {} <input.sam> <output.bam>", args[0]);
        std::process::exit(1);
    }

    let sam_file = File::open(&args[1])?;
    let reader = BufReader::new(sam_file);

    let bam_file = File::create(&args[2])?;
    let mut bgzf_writer = GzEncoder::new(bam_file, Compression::default());

    // Parse header
    let mut header_text = Vec::new();
    let mut ref_names = Vec::new();
    let mut ref_lengths = Vec::new();
    let mut records = Vec::new();

    for line in reader.lines() {
        let line = line?;
        if line.starts_with('@') {
            header_text.extend_from_slice(line.as_bytes());
            header_text.push(b'\n');

            if line.starts_with("@SQ") {
                let mut name = None;
                let mut len = None;
                for field in line.split('\t') {
                    if field.starts_with("SN:") {
                        name = Some(field[3..].to_string());
                    } else if field.starts_with("LN:") {
                        len = field[3..].parse::<u32>().ok();
                    }
                }
                if let (Some(n), Some(l)) = (name, len) {
                    ref_names.push(n);
                    ref_lengths.push(l);
                }
            }
        } else {
            records.push(line);
        }
    }

    // Write BAM magic
    bgzf_writer.write_all(BAM_MAGIC)?;

    // Write header text
    bgzf_writer.write_all(&(header_text.len() as u32).to_le_bytes())?;
    bgzf_writer.write_all(&header_text)?;

    // Write reference sequences
    bgzf_writer.write_all(&(ref_names.len() as u32).to_le_bytes())?;
    for (name, len) in ref_names.iter().zip(ref_lengths.iter()) {
        let name_cstr = format!("{}\0", name);
        bgzf_writer.write_all(&(name_cstr.len() as u32).to_le_bytes())?;
        bgzf_writer.write_all(name_cstr.as_bytes())?;
        bgzf_writer.write_all(&len.to_le_bytes())?;
    }

    // Write dummy alignment records (for demonstration)
    for record in records {
        let fields: Vec<&str> = record.split('\t').collect();
        if fields.len() < 11 {
            continue;
        }

        let qname = fields[0];
        let flag: u16 = fields[1].parse().unwrap_or(0);
        let rname = fields[2];
        let pos: i32 = fields[3].parse().unwrap_or(1) - 1;
        let mapq: u8 = fields[4].parse().unwrap_or(255);
        let cigar = fields[5];
        let rnext = fields[6];
        let pnext: i32 = fields[7].parse().unwrap_or(1) - 1;
        let tlen: i32 = fields[8].parse().unwrap_or(0);
        let seq = fields[9];
        let qual = fields[10];

        let tid = ref_names.iter().position(|r| r == rname).unwrap_or(0) as i32;
        let next_tid = if rnext == "*" {
            -1
        } else if rnext == "=" {
            tid
        } else {
            ref_names.iter().position(|r| r == rnext).unwrap_or(0) as i32
        };

        let l_read_name = qname.len() + 1;
        let n_cigar_op = cigar.matches(|c: char| c.is_ascii_alphabetic()).count() as u16;
        let l_seq = seq.len();
        let bin = 0u16;

        let mut record = Vec::new();
        record.extend_from_slice(&[0u8; 4]); // placeholder for block_size

        record.extend_from_slice(&(tid as i32).to_le_bytes());
        record.extend_from_slice(&pos.to_le_bytes());
        record.push(l_read_name as u8);
        record.push(mapq);
        record.extend_from_slice(&bin.to_le_bytes());
        record.extend_from_slice(&n_cigar_op.to_le_bytes());
        record.extend_from_slice(&flag.to_le_bytes());
        record.extend_from_slice(&(l_seq as u32).to_le_bytes());
        record.extend_from_slice(&(next_tid as i32).to_le_bytes());
        record.extend_from_slice(&pnext.to_le_bytes());
        record.extend_from_slice(&tlen.to_le_bytes());

        record.extend_from_slice(qname.as_bytes());
        record.push(0); // null terminator for read name

        let cigar_bytes = encode_cigar(cigar);
        record.extend_from_slice(&cigar_bytes);

        let seq_bytes = encode_seq(seq);
        record.extend_from_slice(&seq_bytes);

        let qual_bytes = encode_qual(qual);
        record.extend_from_slice(&qual_bytes);

        for tag_field in fields.iter().skip(11) {
            let parts: Vec<&str> = tag_field.splitn(3, ':').collect();
            if parts.len() != 3 {
                eprintln!("⚠️ Skipping malformed tag '{}'", tag_field);
                continue;
            }
        
            let tag = parts[0];
            let type_char = parts[1];
            let value = parts[2];
        
            match type_char {
                "A" => {
                    record.extend_from_slice(tag.as_bytes());
                    record.extend_from_slice(b"A");
                    record.push(value.as_bytes()[0]);
                }
                "i" => {
                    if let Ok(val) = value.parse::<i32>() {
                        record.extend_from_slice(tag.as_bytes());
                        record.extend_from_slice(b"i");
                        record.extend_from_slice(&val.to_le_bytes());
                    } else {
                        eprintln!("⚠️ Could not parse integer tag '{}'", tag_field);
                    }
                }
                "Z" => {
                    record.extend_from_slice(tag.as_bytes());
                    record.extend_from_slice(b"Z");
                    record.extend_from_slice(value.as_bytes());
                    record.push(0);
                }
                "f" | "B" => {
                    eprintln!("⚠️ Skipping unsupported tag type '{}' on record {}: {}", type_char, qname, tag_field);
                }
                _ => {
                    eprintln!("⚠️ Unknown tag type '{}' in record {}: {}", type_char, qname, tag_field);
                }
            }
        }        

        let block_size = (record.len() - 4) as u32;
        record[0..4].copy_from_slice(&block_size.to_le_bytes());
        bgzf_writer.write_all(&record)?;
    }

    // Finish BGZF
    let mut writer = bgzf_writer.finish()?;
    writer.write_all(&BAM_EOF)?;
    writer.flush()?;

    println!("✅ BAM written to {}", args[2]);
    Ok(())
}

fn encode_cigar(cigar: &str) -> Vec<u8> {
    let mut result = Vec::new();
    let mut num = String::new();
    for c in cigar.chars() {
        if c.is_ascii_digit() {
            num.push(c);
        } else {
            let len: u32 = num.parse().unwrap_or(0);
            let op_code = match c {
                'M' => 0, 'I' => 1, 'D' => 2, 'N' => 3, 'S' => 4,
                'H' => 5, 'P' => 6, '=' => 7, 'X' => 8, _ => 15,
            };
            let encoded = (len << 4) | (op_code as u32);
            result.extend_from_slice(&encoded.to_le_bytes());
            num.clear();
        }
    }
    result
}

fn encode_seq(seq: &str) -> Vec<u8> {
    let mut result = Vec::new();
    let mut byte = 0u8;
    for (i, base) in seq.chars().enumerate() {
        let code = match base.to_ascii_uppercase() {
            '=' => 0, 'A' => 1, 'C' => 2, 'M' => 3,
            'G' => 4, 'R' => 5, 'S' => 6, 'V' => 7,
            'T' => 8, 'W' => 9, 'Y' => 10, 'H' => 11,
            'K' => 12, 'D' => 13, 'B' => 14, 'N' => 15,
            _ => 15,
        };
        if i % 2 == 0 {
            byte = code << 4;
        } else {
            byte |= code;
            result.push(byte);
            byte = 0;
        }
    }
    if seq.len() % 2 != 0 {
        result.push(byte);
    }
    result
}

fn encode_qual(qual: &str) -> Vec<u8> {
    if qual == "*" {
        vec![0xFF]
    } else {
        qual.bytes().map(|b| b.saturating_sub(33)).collect()
    }
}
