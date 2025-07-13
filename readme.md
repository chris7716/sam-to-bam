# SAM to BAM Converter (Rust, Streaming Version)

This is a **minimal and memory-efficient** tool that converts a SAM file into a BAM file written in **pure Rust** (without using `htslib`). It streams input line-by-line and compresses output using GZIP.

---

## Features

- Written in **pure Rust**
- Converts SAM to BAM (including header and alignment records)
- Streams input → low memory usage
- Basic support for CIGAR, SEQ, QUAL, and `A`, `i`, `Z` optional tags

---

## Requirements

- Rust toolchain: [Install Rust](https://www.rust-lang.org/tools/install)
- A SAM file (e.g., `example.sam`)

---

## Usage

### 1. Clone or Download

```bash
git clone https://github.com/your-repo/sam-to-bam-rust.git
cd sam-to-bam-rust
```

### 2. Build

```bash
cargo build --release
```

### 3. Run

```bash
./target/release/sam_to_bam <input.sam> <output.bam>
```

Example:

```bash
./target/release/sam_to_bam test_data/example.sam test_data/output.bam
```

---

## Output

- The output BAM file is compressed using standard GZIP (not BGZF).
- Most BAM-compatible tools like `samtools` can still read it, but **random access (e.g. indexing)** is not supported unless re-blocked with `samtools reheader` or `bgzip`.

---

## Limitations

- Only basic SAM optional tags are supported: `A`, `i`, `Z`
- Output uses GZIP compression via `flate2`, not true BGZF (no indexing support)
- No support yet for tag types `B`, `f`, or BAM random access
- No multithreading
