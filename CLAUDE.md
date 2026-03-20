# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

rust-htslib provides safe Rust bindings to HTSlib, the standard C library for reading/writing HTS bioinformatics file formats (SAM/BAM/CRAM, VCF/BCF). The FFI layer comes from the `hts-sys` crate (a separate repo); this crate wraps it in idiomatic Rust APIs.

## Build & Test Commands

```bash
cargo build                          # build with default features (bzip2, lzma, curl)
cargo test -- --test-threads 1       # run all tests (sequential required)
cargo test <test_name> -- --test-threads 1  # run a single test
cargo fmt -- --check                 # check formatting
cargo clippy --all-features --all-targets -- -Dclippy::all -Dunused_imports  # lint (CI uses RUSTFLAGS="-Dwarnings")
cargo test --no-default-features     # test without optional features
cargo test --all-features            # test with all features enabled
```

Tests **must** run with `--test-threads 1` due to shared HTSlib resources.

### System Dependencies (Linux)

```bash
sudo apt-get install zlib1g-dev libbz2-dev musl musl-dev musl-tools clang libc6-dev
```

On macOS, htslib dependencies are typically available via Homebrew.

## Architecture

Each module wraps a distinct HTSlib file format with a consistent Reader/Writer pattern:

- **`bam/`** — SAM/BAM/CRAM reading and writing (~72% of codebase)
  - `mod.rs`: `Reader`, `IndexedReader`, `Writer`; the `Read` trait defines the core interface
  - `record.rs`: `Record` type with CIGAR, aux tags, sequence, quality access
  - `ext.rs`: Extension traits with additional record utilities
  - `pileup.rs`: Pileup iteration over genomic positions
  - `header.rs`: `Header` (mutable) and `HeaderView` (read-only, shared via `Arc`)
  - `buffer.rs`: `RecordBuffer` for bulk processing
  - `record_serde.rs`: Optional serde support (behind `serde_feature` flag)
- **`bcf/`** — VCF/BCF reading and writing
  - Same Reader/Writer/Record/Header pattern as BAM
  - INFO/FORMAT/genotype field access on records
- **`bgzf/`** — BGZIP compression Reader/Writer (implements `std::io::Read`/`Write`)
- **`faidx/`** — Indexed FASTA access
- **`tbx/`** — Tabix-indexed text file access (BED, GFF, VCF)
- **`tpool/`** — HTSlib thread pool wrapper (`Rc<RefCell<>>`)
- **`errors.rs`** — Error types via `thiserror`; crate-wide `Result<T>` alias
- **`htslib.rs`** — Re-exports `hts_sys::*` FFI bindings

### Key Patterns

- Unsafe FFI calls are confined to module internals; public APIs are safe Rust
- Both BAM and BCF define their own `Read` trait (same name, different modules)
- Headers use `Arc<HeaderView>` for shared, thread-safe access
- Test data lives in `test/` directory; tests are inline (`#[cfg(test)]` modules)

## Cargo Features

Default: `bzip2`, `lzma`, `curl`. Optional: `s3`, `gcs`, `libdeflate`, `bindgen`, `static`, `serde_feature`.

## FFI-to-Rust Replacement: Differential Proptest Strategy

When replacing C (hts-sys) FFI calls with pure Rust implementations, use
**differential property-based testing** to prove equivalence:

1. **Write the Rust function as a stub** that initially delegates to the C FFI call
   (with a `// TODO: replace with pure Rust` comment). This establishes the function
   signature and call sites.
2. **Write proptests that compare both implementations**: generate random inputs, call
   both the C function (via FFI) and the Rust function, and assert identical results.
   The C implementation is the oracle — the Rust function must match it exactly on all
   generated inputs.
3. **Replace the stub body** with the real pure Rust implementation.
4. **Run the proptests** — they now verify the Rust code against C on hundreds of
   random inputs, catching edge cases that hand-written tests would miss.
5. Run the full test suite (`cargo test -- --test-threads 1`) to catch regressions.

### Why this approach

- The C bindings remain available throughout the transition, so we can always call
  both implementations side-by-side.
- Proptests explore the input space far more thoroughly than hand-written examples.
- The C implementation serves as a living oracle — no need to manually reverse-engineer
  edge-case behavior.
- When the proptest passes, we have high confidence the replacement is correct.

### Example structure

```rust
/// Pure Rust replacement for htslib::some_function.
/// Initially delegates to C; replace body with Rust implementation.
fn some_function_rs(args) -> Result {
    // TODO: replace with pure Rust
    unsafe { htslib::some_function(args) }
}

#[cfg(test)]
mod proptest_some_function {
    proptest! {
        #[test]
        fn matches_c_implementation(input in arbitrary_input_strategy()) {
            let c_result = unsafe { htslib::some_function(input) };
            let rs_result = some_function_rs(input);
            prop_assert_eq!(c_result, rs_result);
        }
    }
}
```

See `.claude/notes/ffi-tree.md` for the full FFI dependency tree and replacement plan.
