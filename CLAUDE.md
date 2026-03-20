# rust-htslib Development Guidelines

## FFI-to-Rust Replacement: Red-Green Testing

When replacing C (hts-sys) FFI calls with pure Rust implementations, always follow
red-green testing:

1. **Write a failing test first** that exercises the FFI function being replaced,
   covering edge cases and the exact behavior of the C implementation.
2. **Replace the FFI call** with pure Rust code.
3. **Verify the test passes** — the Rust implementation must match C behavior exactly.
4. Run the full test suite (`cargo test`) to catch regressions.

This ensures behavioral equivalence between C and Rust implementations and prevents
silent correctness bugs during incremental replacement.

See `.claude/notes/ffi-tree.md` for the full FFI dependency tree and replacement plan.
