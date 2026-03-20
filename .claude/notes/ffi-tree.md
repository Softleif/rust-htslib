# hts-sys FFI Call Tree

This document maps every `hts-sys` (C htslib) function called from `rust-htslib`,
organized as a dependency tree. **Leaf functions** — those whose C implementations
do not call other hts-sys public API functions used by this crate — are the best
candidates for incremental replacement with pure Rust.

## How to read this document

- **→** means "calls (in its C implementation)"
- **[LEAF]** means the function's C implementation does not call any other function
  in this crate's used set — it can be replaced independently
- **Complexity** ratings reflect the C implementation size and algorithmic difficulty
- **Rust feasibility** considers whether standard Rust crates (flate2, std::collections, etc.)
  can replace the C dependencies

---

## 1. Utility / Bitset (kbitset.h — header-only)

All three are trivial inline functions. Best starting point for pure Rust replacement.

```
kbs_init          [LEAF]  ~2 lines   | Complexity: trivial  | Allocate bitset
kbs_insert        [LEAF]  ~2 lines   | Complexity: trivial  | Set bit (bitwise OR)
kbs_destroy       [LEAF]  ~2 lines   | Complexity: trivial  | Free bitset
```

**Rust replacement**: Use `bit_set` crate or a simple `Vec<u64>` wrapper.
These are used only by `bcf_trim_alleles` / `bcf_remove_allele_set` in `src/bcf/record.rs`.

---

## 2. BGZF Module (bgzf.c)

BGZF (Blocked Gzip Format) is a foundational I/O layer. Most higher-level functions
depend on it, so replacing BGZF unlocks many upstream replacements.

```
bgzf_is_bgzf     [LEAF]  ~11 lines  | Complexity: low    | Check BGZF magic header
bgzf_seek         [LEAF]  ~16 lines  | Complexity: low    | Seek in BGZF stream (delegates to internal bgzf_seek_common)
bgzf_open                 ~21 lines  | Complexity: low    | Open BGZF file handle
  → (internal: hopen, bgzf_read_init, bgzf_write_init — not in our set)
bgzf_read                 ~51 lines  | Complexity: medium | Read from BGZF stream
  → (internal: bgzf_read_block — not in our set)
bgzf_write                ~26 lines  | Complexity: low    | Write to BGZF stream
  → bgzf_flush
bgzf_flush                ~45 lines  | Complexity: medium | Flush BGZF buffers
  → (internal: deflate_block, hwrite — not in our set)
bgzf_close                ~45 lines  | Complexity: medium | Close BGZF handle
  → bgzf_flush
bgzf_thread_pool          ~41 lines  | Complexity: medium | Attach thread pool to BGZF
  → (internal: hts_tpool functions — but hts_tpool_init IS in our set)
  → hts_tpool_init (indirect, via pool creation)
```

**Dependency chain**: `bgzf_close → bgzf_flush`, `bgzf_write → bgzf_flush`

**Rust replacement strategy**: Use `flate2` crate for deflate/inflate. The `noodles`
crate already has a pure-Rust BGZF implementation that could serve as reference.
Start with `bgzf_is_bgzf` and `bgzf_seek` (leaves), then tackle `bgzf_flush` to
unlock `bgzf_close` and `bgzf_write`.

---

## 3. Thread Pool (thread_pool.c)

```
hts_tpool_init            ~90 lines  | Complexity: high   | Create pthread-based thread pool
  → (internal: pthread_create, pthread_mutex_init, etc.)
hts_tpool_destroy         ~35 lines  | Complexity: medium | Destroy thread pool
  → (internal: pthread_join, pthread_mutex_destroy, etc.)
```

Neither is a leaf in the strict sense (they use pthreads, which is a system dependency),
but they don't call other functions in our set.

**Effective leaves**: Yes — no hts-sys cross-dependencies.
**Rust replacement**: Use `rayon` thread pool or `std::thread::scope`. Moderate effort
due to different threading model.

---

## 4. HTS Core (hts.c)

```
hts_get_bgzfp     [LEAF]  ~6 lines   | Complexity: trivial | Get BGZF pointer from htsFile
hts_get_format    [LEAF]  ~3 lines   | Complexity: trivial | Get format struct from htsFile
hts_itr_destroy   [LEAF]  ~13 lines  | Complexity: low     | Destroy iterator
hts_idx_destroy   [LEAF]  ~26 lines  | Complexity: low     | Destroy index
hts_set_fai_filename [LEAF] ~15 lines | Complexity: low    | Set reference FAI path
hts_getline       [LEAF]  ~30 lines  | Complexity: low     | Read a line (format-dispatched)
hts_itr_query     [LEAF]  ~150 lines | Complexity: high    | Query index for region (bin-based lookup)
cram_seek         [LEAF]  ~7 lines   | Complexity: low     | Seek in CRAM file

hts_open                  ~2 lines   | Complexity: low     | Open HTS file (wrapper)
  → (internal: hts_open_format → bgzf_open for BAM/BCF)
  → bgzf_open (indirect)
hts_close                 ~60 lines  | Complexity: medium  | Close HTS file
  → bgzf_close
hts_set_threads           ~10 lines  | Complexity: low     | Set thread count
  → hts_set_opt, bgzf_thread_pool (indirect)
hts_set_thread_pool       ~10 lines  | Complexity: low     | Attach thread pool
  → bgzf_thread_pool
hts_set_opt               ~150 lines | Complexity: high    | Set generic option (variadic)
  → hts_set_threads, hts_set_thread_pool
hts_itr_next              ~60 lines  | Complexity: medium  | Advance iterator
  → bgzf_seek (via seek), readrec callback
```

**Dependency chain**:
```
hts_open → bgzf_open
hts_close → bgzf_close → bgzf_flush
hts_set_opt → hts_set_threads → bgzf_thread_pool → hts_tpool_init
hts_itr_next → bgzf_seek
```

---

## 5. FAIDX Module (faidx.c)

```
faidx_nseq        [LEAF]  ~4 lines   | Complexity: trivial | Count of sequences
faidx_iseq        [LEAF]  ~4 lines   | Complexity: trivial | Sequence name by index
faidx_seq_len64   [LEAF]  ~6 lines   | Complexity: trivial | Sequence length by name

faidx_fetch_seq64         ~128 lines | Complexity: medium  | Fetch subsequence
  → bgzf_seek, bgzf_read (via internal fai_retrieve)
fai_load                  ~126 lines | Complexity: medium  | Load FASTA index
  → bgzf_open, bgzf_close (via internal fai_load3_core)
  → fai_build (conditionally, if FAI_CREATE flag)
fai_build                 ~96 lines  | Complexity: high    | Build FASTA index from scratch
  → bgzf_open, bgzf_close, fai_destroy
fai_destroy               ~10 lines  | Complexity: low     | Destroy FAIDX handle
  → bgzf_close
```

**Dependency chain**:
```
fai_build → bgzf_open, bgzf_close, fai_destroy
fai_load → bgzf_open, bgzf_close, fai_build (conditional)
fai_destroy → bgzf_close
faidx_fetch_seq64 → bgzf_seek, bgzf_read
```

**Rust replacement strategy**: The three leaf accessors (`faidx_nseq`, `faidx_iseq`,
`faidx_seq_len64`) are trivial to replace. The full FAIDX module is a good candidate
for a complete pure-Rust rewrite since the `.fai` format is simple (TSV) and sequence
retrieval is straightforward with a Rust BGZF reader.

---

## 6. Tabix Module (tbx.c)

```
tbx_name2id       [LEAF]  ~4 lines   | Complexity: low     | Sequence name → ID (hash lookup)
tbx_destroy       [LEAF]  ~13 lines  | Complexity: low     | Destroy tabix index
tbx_seqnames      [LEAF]  ~28 lines  | Complexity: low     | Get all sequence names
tbx_readrec       [LEAF]  ~21 lines  | Complexity: medium  | Parse one tabix record
  → (internal: bgzf_getline — not in our used set)

tbx_index_load            ~50 lines  | Complexity: low     | Load tabix index from disk
  → hts_idx_load3 (internal, but uses hts index infrastructure)
tbx_index_build3          ~15 lines  | Complexity: medium  | Build tabix index
  → bgzf_open, bgzf_close, hts_idx_save_as (internal), tbx_destroy
```

**Dependency chain**:
```
tbx_index_build3 → bgzf_open, bgzf_close, tbx_destroy
tbx_index_load → (internal hts_idx functions)
```

---

## 7. SAM/BAM Module (sam.c, header.c)

### 7a. Header functions

```
sam_hdr_tid2name  [LEAF]  ~14 lines  | Complexity: trivial | TID → name (array access)
sam_hdr_name2tid  [LEAF]  ~18 lines  | Complexity: low     | Name → TID (hash lookup)
sam_hdr_str       [LEAF]  ~6 lines   | Complexity: trivial | Get header text
sam_hdr_line_name [LEAF]  ~40 lines  | Complexity: low     | Get line name from header
sam_hdr_destroy   [LEAF]  ~24 lines  | Complexity: low     | Destroy header (refcounted)
sam_hdr_parse     [LEAF]  ~12 lines  | Complexity: low     | Parse header text
  → (internal: sam_hdr_init, sam_hdr_add_lines — not in our used set)

sam_hdr_read              ~45 lines  | Complexity: medium  | Read header from file
  → sam_hdr_dup, sam_hdr_parse (internal path varies by format)
sam_hdr_write             ~75 lines  | Complexity: medium  | Write header to file
  → (format-dispatched: bam_hdr_write, cram_write_SAM_hdr — internal)
sam_hdr_dup               ~65 lines  | Complexity: medium  | Duplicate header
  → (internal: sam_hdr_init, sam_hrecs_rebuild_text)
```

### 7b. Record functions

```
bam_copy1         [LEAF]  ~9 lines   | Complexity: trivial | Copy BAM record
bam_endpos        [LEAF]  ~6 lines   | Complexity: low     | Calculate alignment end position
bam_aux_get       [LEAF]  ~22 lines  | Complexity: low     | Get auxiliary tag value
bam_aux_append    [LEAF]  ~30 lines  | Complexity: low     | Append auxiliary tag
bam_aux_del       [LEAF]  ~5 lines   | Complexity: trivial | Delete auxiliary tag

bam_aux_update_int        ~65 lines  | Complexity: medium  | Update integer aux tag
  → bam_aux_get
bam_aux_update_float      ~35 lines  | Complexity: low     | Update float aux tag
  → bam_aux_get
bam_aux_update_str        ~65 lines  | Complexity: medium  | Update string aux tag
  → bam_aux_get
bam_aux_update_array      ~55 lines  | Complexity: medium  | Update array aux tag
  → bam_aux_get
```

### 7c. I/O functions

```
sam_read1                 ~50 lines  | Complexity: high    | Read one record
  → (format-dispatched: bam_read1, sam_read1_sam, cram — all internal)
sam_write1                ~110 lines | Complexity: high    | Write one record
  → (format-dispatched: bam_write_idx1, sam_format1 — internal)
sam_parse1                ~180 lines | Complexity: high    | Parse SAM text → record
  → (internal: bam_name2id, bam_parse_cigar, etc.)
```

### 7d. Index functions

```
sam_index_load    [LEAF]  ~4 lines   | Complexity: trivial | Load BAM index (wrapper)
sam_index_load2   [LEAF]  ~3 lines   | Complexity: trivial | Load BAM index with custom path
sam_itr_querys            ~7 lines   | Complexity: low     | Query by region string
  → hts_itr_query, sam_hdr_name2tid (indirect via hts_itr_querys)
sam_itr_queryi            ~10 lines  | Complexity: low     | Query by numeric TID/coords
  → hts_itr_query
sam_index_build3          ~40 lines  | Complexity: medium  | Build BAM index
  → hts_open, hts_set_threads, hts_close, hts_idx_destroy
```

### 7e. Pileup functions

```
bam_plp_set_maxcnt [LEAF] ~3 lines   | Complexity: trivial | Set max pileup depth
bam_plp_init      [LEAF]  ~14 lines  | Complexity: low     | Init pileup iterator
bam_plp_destroy   [LEAF]  ~15 lines  | Complexity: low     | Destroy pileup iterator
bam_plp_reset     [LEAF]  ~11 lines  | Complexity: low     | Reset pileup state

bam_plp_auto              ~14 lines  | Complexity: low     | Auto-advance pileup
  → (internal: bam_plp64_auto which calls the read callback → sam_read1)
```

---

## 8. BCF/VCF Module (vcf.c, synced_bcf_reader.c, vcfutils.c)

### 8a. Record lifecycle

```
bcf_init          [LEAF]  ~6 lines   | Complexity: trivial | Allocate BCF record
bcf_destroy       [LEAF]  ~5 lines   | Complexity: trivial | Free BCF record
bcf_clear         [LEAF]  ~30 lines  | Complexity: low     | Reset BCF record fields

bcf_dup                   ~4 lines   | Complexity: trivial | Duplicate BCF record
  → bcf_init, bcf_copy
bcf_copy                  ~30 lines  | Complexity: low     | Copy BCF record data
  → bcf_clear
```

### 8b. Header functions

```
bcf_hdr_id2int    [LEAF]  ~7 lines   | Complexity: trivial | Tag name → ID (hash lookup)
bcf_hdr_get_version [LEAF] ~9 lines  | Complexity: trivial | Get VCF version string
bcf_hdr_add_sample [LEAF] ~8 lines   | Complexity: low     | Add sample to header
bcf_hdr_append    [LEAF]  ~8 lines   | Complexity: low     | Append header line
bcf_hdr_sync      [LEAF]  ~33 lines  | Complexity: low     | Rebuild header ID index
bcf_hdr_destroy   [LEAF]  ~35 lines  | Complexity: low     | Destroy header
bcf_hdr_remove            ~63 lines  | Complexity: medium  | Remove header record
  → (internal: bcf_hdr_get_hrec, bcf_hdr_unregister_hrec)

bcf_hdr_init              ~40 lines  | Complexity: low     | Init header
  → bcf_hdr_append (adds fileformat line)
bcf_hdr_dup               ~18 lines  | Complexity: low     | Duplicate header
  → bcf_hdr_init, bcf_hdr_parse (internal)
bcf_hdr_subset            ~80 lines  | Complexity: high    | Subset header by samples
  → bcf_hdr_init, bcf_hdr_id2int, bcf_hdr_destroy
bcf_hdr_read              ~60 lines  | Complexity: medium  | Read header from file
  → bcf_hdr_init, bcf_hdr_destroy, bgzf_read (via internal)
bcf_hdr_write             ~40 lines  | Complexity: medium  | Write header to file
  → bcf_hdr_sync, bgzf_write, bgzf_flush (via internal)
```

### 8c. Record modification

```
bcf_unpack        [LEAF]  ~68 lines  | Complexity: high    | Unpack binary BCF record
  → (internal: bcf_dec_size, bcf_dec_int1, etc.)

bcf_update_id             ~12 lines  | Complexity: low     | Set record ID
  → bcf_unpack
bcf_add_id                ~30 lines  | Complexity: low     | Append to record ID
  → bcf_unpack
bcf_update_filter         ~11 lines  | Complexity: low     | Set FILTER field
  → bcf_unpack
bcf_add_filter            ~17 lines  | Complexity: low     | Add FILTER value
  → bcf_unpack
bcf_remove_filter         ~12 lines  | Complexity: low     | Remove FILTER value
  → bcf_unpack
bcf_update_alleles        ~80 lines  | Complexity: medium  | Set alleles
  → bcf_unpack
bcf_update_info           ~120 lines | Complexity: high    | Update INFO field
  → bcf_hdr_id2int, bcf_unpack
bcf_update_format         ~120 lines | Complexity: high    | Update FORMAT field
  → bcf_hdr_id2int, bcf_unpack
bcf_update_format_string  ~25 lines  | Complexity: low     | Update string FORMAT
  → bcf_update_format
bcf_translate             ~60 lines  | Complexity: medium  | Translate record between headers
  → bcf_unpack, bcf_hdr_id2int
bcf_subset                ~25 lines  | Complexity: medium  | Subset samples in record
  → (internal: bcf_unpack_fmt_core1)
```

### 8d. Data access

```
bcf_get_fmt               ~7 lines   | Complexity: trivial | Get FORMAT field descriptor
  → bcf_hdr_id2int
bcf_get_info_values       ~7 lines   | Complexity: trivial | Get INFO field data
  → bcf_hdr_id2int
bcf_get_format_values     ~7 lines   | Complexity: trivial | Get FORMAT field data
  → bcf_hdr_id2int
```

### 8e. Allele manipulation (vcfutils.c)

```
bcf_trim_alleles          ~54 lines  | Complexity: high    | Remove unused alleles
  → bcf_get_fmt, kbs_init, kbs_insert, kbs_destroy, bcf_remove_allele_set
bcf_remove_allele_set     ~220 lines | Complexity: high    | Remove alleles by bitset
  → bcf_unpack, bcf_update_info, bcf_update_format, bcf_hdr_id2int, (many)
```

### 8f. VCF formatting

```
vcf_format                ~316 lines | Complexity: high    | BCF record → VCF text
  → bcf_unpack
```

### 8g. I/O

```
bcf_read                  ~10 lines  | Complexity: medium  | Read BCF/VCF record
  → (format-dispatched internally)
bcf_write                 ~40 lines  | Complexity: medium  | Write BCF/VCF record
  → bcf_hdr_sync, bgzf_write (via internal)
bcf_index_build3          ~42 lines  | Complexity: medium  | Build BCF index
  → hts_open, hts_set_threads, hts_close, hts_idx_destroy
```

### 8h. Synced reader

```
bcf_sr_set_opt    [LEAF]  ~35 lines  | Complexity: low     | Set reader option
bcf_sr_remove_reader [LEAF] ~14 lines | Complexity: low    | Remove reader from set

bcf_sr_init               ~9 lines   | Complexity: low     | Init synced reader
  → bcf_sr_set_opt
bcf_sr_destroy            ~20 lines  | Complexity: low     | Destroy synced reader
  → bcf_hdr_destroy
bcf_sr_add_reader         ~180 lines | Complexity: high    | Add file to synced reader
  → hts_open, bcf_hdr_read, hts_close, bgzf_thread_pool
bcf_sr_next_line          ~26 lines  | Complexity: high    | Advance to next line
  → (internal: complex state machine with region matching)
bcf_sr_seek               ~28 lines  | Complexity: high    | Seek in synced reader
  → (internal: reader seek + region overlap)
bcf_sr_set_threads        ~14 lines  | Complexity: low     | Set reader thread count
  → hts_tpool_init
```

---

## Summary: All Leaf Functions

Sorted by recommended replacement priority (considering complexity, impact, and
dependency unlocking potential):

### Priority 1: Trivial leaves (< 10 lines, immediate wins)

| Function | Module | Lines | Rust difficulty | Notes |
|----------|--------|-------|-----------------|-------|
| `kbs_init` | kbitset | 2 | trivial | Replace with `BitSet` |
| `kbs_insert` | kbitset | 2 | trivial | Replace with `BitSet` |
| `kbs_destroy` | kbitset | 2 | trivial | Replace with `BitSet` |
| `faidx_nseq` | faidx | 4 | trivial | Struct field access |
| `faidx_iseq` | faidx | 4 | trivial | Array index access |
| `faidx_seq_len64` | faidx | 6 | trivial | HashMap lookup |
| `bcf_init` | bcf | 6 | trivial | Struct allocation |
| `bcf_destroy` | bcf | 5 | trivial | Struct deallocation |
| `bcf_clear` | bcf | 30 | easy | Field reset |
| `bcf_hdr_id2int` | bcf | 7 | trivial | HashMap lookup |
| `bcf_hdr_get_version` | bcf | 9 | trivial | String accessor |
| `hts_get_bgzfp` | hts | 6 | trivial | Struct field access |
| `hts_get_format` | hts | 3 | trivial | Struct field access |
| `bam_endpos` | bam | 6 | easy | CIGAR arithmetic |
| `bam_copy1` | bam | 9 | easy | Memcpy wrapper |
| `bam_aux_del` | bam | 5 | trivial | Delegation |
| `bam_plp_set_maxcnt` | bam | 3 | trivial | Simple setter |
| `sam_hdr_tid2name` | sam | 14 | trivial | Array access |
| `sam_hdr_str` | sam | 6 | trivial | String accessor |
| `sam_index_load` | sam | 4 | trivial | Wrapper |
| `sam_index_load2` | sam | 3 | trivial | Wrapper |

### Priority 2: Simple leaves (10-40 lines, straightforward logic)

| Function | Module | Lines | Rust difficulty | Notes |
|----------|--------|-------|-----------------|-------|
| `bgzf_is_bgzf` | bgzf | 11 | easy | Read 16 bytes, check magic |
| `bgzf_seek` | bgzf | 16 | easy | File position management |
| `bam_aux_get` | bam | 22 | easy | Linear tag search |
| `bam_aux_append` | bam | 30 | easy | Buffer append |
| `bam_plp_init` | bam | 14 | easy | Struct initialization |
| `bam_plp_destroy` | bam | 15 | easy | Cleanup |
| `bam_plp_reset` | bam | 11 | easy | State reset |
| `sam_hdr_destroy` | sam | 24 | easy | Refcounted cleanup |
| `sam_hdr_parse` | sam | 12 | easy | Text parsing wrapper |
| `sam_hdr_name2tid` | sam | 18 | easy | HashMap lookup |
| `sam_hdr_line_name` | sam | 40 | easy | Switch-based lookup |
| `hts_itr_destroy` | hts | 13 | easy | Cleanup |
| `hts_idx_destroy` | hts | 26 | easy | Cleanup |
| `hts_set_fai_filename` | hts | 15 | easy | String storage |
| `hts_getline` | hts | 30 | easy | Line reader dispatch |
| `bcf_hdr_add_sample` | bcf | 8 | easy | Append to array |
| `bcf_hdr_append` | bcf | 8 | easy | Parse + add record |
| `bcf_hdr_sync` | bcf | 33 | easy | Rebuild hash index |
| `bcf_hdr_destroy` | bcf | 35 | easy | Hash cleanup |
| `bcf_sr_set_opt` | bcf | 35 | easy | Option setter |
| `bcf_sr_remove_reader` | bcf | 14 | easy | Array removal |
| `tbx_name2id` | tbx | 4 | easy | Hash lookup |
| `tbx_destroy` | tbx | 13 | easy | Cleanup |
| `tbx_seqnames` | tbx | 28 | easy | Extract from hash |
| `cram_seek` | hts | 7 | easy | Queue drain + seek |

### Priority 3: Medium leaves (40-150 lines, some algorithmic complexity)

| Function | Module | Lines | Rust difficulty | Notes |
|----------|--------|-------|-----------------|-------|
| `bcf_unpack` | bcf | 68 | medium | Binary BCF unpacking with variable-length encoding |
| `tbx_readrec` | tbx | 21 | medium | Tab-separated record parsing |
| `hts_itr_query` | hts | 150 | hard | Bin-based index query algorithm |
| `bcf_hdr_remove` | bcf | 63 | medium | Search and remove from header |

---

## Dependency Graph (what unlocks what)

Replacing leaf functions is valuable on its own, but some replacements **unlock**
the ability to replace upstream functions. Here are the key chains:

```
                    ┌─────────────┐
                    │  kbs_*      │ ← trivial, unlocks bcf_trim_alleles
                    └─────────────┘

         ┌──────────────────────────────────┐
         │  bgzf_is_bgzf, bgzf_seek [LEAF] │
         └──────────────────────────────────┘
                         │
                    ┌────▼────────┐
                    │ bgzf_flush  │ ← unlocks bgzf_close, bgzf_write
                    └─────────────┘
                         │
              ┌──────────▼──────────┐
              │ bgzf_close          │ ← unlocks hts_close, fai_destroy
              │ bgzf_write          │ ← unlocks bcf_hdr_write, bcf_write
              │ bgzf_open           │ ← unlocks hts_open, fai_load, fai_build
              └─────────────────────┘
                         │
              ┌──────────▼──────────┐
              │ hts_open, hts_close │ ← unlocks sam/bcf readers/writers
              └─────────────────────┘

         ┌──────────────────────────────────┐
         │  bcf_unpack [LEAF]               │ ← unlocks nearly all bcf_update_*,
         │                                  │   vcf_format, bcf_translate, etc.
         └──────────────────────────────────┘

         ┌──────────────────────────────────┐
         │  bam_aux_get [LEAF]              │ ← unlocks all bam_aux_update_* fns
         └──────────────────────────────────┘

         ┌──────────────────────────────────┐
         │  bcf_hdr_id2int [LEAF]           │ ← unlocks bcf_get_fmt,
         │                                  │   bcf_get_info_values,
         │                                  │   bcf_update_info/format, etc.
         └──────────────────────────────────┘
```

### High-value leaf replacements (unlock the most upstream functions)

1. **`bcf_unpack`** — unlocks 10+ BCF record modification functions
2. **`bcf_hdr_id2int`** — unlocks all BCF data access and update functions
3. **`bam_aux_get`** — unlocks all 4 `bam_aux_update_*` functions
4. **`bgzf_is_bgzf` + `bgzf_seek`** — first steps toward full BGZF replacement
5. **`kbs_*`** — unlocks `bcf_trim_alleles` path

---

## Module-level Rust replacement feasibility

| Module | # Functions | # Leaves | Full Rust feasible? | Key blocker |
|--------|------------|----------|--------------------|-|
| kbitset | 3 | 3 | Yes | None |
| faidx | 7 | 3 | Yes | Needs BGZF reader |
| tbx | 6 | 4 | Yes | Needs BGZF + index |
| bgzf | 8 | 2 | Yes | Needs flate2, threading |
| hts core | 14 | 8 | Partially | hts_open/close deeply embedded |
| sam/bam | 29 | 17 | Partially | sam_read1/write1 format complexity |
| bcf/vcf | 44 | 11 | Partially | Binary encoding/decoding |
| tpool | 2 | 2* | Yes | Threading model differs |

*tpool functions are effective leaves (no hts-sys cross-deps) despite using pthreads.

---

## Recommended replacement order

1. **Phase 1 — Utility leaves**: kbs_*, faidx accessors, trivial struct accessors
2. **Phase 2 — Strategic leaves**: bcf_unpack, bcf_hdr_id2int, bam_aux_get (high unlock value)
3. **Phase 3 — BGZF foundation**: bgzf_is_bgzf → bgzf_seek → bgzf_flush → bgzf_close/write/open
4. **Phase 4 — Module rewrites**: Complete faidx, tbx modules (once BGZF is in Rust)
5. **Phase 5 — Record I/O**: sam_read1/write1, bcf_read/write (requires hts_open/close)
