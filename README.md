# stats_function_group

Precise, rule-based functional-group recognition from PDB connectivity.

`stats_function_group` reads atoms and `CONECT` records from PDB files, builds a molecular graph, and reports functional-group counts for individual structures or entire directories. It optionally uses OpenBabel to improve bond-order and aromaticity perception while retaining a geometry-based fallback.

中文简介：本工具用于基于PDB连接关系进行精确、可解释的官能团识别，支持单个结构和批量统计，并输出稳定列顺序的CSV结果。

## Features

- Parses `ATOM`, `HETATM`, element, residue, coordinate, and `CONECT` records.
- Detects 5–8-membered rings directly from the molecular graph.
- Uses mutually exclusive ring classes to avoid duplicate ring counting.
- Recognizes oxygen-, nitrogen-, sulfur-, phosphorus-, halogen-, and carbonyl-centered groups.
- Supports per-residue annotation, per-file batch tables, and aggregate summaries.
- Uses only the Python standard library; OpenBabel is optional but recommended.
- Includes 28 focused regression structures and expected CSV outputs.

## Recognized functional groups

| Category | Labels |
|---|---|
| Rings | `heteroring5`, `aromatic_ring5`, `aliphatic_ring5`, `heteroring6`, `aromatic_ring6`, `aliphatic_ring6`, `ring7`, `ring8` |
| Oxygen | `phenol`, `alcohol`, `ether`, `epoxide` |
| Nitrogen | `nitrile`, `nitro`, `amide`, `amine` |
| Sulfur | `thioether`, `sulfone_or_sulfonyl`, `sulfoxide_or_sulfenyl_oxide`, `thioester` |
| Phosphorus | `phosphate`, `phosphorus_containing` |
| Halogens | `halide` |
| Carbonyl-centered | `carboxylic_acid`, `ester`, `acyl_chloride`, `ketone`, `aldehyde` |

Unclassified element centers may be emitted as `element_<symbol>`.

## Requirements

- Python 3.9 or newer
- Optional: OpenBabel Python bindings for improved bond-order and aromaticity perception

One convenient OpenBabel installation route is:

```bash
conda install -c conda-forge openbabel
```

The program still runs without OpenBabel by using PDB connectivity and bond-length heuristics.

## Usage

Single PDB:

```bash
python fg_from_pdb_connect.py \
  --pdb molecule.pdb \
  --out-annotated annotated.csv \
  --out-summary summary.csv
```

Batch processing:

```bash
python fg_from_pdb_connect.py \
  --pdb-dir structures \
  --recursive \
  --out-table fg_table.csv \
  --out-total fg_total_summary.csv
```

Disable OpenBabel explicitly:

```bash
python fg_from_pdb_connect.py --pdb molecule.pdb --no-use-openbabel
```

Run the included regression set:

```bash
python fg_from_pdb_connect.py \
  --pdb-dir test_stability/test_pdb \
  --out-table test_stability/fg_table.generated.csv \
  --out-total test_stability/fg_total_summary.generated.csv
```

The generated results can be compared with:

- `test_stability/fg_table.csv`
- `test_stability/fg_total_summary.csv`

## Output files

Single-file mode produces:

- `annotated.csv`: one row per recognized residue, with functional-group labels and counts.
- `summary.csv`: total counts for the analyzed PDB.

Batch mode produces:

- `fg_table.csv`: one row per PDB and one column per functional group.
- `fg_total_summary.csv`: aggregate counts across all processed PDB files.

## Accuracy notes

- Reliable PDB element columns and complete `CONECT` records give the best results.
- OpenBabel is recommended when aromaticity and bond orders are not explicit.
- The method is deterministic and interpretable, but it is rule-based; unusual valence states, incomplete structures, organometallic bonding, resonance, protonation, and tautomerism may require manual review.
- Ring detection is graph-based. Five- and six-membered rings are assigned to one mutually exclusive class; seven- and eight-membered rings are reported without finer aromatic/aliphatic subdivision.

## Repository contents

- `fg_from_pdb_connect.py`: functional-group recognition and CSV export.
- `test_stability/test_pdb/`: focused regression PDB files.
- `test_stability/*.csv`: expected batch results and test index.

The large local `all_run_pdb` reference collection is intentionally excluded from the public repository.

## License

Released under the [MIT License](LICENSE).
