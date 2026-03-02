# 🧬 AutoDock-VS: Automated Virtual Screening Pipeline

A fully automated, end-to-end virtual screening pipeline for structure-based drug discovery. Given a receptor PDB and a library of SMILES, it handles protein preparation, ligand 3D conformer generation, active-site detection, parallel Gnina docking, and ranked results export — all in a single command.

---

## Features

- **Automated protein preparation** — strips waters, crystallographic artefacts, and common buffer molecules; adds polar hydrogens at physiological pH via OpenBabel (with graceful fallback)
- **3-D ligand conformer generation** — reads SMILES from Excel, embeds with ETKDGv3, and minimises with MMFF94 force field
- **Three-tier active-site detection** — cascades through co-crystal ligand centroid → fpocket → built-in grid-based cavity detector (pure Python/NumPy, no external tools required)
- **Parallel Gnina docking** — configurable worker count and CPU-per-worker for efficient multi-core throughput
- **Ranked results export** — full CSV summary + formatted Excel workbook of the top-100 compounds by CNN affinity, preserving your original compound IDs
- **Recovery mode** — re-scores from already-docked SDF files without re-running Gnina (`--recover` flag)

---

## Pipeline Overview

```
Step 0 — Protein Preparation     (remove solvent/HETATM, add H via OpenBabel)
Step 1 — Ligand Preparation      (SMILES → 3D SDF via RDKit ETKDGv3 + MMFF94)
Step 2 — Active-Site Detection   (Tier 1: co-crystal → Tier 2: fpocket → Tier 3: grid)
Step 3 — Gnina Docking           (parallel, configurable exhaustiveness)
Step 4 — Results Extraction      (CSV + top-100 Excel with compound IDs)
```

---

## Requirements

### Python dependencies

```bash
pip install rdkit pandas tqdm openpyxl scipy
```

### External tools

| Tool | Required | Notes |
|------|----------|-------|
| [Gnina](https://github.com/gnina/gnina) | ✅ Yes | Must be on `PATH` |
| [OpenBabel](https://openbabel.org) | ⚠️ Recommended | Improves hydrogen placement; falls back gracefully if absent |
| [fpocket](https://github.com/Discngine/fpocket) | ⚠️ Optional | Used as Tier 2 pocket detector; skipped if absent |

Install OpenBabel via conda:
```bash
conda install -c conda-forge openbabel
```

---

## Input Files

| File | Description |
|------|-------------|
| `blind_set.xlsx` | Compound library — must contain a `SMILES` column; an `ID` column is strongly recommended |
| `fixed_version.pdb` | Raw receptor PDB downloaded from the PDB or elsewhere |

Column names and file paths are configured at the top of the script via global constants.

---

## Usage

### Run the full pipeline

```bash
python virtual_screening_fast.py
```

### Re-score from existing docked SDF files (no re-docking)

Useful if docking completed successfully but result extraction failed:

```bash
python virtual_screening_fast.py --recover
```

---

## Configuration

Edit the constants at the top of `virtual_screening_fast.py`:

```python
INPUT_EXCEL            = "blind_set.xlsx"       # compound library
PROTEIN_PDB            = "fixed_version.pdb"    # raw receptor
BOX_SIZE               = 25                     # docking box (Å)
GNINA_EXHAUSTIVENESS   = 8
GNINA_NUM_WORKERS      = 4    # parallel Gnina processes
GNINA_CPU_PER_WORKER   = 2    # --cpu passed to each Gnina process
```

> Total CPU threads ≈ `GNINA_NUM_WORKERS × GNINA_CPU_PER_WORKER`. Keep this ≤ your logical core count.

To retain a metal cofactor during protein preparation (e.g. zinc):
```python
prepare_protein(..., keep_hetatm_residues={"ZN"})
```

---

## Outputs

| File | Description |
|------|-------------|
| `fixed_version_prepared.pdb` | Cleaned, protonated receptor |
| `prepared_ligands.sdf` | 3D-embedded compound library |
| `docking_outputs/` | Per-compound docked SDF files from Gnina |
| `results_summary.csv` | Full results for all compounds |
| `top100_results.xlsx` | Top-100 compounds ranked by CNN affinity, with original IDs |

---

## Active-Site Detection Logic

The pipeline uses a three-tier cascade — stopping as soon as a tier succeeds:

1. **Tier 1 — Co-crystal ligand** (fastest): Searches the *original* (raw) PDB for HETATM records, ignores solvent/ions, and uses the centroid of the largest ligand. Correctly uses the raw PDB rather than the prepared protein (which has HETATM stripped).

2. **Tier 2 — fpocket**: Runs `fpocket` on the prepared receptor and extracts the top-ranked pocket centroid.

3. **Tier 3 — Grid-based cavity detector** (pure Python): Builds a 3D occupancy grid, flood-fills the exterior, and clusters buried void voxels. Scored by `vol^1.5 × burial^0.5 × sphericity` to identify drug-like pockets including partially-open kinase grooves. Includes a manual 3D dilation fallback when SciPy is unavailable.

---

## License

MIT License — see `LICENSE` for details.

---

## Citation

If you use this pipeline in published work, please cite [Gnina](https://github.com/gnina/gnina) and RDKit as appropriate.
