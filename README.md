# EDTA: Effect-Dependent Transcriptome Analysis

## Overview

This repository organizes the analysis code used in the EDTA (Effect-Dependent Transcriptome Analysis) study to support computational reproducibility. The currently documented module identifies ATAC-seq regions overlapping selected zebrafish genes and performs transcription factor binding motif enrichment using HOMER.

**Status:** Only Module 04 (CRE association and TFBS motif enrichment) is documented here. Other modules will be added after verification.

## Repository structure

The following tree reflects **files actually packaged** in the updated ZIP (plus the new scripts). Some original input data and HOMER motif matrices are still missing, as noted below.

```text
Phenotype-Guided-Transcriptome-Analysis-main/
├── data/
│   ├── BET_posi.txt
│   └── ATAC_GRCz11_peaks_434221.bed
├── scripts/
│   └── 04_TFBS/
│       ├── 01_gene_info_from_gtf_by_gene_name.py
│       ├── 02_gene2CREs.py
│       ├── 03_pos4homer.py
│       ├── 04_run_homer.slurm
│       ├── 05_annotate_peaks.slurm              # added
│       └── 06_filter_motif_scores.py           # added
├── results/
│   ├── BET_posi.homer.pos
│   ├── FLU_posi.homer.pos
│   ├── BET_posi_CREs.tsv
│   ├── FLU_posi_CREs.tsv
│   ├── BET_posi_motif/
│   │   ├── knownResults.txt
│   │   └── knownResults.html
│   └── FLU_posi_motif/
│       ├── knownResults.txt
│       └── knownResults.html
├── environment/
├── logs/
├── README.md
└── annotatePeaks_补充位置与文件说明.md
```

Other archived results (including gene tables, unmatched genes, and the ATAC background) are retained in `results/`. The `knownResults/` **matrix directories shown in the Step 05 instructions are not present** in this ZIP.

The `data/` directory lists the required local inputs. **In the supplied repository ZIP, `FLU_posi.txt` and `Danio_rerio.GRCz11.115.gtf` are absent**, so Steps 01–02 cannot be fully rerun from the ZIP alone. Step 05 can use existing `*.homer.pos` files without these inputs. Availability, provenance, and redistribution permissions for reference data should be confirmed before distributing files through GitHub.

## Module 04: CRE association, motif enrichment, and motif annotation

### Input data

| File | Description |
| --- | --- |
| `BET_posi.txt` | Selected betamethasone (BET) gene identifiers, one per line |
| `FLU_posi.txt` | Selected fludrocortisone (FLU) gene identifiers, one per line |
| `Danio_rerio.GRCz11.115.gtf` | Zebrafish gene annotations for the GRCz11 assembly |
| `ATAC_GRCz11_peaks_434221.bed` | ATAC-seq accessible chromatin regions mapped to GRCz11 |

### Software environment

The analyses were run on a Linux high-performance computing (HPC) cluster using Conda-managed software. The `environment/` directory contains two exports of the analysis environment:

- `meme_environment.yml`: Conda environment specification exported from explicitly requested packages (`--from-history`).
- `meme_environment_full.yml`: Full export of the installed Conda environment for recording software dependencies and versions.

The analysis requires Python 3 with `pandas`, `pybedtools`, and `tqdm`; BEDTools; HOMER with the `danRer11` genome installed; and SLURM for the provided batch script.

To create an environment from the portable specification, run:

```bash
conda env create -f environment/meme_environment.yml
conda activate meme
```

The environment name and any cluster-specific paths may need adjustment. A `--from-history` export may omit packages installed with `pip` or otherwise not explicitly recorded by Conda. Installation of HOMER and its genome resources should therefore be confirmed separately. The full environment export is an archival record rather than a guaranteed portable installation file.

### Analysis workflow

**Step 01 — Gene annotation:** Match input identifiers to `gene_id` or `gene_name` in the GRCz11 GTF and export gene coordinates. Coordinates remain 1-based and inclusive in this output.

**Step 02 — Gene–ATAC associations:** Convert gene positions to 0-based, half-open BED coordinates and use BEDTools `intersect -wa -wb` to retain ATAC peaks overlapping each gene body. These are gene-body overlaps, **not** specifically TSS-centered or promoter-only peaks.

**Step 03 — HOMER input preparation:** Export each unique associated ATAC peak as a six-column position file for HOMER. The script also generates `ATAC_background.homer.pos` from all unique ATAC peaks for recordkeeping or alternative analyses; **this background file is not used in the original-parameter HOMER run described below**.

**Step 04 — Motif enrichment:** Run HOMER against the `danRer11` reference genome with `-size given -p 8 -mask`, without a user-specified `-bg` option. Thus, HOMER uses its internally generated/default background. The `-size given` option retains the supplied genomic interval lengths.

**Step 05 — Motif occurrence annotation:** Use `annotatePeaks.pl` with `-size given -mask -m ... -mscore -mbed` to scan the same ATAC peak intervals for motif occurrences and report the highest motif log-odds score for each region. Motif inputs must be the **actual HOMER matrices**, normally `knownResults/known1.motif`, `known2.motif`, etc. The `knownResults.txt` and HTML enrichment tables are **not** motif matrices. By default the script scans every known motif file and requires a complete set; `MOTIF_RANKS` allows explicitly selected ranks for focused tests. The script uses `-noann -nogene` because candidate peak–gene links are defined by Step 02, not HOMER nearest-gene annotation.

**Step 06 — Score filtering and candidate links:** Extract `Best Motif log-odds Score` columns from the Step 05 output; retain peak–motif pairs with **score > 8** (study-defined threshold). Join the retained peaks back to the `*_posi_CREs.tsv` genomic intervals to produce candidate motif–target gene links. These are **computational predictions**, not proof of direct binding or transcriptional activation/repression.

### Run commands

Run the following commands from the repository root:

```bash
cd /path/to/Phenotype-Guided-Transcriptome-Analysis-main

# Step 01: Gene annotations
python scripts/04_TFBS/01_gene_info_from_gtf_by_gene_name.py --compound BET
python scripts/04_TFBS/01_gene_info_from_gtf_by_gene_name.py --compound FLU

# Step 02: Gene–ATAC overlaps
python scripts/04_TFBS/02_gene2CREs.py --compound BET
python scripts/04_TFBS/02_gene2CREs.py --compound FLU

# Step 03: Target positions (and optional ATAC background file)
python scripts/04_TFBS/03_pos4homer.py --compound BET
python scripts/04_TFBS/03_pos4homer.py --compound FLU

# Step 04: HOMER motif enrichment on SLURM
mkdir -p logs
sbatch --export=ALL,COMPOUND=BET scripts/04_TFBS/04_run_homer.slurm
sbatch --export=ALL,COMPOUND=FLU scripts/04_TFBS/04_run_homer.slurm

# BEFORE Step 05: restore the original motif matrix subdirectories:
#   results/04_TFBS/BET_posi_motif/knownResults/known*.motif
#   results/04_TFBS/FLU_posi_motif/knownResults/known*.motif
# or the same motif subdirectories under results/ (supplied ZIP layout).
# Do not attempt to use knownResults.txt/html in place of *.motif.

# Step 05: scan all known motif matrices and output peak-level log-odds scores
sbatch --export=ALL,COMPOUND=BET scripts/04_TFBS/05_annotate_peaks.slurm
sbatch --export=ALL,COMPOUND=FLU scripts/04_TFBS/05_annotate_peaks.slurm

# Optional GRE-only trial: the supplied knownResults.txt ranks GRE at
# rank 7 for BET and rank 4 for FLU (not the top rank overall).
# sbatch --export=ALL,COMPOUND=BET,MOTIF_RANKS=7 scripts/04_TFBS/05_annotate_peaks.slurm
# sbatch --export=ALL,COMPOUND=FLU,MOTIF_RANKS=4 scripts/04_TFBS/05_annotate_peaks.slurm

# Step 06: run after BOTH annotation jobs finish successfully
python scripts/04_TFBS/06_filter_motif_scores.py --compound BET --threshold 8
python scripts/04_TFBS/06_filter_motif_scores.py --compound FLU --threshold 8
```

**Cluster note:** Both batch scripts locate the repository root from their own script paths. The SLURM partition and Conda activation paths in `04_run_homer.slurm` and `05_annotate_peaks.slurm` are specific to the originating cluster and must be adjusted on other systems. The scripts begin with `#!/bin/bash`.

**Path note:** Preprocessing scripts write to `results/04_TFBS/`, while the supplied ZIP stores archived outputs directly under `results/`. The archived `*_posit_gene_table.tsv` filenames also differ from the Step 01/02 scripts' expected `*_posi_gene_table.tsv` names; normalize filenames before rerunning Step 02.  Steps 05–06 accept either layout. Do not mix BET or FLU inputs from different runs. In the supplied ZIP, the HOMER `knownResults/known*.motif` files are absent; Step 05 **cannot complete** until the original motif files are restored or the corresponding HOMER run is regenerated under a matched environment.

### Expected output

New preprocessing outputs are written under `results/04_TFBS/` by the scripts; the supplied ZIP places several existing results directly in `results/`. Step 05 and Step 06 detect either location.

| Output | Description |
| --- | --- |
| `BET_posi_gene_table.tsv` / `FLU_posi_gene_table.tsv` | Output filenames expected from Step 01 (the supplied ZIP uses `BET_posit_gene_table.tsv` and `FLU_posit_gene_table.tsv` instead) |
| `BET_unmatched_genes.txt` / `FLU_unmatched_genes.txt` | Unmatched input identifiers |
| `BET_posi_CREs.tsv` / `FLU_posi_CREs.tsv` | Gene–ATAC peak overlap records |
| `BET_posi.homer.pos` / `FLU_posi.homer.pos` | Unique target ATAC regions for HOMER |
| `ATAC_background.homer.pos` | Optional full ATAC peak background file; not used in the original-parameter HOMER run |
| `BET_posi_motif/` / `FLU_posi_motif/` | HOMER motif enrichment results |
| `*/knownResults/known*.motif` | Required motif position-weight matrices; absent from the supplied ZIP |
| `*/annotation/*_annotated_motif_scores.tsv` | Per-region HOMER best motif scores (Step 05) |
| `*/annotation/*_motif_sites.bed` | HOMER-predicted genomic motif sites (Step 05) |
| `*/annotation/*_motif_files_used.tsv` | Explicit ordered list of motif matrices (Step 05) |
| `*/annotation/*_high_score_peak_motif.tsv` | Motif–peak scores > 8 (Step 06) |
| `*/annotation/*_candidate_TF_gene_links.tsv` | Motif–peak–gene candidate associations (Step 06) |

### Reproducibility considerations

The HOMER command used for the current reproduction follows the archived BET run parameters (`-size given -p 8 -mask`, without `-bg`). The repository preprocessing implementation differs from some descriptions of the earlier analytical workflow in two respects: (i) GTF gene starts are converted correctly to BED coordinates during overlap analysis; and (ii) duplicate ATAC peak intervals are removed when preparing HOMER targets, with peak-specific identifiers assigned in place of gene identifiers. These changes can affect the selected input regions and the resulting motif enrichment statistics. Accordingly, numerical equivalence to the manuscript results should be assessed against the original input files and archived output before claiming exact reproduction. Matching command-line parameters alone does not establish equivalent results.

**Motif-specific interpretation:** In the archived `knownResults.txt`, the RAW264.7 GRE motif ranks 7th in BET (q = 0.0221) and 4th in FLU (q = 0.0027); these are not overall first-place motifs. If reporting "GRE ranked first", specify the comparison set and verify that it refers to an analysis that was actually performed. The *enrichment q-value*, *motif match score*, and *motif rank* are distinct statistics.

**Re-run caveat:** Step 05/06 were reconstructed as a supplementary reproducibility workflow. The current archive has no Step 05 results, so there is no evidence in the ZIP that its outputs exactly match the motif-to-gene network in the manuscript. Retain the original motif matrices and software version, and compare derived predictions with original records before claiming exact numerical reproduction.

## Citation

Citation details will be added after the manuscript metadata are finalized.

## License

A license has not yet been selected.
