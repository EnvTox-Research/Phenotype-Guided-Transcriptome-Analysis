#!/bin/bash
#SBATCH --job-name=EDTA_annotate
#SBATCH --partition=64c512g
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=08:00:00
#SBATCH --output=logs/EDTA_annotate_%j.out
#SBATCH --error=logs/EDTA_annotate_%j.err

# Step 05: HOMER motif occurrence annotation and best log-odds scores.
# Example:
#   sbatch --export=ALL,COMPOUND=BET scripts/04_TFBS/05_annotate_peaks.slurm
# Optional GRE-only test (BET=rank 7; FLU=rank 4 in archived reports):
#   sbatch --export=ALL,COMPOUND=BET,MOTIF_RANKS=7 scripts/04_TFBS/05_annotate_peaks.slurm
#   sbatch --export=ALL,COMPOUND=FLU,MOTIF_RANKS=4 scripts/04_TFBS/05_annotate_peaks.slurm
# For publication-wide scanning, omit MOTIF_RANKS (all known motifs).

set -euo pipefail

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$PROJECT_DIR"

# Same HPC environment as Step 04. Adjust these paths if your environment differs.
source /dssg/home/acct-esezk/esezk-wangpeng6/.conda/etc/profile.d/conda.sh 2>/dev/null || true
module load miniconda3
source activate /dssg/home/acct-esezk/esezk-wangpeng6/.conda/envs/meme/

COMPOUND="${COMPOUND:?Set COMPOUND=BET or FLU in sbatch --export}"
if [[ "$COMPOUND" != "BET" && "$COMPOUND" != "FLU" ]]; then
    echo "Error: COMPOUND must be BET or FLU." >&2
    exit 1
fi
GENOME="danRer11"

# Use the directory written by scripts 01-04, or fall back to
# the results/ layout in the supplied ZIP archive.
if [[ -s "results/04_TFBS/${COMPOUND}_posi.homer.pos" ]]; then
    INPUT_DIR="$PROJECT_DIR/results/04_TFBS"
elif [[ -s "results/${COMPOUND}_posi.homer.pos" ]]; then
    INPUT_DIR="$PROJECT_DIR/results"
else
    echo "Error: ${COMPOUND}_posi.homer.pos was not found." >&2
    exit 1
fi

TARGET="$INPUT_DIR/${COMPOUND}_posi.homer.pos"
MOTIF_OUTPUT="$INPUT_DIR/${COMPOUND}_posi_motif"
MOTIF_DIR="$MOTIF_OUTPUT/knownResults"
ANNOTATION_DIR="$MOTIF_OUTPUT/annotation"

if [[ ! -d "$MOTIF_DIR" ]]; then
    echo "Error: missing $MOTIF_DIR" >&2
    echo "Restore the original knownResults/known*.motif files from HOMER output." >&2
    echo "The supplied ZIP contains only knownResults.txt/html, not .motif matrices." >&2
    exit 1
fi

MOTIFS=()
if [[ -n "${MOTIF_RANKS:-}" ]]; then
    # Explicit subset for focused scanning or a rapid test.
    IFS=',' read -ra RANKS <<< "$MOTIF_RANKS"
    for rank in "${RANKS[@]}"; do
        if [[ ! "$rank" =~ ^[1-9][0-9]*$ ]]; then
            echo "Invalid motif rank: $rank" >&2
            exit 1
        fi
        motif="$MOTIF_DIR/known${rank}.motif"
        if [[ ! -s "$motif" ]]; then
            echo "Error: missing $motif" >&2
            exit 1
        fi
        MOTIFS+=("$motif")
    done
else
    # Full annotation: include all archived known motif matrices.
    shopt -s nullglob
    MOTIFS=("$MOTIF_DIR"/known[0-9]*.motif)
    shopt -u nullglob
    if [[ -s "$MOTIF_OUTPUT/knownResults.txt" ]]; then
        expected=$(( $(wc -l < "$MOTIF_OUTPUT/knownResults.txt") - 1 ))
        if (( ${#MOTIFS[@]} != expected )); then
            echo "Error: expected $expected known motifs, found ${#MOTIFS[@]}." >&2
            echo "Do not run a partial panel as if it were a full scan." >&2
            exit 1
        fi
    fi
fi

if (( ${#MOTIFS[@]} == 0 )); then
    echo "Error: no motif matrix files were selected." >&2
    exit 1
fi

# Reproducible numeric ordering: known2 before known10.
mapfile -t MOTIFS < <(printf '%s\n' "${MOTIFS[@]}" | sort -V -u)
mkdir -p "$ANNOTATION_DIR" "$PROJECT_DIR/logs"

if ! command -v annotatePeaks.pl >/dev/null 2>&1; then
    echo "Error: annotatePeaks.pl is not available after activating the environment." >&2
    exit 1
fi

MANIFEST="$ANNOTATION_DIR/${COMPOUND}_motif_files_used.tsv"
printf 'motif_order\tmotif_rank\tmotif_file\tmotif_header\n' > "$MANIFEST"
index=0
for motif in "${MOTIFS[@]}"; do
    index=$((index + 1))
    filename="$(basename "$motif")"
    rank="${filename#known}"
    rank="${rank%.motif}"
    header="$(head -n 1 "$motif" | tr '\t' ' ' | tr '\r\n' '  ')"
    printf '%s\t%s\t%s\t%s\n' "$index" "$rank" "$filename" "$header" >> "$MANIFEST"
done

echo "Compound: $COMPOUND"
echo "Regions: $TARGET"
echo "Motifs selected: ${#MOTIFS[@]}"
echo "Input motif folder: $MOTIF_DIR"
echo "Scanning exact supplied intervals (-size given) with repeat masking (-mask)."

# -noann -nogene: do not substitute HOMER nearest-gene assignments
# for the gene-peak relationships defined in 02_gene2CREs.py.
annotatePeaks.pl \
    "$TARGET" \
    "$GENOME" \
    -size given \
    -mask \
    -noann \
    -nogene \
    -m "${MOTIFS[@]}" \
    -mscore \
    -mbed "$ANNOTATION_DIR/${COMPOUND}_motif_sites.bed" \
    > "$ANNOTATION_DIR/${COMPOUND}_annotated_motif_scores.tsv"

echo "Annotation complete: $ANNOTATION_DIR"
echo "Next: python scripts/04_TFBS/06_filter_motif_scores.py --compound $COMPOUND --threshold 8"
