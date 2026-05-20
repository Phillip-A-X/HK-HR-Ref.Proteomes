from pathlib import Path
import pandas as pd

# Input files
hrt_file = Path("data/raw/Housekeeping_GenesHuman_HRTatlas.tsv")
unmapped_file = Path("data/processed/hrt_atlas_unmapped_enst_from_biomart.tsv")

# Output
output_file = Path("data/processed/hrt_atlas_unmapped_enst_with_alternative_ids.tsv")

# Original HRT Atlas file has no header
hrt = pd.read_csv(
    hrt_file,
    sep="\t",
    header=None,
    names=["enst", "gene_symbol", "refseq", "ccds_id"]
)

unmapped = pd.read_csv(unmapped_file, sep="\t")

# Falls die Spalte anders heißt, erste Spalte nehmen
unmapped_col = unmapped.columns[0]

result = hrt[hrt["enst"].isin(unmapped[unmapped_col])].copy()

result.to_csv(output_file, sep="\t", index=False)

print(result)
print(f"Written to: {output_file}")