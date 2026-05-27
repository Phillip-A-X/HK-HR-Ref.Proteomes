from pathlib import Path
import pandas as pd

input_file = Path("data/processed/human_all_genes_with_hrt_atlas_hk_status_rescued.tsv")
output_file = Path("data/processed/hrt_atlas_final_hk_ensg_ids_for_biomart.txt")

df = pd.read_csv(input_file, sep="\t")

hk_genes = (
    df.loc[df["hk_status"] == "HK", "Gene stable ID"]
    .dropna()
    .drop_duplicates()
    .sort_values()
)

hk_genes.to_csv(output_file, index=False, header=False)

print(f"HK ENSG IDs written: {len(hk_genes)}")
print(f"Output file: {output_file}")