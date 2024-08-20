import requests
import pandas as pd

from retrying import retry

osa_pathway = pd.read_table("https://rest.kegg.jp/list/pathway/osa", header=None)
osa_pathway.columns = ["Pathway", "Description"]

osa_path_to_gene = pd.DataFrame(columns=["Pathway", "Gene"])

@retry(wait_fixed=2000)
def get_p_to_g(p):
    p_to_g = pd.read_table(f"https://rest.kegg.jp/link/osa/path:{p}", header=None)
    return(p_to_g)

for i, p in enumerate(osa_pathway["Pathway"]):
    print(f"Processing {p}, {int(i/len(osa_pathway['Pathway']) * 100)}%.")
    p_to_g = get_p_to_g(p)
    p_to_g.columns = ["Pathway", "Gene"]
    osa_path_to_gene = pd.concat([osa_path_to_gene, p_to_g], axis=0)

osa_pathway["Description"] = osa_pathway["Description"].str.replace(" - Oryza sativa japonica (Japanese rice)", "")
osa_path_to_gene["Pathway"] = osa_path_to_gene["Pathway"].str.replace("path:", "")
osa_path_to_gene["Gene"] = osa_path_to_gene["Gene"].str.replace("osa:", "")

osa_pathway.to_csv("results/osa_pathway/pathway_to_description.tsv", header=False, index=False, sep="\t")
osa_path_to_gene.to_csv("results/osa_pathway/pathway_to_gene.tsv", header=False, index=False, sep="\t")

