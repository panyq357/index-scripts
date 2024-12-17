import requests
import pandas as pd

from retrying import retry

config = {
    "organism_code": snakemake.params["organism_code"],
    "pathway_to_description": snakemake.output["pathway_to_description"],
    "pathway_to_gene": snakemake.output["pathway_to_gene"],
    "omit_str_in_description": snakemake.params["omit_str_in_description"] if "omit_str_in_description" in snakemake.params.keys() else False
}

def main():

    pathway = pd.read_table(f"https://rest.kegg.jp/list/pathway/{config['organism_code']}", header=None)
    pathway.columns = ["Pathway", "Description"]

    if config["omit_str_in_description"]:
        pathway["Description"] = pathway["Description"].str.replace(config["omit_str_in_description"], "", regex=False)

    path_to_gene = pd.DataFrame(columns=["Pathway", "Gene"])

    @retry(wait_fixed=2000)
    def get_p_to_g(p):
        p_to_g = pd.read_table(f"https://rest.kegg.jp/link/{config['organism_code']}/path:{p}", header=None)
        return(p_to_g)

    for i, p in enumerate(pathway["Pathway"]):
        print(f"Processing {p}, {int(i/len(pathway['Pathway']) * 100)}%.")
        p_to_g = get_p_to_g(p)
        p_to_g.columns = ["Pathway", "Gene"]
        path_to_gene = pd.concat([path_to_gene, p_to_g], axis=0)

    path_to_gene["Pathway"] = path_to_gene["Pathway"].str.replace("path:", "")
    path_to_gene["Gene"] = path_to_gene["Gene"].str.replace(f"{config['organism_code']}:", "")

    pathway.to_csv(config["pathway_to_description"], header=False, index=False, sep="\t")
    path_to_gene.to_csv(config["pathway_to_gene"], header=False, index=False, sep="\t")


if __name__ == "__main__":
    main()
