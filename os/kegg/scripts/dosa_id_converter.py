from io import StringIO

import requests
import pandas as pd


id_converter = pd.DataFrame()

# Fetch one chromosome at a time.
for i in range(1, 13):
    response = requests.get("https://rest.kegg.jp/find/genes/dosa:Os%02dt" % i)
    df = pd.read_table(StringIO(response.text), header=None)
    tx_id = df[0].str.replace("dosa:", "")
    gene_id = tx_id.str.replace("t", "g").str.extract(".*(Os[0-9]{2}g[0-9]{7}).*")
    id_converter = pd.concat([id_converter, pd.concat([tx_id, gene_id], axis=1)], axis=0)

id_converter.to_csv(snakemake.output[0], header=False, index=False, sep="\t")

