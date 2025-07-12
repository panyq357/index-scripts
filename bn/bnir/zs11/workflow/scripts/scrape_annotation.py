from pathlib import Path

import pandas as pd

from selenium import webdriver
from selenium.webdriver.chrome.service import Service
from selenium.webdriver.chrome.options import Options

complete = Path("resources/scrape_annotation/complete")
# complete = Path(snakemake.output[0])

service = Service("/home/panyq/Tools/chromedriver-linux64/chromedriver")

options = Options()
prefs = {
   "download.default_directory": str(complete.parent),
   "savefile.default_directory": str(complete.parent)
}

options.add_experimental_option("prefs", prefs)

driver = webdriver.Chrome(service=service, options=options)

chromosomes = [
    "A01",
    "A02",
    "A03",
    "A04",
    "A05",
    "A06",
    "A07",
    "A08",
    "A09",
    "A10",
    "C01",
    "C02",
    "C03",
    "C04",
    "C05",
    "C06",
    "C07",
    "C08",
    "C09"
]

for chromosome in chromosomes:
    driver.get(f"https://yanglab.hzau.edu.cn/BnIR/search_result?id={chromosome}:0..999999999")
    driver.find_element(by="xpath", value="/html/body/div[1]/div[2]/div[4]/div/div[1]/button").click()

df_list = []
for xlsx_path in Path("/home/panyq/Downloads").glob("gene basic information*"):
    df_list.append(pd.read_excel(xlsx_path, skiprows=1))

all_anno = pd.concat(df_list, axis=0).sort_values("ZS11 Gene ID")

all_anno.to_csv("results/bnir_scrape_annotation.csv", index=False)
