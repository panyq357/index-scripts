import gzip
import re

RAP = re.compile(r"Os[0-9]{2}g[0-9]{7}")
MSU = re.compile(r"LOC_Os[0-9]{2}g[0-9]{5}")

with gzip.open("rawdata/RAP-MSU_2023-09-07.txt.gz", "rt") as f:
    lines = f.readlines()

out_lines = []
for line in lines:
    rap_list = RAP.findall(line)
    if len(rap_list) == 1:
        rap = rap_list[0]
    elif len(rap_list) == 0:
        rap = None
    else:
        raise Exception(f"Multiple RAP ID found in line: {line}")
    msu_list = MSU.findall(line)
    for msu in msu_list:
        out_lines.append(f"{rap}\t{msu}\n")
out_lines = sorted(list(set(out_lines)))

with open("results/id_converter.txt", "wt") as f:
    f.writelines(out_lines)

