import json
import shutil

from pathlib import Path

fa = Path(snakemake.input.fa)
fai = Path(snakemake.input.fai)
gff = Path(snakemake.input.gff)

outdir = Path(snakemake.output.outdir)

if not outdir.exists():
    outdir.mkdir(parents=True)

shutil.copyfile(fa, outdir / fa.name)
shutil.copyfile(fai, outdir / fai.name)
shutil.copyfile(gff, outdir / gff.name)

out = {
  "id": snakemake.params.id,
  "name": snakemake.params.name,
  "fastaURL": fa.name,
  "indexURL": fai.name,
#   "chromosomeOrder": snakemake.params.chromosome_order,
  "tracks": [
      {
          "name": "GFF",
          "format": "gff",
          "url": gff.name
      }
  ]
}

with open(snakemake.output.json, "wt") as f:
    json.dump(out, f)
