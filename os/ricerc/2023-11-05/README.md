Build indices.

```bash
snakemake all --cores 24
```

View pipeline.

```bash
snakemake all --dag | dot -Tsvg > dag.svg
```

![dag.svg](dag.svg)
