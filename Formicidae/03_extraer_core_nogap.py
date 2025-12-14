#!/usr/bin/env python3

# Convierte formicidae_core_aln.fasta en una versión sin gaps
# para permitir sliding window real sobre el core.

IN_FILE = "formicidae_core_aln.fasta"
OUT_FILE = "formicidae_core_nogap.fasta"

with open(IN_FILE) as f:
    ids = []
    seqs = []
    current_id = None
    current_seq = []

    for line in f:
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if current_id is not None:
                seqs.append("".join(current_seq).replace("-", ""))
            current_id = line[1:].split()[0]
            ids.append(current_id)
            current_seq = []
        else:
            current_seq.append(line)
    if current_id is not None:
        seqs.append("".join(current_seq).replace("-", ""))

with open(OUT_FILE, "w") as out:
    for sid, s in zip(ids, seqs):
        out.write(f">{sid}\n{s}\n")

print(f"Archivo generado: {OUT_FILE}")

