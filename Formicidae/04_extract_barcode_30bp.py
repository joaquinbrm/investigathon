#!/usr/bin/env python3

"""
Extrae la ventana de 30 bp (columnas 49–78) del core de 82 bp
y genera un FASTA nuevo con el barcode.

Entrada: formicidae_trimmed_ref.fasta  (longitud alineada = 82)
Salida:  formicidae_barcode_30bp.fasta (longitud alineada = 30)
"""

IN_FILE = "formicidae_trimmed_ref.fasta"
OUT_FILE = "formicidae_barcode_30bp.fasta"

START = 49     # columna inicial (0-based, incluida)
END = 79       # columna final exclusiva (49..78 => 30 bp)

ids = []
seqs = []

with open(IN_FILE) as f:
    sid = None
    seq = []
    for line in f:
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if sid is not None:
                seqs.append("".join(seq))
            sid = line[1:].split()[0]
            ids.append(sid)
            seq = []
        else:
            seq.append(line)
    if sid is not None:
        seqs.append("".join(seq))

L = len(seqs[0])
print(f"Longitud del core: {L} columnas")

if END > L:
    raise SystemExit(f"ERROR: END={END} > longitud {L}")

with open(OUT_FILE, "w") as out:
    for sid, s in zip(ids, seqs):
        sub = s[START:END]
        out.write(f">{sid}\n{sub}\n")

print(f"Barcode extraído: columnas {START}–{END-1} (len={END-START})")
print(f"FASTA escrito en: {OUT_FILE}")

