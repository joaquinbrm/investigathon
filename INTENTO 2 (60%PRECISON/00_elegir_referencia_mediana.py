#!/usr/bin/env python3

"""
Elige una secuencia de referencia desde formicidae_ge600.fasta.
Criterio: la secuencia cuya longitud esté más cerca de la mediana.
Escribe un archivo COI_ref.fasta con ID 'COI_REF'.
"""

IN_FILE = "formicidae_ge600.fasta"
OUT_FILE = "COI_ref.fasta"

# 1) Leer todas las secuencias y sus longitudes
seqs = []  # lista de (id, seq)
lengths = []

with open(IN_FILE) as f:
    sid = None
    seq = []
    for line in f:
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if sid is not None:
                full = "".join(seq)
                seqs.append((sid, full))
                lengths.append(len(full))
            sid = line[1:].split()[0]
            seq = []
        else:
            seq.append(line)
    if sid is not None:
        full = "".join(seq)
        seqs.append((sid, full))
        lengths.append(len(full))

if not seqs:
    raise SystemExit("No se leyeron secuencias en formicidae_ge600.fasta")

# 2) Calcular la mediana de longitudes
lengths_sorted = sorted(lengths)
n = len(lengths_sorted)
if n % 2 == 1:
    median_len = lengths_sorted[n // 2]
else:
    median_len = (lengths_sorted[n // 2 - 1] + lengths_sorted[n // 2]) / 2

print(f"Mediana de longitudes: {median_len:.1f} bp")

# 3) Elegir la secuencia cuya longitud esté más cerca de la mediana
best_sid = None
best_seq = None
best_diff = None

for (sid, seq) in seqs:
    L = len(seq)
    diff = abs(L - median_len)
    if best_diff is None or diff < best_diff:
        best_diff = diff
        best_sid = sid
        best_seq = seq

print(f"Secuencia elegida como referencia: {best_sid}")
print(f"Longitud: {len(best_seq)} bp (diferencia con la mediana: {best_diff:.1f} bp)")

# 4) Escribir como COI_REF
with open(OUT_FILE, "w") as out:
    out.write(">COI_REF\n")
    for i in range(0, len(best_seq), 60):
        out.write(best_seq[i:i+60] + "\n")

print(f"Referencia escrita en: {OUT_FILE} con ID 'COI_REF'")

