#!/usr/bin/env python3

from collections import Counter
import math

IN_FILE = "formicidae_trimmed_ref.fasta"
OUT_FILE = "col_stats.csv"

# Leer alineamiento recortado
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

n = len(seqs)
L = len(seqs[0])

print(f"Secuencias: {n}, Longitud: {L}")

def entropia(lista_bases):
    if not lista_bases:
        return 0
    c = Counter(lista_bases)
    total = sum(c.values())
    H = 0
    for v in c.values():
        p = v / total
        H -= p * math.log2(p)
    return H

# Abrimos salida
out = open(OUT_FILE, "w")
out.write("columna,cobertura,identidad,entropia\n")

# Analizamos columna por columna
for col in range(L):
    bases = [s[col] for s in seqs if s[col] != "-"]
    if len(bases) == 0:
        cobertura = 0
        identidad = 0
        H = 0
    else:
        cobertura = len(bases) / n
        freqs = Counter(bases)
        mayor = freqs.most_common(1)[0][1]
        identidad = mayor / len(bases)
        H = entropia(bases)

    out.write(f"{col},{cobertura:.4f},{identidad:.4f},{H:.4f}\n")

out.close()

print(f"Estadísticas escritas en {OUT_FILE}")

