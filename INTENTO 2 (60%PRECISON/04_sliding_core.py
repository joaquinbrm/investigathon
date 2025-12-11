#!/usr/bin/env python3

import math
from collections import Counter

# Archivo sin gaps del core
CORE_FILE = "formicidae_core_nogap.fasta"

# Tamaños de ventana a evaluar
WINDOW_SIZES = [30, 40, 50, 60]

# Leer secuencias
ids = []
seqs = []

with open(CORE_FILE) as f:
    seq = ""
    sid = None
    for line in f:
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if sid is not None:
                seqs.append(seq)
            sid = line[1:].split()[0]
            ids.append(sid)
            seq = ""
        else:
            seq += line
    if sid is not None:
        seqs.append(seq)

if not seqs:
    raise SystemExit("No se leyeron secuencias del core.")

n = len(seqs)
L = len(seqs[0])
print(f"Secuencias: {n}, longitud del core sin gaps: {L}")

# Función de entropía
def H_shannon(chars):
    counts = Counter(chars)
    total = sum(counts.values())
    H = 0
    for c in counts.values():
        p = c / total
        H -= p * math.log2(p)
    return H

# Sliding window por cada tamaño
for WIN in WINDOW_SIZES:
    print("\n=====================================")
    print(f"Analizando ventana de {WIN} bp")
    print("=====================================")

    results = []

    for start in range(0, L - WIN + 1):
        end = start + WIN
        window_entropy = 0

        # calcular entropía promedio por columna
        for col in range(start, end):
            bases = [s[col] for s in seqs]
            H = H_shannon(bases)
            window_entropy += H

        window_entropy /= WIN

        results.append((start, end, window_entropy))

    # ordenar por entropía descendente
    results.sort(key=lambda x: x[2], reverse=True)

    # mostrar top 5
    print("Top 5 ventanas:")
    for s, e, Hm in results[:5]:
        print(f"  {s}-{e}  Hmedio={Hm:.4f}")

