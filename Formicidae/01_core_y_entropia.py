#!/usr/bin/env python3

"""
Script para:
1) Leer un FASTA alineado de COI (Formicidae).
2) Calcular cobertura y entropía de Shannon por columna.
3) Encontrar el bloque "core" con cobertura >= coverage_threshold.
4) Escribir:
   - el alineamiento recortado al core,
   - una tabla TSV con cobertura y entropía por columna.
"""

import math
from collections import Counter

# ==========================
# CONFIGURACIÓN
# ==========================

# Nombre del archivo de alineamiento (FASTA) que ya tenés
ALN_FILE = "formicidae_ge600_aln.fasta"

# Umbral de cobertura mínimo (0.85 = 85%)
COVERAGE_THRESHOLD = 0.85

# Archivos de salida
CORE_OUT_FASTA = "formicidae_core_aln.fasta"
COL_STATS_TSV = "formicidae_col_stats.tsv"

print(f"Usando alineamiento: {ALN_FILE}")
print(f"Umbral de cobertura: {COVERAGE_THRESHOLD*100:.1f}%")

# ==========================
# 1) LEER EL ALINEAMIENTO
# ==========================

ids = []
seqs = []

with open(ALN_FILE) as f:
    current_id = None
    current_seq = []

    for line in f:
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            # Si estábamos leyendo una secuencia anterior, la guardamos
            if current_id is not None:
                seqs.append("".join(current_seq))
            # Nuevo ID
            current_id = line[1:].split()[0]
            ids.append(current_id)
            current_seq = []
        else:
            current_seq.append(line)

    # Ultima secuencia
    if current_id is not None:
        seqs.append("".join(current_seq))

if not seqs:
    raise SystemExit("ERROR: No se leyeron secuencias del archivo de alineamiento.")

n_seq = len(seqs)
L = len(seqs[0])

print(f"Secuencias leídas: {n_seq}")
print(f"Longitud del alineamiento: {L} columnas")

# Verificar que todas las secuencias tengan igual longitud
for sid, s in zip(ids, seqs):
    if len(s) != L:
        raise SystemExit(f"ERROR: La secuencia {sid} tiene longitud {len(s)} distinta de {L}.")

# ==========================
# 2) FUNCIONES AUXILIARES
# ==========================

def shannon_entropy(bases):
    """
    Calcula la entropía de Shannon (en bits) de una lista de bases (A/C/G/T).
    Ignora gaps: se espera que 'bases' ya no contenga '-'.
    """
    if not bases:
        return 0.0
    counts = Counter(bases)
    total = sum(counts.values())
    H = 0.0
    for base, c in counts.items():
        p = c / total
        H -= p * math.log2(p)
    return H

# ==========================
# 3) COBERTURA Y ENTROPÍA POR COLUMNA
# ==========================

coverages = []
entropies = []

for pos in range(L):
    column = [s[pos] for s in seqs]
    non_gaps = [b for b in column if b != "-"]
    cov = len(non_gaps) / n_seq
    H = shannon_entropy(non_gaps)
    coverages.append(cov)
    entropies.append(H)

max_H = max(entropies)
min_H = min(entropies)
avg_H = sum(entropies) / len(entropies)

print(f"Entropía mínima por columna: {min_H:.3f}")
print(f"Entropía máxima por columna: {max_H:.3f}")
print(f"Entropía promedio por columna: {avg_H:.3f}")

# ==========================
# 4) ENCONTRAR EL BLOQUE CORE POR COBERTURA
# ==========================

good = [cov >= COVERAGE_THRESHOLD for cov in coverages]

best_start = 0
best_end = 0  # rango [best_start, best_end)
current_start = None

for i, is_good in enumerate(good):
    if is_good and current_start is None:
        # Comienza un bloque bueno
        current_start = i
    if (not is_good or i == L - 1) and current_start is not None:
        # Terminó un bloque bueno
        end = i if not is_good else i + 1
        if end - current_start > best_end - best_start:
            best_start, best_end = current_start, end
        current_start = None

core_length = best_end - best_start

if core_length <= 0:
    raise SystemExit(
        "ERROR: No se encontró ningún bloque 'core' con la cobertura indicada. "
        "Probá bajar COVERAGE_THRESHOLD (por ejemplo a 0.80)."
    )

print(f"Core encontrado: columnas {best_start} - {best_end} (longitud: {core_length})")
print(f"Proporción de alineamiento conservada como core: {core_length / L:.2%}")

# ==========================
# 5) ESCRIBIR ALINEAMIENTO CORE
# ==========================

with open(CORE_OUT_FASTA, "w") as out:
    for sid, s in zip(ids, seqs):
        sub = s[best_start:best_end]
        out.write(f">{sid}\n{sub}\n")

print(f"Alineamiento core escrito en: {CORE_OUT_FASTA}")

# ==========================
# 6) ESCRIBIR TABLA POR COLUMNA
# ==========================

with open(COL_STATS_TSV, "w") as out:
    out.write("columna\tcoverage\tentropy\n")
    for i, (cov, H) in enumerate(zip(coverages, entropies)):
        out.write(f"{i}\t{cov:.5f}\t{H:.5f}\n")

print(f"Tabla de cobertura y entropía por columna escrita en: {COL_STATS_TSV}")
print("Listo. Ahora podés mirar el core y analizar las columnas/ventanas más informativas.")

