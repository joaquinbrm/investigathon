#!/usr/bin/env python3

"""
Trunca un alineamiento global usando una secuencia de referencia.

1) Lee formicidae_ref_aln.fasta (alineado con MAFFT, donde una secuencia es COI_REF).
2) Calcula cobertura por columna.
3) Encuentra el bloque contiguo más largo donde:
   - la referencia no tiene gaps
   - la cobertura es >= COV_THRESHOLD
4) Recorta TODAS las secuencias a ese bloque y escribe un FASTA nuevo.
"""

# Configuración
ALN_FILE = "formicidae_ref_aln.fasta"
REF_ID = "COI_REF"
COV_THRESHOLD = 0.80
OUT_FILE = "formicidae_trimmed_ref.fasta"

print(f"Usando alineamiento: {ALN_FILE}")
print(f"ID referencia: {REF_ID}")
print(f"Umbral de cobertura: {COV_THRESHOLD*100:.1f}%")

# Leer alineamiento
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
            if current_id is not None:
                seqs.append("".join(current_seq))
            current_id = line[1:].split()[0]
            ids.append(current_id)
            current_seq = []
        else:
            current_seq.append(line)
    if current_id is not None:
        seqs.append("".join(current_seq))

if not seqs:
    raise SystemExit("ERROR: No se leyeron secuencias del alineamiento.")

n_seq = len(seqs)
L = len(seqs[0])

print(f"Secuencias leídas: {n_seq}")
print(f"Longitud del alineamiento: {L} columnas")

for sid, s in zip(ids, seqs):
    if len(s) != L:
        raise SystemExit(f"ERROR: La secuencia {sid} tiene longitud distinta ({len(s)}) de {L}.")

# Encontrar la referencia
try:
    ref_idx = ids.index(REF_ID)
except ValueError:
    raise SystemExit(f"ERROR: No se encontró la referencia con ID '{REF_ID}'.")

ref_seq = seqs[ref_idx]
print("Referencia encontrada en el índice:", ref_idx)

# Cobertura por columna
coverages = []
for pos in range(L):
    column = [s[pos] for s in seqs]
    non_gaps = [b for b in column if b != "-"]
    cov = len(non_gaps) / n_seq
    coverages.append(cov)

# Encontrar bloque bueno (ref != '-' y cobertura >= umbral)
good = []
for pos in range(L):
    is_good = (ref_seq[pos] != "-") and (coverages[pos] >= COV_THRESHOLD)
    good.append(is_good)

best_start = 0
best_end = 0
current_start = None

for i, g in enumerate(good):
    if g and current_start is None:
        current_start = i
    if (not g or i == L - 1) and current_start is not None:
        end = i if not g else i + 1
        if end - current_start > best_end - best_start:
            best_start, best_end = current_start, end
        current_start = None

core_len = best_end - best_start
if core_len <= 0:
    raise SystemExit(
        "ERROR: No se encontró ningún bloque con ref!=gap y cobertura suficiente. "
        "Probá bajar COV_THRESHOLD."
    )

print(f"Bloque recortado: columnas {best_start} - {best_end} (longitud {core_len})")
print(f"Proporción del alineamiento conservada: {core_len / L:.2%}")

core_covs = coverages[best_start:best_end]
print(f"Cobertura mínima en el bloque: {min(core_covs):.3f}")
print(f"Cobertura media en el bloque: {sum(core_covs)/len(core_covs):.3f}")

# Escribir FASTA recortado
with open(OUT_FILE, "w") as out:
    for sid, s in zip(ids, seqs):
        sub = s[best_start:best_end]
        out.write(f">{sid}\n{sub}\n")

print(f"FASTA recortado escrito en: {OUT_FILE}")

