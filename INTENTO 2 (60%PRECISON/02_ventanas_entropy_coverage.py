#!/usr/bin/env python3

"""
Lee formicidae_col_stats.tsv (columna, coverage, entropy)
y calcula ventanas deslizantes con:
- cobertura media
- entropía media

Devuelve las mejores ventanas según entropía media,
filtrando por una cobertura media mínima.
"""

STATS_FILE = "formicidae_col_stats.tsv"

WIN_SIZE = 100       # tamaño de ventana en columnas
STEP = 20            # paso entre ventanas
MIN_MEAN_COV = 0.70  # cobertura media mínima para considerar una ventana

OUT_FILE = "ventanas_entropy_coverage.tsv"

print(f"Usando archivo de stats: {STATS_FILE}")
print(f"Tamaño de ventana: {WIN_SIZE}, paso: {STEP}")
print(f"Cobertura media mínima: {MIN_MEAN_COV}")

# 1) Leer stats por columna
cols = []
covs = []
ents = []

with open(STATS_FILE) as f:
    header = f.readline().strip().split("\t")
    # esperamos: columna\tcoverage\tentropy
    for line in f:
        if not line.strip():
            continue
        c, cov, ent = line.strip().split("\t")
        cols.append(int(c))
        covs.append(float(cov))
        ents.append(float(ent))

if not cols:
    raise SystemExit("No se pudieron leer datos de formicidae_col_stats.tsv")

n = len(cols)
print(f"Columnas leídas: {n}")

# 2) Sliding window
windows = []

for start_idx in range(0, n - WIN_SIZE + 1, STEP):
    end_idx = start_idx + WIN_SIZE  # índice exclusivo
    win_cov = covs[start_idx:end_idx]
    win_ent = ents[start_idx:end_idx]
    mean_cov = sum(win_cov) / len(win_cov)
    mean_ent = sum(win_ent) / len(win_ent)

    if mean_cov >= MIN_MEAN_COV:
        col_start = cols[start_idx]
        col_end = cols[end_idx - 1]
        windows.append((col_start, col_end, mean_cov, mean_ent))

# 3) Ordenar por entropía media (descendente)
windows.sort(key=lambda x: x[3], reverse=True)

if not windows:
    raise SystemExit("No hay ninguna ventana que cumpla el criterio de cobertura. Probá bajar MIN_MEAN_COV.")

# 4) Escribir resultados
with open(OUT_FILE, "w") as out:
    out.write("col_start\tcol_end\tmean_coverage\tmean_entropy\n")
    for (cs, ce, mc, me) in windows:
        out.write(f"{cs}\t{ce}\t{mc:.5f}\t{me:.5f}\n")

print(f"Resultados escritos en: {OUT_FILE}")

print("Top 10 ventanas:")
for (cs, ce, mc, me) in windows[:10]:
    print(f"  {cs}-{ce}  cov={mc:.3f}  H={me:.3f}")
