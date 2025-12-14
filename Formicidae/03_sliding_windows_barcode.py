#!/usr/bin/env python3

import csv

# Archivo con las estadísticas por columna
COL_FILE = "col_stats.csv"

# Tamaños de ventana a evaluar (en columnas)
WINDOW_SIZES = [20, 30, 40]

# Leer estadísticas por columna
cols = []  # lista de dicts: {idx, cobertura, identidad, entropia}

with open(COL_FILE) as f:
    reader = csv.DictReader(f)
    for row in reader:
        idx = int(row["columna"])
        cov = float(row["cobertura"])
        ident = float(row["identidad"])
        entro = float(row["entropia"])
        cols.append({
            "idx": idx,
            "cov": cov,
            "ident": ident,
            "entro": entro
        })

L = len(cols)
print(f"Columnas totales: {L}")

# Para cada tamaño de ventana hacemos sliding
for WIN in WINDOW_SIZES:
    if WIN > L:
        print(f"\nVentana de {WIN} columnas es mayor que la longitud ({L}), se omite.")
        continue

    print("\n" + "="*50)
    print(f"Analizando ventanas de {WIN} columnas")
    print("="*50)

    ventanas = []

    for start in range(0, L - WIN + 1):
        end = start + WIN  # exclusivo

        sub = cols[start:end]

        # medidas agregadas
        mean_cov = sum(c["cov"] for c in sub) / WIN
        mean_ident = sum(c["ident"] for c in sub) / WIN
        mean_entro = sum(c["entro"] for c in sub) / WIN

        # Podés definir un puntaje combinado si querés:
        # por ejemplo favorecer alta entropía y no demasiada identidad
        # aquí sólo reportamos y ordenamos por entropía
        ventanas.append({
            "start": start,
            "end": end,
            "mean_cov": mean_cov,
            "mean_ident": mean_ident,
            "mean_entro": mean_entro
        })

    # Ordenamos por entropía media descendente (más variabilidad primero)
    ventanas.sort(key=lambda x: x["mean_entro"], reverse=True)

    # Mostramos las top 5
    print("Top 5 ventanas por entropía media:")
    for v in ventanas[:5]:
        s = v["start"]
        e = v["end"]
        print(
            f"  columnas {s}-{e-1} (len={e-s})  "
            f"Hmean={v['mean_entro']:.4f}  "
            f"Ident_mean={v['mean_ident']:.4f}  "
            f"Cov_mean={v['mean_cov']:.4f}"
        )

