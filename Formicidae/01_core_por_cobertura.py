#!/usr/bin/env python3

from collections import defaultdict

# Configuración
aln_file = "formicidae_core_aln.fasta"     # alineamiento recortado
meta_file = "metadata_formicidae.tsv"      # metadata con seq_id y species
win_size = 60                              # tamaño de ventana (en columnas del alineamiento)
step = 10                                  # avance entre ventanas
max_per_species = 20                       # máximo de secuencias por especie para el cálculo
windows_scores_out = "windows_scores.tsv"  # tabla con métricas por ventana
best_barcode_aln_out = "formicidae_best_barcode_aln.fasta"
best_barcode_ungapped_out = "formicidae_best_barcode.fasta"

print(f"Usando alineamiento core: {aln_file}")
print(f"Usando metadata: {meta_file}")
print(f"Tamaño de ventana: {win_size}, paso: {step}")
print(f"Máximo de secuencias por especie: {max_per_species}")

# 1) Cargar especie por ID
species_by_id = {}
with open(meta_file) as meta:
    header = meta.readline().rstrip("\n").split("\t")
    try:
        id_idx = header.index("seq_id")
        sp_idx = header.index("species")
    except ValueError:
        raise SystemExit("La metadata debe tener columnas 'seq_id' y 'species' separadas por TAB.")

    for line in meta:
        if not line.strip():
            continue
        cols = line.rstrip("\n").split("\t")
        if len(cols) <= max(id_idx, sp_idx):
            continue
        species_by_id[cols[id_idx]] = cols[sp_idx]

if not species_by_id:
    raise SystemExit("No se pudo cargar ninguna especie desde la metadata.")

print(f"Entradas de metadata leídas: {len(species_by_id)}")

# 2) Cargar alineamiento core
ids = []
seqs = []

with open(aln_file) as f:
    seq = ""
    sid = None
    for line in f:
        line = line.rstrip("\n")
        if not line:
            continue
        if line.startswith(">"):
            if sid is not None:
                seqs.append(seq)
            sid = line[1:].strip().split()[0]
            ids.append(sid)
            seq = ""
        else:
            seq += line.strip()
    if sid is not None:
        seqs.append(seq)

if not seqs:
    raise SystemExit("No se leyeron secuencias del alineamiento core.")

n_seq = len(seqs)
L = len(seqs[0])
print(f"Secuencias en el core: {n_seq}")
print(f"Longitud del core: {L} columnas")

for i, s in enumerate(seqs):
    if len(s) != L:
        raise SystemExit(f"La secuencia {ids[i]} tiene longitud distinta en el core.")

# 3) Seleccionar un subconjunto de secuencias por especie para que sea manejable
species_to_indices = defaultdict(list)
for idx, sid in enumerate(ids):
    sp = species_by_id.get(sid)
    if sp is not None:
        species_to_indices[sp].append(idx)

selected_indices = []
for sp, idxs in species_to_indices.items():
    if len(idxs) > max_per_species:
        selected_indices.extend(idxs[:max_per_species])
    else:
        selected_indices.extend(idxs)

selected_indices = sorted(set(selected_indices))

print(f"Especies totales en metadata/alineamiento: {len(species_to_indices)}")
print(f"Secuencias seleccionadas para cálculo: {len(selected_indices)} (máx {max_per_species} por especie)")

# 4) Función de distancia p (ignorando gaps)
def p_distance(a, b):
    dif, comp = 0, 0
    for x, y in zip(a, b):
        if x == "-" or y == "-":
            continue
        comp += 1
        if x != y:
            dif += 1
    if comp == 0:
        return 0.0
    return dif / comp

# 5) Sliding window y cálculo intra / inter
results = []

for start in range(0, L - win_size + 1, step):
    end = start + win_size
    intra_dists = []
    inter_dists = []

    for i_idx in range(len(selected_indices)):
        i = selected_indices[i_idx]
        sid_i = ids[i]
        sp_i = species_by_id.get(sid_i)
        if sp_i is None:
            continue
        frag_i = seqs[i][start:end]

        for j_idx in range(i_idx + 1, len(selected_indices)):
            j = selected_indices[j_idx]
            sid_j = ids[j]
            sp_j = species_by_id.get(sid_j)
            if sp_j is None:
                continue
            frag_j = seqs[j][start:end]
            d = p_distance(frag_i, frag_j)
            if sp_i == sp_j:
                intra_dists.append(d)
            else:
                inter_dists.append(d)

    if intra_dists and inter_dists:
        mean_intra = sum(intra_dists) / len(intra_dists)
        mean_inter = sum(inter_dists) / len(inter_dists)
        ratio = mean_inter / (mean_intra + 1e-6)
        results.append((start, end, mean_intra, mean_inter, ratio))
        print(f"Ventana {start}-{end}: intra={mean_intra:.4f} inter={mean_inter:.4f} ratio={ratio:.2f}")
    else:
        print(f"Ventana {start}-{end}: sin suficientes datos (intra o inter vacíos)")

if not results:
    raise SystemExit("No se obtuvo ninguna ventana con datos intra e inter suficientes.")

# 6) Guardar resultados en tabla
results_sorted = sorted(results, key=lambda x: x[4], reverse=True)

with open(windows_scores_out, "w") as out:
    out.write("start\tend\tmean_intra\tmean_inter\tratio\n")
    for (start, end, mi, me, r) in results_sorted:
        out.write(f"{start}\t{end}\t{mi:.6f}\t{me:.6f}\t{r:.4f}\n")

print(f"Resultados por ventana guardados en: {windows_scores_out}")

# 7) Elegir la mejor ventana (mayor ratio inter/intra)
best_start, best_end, best_mi, best_me, best_ratio = results_sorted[0]
print("MEJOR VENTANA:")
print(f"  columnas {best_start}-{best_end}")
print(f"  mean_intra = {best_mi:.4f}")
print(f"  mean_inter = {best_me:.4f}")
print(f"  ratio = {best_ratio:.2f}")

# 8) Exportar esa región como alineamiento recortado (con gaps)
with open(best_barcode_aln_out, "w") as out:
    for sid, s in zip(ids, seqs):
        sub = s[best_start:best_end]
        out.write(f">{sid}\n{sub}\n")

print(f"Alineamiento recortado de la mejor ventana escrito en: {best_barcode_aln_out}")

# 9) Exportar versión sin gaps (por si querés usarla sin alineamiento)
with open(best_barcode_ungapped_out, "w") as out:
    for sid, s in zip(ids, seqs):
        sub = s[best_start:best_end]
        sub_nogap = sub.replace("-", "")
        out.write(f">{sid}\n{sub_nogap}\n")

print(f"Secuencias recortadas sin gaps escritas en: {best_barcode_ungapped_out}")

