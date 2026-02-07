import os
import ast
import numpy as np
import pandas as pd
from neuromaps.datasets import fetch_annotation
from neuromaps import transforms
from enigmatoolbox.utils.parcellation import surface_to_parcel
from enigmatoolbox.datasets import load_summary_stats

# ==============================
# Config (edit these if needed)
# ==============================
PET_INFO = r"E:\xhmhc\Cushing_project\code\neuromaps_PET_info.csv"
OUT_DIR = r"F:\Google Drive\post-doc\Structural_subtype_new\Cue_res\PET_analysis"
APPLY_MINMAX = True  # matches PET_analysis.py (min-max after parcellation)
SAVE_RAW = True


def min_max_scale(x):
    x = np.asarray(x, dtype=float)
    x = x - np.min(x)
    mx = np.max(x)
    return x / mx if mx != 0 else x


def get_region_names():
    try:
        sum_stats = load_summary_stats('22q')
        region_names = sum_stats['CortThick_case_vs_controls']['Structure']
        region_names = [str(r) for r in region_names]
        if len(region_names) == 68:
            return region_names
    except Exception as e:
        print(f"[WARN] Could not load region names from ENIGMA: {e}")
    return [f"Region_{i:02d}" for i in range(1, 69)]


if not os.path.exists(OUT_DIR):
    os.makedirs(OUT_DIR)

# Load PET metadata
pet_df = pd.read_csv(PET_INFO)
region_names = get_region_names()

rows_scaled = []
rows_raw = []
meta_rows = []

for _, row in pet_df.iterrows():
    source, desc, space, res = ast.literal_eval(row['annotation'])
    transmitter = row['transmitter']
    total_n = row['total_n']

    print(f"Fetching: {source}, {desc}, {space}, {res}")

    # Fetch annotation
    annot = fetch_annotation(source=source, desc=desc, space=space)

    # Convert to fsaverage 10k (fsa5)
    if space == 'MNI152':
        fsa5 = transforms.mni152_to_fsaverage(annot, '10k')
    elif space == 'fsaverage':
        if res != '10k':
            fsa5 = transforms.fsaverage_to_fsaverage(annot, '10k')
        else:
            fsa5 = annot
    else:
        raise ValueError(f"Unsupported space: {space}")

    # Vertex data (lh/rh)
    fsa5_l, fsa5_r = fsa5
    v_l = fsa5_l.agg_data()
    v_r = fsa5_r.agg_data()

    # Optional: clip negatives at vertex level (kept to match PET_analysis.py)
    v_l = np.where(v_l < 0, 0, v_l)
    v_r = np.where(v_r < 0, 0, v_r)

    # Concatenate (lh + rh)
    fsa5_array = np.append(v_l, v_r)

    # Parcellate to DK (aparc_fsa5)
    dk_raw = surface_to_parcel(fsa5_array, 'aparc_fsa5')

    # remove 0: unknown, 4: left corpus callosum, 39: right corpus callosum
    dk_raw = np.delete(dk_raw, [0, 4, 39], axis=0)

    if APPLY_MINMAX:
        dk_scaled = min_max_scale(dk_raw)
    else:
        dk_scaled = dk_raw

    meta_rows.append({
        'Transmitter': transmitter,
        'Source': source,
        'Desc': desc,
        'Space': space,
        'Res': res,
        'total_n': total_n
    })

    rows_scaled.append(dk_scaled)
    rows_raw.append(dk_raw)

# Build long tables
meta_df = pd.DataFrame(meta_rows)
scaled_df = pd.DataFrame(rows_scaled, columns=region_names)
long_scaled = pd.concat([meta_df, scaled_df], axis=1)

out_scaled = os.path.join(OUT_DIR, "PET_DK68_parcellated_long.csv")
long_scaled.to_csv(out_scaled, index=False)
print(f"Saved: {out_scaled}")

if SAVE_RAW:
    raw_df = pd.DataFrame(rows_raw, columns=region_names)
    long_raw = pd.concat([meta_df, raw_df], axis=1)
    out_raw = os.path.join(OUT_DIR, "PET_DK68_parcellated_long_raw.csv")
    long_raw.to_csv(out_raw, index=False)
    print(f"Saved: {out_raw}")

# Weighted mean per transmitter (scaled)
weighted = long_scaled.groupby('Transmitter', sort=False).apply(
    lambda x: np.average(x[region_names].to_numpy(), weights=x['total_n'].to_numpy(), axis=0),
    include_groups=False
)

pet_wide = pd.DataFrame(np.column_stack(weighted.to_list()), columns=weighted.index)
pet_wide.insert(0, 'RegionName', region_names)

out_wide = os.path.join(OUT_DIR, "PET_DK68_weighted_mean.csv")
pet_wide.to_csv(out_wide, index=False)
print(f"Saved: {out_wide}")
