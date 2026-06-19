#!/usr/bin/env python3
"""
Lung CT GMM Batch Analysis v3
  - Multiple study dates per folder supported
  - Verification images output
  - Auto-selects: thin slice (1mm) + lung kernel + non-contrast per date
  - Lung segmentation: Lungmask AI (U-net R231)
  - Tissue classification: Zaffino et al. 2021 GMM (5 components)
"""

import os, sys, json, csv, time, argparse
import numpy as np
from datetime import datetime

INPUT_DIR = "/tmp/ILD_Data"
OUTPUT_DIR = "/tmp/CT_results/gmm_v3"

GMM_N_COMPONENTS = 5
GMM_COVARIANCE_TYPE = "full"
GMM_N_INIT = 6
GMM_INIT_PARAMS = "kmeans"
GMM_RANDOM_STATE = 42
GMM_MAX_ITER = 200
GMM_TOL = 1e-3
MEDIAN_FILTER_SIZE = 4
TISSUE_NAMES = ["Air", "Healthy_Lung", "GGO", "Consolidation", "Dense_Tissue"]
LUNG_KERNELS = ["FC51", "FC52", "I70f", "LUNG"]
N_VERIFY_SLICES = 5


def find_dicom_patients(input_dir):
    patients = []
    for item in sorted(os.listdir(input_dir)):
        item_path = os.path.join(input_dir, item)
        if not os.path.isdir(item_path): continue
        dcmdt_path = os.path.join(item_path, "DCMDT")
        if os.path.isdir(dcmdt_path):
            dcm_files = [f for f in os.listdir(dcmdt_path)
                         if os.path.isfile(os.path.join(dcmdt_path, f))
                         and not f.startswith('.')]
            if len(dcm_files) > 0:
                patients.append({"name": item, "path": dcmdt_path, "n_slices": len(dcm_files)})
    return patients


def is_lung_kernel(kernel_str):
    kernel_str = str(kernel_str).upper()
    for lk in LUNG_KERNELS:
        if lk.upper() in kernel_str: return True
    return False


def is_thin_slice(thickness_str):
    try: return float(str(thickness_str).strip()) <= 1.5
    except: return False


def is_contrast(desc_str):
    desc = str(desc_str).upper()
    return ",CE," in f",{desc}," or desc.startswith("CE,")


def find_study_dates(dicom_dir):
    import pydicom
    dates = set()
    for f in os.listdir(dicom_dir):
        if f.startswith('.'): continue
        try:
            ds = pydicom.dcmread(os.path.join(dicom_dir, f), stop_before_pixels=True,
                                 specific_tags=[(0x0008, 0x0020)])
            d = str(ds.get((0x0008, 0x0020), {}).value) if ds.get((0x0008, 0x0020)) else None
            if d: dates.add(d)
        except: pass
    return sorted(dates)


def select_best_series_for_date(dicom_dir, target_date):
    import pydicom, SimpleITK as sitk
    reader = sitk.ImageSeriesReader()
    series_ids = reader.GetGDCMSeriesIDs(dicom_dir)
    if not series_ids: return None, None, None
    candidates = []
    for sid in series_ids:
        files = reader.GetGDCMSeriesFileNames(dicom_dir, sid)
        if len(files) < 50: continue
        try:
            ds = pydicom.dcmread(files[0], stop_before_pixels=True,
                                 specific_tags=[(0x0008,0x0020),(0x0018,0x0050),
                                                (0x0018,0x1210),(0x0008,0x103E)])
            study_date = str(ds.get((0x0008,0x0020), {}).value) if ds.get((0x0008,0x0020)) else ""
            if study_date != target_date: continue
            thick = str(ds.get((0x0018,0x0050), {}).value) if ds.get((0x0018,0x0050)) else ""
            kernel = str(ds.get((0x0018,0x1210), {}).value) if ds.get((0x0018,0x1210)) else ""
            desc = str(ds.get((0x0008,0x103E), {}).value) if ds.get((0x0008,0x103E)) else ""
            candidates.append({"sid": sid, "files": files, "n": len(files),
                               "thick": thick, "kernel": kernel, "desc": desc,
                               "is_thin": is_thin_slice(thick),
                               "is_lung": is_lung_kernel(kernel),
                               "is_ce": is_contrast(desc)})
        except: pass
    if not candidates: return None, None, None
    def score(c):
        s = 0
        if c["is_thin"]: s += 1000
        if c["is_lung"]: s += 500
        if not c["is_ce"]: s += 100
        s += c["n"]
        return s
    candidates.sort(key=score, reverse=True)
    best = candidates[0]
    print(f"  Selected: {best['n']} slices | Thick={best['thick']} | "
          f"Kernel={best['kernel']} | CE={'Yes' if best['is_ce'] else 'No'} | {best['desc'][:50]}")
    reader2 = sitk.ImageSeriesReader()
    reader2.SetFileNames(best["files"])
    reader2.MetaDataDictionaryArrayUpdateOn()
    reader2.LoadPrivateTagsOn()
    image = reader2.Execute()
    volume_array = sitk.GetArrayFromImage(image)
    print(f"  Volume: {volume_array.shape}, spacing={image.GetSpacing()}")
    series_info = {"series_uid": best["sid"], "n_slices": best["n"],
                   "slice_thickness": best["thick"], "kernel": best["kernel"],
                   "description": best["desc"], "contrast": best["is_ce"],
                   "study_date": target_date, "n_candidate_series": len(candidates)}
    return image, volume_array, series_info


def segment_lungs_ai(image):
    from lungmask import LMInferer
    print("  Running Lungmask U-net(R231)...")
    inferer = LMInferer()
    segmentation = inferer.apply(image)
    lung_mask = segmentation > 0
    n_total = int(np.sum(lung_mask))
    n_right = int(np.sum(segmentation == 1))
    n_left = int(np.sum(segmentation == 2))
    print(f"  Lung voxels: {n_total:,} (R: {n_right:,}, L: {n_left:,})")
    return lung_mask, segmentation, {"total": n_total, "right": n_right, "left": n_left}


def run_gmm(volume_array, lung_mask):
    from sklearn.mixture import GaussianMixture
    from scipy.ndimage import median_filter
    lung_voxels = np.sum(lung_mask)
    if lung_voxels == 0: return None, None, None
    lung_intensities = volume_array[lung_mask].astype(np.float64).reshape(-1, 1)
    print(f"  GMM fitting ({GMM_N_COMPONENTS} components, {lung_voxels:,} voxels)...")
    t0 = time.time()
    gmm = GaussianMixture(n_components=GMM_N_COMPONENTS, covariance_type=GMM_COVARIANCE_TYPE,
                           n_init=GMM_N_INIT, init_params=GMM_INIT_PARAMS,
                           random_state=GMM_RANDOM_STATE, max_iter=GMM_MAX_ITER, tol=GMM_TOL)
    gmm.fit(lung_intensities)
    fit_time = time.time() - t0
    print(f"  Converged: {gmm.converged_} ({fit_time:.1f}s)")
    raw_labels = gmm.predict(lung_intensities)
    means = gmm.means_.flatten()
    weights = gmm.weights_.flatten()
    covariances = gmm.covariances_.flatten()
    sorted_indices = np.argsort(means)
    label_mapping = {old: new for new, old in enumerate(sorted_indices)}
    mapped_labels = np.array([label_mapping[l] for l in raw_labels])
    label_map = np.zeros(volume_array.shape, dtype=np.int8)
    label_map[lung_mask] = mapped_labels + 1
    sorted_means = means[sorted_indices]
    sorted_weights = weights[sorted_indices]
    sorted_covariances = covariances[sorted_indices]
    for i, (name, mean) in enumerate(zip(TISSUE_NAMES, sorted_means)):
        print(f"    {i+1}: {name} = {mean:.1f} HU (weight={sorted_weights[i]:.3f})")
    mask_nz = label_map > 0
    filtered = median_filter(label_map, size=MEDIAN_FILTER_SIZE)
    label_map[mask_nz] = filtered[mask_nz]
    gmm_info = {"means_HU": sorted_means.tolist(), "weights": sorted_weights.tolist(),
                "covariances": sorted_covariances.tolist(), "converged": bool(gmm.converged_),
                "n_iter": int(gmm.n_iter_), "fit_time_s": round(fit_time, 1)}
    return label_map, TISSUE_NAMES, gmm_info


def compute_statistics(volume_array, label_map, lung_mask, tissue_names, spacing, lung_voxel_info):
    voxel_vol_cm3 = (spacing[0] * spacing[1] * spacing[2]) / 1000.0
    result = {}; total_vol = 0
    for i, name in enumerate(tissue_names):
        n_voxels = int(np.sum(label_map == (i + 1)))
        vol_cm3 = n_voxels * voxel_vol_cm3
        result[f"{name}_voxels"] = n_voxels
        result[f"{name}_cm3"] = round(vol_cm3, 2)
        total_vol += vol_cm3
    result["total_lung_cm3"] = round(total_vol, 2)
    for name in tissue_names:
        vol = result[f"{name}_cm3"]
        result[f"{name}_pct"] = round(vol / total_vol * 100, 2) if total_vol > 0 else 0
    result["lung_involvement_pct"] = round(
        result.get("GGO_pct",0) + result.get("Consolidation_pct",0) + result.get("Dense_Tissue_pct",0), 2)
    lung_hu = volume_array[lung_mask]
    result["lung_mean_HU"] = round(float(np.mean(lung_hu)), 1)
    result["lung_median_HU"] = round(float(np.median(lung_hu)), 1)
    result["lung_std_HU"] = round(float(np.std(lung_hu)), 1)
    result["right_lung_voxels"] = lung_voxel_info["right"]
    result["left_lung_voxels"] = lung_voxel_info["left"]
    result["spacing_x_mm"] = round(spacing[0], 4)
    result["spacing_y_mm"] = round(spacing[1], 4)
    result["spacing_z_mm"] = round(spacing[2], 4)
    result["voxel_vol_mm3"] = round(spacing[0] * spacing[1] * spacing[2], 4)
    return result


def save_verification_images(volume_array, label_map, lung_mask, output_name, output_dir):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap

    pt_dir = os.path.join(output_dir, output_name)
    os.makedirs(pt_dir, exist_ok=True)

    # Identify the slice range containing lung
    lung_slices = np.where(np.any(lung_mask, axis=(1, 2)))[0]
    if len(lung_slices) == 0: return

    z_min, z_max = lung_slices[0], lung_slices[-1]
    z_range = z_max - z_min

    # Select slices at 10%, 30%, 50%, 70%, 90% positions of the lung extent
    percentiles = [0.10, 0.30, 0.50, 0.70, 0.90]
    slice_indices = [int(z_min + z_range * p) for p in percentiles]
    slice_labels = ["upper", "upper_mid", "middle", "lower_mid", "lower"]

    # GMM label colormap: 0=transparent, 1=Air, 2=Healthy, 3=GGO, 4=Consolidation, 5=Dense
    colors_rgba = [
        (0, 0, 0, 0),        # 0: background (transparent)
        (0, 0, 0, 0.4),      # 1: Air (dark, semi-transparent)
        (0, 0.5, 1, 0.4),    # 2: Healthy (blue)
        (1, 1, 0, 0.6),      # 3: GGO (yellow)
        (1, 0.5, 0, 0.6),    # 4: Consolidation (orange)
        (1, 0, 0, 0.6),      # 5: Dense (red)
    ]
    gmm_cmap = ListedColormap(colors_rgba)

    for idx, (z, label) in enumerate(zip(slice_indices, slice_labels)):
        if z < 0 or z >= volume_array.shape[0]: continue

        ct_slice = volume_array[z, :, :]
        lbl_slice = label_map[z, :, :]

        fig, axes = plt.subplots(1, 3, figsize=(18, 6))

        # Left: CT only
        axes[0].imshow(ct_slice, cmap='gray', vmin=-1000, vmax=400)
        axes[0].set_title(f'CT (slice {z})', fontsize=12)
        axes[0].axis('off')

        # Middle: CT + GMM overlay
        axes[1].imshow(ct_slice, cmap='gray', vmin=-1000, vmax=400)
        masked_lbl = np.ma.masked_where(lbl_slice == 0, lbl_slice)
        axes[1].imshow(masked_lbl, cmap=gmm_cmap, vmin=0, vmax=5, interpolation='nearest')
        axes[1].set_title(f'CT + GMM overlay ({label})', fontsize=12)
        axes[1].axis('off')

        # Right: GMM labels only
        label_colors_solid = [
            (0, 0, 0),        # 0: background
            (0.2, 0.2, 0.2),  # 1: Air
            (0, 0.5, 1),      # 2: Healthy
            (1, 1, 0),        # 3: GGO
            (1, 0.5, 0),      # 4: Consolidation
            (1, 0, 0),        # 5: Dense
        ]
        gmm_cmap_solid = ListedColormap(label_colors_solid)
        axes[2].imshow(lbl_slice, cmap=gmm_cmap_solid, vmin=0, vmax=5, interpolation='nearest')
        axes[2].set_title(f'GMM labels ({label})', fontsize=12)
        axes[2].axis('off')

        # Legend
        legend_text = "  ".join([f"{i}: {n}" for i, n in enumerate(TISSUE_NAMES, 1)])
        fig.text(0.5, 0.02, legend_text, ha='center', fontsize=10,
                 bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

        fig.suptitle(f'{output_name}  —  {label} (z={z}/{volume_array.shape[0]})', fontsize=14)
        plt.tight_layout(rect=[0, 0.05, 1, 0.95])

        png_path = os.path.join(pt_dir, f"verify_{idx+1}_{label}.png")
        plt.savefig(png_path, dpi=150, bbox_inches='tight')
        plt.close(fig)

    print(f"  Verification images: {len(slice_indices)} PNGs saved")


def save_patient_results(result, output_name, output_dir, gmm_info=None,
                         label_map=None, image=None, series_info=None):
    import SimpleITK as sitk
    pt_dir = os.path.join(output_dir, output_name); os.makedirs(pt_dir, exist_ok=True)
    json_data = {"method": "Lungmask_AI + Zaffino_GMM",
                 "lung_segmentation": "Lungmask U-net(R231)",
                 "tissue_classification": "GMM 5-component (Zaffino et al. 2021)",
                 "series_selection": series_info,
                 "gmm_parameters": {"n_components": GMM_N_COMPONENTS, "covariance_type": GMM_COVARIANCE_TYPE,
                                    "n_init": GMM_N_INIT, "init_params": GMM_INIT_PARAMS,
                                    "random_state": GMM_RANDOM_STATE, "max_iter": GMM_MAX_ITER,
                                    "tol": GMM_TOL, "median_filter_size": MEDIAN_FILTER_SIZE},
                 "gmm_fit_results": gmm_info,
                 "output_name": output_name,
                 "analysis_date": datetime.now().isoformat(),
                 "results": {k: v for k, v in result.items()
                             if k not in ("patient","folder_name","study_date","output_name")}}
    with open(os.path.join(pt_dir, "gmm_results.json"), 'w', encoding='utf-8') as f:
        json.dump(json_data, f, indent=2, ensure_ascii=False)
    with open(os.path.join(pt_dir, "segment_statistics.csv"), 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)
        writer.writerow(["Segment", "Voxels", "Volume_cm3", "Percentage"])
        for name in TISSUE_NAMES:
            writer.writerow([name, result.get(f"{name}_voxels",0),
                             result.get(f"{name}_cm3",0), result.get(f"{name}_pct",0)])
    if label_map is not None and image is not None:
        label_image = sitk.GetImageFromArray(label_map.astype(np.int16))
        label_image.CopyInformation(image)
        sitk.WriteImage(label_image, os.path.join(pt_dir, "gmm_labelmap.nii.gz"))
    print(f"  Saved to: {pt_dir}")


def collect_dicom_metadata(patients, output_dir):
    import pydicom
    print("\nCollecting DICOM metadata...")
    tags = [("PatientID",0x0010,0x0020),("PatientName",0x0010,0x0010),
            ("PatientBirthDate",0x0010,0x0030),("PatientSex",0x0010,0x0040),
            ("PatientAge",0x0010,0x1010),("StudyDate",0x0008,0x0020),
            ("StudyDescription",0x0008,0x1030),("Modality",0x0008,0x0060),
            ("Manufacturer",0x0008,0x0070),("ManufacturerModelName",0x0008,0x1090),
            ("InstitutionName",0x0008,0x0080)]
    rows = []
    for p in patients:
        dcm_files = [f for f in os.listdir(p["path"]) if not f.startswith('.')]
        if not dcm_files: continue
        dates = find_study_dates(p["path"])
        for d in dates:
            for f in dcm_files:
                try:
                    ds = pydicom.dcmread(os.path.join(p["path"], f), stop_before_pixels=True)
                    charset = ds.get((0x0008,0x0005), None)
                    if charset: ds.SpecificCharacterSet = charset.value
                    sd = str(ds.get((0x0008,0x0020), {}).value) if ds.get((0x0008,0x0020)) else ""
                    if sd != d: continue
                    row = {"folder_name": p["name"], "study_date": d,
                           "output_name": f"{p['name']}_{d}"}
                    for tag_name, g, e in tags:
                        elem = ds.get((g,e), None)
                        if elem is not None:
                            v = elem.value
                            if tag_name == "PatientName":
                                pn = str(v); parts = pn.split('=')
                                row[tag_name] = parts[-1] if len(parts) > 1 else parts[0]
                            else: row[tag_name] = str(v)
                        else: row[tag_name] = ""
                    rows.append(row)
                    print(f"  {row['output_name']}: {row.get('PatientName','')} / {d}")
                    break
                except: pass
    fieldnames = ["folder_name", "study_date", "output_name"] + [t[0] for t in tags]
    metadata_csv = os.path.join(output_dir, "dicom_metadata.csv")
    with open(metadata_csv, 'w', newline='', encoding='utf-8-sig') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction='ignore')
        writer.writeheader(); writer.writerows(rows)
    print(f"  dicom_metadata.csv: {metadata_csv}")


def process_one_study(folder_name, dicom_dir, study_date, output_dir):
    output_name = f"{folder_name}_{study_date}"
    print(f"\n{'='*60}")
    print(f"  {output_name} (folder={folder_name}, date={study_date})")
    print(f"{'='*60}")
    t0 = time.time()

    print("\n[1/5] Loading DICOM (smart series selection)...")
    image, volume_array, series_info = select_best_series_for_date(dicom_dir, study_date)
    if volume_array is None:
        print("  SKIP: No suitable series found"); return None

    print("\n[2/5] AI lung segmentation (Lungmask)...")
    lung_mask, _, lung_voxel_info = segment_lungs_ai(image)
    if np.sum(lung_mask) == 0:
        print("  SKIP: No lung detected"); return None

    print("\n[3/5] GMM tissue classification...")
    label_map, tissue_names, gmm_info = run_gmm(volume_array, lung_mask)
    if label_map is None: return None

    print("\n[4/5] Statistics & save...")
    spacing = image.GetSpacing()
    result = compute_statistics(volume_array, label_map, lung_mask, tissue_names, spacing, lung_voxel_info)
    result["folder_name"] = folder_name
    result["study_date"] = study_date
    result["output_name"] = output_name
    result["series_kernel"] = series_info["kernel"]
    result["series_thickness"] = series_info["slice_thickness"]
    result["series_contrast"] = "CE" if series_info["contrast"] else "plain"
    result["series_description"] = series_info["description"][:50]
    save_patient_results(result, output_name, output_dir, gmm_info, label_map, image, series_info)

    print("\n[5/5] Verification images...")
    save_verification_images(volume_array, label_map, lung_mask, output_name, output_dir)

    elapsed = time.time() - t0
    result["processing_time_s"] = round(elapsed, 1)

    print(f"\n  --- Results: {output_name} ---")
    for tname in tissue_names:
        vol = result[f"{tname}_cm3"]; pct = result[f"{tname}_pct"]
        print(f"  {tname:20s}: {vol:8.1f} cm3 ({pct:5.1f}%) {'#' * int(pct/2)}")
    print(f"  {'Total':20s}: {result['total_lung_cm3']:8.1f} cm3")
    print(f"  Lung involvement   = {result['lung_involvement_pct']:.1f}%")
    print(f"  Time: {elapsed:.1f}s")
    return result


def main():
    parser = argparse.ArgumentParser(description="Lung GMM v3 (Multi-date + Verification)")
    parser.add_argument("--test", action="store_true", help="Process only the first folder")
    parser.add_argument("--patient", type=str, help="Specific folder name")
    parser.add_argument("--input", type=str, default=INPUT_DIR)
    parser.add_argument("--output", type=str, default=OUTPUT_DIR)
    args = parser.parse_args()
    input_dir = args.input; output_dir = args.output; os.makedirs(output_dir, exist_ok=True)

    print("="*60)
    print("  Lung CT GMM Batch Analysis v3 (Multi-date + Verify)")
    print("  Series: 1mm + Lung kernel + non-CE preferred")
    print("  Segmentation: Lungmask AI (U-net R231)")
    print("  Classification: GMM 5-component (Zaffino 2021)")
    print("="*60)

    patients = find_dicom_patients(input_dir)
    print(f"\nFound {len(patients)} folders\n")

    if args.test: patients = patients[:1]
    elif args.patient:
        patients = [p for p in patients if p["name"] == args.patient]
        if not patients: print("Not found."); return

    tasks = []
    for p in patients:
        dates = find_study_dates(p["path"])
        for d in dates:
            tasks.append({"folder_name": p["name"], "path": p["path"], "date": d})
        print(f"  {p['name']}: {len(dates)} date(s) -> {', '.join(dates)}")

    print(f"\nTotal tasks: {len(tasks)} studies across {len(patients)} folders\n")

    all_results = []; skipped = []
    total_start = time.time()

    for i, task in enumerate(tasks):
        print(f"\n>>> [{i+1}/{len(tasks)}]")
        try:
            result = process_one_study(task["folder_name"], task["path"], task["date"], output_dir)
            if result: all_results.append(result)
            else: skipped.append(f"{task['folder_name']}_{task['date']}")
        except Exception as e:
            print(f"\n  ERROR: {e}"); import traceback; traceback.print_exc()
            skipped.append(f"{task['folder_name']}_{task['date']}")

    total_elapsed = time.time() - total_start

    print(f"\n{'='*60}\n  BATCH SUMMARY\n{'='*60}")
    print(f"  Studies   : {len(all_results)}/{len(tasks)} completed")
    print(f"  Total time: {total_elapsed/60:.1f} min")
    if skipped: print(f"  Skipped   : {len(skipped)} ({', '.join(skipped)})")

    if all_results:
        summary_path = os.path.join(output_dir, "batch_summary.csv")
        columns = ["output_name", "folder_name", "study_date",
                    "Air_cm3","Healthy_Lung_cm3","GGO_cm3","Consolidation_cm3","Dense_Tissue_cm3","total_lung_cm3",
                    "Air_pct","Healthy_Lung_pct","GGO_pct","Consolidation_pct","Dense_Tissue_pct","lung_involvement_pct",
                    "lung_mean_HU","lung_median_HU","lung_std_HU",
                    "right_lung_voxels","left_lung_voxels",
                    "spacing_x_mm","spacing_y_mm","spacing_z_mm","voxel_vol_mm3",
                    "series_kernel","series_thickness","series_contrast","series_description",
                    "processing_time_s"]
        with open(summary_path, 'w', newline='', encoding='utf-8') as f:
            writer = csv.DictWriter(f, fieldnames=columns, extrasaction='ignore')
            writer.writeheader(); writer.writerows(all_results)
        print(f"\n  batch_summary.csv: {summary_path}")

        print(f"\n  {'Output':25s} {'Date':>10s} {'Involve%':>9s} {'Healthy%':>9s} "
              f"{'GGO%':>7s} {'Kernel':>10s} {'CE':>5s}")
        print(f"  {'-'*80}")
        for r in all_results:
            print(f"  {r['output_name']:25s} {r.get('study_date',''):>10s} "
                  f"{r.get('lung_involvement_pct',0):8.1f}% "
                  f"{r.get('Healthy_Lung_pct',0):8.1f}% "
                  f"{r.get('GGO_pct',0):6.1f}% "
                  f"{r.get('series_kernel',''):>10s} "
                  f"{r.get('series_contrast',''):>5s}")

    collect_dicom_metadata(patients, output_dir)

    config = {"analysis_date": datetime.now().isoformat(),
              "n_folders": len(patients), "n_studies": len(tasks),
              "n_completed": len(all_results), "n_skipped": len(skipped),
              "total_time_min": round(total_elapsed/60, 1),
              "series_selection": {"priority": "1mm + lung kernel + non-CE", "min_slices": 50},
              "lung_segmentation": "Lungmask U-net(R231)",
              "gmm_parameters": {"n_components": GMM_N_COMPONENTS, "random_state": GMM_RANDOM_STATE},
              "verification_images": {"n_slices": N_VERIFY_SLICES, "positions": "10%, 30%, 50%, 70%, 90% of lung"},
              "software_versions": {}}
    try: import torch; config["software_versions"]["torch"] = torch.__version__
    except: pass
    try: import sklearn; config["software_versions"]["scikit-learn"] = sklearn.__version__
    except: pass
    config["software_versions"]["python"] = sys.version
    with open(os.path.join(output_dir, "analysis_config.json"), 'w', encoding='utf-8') as f:
        json.dump(config, f, indent=2, ensure_ascii=False)
    print(f"\n{'='*60}\n  Done.\n{'='*60}")

if __name__ == "__main__": main()
