"""P5 Vina 전체 검증 스크립트 (v3 - 파싱 버그 수정)
200개 생성 분자 전체 + 기준 약물 7종을 AutoDock Vina로 도킹.
GNINA(CNN scoring) vs Vina(physics-based scoring) 종합 비교.

실행:
  python scripts/vina_validation.py \
      --candidates p5_200_candidates.csv \
      --top_n 200 \
      --exhaustiveness 16
"""
import argparse
import os
import subprocess
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem

# ================================================================
# PDB 구조
# ================================================================
PDB_STRUCTURES = {
    "2ITY_WT": {
        "protein": "2ITY_protein.pdb",
        "autobox_ligand": "2ITY_IRE.pdb",
        "target_type": "WT",
    },
    "4WKQ_WT": {
        "protein": "4WKQ_protein.pdb",
        "autobox_ligand": "4WKQ_IRE.pdb",
        "target_type": "WT",
    },
    "4HJO_inactive": {
        "protein": "4HJO_protein.pdb",
        "autobox_ligand": "4HJO_AQ4.pdb",
        "target_type": "WT",
    },
    "4I23_WT": {
        "protein": "4I23_protein.pdb",
        "autobox_ligand": "4I23_1C9.pdb",
        "target_type": "WT",
    },
    "4G5J_T790M": {
        "protein": "4G5J_protein.pdb",
        "autobox_ligand": "4G5J_0WM.pdb",
        "target_type": "T790M",
    },
    "4ZAU_T790M": {
        "protein": "4ZAU_protein.pdb",
        "autobox_ligand": "4ZAU_YY3.pdb",
        "target_type": "T790M",
    },
    "5D41_allosteric": {
        "protein": "5D41_protein.pdb",
        "autobox_ligand": "5D41_57N.pdb",
        "target_type": "T790M_allosteric",
    },
}

# 기준 약물
REFERENCE_DRUGS = {
    "Gefitinib_1G": {
        "smiles": "COc1cc2ncnc(Nc3ccc(F)c(Cl)c3)c2cc1OCCCN1CCOCC1",
        "generation": "1G",
        "mechanism": "reversible",
    },
    "Erlotinib_1G": {
        "smiles": "COc1cc2ncnc(Nc3cccc(C#C)c3)c2cc1OCCOC",
        "generation": "1G",
        "mechanism": "reversible",
    },
    "Afatinib_2G": {
        "smiles": "CN(C)C/C=C/C(=O)Nc1cc2c(Nc3ccc(F)c(Cl)c3)ncnc2cc1OC1CCOC1",
        "generation": "2G",
        "mechanism": "covalent (pan-ErbB)",
    },
    "Dacomitinib_2G": {
        "smiles": "COc1cc2ncnc(Nc3ccc(F)c(Cl)c3)c2cc1NC(=O)/C=C/CN1CCCCC1",
        "generation": "2G",
        "mechanism": "covalent (pan-ErbB)",
    },
    "Osimertinib_3G": {
        "smiles": "COc1cc(N(C)CCN(C)C)c(NC(=O)/C=C/CN(C)C)cc1Nc1nccc(-c2cn(C)c3ccccc23)n1",
        "generation": "3G",
        "mechanism": "covalent (T790M selective)",
    },
    "EAI045_4G": {
        "smiles": "COc1cc(C2CC(=O)N(c3cccc(C(F)(F)F)c3)C2=O)cc(OC)c1O",
        "generation": "4G",
        "mechanism": "allosteric (C797S bypass)",
    },
    "BDTX1535_4G": {
        "smiles": "CN1CCN(c2ccc(-c3cc4cnc(Nc5ccc(N6CCN(C)CC6)c(C)c5)nc4s3)cn2)CC1",
        "generation": "4G",
        "mechanism": "non-covalent (C797S bypass)",
    },
}


def parse_args():
    parser = argparse.ArgumentParser(
        description="P5 Vina: 200개 전체 + 기준 약물 + GNINA cross-validation")
    parser.add_argument("--candidates", type=str, required=True)
    parser.add_argument("--top_n", type=int, default=200)
    parser.add_argument("--docking_dir", type=str,
                        default="/home/nudge/Project/elix/docking")
    parser.add_argument("--output_dir", type=str,
                        default="/home/nudge/Project/elix/docking/results/vina_validation")
    parser.add_argument("--exhaustiveness", type=int, default=16)
    parser.add_argument("--refine_top_n", type=int, default=30)
    parser.add_argument("--refine_exhaustiveness", type=int, default=32)
    parser.add_argument("--refine_seeds", type=int, nargs="+",
                        default=[0, 42, 123])
    parser.add_argument("--skip_refine", action="store_true")
    parser.add_argument("--skip_reference", action="store_true")
    parser.add_argument("--pdb_targets", type=str, nargs="+",
                        default=list(PDB_STRUCTURES.keys()))
    parser.add_argument("--gnina_results", type=str,
                        default="results/generated_molecules/docking_results.csv")
    parser.add_argument("--gnina_refined", type=str,
                        default="results/generated_molecules/docking_results_refined.csv")
    parser.add_argument("--gnina_reference", type=str,
                        default="results/generated_molecules/reference_drugs_docking.csv")
    return parser.parse_args()


def smiles_to_sdf(smiles, output_path):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = 42
    result = AllChem.EmbedMolecule(mol, params)
    if result == -1:
        params.useRandomCoords = True
        result = AllChem.EmbedMolecule(mol, params)
        if result == -1:
            return False
    try:
        AllChem.MMFFOptimizeMolecule(mol, maxIters=500)
    except:
        pass
    writer = Chem.SDWriter(output_path)
    writer.write(mol)
    writer.close()
    return True


def sdf_to_pdbqt(sdf_path, pdbqt_path):
    cmd = ["obabel", sdf_path, "-O", pdbqt_path]
    try:
        subprocess.run(cmd, capture_output=True, text=True, timeout=60)
        return os.path.exists(pdbqt_path)
    except:
        return False


def pdb_to_pdbqt(pdb_path, pdbqt_path):
    cmd = ["obabel", pdb_path, "-O", pdbqt_path, "-xr"]
    try:
        subprocess.run(cmd, capture_output=True, text=True, timeout=60)
        return os.path.exists(pdbqt_path)
    except:
        return False


def get_box_from_ligand(ligand_pdb_path, padding=4.0):
    coords = []
    with open(ligand_pdb_path, "r") as f:
        for line in f:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                coords.append([x, y, z])
    if not coords:
        return None
    coords = np.array(coords)
    center = coords.mean(axis=0)
    size = (coords.max(axis=0) - coords.min(axis=0)) + 2 * padding
    return {
        "center_x": center[0], "center_y": center[1], "center_z": center[2],
        "size_x": size[0], "size_y": size[1], "size_z": size[2],
    }


def run_vina(receptor_pdbqt, ligand_pdbqt, box, output_pdbqt,
             exhaustiveness=16, seed=0):
    """AutoDock Vina 실행.
    Vina 출력 포맷:
       1       -8.355          0          0
       2       -8.283      3.464      6.894
    """
    cmd = [
        "vina",
        "--receptor", receptor_pdbqt,
        "--ligand", ligand_pdbqt,
        "--center_x", str(box["center_x"]),
        "--center_y", str(box["center_y"]),
        "--center_z", str(box["center_z"]),
        "--size_x", str(box["size_x"]),
        "--size_y", str(box["size_y"]),
        "--size_z", str(box["size_z"]),
        "--exhaustiveness", str(exhaustiveness),
        "--seed", str(seed),
        "--out", output_pdbqt,
    ]
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
        output = result.stdout + result.stderr
        for line in output.split("\n"):
            parts = line.split()
            if len(parts) >= 4 and parts[0] == "1":
                try:
                    aff = float(parts[1])
                    return {"vina_affinity": aff}
                except ValueError:
                    continue
    except Exception as e:
        print(f"    Vina ERROR: {e}")
    return None


def dock_molecule_all_targets_vina(name, pdbqt_path, output_dir,
                                    receptor_pdbqts, boxes, pdb_targets,
                                    exhaustiveness, seeds=None):
    if seeds is None:
        seeds = [0]

    results = {}
    for pdb_key in pdb_targets:
        affs = []
        for s in seeds:
            seed_tag = f"_s{s}" if len(seeds) > 1 else ""
            out_path = os.path.join(output_dir, f"{name}_{pdb_key}{seed_tag}_vina.pdbqt")
            res = run_vina(receptor_pdbqts[pdb_key], pdbqt_path,
                           boxes[pdb_key], out_path, exhaustiveness, s)
            if res:
                affs.append(res["vina_affinity"])

        if affs:
            results[pdb_key] = {
                "vina_aff_mean": np.mean(affs),
                "vina_aff_std": np.std(affs) if len(seeds) > 1 else 0.0,
            }
        else:
            results[pdb_key] = None

    return results


def flatten_vina_results(results, prefix=""):
    flat = {}
    for pdb_key, res in results.items():
        p = f"{prefix}{pdb_key}" if prefix else pdb_key
        if res:
            flat[f"{p}_vina_aff"] = res["vina_aff_mean"]
            flat[f"{p}_vina_aff_std"] = res["vina_aff_std"]
        else:
            flat[f"{p}_vina_aff"] = np.nan
            flat[f"{p}_vina_aff_std"] = np.nan
    return flat


def compute_vina_selectivity(flat, pdb_targets):
    metrics = {}
    t790m_keys = [k for k in pdb_targets
                  if PDB_STRUCTURES[k]["target_type"] == "T790M"]
    wt_keys = [k for k in pdb_targets
               if PDB_STRUCTURES[k]["target_type"] == "WT"]

    t790m_vals = [flat[f"{k}_vina_aff"] for k in t790m_keys
                  if not np.isnan(flat.get(f"{k}_vina_aff", np.nan))]
    wt_vals = [flat[f"{k}_vina_aff"] for k in wt_keys
               if not np.isnan(flat.get(f"{k}_vina_aff", np.nan))]

    if t790m_vals:
        metrics["T790M_vina_avg"] = np.mean(t790m_vals)
    if wt_vals:
        metrics["WT_vina_avg"] = np.mean(wt_vals)
    if t790m_vals and wt_vals:
        # Vina: 더 음수 = 더 좋음. selectivity = WT - T790M (양수면 T790M 선택적)
        metrics["vina_selectivity"] = np.mean(wt_vals) - np.mean(t790m_vals)

    zau = flat.get("4ZAU_T790M_vina_aff", np.nan)
    g5j = flat.get("4G5J_T790M_vina_aff", np.nan)
    if not np.isnan(zau) and not np.isnan(g5j):
        metrics["T790M_vina_cross_pdb_diff"] = abs(zau - g5j)
    ity = flat.get("2ITY_WT_vina_aff", np.nan)
    wkq = flat.get("4WKQ_WT_vina_aff", np.nan)
    if not np.isnan(ity) and not np.isnan(wkq):
        metrics["WT_vina_cross_pdb_diff"] = abs(ity - wkq)

    allo = flat.get("5D41_allosteric_vina_aff", np.nan)
    if not np.isnan(allo):
        metrics["allosteric_vina_aff"] = allo

    return metrics


def main():
    args = parse_args()
    docking_dir = args.docking_dir
    os.makedirs(args.output_dir, exist_ok=True)
    sdf_dir = os.path.join(args.output_dir, "sdf_files")
    pdbqt_dir = os.path.join(args.output_dir, "pdbqt_files")
    receptor_dir = os.path.join(args.output_dir, "receptor_pdbqt")
    os.makedirs(sdf_dir, exist_ok=True)
    os.makedirs(pdbqt_dir, exist_ok=True)
    os.makedirs(receptor_dir, exist_ok=True)

    # ================================================================
    # Step 0: PDB 준비
    # ================================================================
    print("=" * 70)
    print("Vina 준비: receptor PDBQT 변환 + 도킹 박스 계산")
    print("=" * 70)

    valid_targets = []
    receptor_pdbqts = {}
    boxes = {}

    for pdb_key in args.pdb_targets:
        info = PDB_STRUCTURES[pdb_key]
        protein_pdb = os.path.join(docking_dir, "pdb", info["protein"])
        autobox_pdb = os.path.join(docking_dir, "ligands", info["autobox_ligand"])

        if not os.path.exists(protein_pdb) or not os.path.exists(autobox_pdb):
            print(f"  {pdb_key:20s}: MISSING")
            continue

        rec_pdbqt = os.path.join(receptor_dir, f"{pdb_key}_receptor.pdbqt")
        if not os.path.exists(rec_pdbqt):
            if not pdb_to_pdbqt(protein_pdb, rec_pdbqt):
                print(f"  {pdb_key:20s}: receptor PDBQT 변환 실패")
                continue

        box = get_box_from_ligand(autobox_pdb)
        if box is None:
            print(f"  {pdb_key:20s}: 도킹 박스 계산 실패")
            continue

        receptor_pdbqts[pdb_key] = rec_pdbqt
        boxes[pdb_key] = box
        valid_targets.append(pdb_key)
        print(f"  {pdb_key:20s}: OK (center=[{box['center_x']:.1f}, {box['center_y']:.1f}, {box['center_z']:.1f}])")

    if not valid_targets:
        print("\n사용 가능한 PDB 없음 — 종료")
        return

    pdb_targets = valid_targets
    total_docks = (
        (7 * len(args.refine_seeds) * len(pdb_targets) if not args.skip_reference else 0) +
        args.top_n * len(pdb_targets) +
        (args.refine_top_n * len(args.refine_seeds) * len(pdb_targets) if not args.skip_refine else 0)
    )
    print(f"\n사용할 PDB: {len(pdb_targets)}개, 예상 도킹 횟수: ~{total_docks}회")

    # ================================================================
    # Step 1: 기준 약물 Vina 도킹
    # ================================================================
    ref_results = []
    if not args.skip_reference:
        print("\n" + "=" * 70)
        print("[Step 1] 기준 약물 Vina 도킹")
        print(f"  exhaustiveness={args.exhaustiveness}, seeds={args.refine_seeds}")
        print("=" * 70)

        for drug_name, drug_info in REFERENCE_DRUGS.items():
            print(f"\n  {drug_name} ({drug_info['generation']})")

            sdf_path = os.path.join(sdf_dir, f"{drug_name}.sdf")
            pdbqt_path = os.path.join(pdbqt_dir, f"{drug_name}.pdbqt")
            if not smiles_to_sdf(drug_info["smiles"], sdf_path):
                print(f"    SDF 변환 실패"); continue
            if not sdf_to_pdbqt(sdf_path, pdbqt_path):
                print(f"    PDBQT 변환 실패"); continue

            dock_res = dock_molecule_all_targets_vina(
                drug_name, pdbqt_path, args.output_dir,
                receptor_pdbqts, boxes, pdb_targets,
                args.exhaustiveness, args.refine_seeds,
            )
            flat = flatten_vina_results(dock_res)
            sel = compute_vina_selectivity(flat, pdb_targets)

            entry = {
                "drug_name": drug_name,
                "generation": drug_info["generation"],
                "mechanism": drug_info["mechanism"],
                "smiles": drug_info["smiles"],
                **flat, **sel,
            }
            ref_results.append(entry)

            for pk in pdb_targets:
                v = flat.get(f"{pk}_vina_aff", np.nan)
                std = flat.get(f"{pk}_vina_aff_std", 0)
                if not np.isnan(v):
                    print(f"    {pk:20s}: {v:.2f} +/- {std:.3f} kcal/mol")
            if "vina_selectivity" in sel:
                print(f"    선택성 (WT-T790M): {sel['vina_selectivity']:+.2f}")

        df_ref = pd.DataFrame(ref_results)
        ref_path = os.path.join(args.output_dir, "reference_drugs_vina.csv")
        df_ref.to_csv(ref_path, index=False)
        print(f"\n기준 약물 Vina 결과 저장: {ref_path}")

        # 요약
        print("\n" + "=" * 70)
        print("기준 약물 Vina 요약 (kcal/mol, 음수가 좋음)")
        print("=" * 70)
        print(f"\n{'Drug':25s}", end="")
        for pk in pdb_targets:
            print(f" | {pk.split('_')[0]:>8s}", end="")
        print(" | sel(W-T)")
        print("-" * (25 + len(pdb_targets) * 11 + 12))
        for _, row in df_ref.iterrows():
            print(f"{row['drug_name']:25s}", end="")
            for pk in pdb_targets:
                v = row.get(f"{pk}_vina_aff", np.nan)
                print(f" | {v:8.2f}" if not np.isnan(v) else " |      N/A", end="")
            s = row.get("vina_selectivity", np.nan)
            print(f" | {s:+.2f}" if not np.isnan(s) else " |    N/A")
        print("=" * 70)

    # ================================================================
    # Step 2: 생성 분자 전체 1차 Vina 스크리닝
    # ================================================================
    df = pd.read_csv(args.candidates)
    df = df.head(args.top_n)
    print(f"\n[Step 2] 생성 분자 Vina 1차 스크리닝 — {len(df)} molecules x {len(pdb_targets)} PDB")

    print("\n  SMILES -> SDF -> PDBQT 변환 중...")
    pdbqt_paths = {}
    for idx, row in df.iterrows():
        sdf_path = os.path.join(sdf_dir, f"mol_{idx:04d}.sdf")
        pdbqt_path = os.path.join(pdbqt_dir, f"mol_{idx:04d}.pdbqt")
        if smiles_to_sdf(row["smiles"], sdf_path) and sdf_to_pdbqt(sdf_path, pdbqt_path):
            pdbqt_paths[idx] = pdbqt_path
        else:
            print(f"    FAILED: mol_{idx:04d}")
    print(f"  변환 성공: {len(pdbqt_paths)} / {len(df)}")

    results = []
    for i, (idx, pdbqt_path) in enumerate(pdbqt_paths.items()):
        row = df.loc[idx]
        smiles = row["smiles"]
        print(f"\n  [{i+1}/{len(pdbqt_paths)}] {smiles[:50]}...")

        dock_res = dock_molecule_all_targets_vina(
            f"mol_{idx:04d}", pdbqt_path, args.output_dir,
            receptor_pdbqts, boxes, pdb_targets,
            args.exhaustiveness, seeds=None,
        )
        flat = flatten_vina_results(dock_res)
        sel = compute_vina_selectivity(flat, pdb_targets)

        entry = {"mol_idx": idx, "smiles": smiles, "reward": row.get("reward", 0)}
        for col in ["EGFR", "Stab", "Perm", "EGFR_sim", "MW", "LogP", "TPSA",
                     "SAscore", "BBB_score", "max_Tanimoto_vs_drugs", "sampling_group"]:
            if col in row:
                entry[col] = row[col]
        entry.update(flat)
        entry.update(sel)
        results.append(entry)

        t = sel.get("T790M_vina_avg", np.nan)
        w = sel.get("WT_vina_avg", np.nan)
        s = sel.get("vina_selectivity", np.nan)
        print(f"    T790M={t:.2f} | WT={w:.2f} | sel={s:+.2f}" if not np.isnan(s) else "    일부 실패")

    df_results = pd.DataFrame(results)
    results_path = os.path.join(args.output_dir, "vina_docking_results.csv")
    df_results.to_csv(results_path, index=False)
    print(f"\n1차 Vina 결과 저장: {results_path}")

    # ================================================================
    # 1차 Vina 결과 보고서
    # ================================================================
    print("\n" + "=" * 70)
    print("Vina 1차 스크리닝 결과 보고서")
    print("=" * 70)

    docked = df_results.dropna(subset=["T790M_vina_avg", "WT_vina_avg"])
    print(f"\nDocking 성공: {len(docked)} / {len(df_results)}")

    if len(docked) > 0:
        print(f"\n--- Vina 통계 (kcal/mol, 음수가 좋음) ---")
        print(f"  T790M avg: {docked['T790M_vina_avg'].mean():.2f} +/- {docked['T790M_vina_avg'].std():.2f}")
        print(f"  WT avg:    {docked['WT_vina_avg'].mean():.2f} +/- {docked['WT_vina_avg'].std():.2f}")
        print(f"  선택성:    {docked['vina_selectivity'].mean():.2f} +/- {docked['vina_selectivity'].std():.2f}")
        selective = (docked["vina_selectivity"] > 0).sum()
        print(f"  T790M 선택적: {selective} / {len(docked)} ({100*selective/len(docked):.1f}%)")

        if ref_results:
            print(f"\n--- vs 기존 약물 (4ZAU Vina) ---")
            col_4zau = "4ZAU_T790M_vina_aff"
            if col_4zau in docked.columns:
                for ref in ref_results:
                    ref_val = ref.get("4ZAU_T790M_vina_aff", np.nan)
                    if np.isnan(ref_val):
                        continue
                    better = (docked[col_4zau] < ref_val).sum()
                    print(f"  vs {ref['drug_name']:25s} ({ref_val:.2f}): {better:3d} / {len(docked)} 초과 ({100*better/len(docked):.1f}%)")

        print(f"\n--- Top 5 (4ZAU Vina, 가장 음수) ---")
        col_4zau = "4ZAU_T790M_vina_aff"
        if col_4zau in docked.columns:
            for i, (_, row) in enumerate(docked.nsmallest(5, col_4zau).iterrows()):
                print(f"  #{i+1} mol_{row['mol_idx']:04d} | 4ZAU_vina={row[col_4zau]:.2f} | sel={row.get('vina_selectivity', np.nan):+.2f} | EGFR(ML)={row.get('EGFR', 'N/A')}")

    print("=" * 70)

    # ================================================================
    # Step 3: 2차 정밀 Vina 도킹
    # ================================================================
    if args.skip_refine or len(docked) == 0:
        if args.skip_refine:
            print("\n2차 정밀 도킹 건너뜀")
        else:
            print("\n1차 도킹 성공 분자 없음")

        # Step 4로 진행
        self_crossval(args, df_results, ref_results, pdb_targets)
        return

    refine_n = min(args.refine_top_n, len(docked))
    seeds = args.refine_seeds
    print(f"\n[Step 3] 2차 Vina 정밀 도킹 — Top {refine_n}, exhaustiveness={args.refine_exhaustiveness}, seeds={seeds}")

    col_4zau = "4ZAU_T790M_vina_aff"
    top_candidates = docked.nsmallest(refine_n, col_4zau) if col_4zau in docked.columns else docked.nsmallest(refine_n, "T790M_vina_avg")

    refined_results = []
    for i, (_, top_row) in enumerate(top_candidates.iterrows()):
        mol_idx = top_row["mol_idx"]
        pdbqt_path = pdbqt_paths.get(mol_idx)
        if pdbqt_path is None:
            continue
        screen_val = top_row.get(col_4zau, "N/A")
        print(f"\n  [{i+1}/{refine_n}] mol_{mol_idx:04d} (1차 4ZAU={screen_val})")

        dock_res = dock_molecule_all_targets_vina(
            f"mol_{mol_idx:04d}_refined", pdbqt_path, args.output_dir,
            receptor_pdbqts, boxes, pdb_targets,
            args.refine_exhaustiveness, seeds,
        )
        flat = flatten_vina_results(dock_res, prefix="refined_")
        sel_raw = flatten_vina_results(dock_res)
        sel = compute_vina_selectivity(sel_raw, pdb_targets)

        entry = {
            "mol_idx": mol_idx,
            "smiles": top_row["smiles"],
            "reward": top_row.get("reward", 0),
            "EGFR": top_row.get("EGFR", None),
            "sampling_group": top_row.get("sampling_group", None),
        }
        for pk in pdb_targets:
            entry[f"screen_{pk}_vina_aff"] = top_row.get(f"{pk}_vina_aff", np.nan)
        entry.update(flat)
        entry.update({f"refined_{k}": v for k, v in sel.items()})
        refined_results.append(entry)

        t = sel.get("T790M_vina_avg", np.nan)
        w = sel.get("WT_vina_avg", np.nan)
        s = sel.get("vina_selectivity", np.nan)
        print(f"    T790M={t:.2f} | WT={w:.2f} | sel={s:+.2f}" if not np.isnan(s) else "    일부 실패")

    if refined_results:
        df_refined = pd.DataFrame(refined_results)
        refined_path = os.path.join(args.output_dir, "vina_docking_refined.csv")
        df_refined.to_csv(refined_path, index=False)
        print(f"\n2차 Vina 결과 저장: {refined_path}")

        print("\n" + "=" * 70)
        print("Vina 2차 정밀 결과 보고서")
        print("=" * 70)
        print(f"\n재현성 (seed 간 변동)")
        for pk in pdb_targets:
            std_col = f"refined_{pk}_vina_aff_std"
            if std_col in df_refined.columns:
                ms = df_refined[std_col].mean()
                print(f"  {pk:20s}: 평균 std={ms:.3f} kcal/mol")

    # Step 4
    self_crossval(args, df_results, ref_results, pdb_targets)


def self_crossval(args, df_vina_results, ref_results, pdb_targets):
    """GNINA vs Vina cross-validation."""
    print("\n" + "=" * 70)
    print("GNINA vs Vina Cross-Validation")
    print("=" * 70)

    # 기준 약물 비교
    if ref_results and os.path.exists(args.gnina_reference):
        print("\n--- 기준 약물: GNINA vs Vina ---")
        df_gnina_ref = pd.read_csv(args.gnina_reference)
        df_vina_ref = pd.DataFrame(ref_results)
        print(f"\n{'Drug':25s} | {'GNINA 4ZAU':>10s} | {'Vina 4ZAU':>10s} | {'GNINA sel':>10s} | {'Vina sel':>10s}")
        print("-" * 75)
        for _, vrow in df_vina_ref.iterrows():
            drug_name = vrow["drug_name"]
            gm = df_gnina_ref[df_gnina_ref["drug_name"] == drug_name]
            g_4zau = gm["4ZAU_T790M_cnn_aff"].values[0] if len(gm) > 0 and "4ZAU_T790M_cnn_aff" in gm.columns else np.nan
            g_sel = gm["selectivity_T790M_vs_WT"].values[0] if len(gm) > 0 and "selectivity_T790M_vs_WT" in gm.columns else np.nan
            v_4zau = vrow.get("4ZAU_T790M_vina_aff", np.nan)
            v_sel = vrow.get("vina_selectivity", np.nan)

            g4 = f"{g_4zau:10.2f}" if not np.isnan(g_4zau) else "       N/A"
            v4 = f"{v_4zau:10.2f}" if not np.isnan(v_4zau) else "       N/A"
            gs = f"{g_sel:+10.2f}" if not np.isnan(g_sel) else "       N/A"
            vs = f"{v_sel:+10.2f}" if not np.isnan(v_sel) else "       N/A"
            print(f"{drug_name:25s} | {g4} | {v4} | {gs} | {vs}")
    else:
        if not os.path.exists(args.gnina_reference):
            print("\nGNINA 기준 약물 결과 없음 — GNINA 도킹 완료 후 다시 실행하세요")

    # 생성 분자 상관관계
    if os.path.exists(args.gnina_results) and len(df_vina_results) > 0:
        print(f"\n--- 생성 분자: GNINA vs Vina 상관관계 ---")
        df_gnina_all = pd.read_csv(args.gnina_results)

        df_merged = pd.merge(
            df_gnina_all[["mol_idx"] + [c for c in df_gnina_all.columns if "cnn_aff" in c and "std" not in c]],
            df_vina_results[["mol_idx"] + [c for c in df_vina_results.columns if "vina_aff" in c and "std" not in c]],
            on="mol_idx", how="inner",
        )
        print(f"  Matched molecules: {len(df_merged)}")

        for pk in pdb_targets:
            gnina_col = f"{pk}_cnn_aff"
            vina_col = f"{pk}_vina_aff"
            if gnina_col in df_merged.columns and vina_col in df_merged.columns:
                valid = df_merged[[gnina_col, vina_col]].dropna()
                if len(valid) >= 5:
                    pearson = valid[gnina_col].corr(valid[vina_col])
                    spearman = valid[gnina_col].rank().corr(valid[vina_col].rank())
                    print(f"  {pk:20s}: Pearson r={pearson:.3f}, Spearman rho={spearman:.3f} (n={len(valid)})")

        # GNINA Top 10 vs Vina 순위
        g4zau = "4ZAU_T790M_cnn_aff"
        v4zau = "4ZAU_T790M_vina_aff"
        if g4zau in df_merged.columns and v4zau in df_merged.columns:
            df_rank = df_merged[["mol_idx", g4zau, v4zau]].dropna()
            if len(df_rank) > 0:
                df_rank["gnina_rank"] = df_rank[g4zau].rank(ascending=False)
                df_rank["vina_rank"] = df_rank[v4zau].rank(ascending=True)
                print(f"\n  GNINA Top 10 -> Vina 순위:")
                for _, row in df_rank.nlargest(10, g4zau).iterrows():
                    print(f"    mol_{int(row['mol_idx']):04d}: GNINA #{int(row['gnina_rank'])}, Vina #{int(row['vina_rank'])}")
    else:
        if not os.path.exists(args.gnina_results):
            print("\nGNINA 1차 결과 없음 — GNINA 도킹 완료 후 다시 실행하세요")

    print("\n" + "=" * 70)
    print("완료. 결과 파일:")
    print(f"  기준 약물 Vina:   {os.path.join(args.output_dir, 'reference_drugs_vina.csv')}")
    vr = os.path.join(args.output_dir, "vina_docking_results.csv")
    if os.path.exists(vr):
        print(f"  1차 Vina 전체:    {vr}")
    vrr = os.path.join(args.output_dir, "vina_docking_refined.csv")
    if os.path.exists(vrr):
        print(f"  2차 Vina 정밀:    {vrr}")
    print("=" * 70)


if __name__ == "__main__":
    main()