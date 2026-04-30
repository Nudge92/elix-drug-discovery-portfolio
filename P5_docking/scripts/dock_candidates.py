"""P5 생성 분자 Docking 자동화 스크립트
v4: 7개 PDB 전체 도킹 + 1~4세대 re-docking + multi-seed 정밀 도킹 + 종합 보고서

PDB 구조:
  - 2ITY: Gefitinib + WT EGFR (1세대 reference)
  - 4WKQ: Gefitinib + WT EGFR (1.85Å 고해상도)
  - 4HJO: Erlotinib + EGFR inactive conformation
  - 4I23: Dacomitinib + WT (covalent 2세대)
  - 4G5J: Afatinib + T790M (2세대 covalent)
  - 4ZAU: Osimertinib + T790M (3세대, 핵심)
  - 5D41: EAI045 + T790M allosteric (4세대)

분석 관점:
  - T790M 활성 (4ZAU, 4G5J): 변이체 결합력
  - WT 선택성 (2ITY, 4WKQ): wild-type 회피
  - Inactive conformation (4HJO): conformation 선택성
  - Allosteric (5D41): 대안적 결합 모드
  - Cross-PDB robustness: 같은 타겟, 다른 PDB에서 결과 일관성
"""
import argparse
import os
import subprocess
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem

# ================================================================
# PDB 구조 정의 (7개 전체)
# ================================================================
PDB_STRUCTURES = {
    "2ITY_WT": {
        "protein": "2ITY_protein.pdb",
        "autobox_ligand": "2ITY_IRE.pdb",
        "description": "WT EGFR active (Gefitinib, 1G reference)",
        "target_type": "WT",
        "conformation": "active",
    },
    "4WKQ_WT": {
        "protein": "4WKQ_protein.pdb",
        "autobox_ligand": "4WKQ_IRE.pdb",
        "description": "WT EGFR active (Gefitinib, 1.85Å high-res)",
        "target_type": "WT",
        "conformation": "active",
    },
    "4HJO_inactive": {
        "protein": "4HJO_protein.pdb",
        "autobox_ligand": "4HJO_AQ4.pdb",
        "description": "EGFR inactive (Erlotinib)",
        "target_type": "WT",
        "conformation": "inactive",
    },
    "4I23_WT": {
        "protein": "4I23_protein.pdb",
        "autobox_ligand": "4I23_1C9.pdb",
        "description": "WT EGFR (Dacomitinib, 2G covalent)",
        "target_type": "WT",
        "conformation": "active",
    },
    "4G5J_T790M": {
        "protein": "4G5J_protein.pdb",
        "autobox_ligand": "4G5J_0WM.pdb",
        "description": "T790M mutant (Afatinib, 2G covalent)",
        "target_type": "T790M",
        "conformation": "active",
    },
    "4ZAU_T790M": {
        "protein": "4ZAU_protein.pdb",
        "autobox_ligand": "4ZAU_YY3.pdb",
        "description": "T790M mutant (Osimertinib, 3G, 핵심)",
        "target_type": "T790M",
        "conformation": "active",
    },
    "5D41_allosteric": {
        "protein": "5D41_protein.pdb",
        "autobox_ligand": "5D41_57N.pdb",
        "description": "T790M allosteric site (EAI045, 4G)",
        "target_type": "T790M_allosteric",
        "conformation": "allosteric",
    },
}

# ================================================================
# 1~4세대 EGFR 억제제 기준 약물
# ================================================================
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
        description="P5 Docking: 7 PDB × 1~4세대 re-docking + 생성 분자 도킹 + 정밀 검증")
    parser.add_argument("--candidates", type=str, required=True)
    parser.add_argument("--top_n", type=int, default=50)
    parser.add_argument("--gnina_path", type=str,
                        default="/home/nudge/Project/elix/docking/bin/gnina")
    parser.add_argument("--docking_dir", type=str,
                        default="/home/nudge/Project/elix/docking")
    parser.add_argument("--output_dir", type=str,
                        default="/home/nudge/Project/elix/docking/results/generated_molecules")
    parser.add_argument("--exhaustiveness", type=int, default=16)
    parser.add_argument("--refine_top_n", type=int, default=30,
                        help="2차 정밀 도킹 대상 수")
    parser.add_argument("--refine_exhaustiveness", type=int, default=32)
    parser.add_argument("--refine_seeds", type=int, nargs="+",
                        default=[0, 42, 123])
    parser.add_argument("--skip_refine", action="store_true",
                        help="2차 정밀 도킹 건너뛰기")
    parser.add_argument("--skip_reference", action="store_true",
                        help="1~4세대 기준 약물 도킹 건너뛰기")
    parser.add_argument("--pdb_targets", type=str, nargs="+",
                        default=list(PDB_STRUCTURES.keys()),
                        help="사용할 PDB 구조 (기본: 전체 7개)")
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


def run_gnina(gnina_path, receptor_pdb, ligand_sdf, autobox_ligand,
              output_sdf, exhaustiveness=16, seed=0):
    cmd = [gnina_path, "-r", receptor_pdb, "-l", ligand_sdf,
           "--autobox_ligand", autobox_ligand, "--autobox_add", "4",
           "--exhaustiveness", str(exhaustiveness), "--seed", str(seed),
           "-o", output_sdf]
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
        output = result.stdout + result.stderr
        for line in output.split("\n"):
            line = line.strip()
            if line and line[0:1].isdigit():
                parts = line.split()
                if len(parts) >= 5 and parts[0] == "1":
                    return {
                        "affinity": float(parts[1]),
                        "cnn_pose_score": float(parts[3]),
                        "cnn_affinity": float(parts[4]),
                    }
    except Exception as e:
        print(f"    ERROR: {e}")
    return None


def resolve_pdb_paths(docking_dir, pdb_key):
    """PDB 키로부터 protein, autobox_ligand 경로를 반환."""
    info = PDB_STRUCTURES[pdb_key]
    protein = os.path.join(docking_dir, "pdb", info["protein"])
    autobox = os.path.join(docking_dir, "ligands", info["autobox_ligand"])
    return protein, autobox


def dock_molecule_all_targets(name, sdf_path, output_dir, gnina_path,
                              docking_dir, pdb_targets, exhaustiveness,
                              seeds=None):
    """단일 분자를 모든 PDB 타겟에 도킹. seeds가 주어지면 multi-seed."""
    if seeds is None:
        seeds = [0]

    results = {}
    for pdb_key in pdb_targets:
        protein, autobox = resolve_pdb_paths(docking_dir, pdb_key)
        affs, poses = [], []

        for s in seeds:
            seed_tag = f"_s{s}" if len(seeds) > 1 else ""
            out_path = os.path.join(output_dir, f"{name}_{pdb_key}{seed_tag}.sdf.gz")
            res = run_gnina(gnina_path, protein, sdf_path, autobox, out_path,
                            exhaustiveness, s)
            if res:
                affs.append(res["cnn_affinity"])
                poses.append(res["cnn_pose_score"])

        if affs:
            results[pdb_key] = {
                "cnn_aff_mean": np.mean(affs),
                "cnn_aff_std": np.std(affs) if len(seeds) > 1 else 0.0,
                "cnn_pose_mean": np.mean(poses),
                "n_seeds_ok": len(affs),
            }
        else:
            results[pdb_key] = None

    return results


def flatten_dock_results(results, prefix=""):
    """도킹 결과 dict를 flat dict로 변환."""
    flat = {}
    for pdb_key, res in results.items():
        p = f"{prefix}{pdb_key}" if prefix else pdb_key
        if res:
            flat[f"{p}_cnn_aff"] = res["cnn_aff_mean"]
            flat[f"{p}_cnn_aff_std"] = res["cnn_aff_std"]
            flat[f"{p}_cnn_pose"] = res["cnn_pose_mean"]
        else:
            flat[f"{p}_cnn_aff"] = np.nan
            flat[f"{p}_cnn_aff_std"] = np.nan
            flat[f"{p}_cnn_pose"] = np.nan
    return flat


def compute_selectivity_metrics(flat, pdb_targets):
    """T790M vs WT 선택성 등 파생 메트릭 계산."""
    metrics = {}

    # T790M 구조들의 평균
    t790m_keys = [k for k in pdb_targets
                  if PDB_STRUCTURES[k]["target_type"] == "T790M"]
    wt_keys = [k for k in pdb_targets
               if PDB_STRUCTURES[k]["target_type"] == "WT"]

    t790m_vals = [flat[f"{k}_cnn_aff"] for k in t790m_keys
                  if not np.isnan(flat.get(f"{k}_cnn_aff", np.nan))]
    wt_vals = [flat[f"{k}_cnn_aff"] for k in wt_keys
               if not np.isnan(flat.get(f"{k}_cnn_aff", np.nan))]

    if t790m_vals:
        metrics["T790M_avg_aff"] = np.mean(t790m_vals)
    if wt_vals:
        metrics["WT_avg_aff"] = np.mean(wt_vals)
    if t790m_vals and wt_vals:
        metrics["selectivity_T790M_vs_WT"] = np.mean(t790m_vals) - np.mean(wt_vals)

    # 4ZAU vs 4G5J cross-PDB consistency
    zau = flat.get("4ZAU_T790M_cnn_aff", np.nan)
    g5j = flat.get("4G5J_T790M_cnn_aff", np.nan)
    if not np.isnan(zau) and not np.isnan(g5j):
        metrics["T790M_cross_pdb_diff"] = abs(zau - g5j)

    # 2ITY vs 4WKQ cross-PDB consistency
    ity = flat.get("2ITY_WT_cnn_aff", np.nan)
    wkq = flat.get("4WKQ_WT_cnn_aff", np.nan)
    if not np.isnan(ity) and not np.isnan(wkq):
        metrics["WT_cross_pdb_diff"] = abs(ity - wkq)

    # Allosteric binding
    allo = flat.get("5D41_allosteric_cnn_aff", np.nan)
    if not np.isnan(allo):
        metrics["allosteric_aff"] = allo

    return metrics


def main():
    args = parse_args()
    docking_dir = args.docking_dir
    pdb_targets = args.pdb_targets
    os.makedirs(args.output_dir, exist_ok=True)
    sdf_dir = os.path.join(args.output_dir, "sdf_files")
    os.makedirs(sdf_dir, exist_ok=True)

    # PDB 파일 존재 확인
    print("=" * 70)
    print("PDB 구조 확인")
    print("=" * 70)
    valid_targets = []
    for pdb_key in pdb_targets:
        protein, autobox = resolve_pdb_paths(docking_dir, pdb_key)
        p_ok = os.path.exists(protein)
        a_ok = os.path.exists(autobox)
        status = "OK" if (p_ok and a_ok) else f"MISSING (protein={p_ok}, ligand={a_ok})"
        print(f"  {pdb_key:20s}: {status}")
        if p_ok and a_ok:
            valid_targets.append(pdb_key)

    if not valid_targets:
        print("\n사용 가능한 PDB 구조 없음 — 종료")
        return

    pdb_targets = valid_targets
    print(f"\n사용할 PDB 구조: {len(pdb_targets)}개")
    total_docks_estimate = (
        (7 * len(args.refine_seeds) * len(pdb_targets) if not args.skip_reference else 0) +
        args.top_n * len(pdb_targets) +
        (args.refine_top_n * len(args.refine_seeds) * len(pdb_targets) if not args.skip_refine else 0)
    )
    print(f"예상 총 도킹 횟수: ~{total_docks_estimate}회")

    # ================================================================
    # Step 0: 1~4세대 기준 약물 Re-docking
    # ================================================================
    ref_results = []
    if not args.skip_reference:
        print("\n" + "=" * 70)
        print("[Step 0] 1~4세대 기준 약물 Re-docking")
        print(f"  PDB: {len(pdb_targets)}개, exhaustiveness={args.exhaustiveness}, seeds={args.refine_seeds}")
        print("=" * 70)

        ref_sdf_dir = os.path.join(args.output_dir, "reference_sdf")
        os.makedirs(ref_sdf_dir, exist_ok=True)

        for drug_name, drug_info in REFERENCE_DRUGS.items():
            print(f"\n  {drug_name} ({drug_info['generation']}, {drug_info['mechanism']})")

            # SDF 변환
            sdf_path = os.path.join(ref_sdf_dir, f"{drug_name}.sdf")
            if not smiles_to_sdf(drug_info["smiles"], sdf_path):
                print(f"    SDF 변환 실패 — 건너뜀")
                continue

            # 모든 PDB에 도킹
            dock_res = dock_molecule_all_targets(
                name=drug_name, sdf_path=sdf_path, output_dir=args.output_dir,
                gnina_path=args.gnina_path, docking_dir=docking_dir,
                pdb_targets=pdb_targets, exhaustiveness=args.exhaustiveness,
                seeds=args.refine_seeds,
            )

            flat = flatten_dock_results(dock_res)
            sel_metrics = compute_selectivity_metrics(flat, pdb_targets)

            entry = {
                "drug_name": drug_name,
                "generation": drug_info["generation"],
                "mechanism": drug_info["mechanism"],
                "smiles": drug_info["smiles"],
                **flat,
                **sel_metrics,
            }
            ref_results.append(entry)

            # 요약 출력
            for pdb_key in pdb_targets:
                aff = flat.get(f"{pdb_key}_cnn_aff", np.nan)
                std = flat.get(f"{pdb_key}_cnn_aff_std", 0)
                if not np.isnan(aff):
                    print(f"    {pdb_key:20s}: {aff:.2f} ± {std:.3f}")
            if "selectivity_T790M_vs_WT" in sel_metrics:
                print(f"    선택성 (T790M-WT): {sel_metrics['selectivity_T790M_vs_WT']:+.2f}")

        df_ref = pd.DataFrame(ref_results)
        ref_path = os.path.join(args.output_dir, "reference_drugs_docking.csv")
        df_ref.to_csv(ref_path, index=False)
        print(f"\n기준 약물 결과 저장: {ref_path}")

        # 기준 약물 요약 보고서
        print("\n" + "=" * 70)
        print("1~4세대 기준 약물 Docking 요약")
        print("=" * 70)
        print(f"\n{'Drug':25s}", end="")
        for pk in pdb_targets:
            label = pk.split("_")[0]
            print(f" | {label:>8s}", end="")
        print(" | sel(T-W)")
        print("-" * (25 + len(pdb_targets) * 11 + 12))
        for _, row in df_ref.iterrows():
            print(f"{row['drug_name']:25s}", end="")
            for pk in pdb_targets:
                v = row.get(f"{pk}_cnn_aff", np.nan)
                print(f" | {v:8.2f}" if not np.isnan(v) else " |      N/A", end="")
            sel = row.get("selectivity_T790M_vs_WT", np.nan)
            print(f" | {sel:+.2f}" if not np.isnan(sel) else " |    N/A")
        print("=" * 70)
    else:
        print("\n1~4세대 기준 약물 도킹 건너뜀 (--skip_reference)")

    # ================================================================
    # Step 1: 생성 분자 SMILES → 3D SDF 변환
    # ================================================================
    df = pd.read_csv(args.candidates)
    df = df.head(args.top_n)
    print(f"\nLoaded {len(df)} candidate molecules")

    print("\n[Step 1] Converting SMILES to 3D SDF...")
    sdf_paths = {}
    for idx, row in df.iterrows():
        sdf_path = os.path.join(sdf_dir, f"mol_{idx:04d}.sdf")
        if smiles_to_sdf(row["smiles"], sdf_path):
            sdf_paths[idx] = sdf_path
        else:
            print(f"  FAILED: {row['smiles']}")
    print(f"  Converted: {len(sdf_paths)} / {len(df)}")

    # ================================================================
    # Step 2: 1차 스크리닝 도킹 (전체, single seed, 7 PDB)
    # ================================================================
    results = []
    print(f"\n[Step 2] 1차 스크리닝 — {len(sdf_paths)} molecules × {len(pdb_targets)} PDB, exhaustiveness={args.exhaustiveness}")

    for i, (idx, sdf_path) in enumerate(sdf_paths.items()):
        row = df.loc[idx]
        smiles = row["smiles"]
        print(f"\n  [{i+1}/{len(sdf_paths)}] {smiles[:50]}...")

        dock_res = dock_molecule_all_targets(
            name=f"mol_{idx:04d}", sdf_path=sdf_path, output_dir=args.output_dir,
            gnina_path=args.gnina_path, docking_dir=docking_dir,
            pdb_targets=pdb_targets, exhaustiveness=args.exhaustiveness,
            seeds=None,  # single seed
        )

        flat = flatten_dock_results(dock_res)
        sel_metrics = compute_selectivity_metrics(flat, pdb_targets)

        entry = {"mol_idx": idx, "smiles": smiles, "reward": row.get("reward", 0)}
        for col in ["EGFR", "Stab", "Perm", "EGFR_sim", "MW", "LogP", "TPSA",
                     "SAscore", "BBB_score", "max_Tanimoto_vs_drugs", "sampling_group"]:
            if col in row:
                entry[col] = row[col]
        entry.update(flat)
        entry.update(sel_metrics)
        results.append(entry)

        # 간단 요약
        t790m = sel_metrics.get("T790M_avg_aff", np.nan)
        wt = sel_metrics.get("WT_avg_aff", np.nan)
        sel = sel_metrics.get("selectivity_T790M_vs_WT", np.nan)
        t_str = f"{t790m:.2f}" if not np.isnan(t790m) else "N/A"
        w_str = f"{wt:.2f}" if not np.isnan(wt) else "N/A"
        s_str = f"{sel:+.2f}" if not np.isnan(sel) else "N/A"
        print(f"    T790M_avg={t_str} | WT_avg={w_str} | sel={s_str}")

    df_results = pd.DataFrame(results)
    results_path = os.path.join(args.output_dir, "docking_results.csv")
    df_results.to_csv(results_path, index=False)
    print(f"\n1차 결과 저장: {results_path}")

    # ================================================================
    # 1차 결과 보고서
    # ================================================================
    print("\n" + "=" * 70)
    print("1차 스크리닝 결과 보고서")
    print("=" * 70)

    docked = df_results.dropna(subset=["T790M_avg_aff", "WT_avg_aff"])
    print(f"\nDocking 성공 (T790M+WT 모두): {len(docked)} / {len(df_results)}")

    if len(docked) > 0:
        print(f"\n--- 생성 분자 통계 ---")
        print(f"  T790M avg aff: {docked['T790M_avg_aff'].mean():.2f} ± {docked['T790M_avg_aff'].std():.2f}")
        print(f"  WT avg aff:    {docked['WT_avg_aff'].mean():.2f} ± {docked['WT_avg_aff'].std():.2f}")
        print(f"  선택성 (T-W):  {docked['selectivity_T790M_vs_WT'].mean():.2f} ± {docked['selectivity_T790M_vs_WT'].std():.2f}")
        selective = (docked["selectivity_T790M_vs_WT"] > 0).sum()
        print(f"  T790M 선택적:  {selective} / {len(docked)} ({100*selective/len(docked):.1f}%)")

        # Cross-PDB robustness
        if "T790M_cross_pdb_diff" in docked.columns:
            print(f"\n--- Cross-PDB Robustness ---")
            print(f"  T790M (4ZAU vs 4G5J) 평균 차이: {docked['T790M_cross_pdb_diff'].mean():.3f}")
        if "WT_cross_pdb_diff" in docked.columns:
            print(f"  WT (2ITY vs 4WKQ) 평균 차이:    {docked['WT_cross_pdb_diff'].mean():.3f}")

        # 기준 약물 대비 비교
        if ref_results:
            print(f"\n--- 생성 분자 vs 기존 약물 (T790M avg 기준) ---")
            for ref in ref_results:
                ref_val = ref.get("T790M_avg_aff", np.nan)
                if np.isnan(ref_val):
                    continue
                ref_name = ref["drug_name"]
                better = (docked["T790M_avg_aff"] > ref_val).sum()
                print(f"  vs {ref_name:25s} ({ref_val:.2f}): {better:3d} / {len(docked)} 초과 ({100*better/len(docked):.1f}%)")

        # Top 5
        print(f"\n--- Top 5 (4ZAU T790M CNN Affinity) ---")
        col_4zau = "4ZAU_T790M_cnn_aff"
        if col_4zau in docked.columns:
            for i, (_, row) in enumerate(docked.nlargest(5, col_4zau).iterrows()):
                print(f"  #{i+1} mol_{row['mol_idx']:04d} | 4ZAU={row[col_4zau]:.2f} | sel={row.get('selectivity_T790M_vs_WT', np.nan):+.2f} | EGFR(ML)={row.get('EGFR', 'N/A')}")

    print("=" * 70)

    # ================================================================
    # Step 3: 2차 정밀 도킹 (Top N, multi-seed, high exhaustiveness)
    # ================================================================
    if args.skip_refine:
        print("\n2차 정밀 도킹 건너뜀 (--skip_refine)")
        return

    if len(docked) == 0:
        print("\n1차 도킹 성공 분자 없음 — 2차 건너뜀")
        return

    refine_n = min(args.refine_top_n, len(docked))
    seeds = args.refine_seeds
    print(f"\n[Step 3] 2차 정밀 도킹 — Top {refine_n} (4ZAU 기준), exhaustiveness={args.refine_exhaustiveness}, seeds={seeds}")

    top_candidates = docked.nlargest(refine_n, "4ZAU_T790M_cnn_aff") if "4ZAU_T790M_cnn_aff" in docked.columns else docked.nlargest(refine_n, "T790M_avg_aff")
    refined_results = []

    for i, (_, top_row) in enumerate(top_candidates.iterrows()):
        mol_idx = top_row["mol_idx"]
        sdf_path = sdf_paths.get(mol_idx)
        if sdf_path is None:
            continue
        print(f"\n  [{i+1}/{refine_n}] mol_{mol_idx:04d}")

        dock_res = dock_molecule_all_targets(
            name=f"mol_{mol_idx:04d}_refined", sdf_path=sdf_path,
            output_dir=args.output_dir, gnina_path=args.gnina_path,
            docking_dir=docking_dir, pdb_targets=pdb_targets,
            exhaustiveness=args.refine_exhaustiveness, seeds=seeds,
        )

        flat = flatten_dock_results(dock_res, prefix="refined_")
        sel_metrics_raw = flatten_dock_results(dock_res)
        sel_metrics = compute_selectivity_metrics(sel_metrics_raw, pdb_targets)

        entry = {
            "mol_idx": mol_idx,
            "smiles": top_row["smiles"],
            "reward": top_row.get("reward", 0),
            "EGFR": top_row.get("EGFR", None),
            "sampling_group": top_row.get("sampling_group", None),
        }
        # 1차 결과 참조
        for pk in pdb_targets:
            entry[f"screen_{pk}_cnn_aff"] = top_row.get(f"{pk}_cnn_aff", np.nan)

        entry.update(flat)
        entry.update({f"refined_{k}": v for k, v in sel_metrics.items()})
        refined_results.append(entry)

        t790m = sel_metrics.get("T790M_avg_aff", np.nan)
        wt = sel_metrics.get("WT_avg_aff", np.nan)
        sel = sel_metrics.get("selectivity_T790M_vs_WT", np.nan)
        print(f"    T790M_avg={t790m:.2f} | WT_avg={wt:.2f} | sel={sel:+.2f}" if not np.isnan(sel) else "    일부 실패")

    # ================================================================
    # 2차 결과 저장 및 최종 보고서
    # ================================================================
    if refined_results:
        df_refined = pd.DataFrame(refined_results)
        refined_path = os.path.join(args.output_dir, "docking_results_refined.csv")
        df_refined.to_csv(refined_path, index=False)
        print(f"\n정밀 결과 저장: {refined_path}")

        print("\n" + "=" * 70)
        print("2차 정밀 도킹 결과 보고서")
        print("=" * 70)

        # 재현성 — 각 PDB별 std 확인
        print(f"\n재현성 (seed 간 변동, std)")
        for pk in pdb_targets:
            std_col = f"refined_{pk}_cnn_aff_std"
            if std_col in df_refined.columns:
                mean_std = df_refined[std_col].mean()
                reliable = (df_refined[std_col] < 0.3).sum()
                print(f"  {pk:20s}: 평균 std={mean_std:.3f}, 신뢰(std<0.3)={reliable}/{len(df_refined)}")

        # 선택성
        sel_col = "refined_selectivity_T790M_vs_WT"
        if sel_col in df_refined.columns:
            selective = (df_refined[sel_col] > 0).sum()
            print(f"\n선택성 (T790M - WT > 0): {selective} / {len(df_refined)}")

        # 기준 약물 대비 최종 비교
        if ref_results:
            print(f"\n--- 정밀 도킹 분자 vs 기존 약물 ---")
            t790m_col = "refined_T790M_avg_aff"
            if t790m_col in df_refined.columns:
                for ref in ref_results:
                    ref_val = ref.get("T790M_avg_aff", np.nan)
                    if np.isnan(ref_val):
                        continue
                    ref_name = ref["drug_name"]
                    better = (df_refined[t790m_col] > ref_val).sum()
                    print(f"  vs {ref_name:25s} ({ref_val:.2f}): {better} / {len(df_refined)} 초과")

        # Top 10
        print(f"\n--- Top 10 (T790M avg, 정밀 도킹) ---")
        sort_col = "refined_T790M_avg_aff"
        if sort_col in df_refined.columns:
            for i, (_, row) in enumerate(df_refined.nlargest(10, sort_col).iterrows()):
                sel = row.get(sel_col, np.nan)
                print(f"  #{i+1} mol_{row['mol_idx']:04d}")
                print(f"    T790M_avg={row[sort_col]:.2f} | sel={sel:+.2f}" if not np.isnan(sel) else f"    T790M_avg={row[sort_col]:.2f}")
                print(f"    EGFR(ML)={row.get('EGFR', 'N/A')} | reward={row.get('reward', 'N/A')}")

        print("=" * 70)

    print("\n완료. 결과 파일:")
    print(f"  기준 약물: {os.path.join(args.output_dir, 'reference_drugs_docking.csv')}")
    print(f"  1차 스크리닝: {results_path}")
    if refined_results:
        print(f"  2차 정밀: {os.path.join(args.output_dir, 'docking_results_refined.csv')}")


if __name__ == "__main__":
    main()