#!/usr/bin/env python3
"""Run selected assertions and regressions without rewriting recorded outputs.

This suite excludes long searches, optional SDP checks, manuscript
compilation, and fresh Lean builds. --upstream adds text/record consistency
checks against the exact recorded Chain Theory commit. --recent includes
selected checks of recent research results.
"""
import argparse
import hashlib
import importlib.metadata
import os
from pathlib import Path
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
SOURCE_COMMIT = "a28db8998bce5bd5ef11195c05f13c86cd9dac5c"
DEPENDENCIES = {"numpy": "2.5.3", "sympy": "1.14.0", "mpmath": "1.3.0",
                "scipy": "1.18.1"}
SELECTED = [
    "chain-theory/alpha/records_independent.py",
    "chain-theory/alpha/config_enumeration.py",
    "chain-theory/alpha/algebraic.py",
    "chain-theory/alpha/mutations_classified.py",
    "chain-theory/alpha/freedom.py",
    "chain-theory/alpha/xcheck_operator_form.py",
    "chain-theory/alpha/conjugation.py",
    "geometric-unity/A_spin_obstruction.py",
    "geometric-unity/xcheck_round5.py",
    "geometric-unity/projection_regression.py",
]


RECENT = [
    ("chain-theory/alpha/records_uniform.py", []),
    ("chain-theory/alpha/xcheck_round_ct5.py", []),
    ("geometric-unity/xcheck_round14.py", []),
    ("geometric-unity/xcheck_round15.py", ["--geometry"]),
    ("geometric-unity/T_singular_residual.py", []),
    ("geometric-unity/U_rank_one_frontier.py", []),
    ("geometric-unity/V_hook_intersections.py", []),
    ("geometric-unity/W_quartic_rank_one.py", []),
    ("geometric-unity/xcheck_W_quartic.py", []),
    ("geometric-unity/X_symmetric_hook.py", []),
    ("geometric-unity/Y_zero_curvature_quartic.py", []),
    ("geometric-unity/Z_cubic_rank_one.py", []),
    ("geometric-unity/AA_conic_traces.py", []),
    ("geometric-unity/AA_exterior_derivative.py", []),
    ("geometric-unity/AA_lambda3_hook.py", []),
    ("geometric-unity/AA_three_form_simple.py", []),
    ("geometric-unity/AF_symmetric_gamma.py", []),
    ("geometric-unity/AG_other_lambda3.py", []),
    ("geometric-unity/AH_lambda5_gram.py", []),
    ("geometric-unity/AH_cyclic_obstruction.py", []),
    ("geometric-unity/AI_flat_lambda3.py", []),
    ("geometric-unity/AJ_curved_endpoint.py", []),
    ("geometric-unity/AK_lambda3_classification.py", []),
    ("geometric-unity/AL_opposite_lambda5_system.py", []),
    ("geometric-unity/AM_global_flat_models.py", []),
    ("geometric-unity/AP_first_lambda5_obstruction.py", []),
    ("geometric-unity/AQ_opposite_lambda5_conservation.py", []),
    ("geometric-unity/AR_opposite_lambda5_geometry.py", []),
    ("geometric-unity/AS_opposite_lambda5_system.py", []),
    ("geometric-unity/AT_opposite_profile_crossing.py", []),
    ("geometric-unity/AU_opposite_ricci_spectrum.py", []),
    ("geometric-unity/AV_opposite_gram_transport.py", []),
    ("geometric-unity/AW_lambda5_conic_classification.py", []),
    ("geometric-unity/AX_hook_derivative.py", []),
    ("geometric-unity/AX_hook_curvature.py", []),
    ("geometric-unity/AX_hook_curvature_xcheck.py", []),
    ("geometric-unity/AZ_flat_lambda.py", []),
    ("geometric-unity/BA_hook_closure.py", []),
    ("geometric-unity/BC_affine_hooks.py", []),
    ("geometric-unity/BC_quadratic_xcheck.py", []),
    ("geometric-unity/BD_h4_affine.py", []),
    ("geometric-unity/BE_h2_affine.py", []),
    ("geometric-unity/BF_h6_casimir.py", []),
    ("geometric-unity/BG_h6_curvature.py", []),
    ("geometric-unity/BG_h6_moments.py", []),
    ("geometric-unity/BH_h6_degree2.py", []),
    ("geometric-unity/BI_h3_stress.py", []),
    ("geometric-unity/BI_mixed_h6.py", []),
    ("geometric-unity/BJ_sphere_boundary.py", []),
    ("geometric-unity/BK_quintic_gram.py", []),
    ("geometric-unity/BK_symmetric_quintics.py", []),
    ("geometric-unity/BL_mixed_stress.py", []),
    ("geometric-unity/BL_symmetric_mixed.py", []),
    ("geometric-unity/BL_real_certificate.py", []),
    ("geometric-unity/BM_mixed_grade9.py", []),
    ("geometric-unity/BN_h6_invariant_grade9.py", []),
    ("geometric-unity/BO_h6_full_kernel.py", []),
    ("geometric-unity/BP_h6_cartan_grade9.py", []),
    ("geometric-unity/BQ_h6_invariant_forms.py", []),
    ("geometric-unity/BR_h6_an_closed.py", ["--parabolic-characters"]),
    ("geometric-unity/BS_h6_partial_value.py", []),
    ("geometric-unity/BT_h6_an73_quadratic.py", []),
    ("geometric-unity/BU_h6_h2_prolongation.py", []),
    ("geometric-unity/BV_h6_an_positive.py", []),
    ("geometric-unity/BX_h4_kernels.py", []),
    ("geometric-unity/BY_data/verify_h4_low.py", ["--exhaustive"]),
    ("geometric-unity/BX_data/verify_h4_lift.py", ["--exhaustive"]),
    ("geometric-unity/CF_h6_fourth_jet.py", []),
    ("geometric-unity/CB_h4_partial_value.py", ["--full"]),
    ("geometric-unity/CD_clifford_hook_identities.py", []),
    ("geometric-unity/CG_hook_no_separator.py", ["--full-lie", "--all-coefficients"]),
    ("geometric-unity/CH_h4_second_order_stress.py", []),
    ("geometric-unity/CI_h4_cartan_support.py", ["--full-independent"]),
    ("geometric-unity/CJ_h6_stress_gradient.py", []),
    ("geometric-unity/CK_h4_invariant_forms.py", []),
    ("geometric-unity/CL_h6_oscillating_carriers.py", ["--full"]),
    ("geometric-unity/CO_h6_twisted_receiver.py", []),
    ("geometric-unity/CN_h4_s5_cubic.py", []),
    ("geometric-unity/CM_h6_sixth_jet.py", []),
    ("geometric-unity/CU_h6_local_field.py", []),
    ("geometric-unity/CU_h6_global_field.py", []),
    ("geometric-unity/CT_h6_potential_obstruction.py", []),
    ("geometric-unity/CR_h6_orbit_closure.py", []),
    ("geometric-unity/CS_h6_highpass.py", []),
    ("geometric-unity/CV_h6_axial.py", []),
    ("geometric-unity/CZ_h6_axial_local.py", []),
    ("geometric-unity/DA_h6_mixed_growth.py", []),
    ("geometric-unity/DB_h6_bounded_average.py", []),
    ("geometric-unity/CW_h4_s5_moment.py", []),
    ("geometric-unity/CX_sparse_hook_support.py", []),
    ("geometric-unity/CX_h4_minimum_external_rank.py", []),
    ("geometric-unity/DE_h6_elliptic_periodic.py", []),
    ("geometric-unity/DF_h6_l2_vanishing.py", []),
    ("geometric-unity/DC_h4_rank4_moment.py", []),
    ("geometric-unity/DS_h6_axial_multisign.py", ["--character", "15", "--geometry-only"]),
    ("geometric-unity/DI_h6_axial_gaussian_lifts.py", []),
    ("geometric-unity/DH_h6_bundle_gap.py", []),
    ("geometric-unity/DG_h4_rank3_moment.py", []),
    ("geometric-unity/DG_h6_rank6_sharpness.py", []),
    ("geometric-unity/DJ_hook_bounds.py", []),
    ("geometric-unity/DK_h6_mixed_gradient.py", []),
    ("geometric-unity/DK_h6_characteristic_cover.py", []),
    ("geometric-unity/DK_h6_characteristic_global.py", []),
    ("geometric-unity/DN_high_cartan.py", []),
    ("geometric-unity/DM_h6_first_tableau.py", []),
    ("geometric-unity/DO_h6_stress_spectrum.py", []),
    ("geometric-unity/DS_h6_axial_multisign.py", ["--eight-signs", "--character", "6", "--geometry-only"]),
    ("geometric-unity/DR_h6_rank_two.py", []),
    ("geometric-unity/BB_weyl_direction.py", []),
    ("geometric-unity/BB_cubic_extension.py", []),
    ("geometric-unity/CC_diagonal_weyl.py", []),
    ("geometric-unity/DD_mixed_weyl.py", []),
    ("geometric-unity/EE_zero_electric_weyl.py", []),
    ("geometric-unity/FF_joint_weyl.py", []),
    ("geometric-unity/GG_ricci_zero_weyl.py", []),
    ("geometric-unity/HH_ricci_zero.py", []),
]
LEAN_SOURCE_SHA256 = "522844e911cf5f6091b603db0d1d042f50544d2ec8d10ed9336b722945c331da"


def check_alpha_source_bytes(path):
    """Preserve the raw-byte pin used by the recent alpha source readers."""
    lean_source = path / "papers/alpha/lean/ChainAlpha.lean"
    if not lean_source.is_file() or hashlib.sha256(lean_source.read_bytes()).hexdigest() != LEAN_SOURCE_SHA256:
        raise ValueError("recent alpha checks require the exact pinned LF source bytes; use a fresh clone with "
                         "'git clone --config core.autocrlf=false ...' and the recorded commit")


def source_commands(path, recent=False):
    compendium = path / "compendium"
    if not (compendium / "index/statements.tsv").is_file():
        raise ValueError("--upstream must point to a Chain Theory checkout")
    head = subprocess.run(["git", "-C", str(path), "rev-parse", "HEAD"],
                          capture_output=True, text=True, encoding="utf-8", timeout=30)
    if head.returncode or head.stdout.strip() != SOURCE_COMMIT:
        raise ValueError(f"--upstream must be checked out at {SOURCE_COMMIT}")
    dirty = subprocess.run(["git", "-C", str(path), "status", "--porcelain",
                            "--untracked-files=all", "--ignored"],
                           capture_output=True, text=True, encoding="utf-8", timeout=30)
    if dirty.returncode or dirty.stdout.strip():
        raise ValueError("--upstream must be pristine, with no tracked changes, untracked files, or ignored files; "
                         "use a fresh clone checked out at the recorded commit")
    prefix = ROOT / "computation/chain-theory/compendium"
    logs = prefix / "outputs"
    commands = [
        ("source census", [str(prefix / "census.py"), str(compendium / "index/statements.tsv")]),
        ("source premises", [str(prefix / "premises.py"), str(compendium)]),
        ("source hydrogen registers", [str(prefix / "hydrogen_registers.py"), str(compendium)]),
        ("Lean coverage using recorded logs", [
            str(prefix / "lean_tables.py"), str(compendium),
            "--axioms", str(logs / "lean_axioms.txt"), "--checks", str(logs / "lean_table_checks.txt"),
            "--round6-tsv", str(logs / "lean_resolve.tsv")]),
    ]
    if recent:
        check_alpha_source_bytes(path)
        alpha = ROOT / "computation/chain-theory/alpha"
        commands += [
            ("recent conjugation source anchors", [str(alpha / "conjugation.py"), "--upstream", str(path)]),
            ("recent record source anchors", [str(alpha / "records_uniform.py"), "--upstream", str(path)]),
        ]
    return commands


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--upstream", type=Path, help="optional pristine Chain Theory checkout (no extra files)")
    parser.add_argument("--recent", action="store_true", help="include selected recent research checks")
    parser.add_argument("--timeout", type=float, default=180, help="seconds per check (default: 180)")
    parser.add_argument("--list", action="store_true", help="list selected checks without running them")
    args = parser.parse_args()
    if args.timeout <= 0:
        parser.error("--timeout must be positive")
    checks = [("focused regressions", ["-m", "unittest", "discover", "-s", "tests", "-v"])]
    checks += [(item, [str(ROOT / "computation" / item)]) for item in SELECTED]
    if args.recent:
        checks += [(item, [str(ROOT / "computation" / item), *options]) for item, options in RECENT]
    if args.upstream:
        try:
            checks += source_commands(args.upstream.resolve(), recent=args.recent)
        except (ValueError, OSError, subprocess.TimeoutExpired) as error:
            parser.error(str(error))
    if args.list:
        for label, _ in checks:
            print(label)
        return 0
    missing = []
    for package, expected in DEPENDENCIES.items():
        try:
            actual = importlib.metadata.version(package)
        except importlib.metadata.PackageNotFoundError:
            actual = "missing"
        if actual != expected:
            missing.append(f"{package}: expected {expected}, found {actual}")
    if missing:
        print("\n".join(missing), file=sys.stderr)
        print("Install dependencies: python -m pip install -r requirements.txt", file=sys.stderr)
        return 2
    environment = dict(os.environ, PYTHONHASHSEED="0", PYTHONIOENCODING="utf-8")
    failed = 0
    for label, arguments in checks:
        started = time.monotonic()
        try:
            result = subprocess.run([sys.executable, "-X", "utf8", *arguments],
                                    cwd=ROOT, env=environment, capture_output=True,
                                    text=True, encoding="utf-8", errors="replace", timeout=args.timeout)
        except (OSError, subprocess.TimeoutExpired) as error:
            print(f"FAIL {label}: {error}", flush=True)
            failed += 1
            continue
        seconds = time.monotonic() - started
        if result.returncode:
            failed += 1
            print(f"FAIL {label} ({seconds:.1f}s; exit {result.returncode})", flush=True)
            print(result.stdout)
            print(result.stderr)
        else:
            print(f"PASS {label} ({seconds:.1f}s)", flush=True)
    print(f"\n{len(checks) - failed}/{len(checks)} selected checks passed.")
    print("Long searches, optional checks, manuscript compilation, and fresh Lean builds are outside this suite.")
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
