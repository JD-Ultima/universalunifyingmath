#!/usr/bin/env python3
"""Build the compendium's Lean libraries with a given core Lean binary, in a fresh scratch directory.

Usage:
  python3 -I lean_build.py --lean LEAN --src SRC/compendium/lean --work WORK
                           [--targets ChainTheory AlphaDerivation ...] [--jobs 3] [--mem-mb 12000]
                           [--record outputs/lean_build_record.tsv]

Requires Linux (including WSL): child processes are reaped with os.wait4, and
ru_maxrss is interpreted using Linux's KiB units.

What it does (nothing is executed from inside SRC; no lake, no lakefile, no script of the source is run):
  1. Discovers the library modules of the six packages from their `lean_lib` names (lakefile.toml read as
     text): module M of library N in package P is the file P/N.lean or P/N/.../M.lean.
  2. Takes the import closure of the targets (default: the `ChainTheory` and `AlphaDerivation` libraries,
     i.e. the 493 modules of the source's `lake build`; `--targets all` adds `ActClosureCapFree`).
  3. Safety gate: refuses to build if any module of the closure contains a command that runs code at
     elaboration time or touches the system (#eval, run_cmd/run_elab/run_meta, initialize, IO., import Lean,
     macro/elab/syntax definitions, @[extern], implemented_by, unsafe, native_decide, #exit).
  4. Copies each source file into WORK/src (fresh; WORK must not exist) and records its sha256.
  5. Compiles in dependency order:  lean -R WORK/src -o WORK/out/M.olean -i WORK/out/M.ilean WORK/src/M.lean
     with LEAN_PATH=WORK/out, cwd WORK, a pool of --jobs processes and a memory budget (per-module estimate =
     the source's own build_record.tsv peak, used for scheduling only).
  6. Writes a per-module record (seconds, peak RSS, exit code, number of warnings, `sorry` warnings,
     error text) and prints a summary.

  --patch FILE (optional): a TSV (module, line, old, new) of compatibility edits applied to the COPIES before
  compiling. Each `old` must occur exactly once on that line of the original file; the patched modules are
  listed in the record. Without --patch the sources are compiled exactly as published.
"""
import argparse, csv, hashlib, os, pathlib, re, shutil, subprocess, sys, time

LIBS = [  # (package directory, lean_lib name); checked against the lakefiles below
    ("AxiomDeclaration", "AxiomDeclaration"), ("ActClosure", "ActClosure"),
    ("AxiomDeclarationTwo", "AxiomDeclarationTwo"), ("ActClosureCapFree", "ActClosureCapFree"),
    ("KernelDerivationC", "KernelDerivationC"), ("ChainTheory", "AlphaDerivation"), ("ChainTheory", "ChainTheory"),
]
FORBIDDEN = [r"#eval", r"\brun_cmd\b", r"\brun_elab\b", r"\brun_meta\b", r"\binitialize\b", r"\bIO\.",
             r"^\s*import\s+Lean\b", r"^\s*(?:local\s+|scoped\s+)?(?:macro|elab|syntax|macro_rules|elab_rules)\b",
             r"@\[\s*extern", r"implemented_by", r"\bunsafe\b", r"native_decide", r"#exit"]


def discover(src):
    libs_in_lakefiles = set()
    for pkg in {p for p, _ in LIBS}:
        txt = (src / pkg / "lakefile.toml").read_text(encoding="utf-8")
        for m in re.finditer(r'\[\[lean_lib\]\]\s*\nname\s*=\s*"([^"]+)"', txt):
            libs_in_lakefiles.add((pkg, m.group(1)))
    assert libs_in_lakefiles == set(LIBS), (sorted(libs_in_lakefiles), LIBS)
    mods = {}
    for pkg, lib in LIBS:
        d = src / pkg
        mods[lib] = d / (lib + ".lean")
        for p in sorted((d / lib).rglob("*.lean")):
            m = ".".join(p.relative_to(d).with_suffix("").parts)
            assert m not in mods, "duplicate module " + m
            mods[m] = p
    imports = {m: re.findall(r"^import\s+(\S+)", p.read_text(encoding="utf-8"), re.M) for m, p in mods.items()}
    return mods, imports


def closure(roots, imports):
    seen, stack = set(), list(roots)
    while stack:
        x = stack.pop()
        if x in seen:
            continue
        if x not in imports:
            raise SystemExit(f"unresolved import {x}")
        seen.add(x)
        stack.extend(imports[x])
    return seen


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--lean", required=True)
    ap.add_argument("--src", required=True)
    ap.add_argument("--work", required=True)
    ap.add_argument("--targets", nargs="*", default=["ChainTheory", "AlphaDerivation"])
    ap.add_argument("--jobs", type=int, default=3)
    ap.add_argument("--mem-mb", type=int, default=12000)
    ap.add_argument("--mem-scale", type=float, default=1.35,
                    help="factor on the source's recorded peak memory (Lean 4.35.0 peaks were up to 1.3x higher)")
    ap.add_argument("--record", default=None)
    ap.add_argument("--patch", default=None)
    ap.add_argument("--resume", action="store_true",
                    help="continue in an existing WORK: modules recorded in WORK/status.tsv with exit 0 are kept")
    a = ap.parse_args()
    if sys.platform != "linux" or not hasattr(os, "wait4"):
        raise SystemExit(
            "lean_build.py requires Linux/WSL (os.wait4 and Linux peak-RSS units). "
            "Use Linux Python and a Linux Lean binary; native Windows is unsupported."
        )
    lean = os.path.abspath(a.lean)
    src = pathlib.Path(a.src).resolve()
    work = pathlib.Path(a.work).resolve()
    ver = subprocess.run([lean, "--version"], capture_output=True, text=True).stdout.strip()
    print("lean:", ver)
    mods, imports = discover(src)
    targets = [t for t, _ in [(l, 0) for _, l in LIBS]] if a.targets == ["all"] else a.targets
    need = closure(targets, imports)
    print(f"library modules discovered: {len(mods)}; in the import closure of {targets}: {len(need)}")
    # safety gate
    bad = []
    for m in sorted(need):
        txt = mods[m].read_text(encoding="utf-8")
        for pat in FORBIDDEN:
            if re.search(pat, txt, re.M):
                bad.append((m, pat))
    if bad:
        for m, pat in bad[:20]:
            print("REFUSED", m, pat)
        raise SystemExit("safety gate: forbidden construct in a module to be built")
    print("safety gate: no forbidden construct in the", len(need), "modules")
    if work.exists() and not a.resume:
        raise SystemExit(f"{work} exists; give a fresh directory (or --resume)")
    fresh = not work.exists()
    (work / "src").mkdir(parents=True, exist_ok=True)
    (work / "out").mkdir(exist_ok=True)
    shas = {m: hashlib.sha256(mods[m].read_bytes()).hexdigest() for m in need}
    if fresh:
        for m in need:
            dst = work / "src" / (m.replace(".", "/") + ".lean")
            dst.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(mods[m], dst)
    agg = hashlib.sha256("".join(f"{m} {shas[m]}\n" for m in sorted(shas)).encode()).hexdigest()
    print("sources " + ("copied" if fresh else "already copied (resume)") + "; sha256 over sorted 'module sha256' lines:", agg)
    patched = {}
    if a.patch and not fresh:  # resume: apply the rows not yet applied to the copies
        for r in csv.DictReader(open(a.patch, encoding="utf-8"), delimiter="\t"):
            m, ln = r["module"], int(r["line"])
            if m not in need:
                continue
            dst = work / "src" / (m.replace(".", "/") + ".lean")
            lines = dst.read_text(encoding="utf-8").split("\n")
            orig = mods[m].read_text(encoding="utf-8").split("\n")
            if lines[ln - 1] == orig[ln - 1]:
                if lines[ln - 1].count(r["old"]) != 1:
                    raise SystemExit(f"patch does not apply: {m}:{ln}")
                lines[ln - 1] = lines[ln - 1].replace(r["old"], r["new"])
                dst.write_text("\n".join(lines), encoding="utf-8")
                print(f"PATCHED (resume) {m}:{ln}: {r['old']!r} -> {r['new']!r}")
            elif lines[ln - 1] != orig[ln - 1].replace(r["old"], r["new"]):
                raise SystemExit(f"copy of {m}:{ln} is neither the source nor the patched line")
            patched.setdefault(m, []).append(ln)
        # a module whose copy changed since it was built must be rebuilt
        st = work / "status.tsv"
        if st.exists():
            keep = [l for l in st.read_text(encoding="utf-8").splitlines() if l.split("\t")[0] not in patched or l.split("\t")[3] == "0"]
            st.write_text("".join(l + "\n" for l in keep), encoding="utf-8")
    if a.patch and fresh:
        for r in csv.DictReader(open(a.patch, encoding="utf-8"), delimiter="\t"):
            m, ln = r["module"], int(r["line"])
            if m not in need:
                continue
            dst = work / "src" / (m.replace(".", "/") + ".lean")
            lines = dst.read_text(encoding="utf-8").split("\n")
            if lines[ln - 1].count(r["old"]) != 1:
                raise SystemExit(f"patch does not apply: {m}:{ln}")
            lines[ln - 1] = lines[ln - 1].replace(r["old"], r["new"])
            dst.write_text("\n".join(lines), encoding="utf-8")
            patched.setdefault(m, []).append(ln)
            print(f"PATCHED {m}:{ln}: {r['old']!r} -> {r['new']!r}")
    # on resume, every copy must equal its source except the patched modules
    for m in need:
        if m not in patched:
            dst = work / "src" / (m.replace(".", "/") + ".lean")
            if hashlib.sha256(dst.read_bytes()).hexdigest() != shas[m]:
                raise SystemExit(f"copy of {m} differs from the source")
    # scheduling estimate from the source's own record (scheduling only)
    est = {}
    rec_path = src / "build_record.tsv"
    if rec_path.exists():
        for r in csv.DictReader(open(rec_path, encoding="utf-8"), delimiter="\t"):
            est[r["module"]] = float(r["peak_memory_MB"]) * a.mem_scale
    secs_est = {}
    if rec_path.exists():
        for r in csv.DictReader(open(rec_path, encoding="utf-8"), delimiter="\t"):
            secs_est[r["module"]] = float(r["seconds"])
    users = {m: [u for u in need if m in imports[u]] for m in need}
    tail = {}
    def chain(m):  # longest remaining chain of estimated seconds from m to a target (scheduling priority)
        if m not in tail:
            tail[m] = secs_est.get(m, 5.0) + max([chain(u) for u in users[m]] or [0.0])
        return tail[m]
    sys.setrecursionlimit(10000)
    for m in need:
        chain(m)
    env = dict(os.environ, LEAN_PATH=str(work / "out"))
    done, failed, running = {}, set(), {}
    status_path = work / "status.tsv"
    if a.resume and status_path.exists():
        for line in status_path.read_text(encoding="utf-8").splitlines():
            f = line.split("\t")
            if len(f) >= 6 and f[0] in need and int(f[3]) == 0 and (work / "out" / (f[0].replace(".", "/") + ".olean")).exists():
                done[f[0]] = (float(f[1]), int(f[2]), 0, int(f[4]), int(f[5]), "")
        print(f"resume: {len(done)} modules already built (status.tsv)")
    status_fh = open(status_path, "a", encoding="utf-8")
    pending = set(need) - set(done)
    t0 = time.time()
    order = sorted(need, key=lambda m: -tail[m])  # critical path first
    while pending or running:
        # launch
        launched = True
        while launched:
            launched = False
            budget = sum(est.get(m, 1500) for m in running.values())
            for m in order:
                if m not in pending:
                    continue
                if any(i in failed for i in imports[m]):
                    pending.discard(m); failed.add(m); done[m] = (0.0, 0, -1, 0, 0, "dependency failed")
                    launched = True
                    break
                if not all(i in done and i not in failed for i in imports[m]):
                    continue
                if len(running) >= a.jobs or (running and budget + est.get(m, 1500) > a.mem_mb):
                    continue
                f = work / "src" / (m.replace(".", "/") + ".lean")
                o = work / "out" / (m.replace(".", "/") + ".olean")
                o.parent.mkdir(parents=True, exist_ok=True)
                logf = open(work / (m + ".log"), "w")
                p = subprocess.Popen([lean, "-R", str(work / "src"), "-o", str(o), "-i", str(o.with_suffix(".ilean")), str(f)],
                                     cwd=work, env=env, stdout=logf, stderr=subprocess.STDOUT)
                running[p.pid] = m
                p._t, p._log = time.time(), logf
                running_procs[p.pid] = p
                pending.discard(m)
                launched = True
                break
        if not running:
            if pending:
                raise SystemExit(f"stuck with {len(pending)} pending modules")
            break
        pid, status, ru = os.wait4(-1, 0)
        if pid not in running:
            continue
        m = running.pop(pid)
        p = running_procs.pop(pid)
        p._log.close()
        secs = time.time() - p._t
        code = os.waitstatus_to_exitcode(status)
        log = (work / (m + ".log")).read_text(encoding="utf-8", errors="replace")
        nwarn = len(re.findall(r": warning:", log))
        nsorry = len(re.findall(r"declaration uses 'sorry'", log))
        err = " | ".join(l for l in log.splitlines() if ": error:" in l).replace(str(work) + "/", "")[:400]
        done[m] = (secs, ru.ru_maxrss // 1024, code, nwarn, nsorry, err)
        status_fh.write(f"{m}\t{secs:.1f}\t{ru.ru_maxrss // 1024}\t{code}\t{nwarn}\t{nsorry}\t{err}\n")
        status_fh.flush()
        if code != 0:
            failed.add(m)
        print(f"[{len(done):>3}/{len(need)}] {m}  {secs:.1f}s  {ru.ru_maxrss // 1024} MB  exit {code}"
              + (f"  warnings {nwarn}" if nwarn else "") + (f"  ERROR {err}" if err else ""), flush=True)
    wall = time.time() - t0
    ok = [m for m in done if done[m][2] == 0]
    print(f"\nbuilt {len(ok)} of {len(need)} modules; failed {len(failed)}; wall {wall:.0f} s; "
          f"sum of module times {sum(v[0] for v in done.values()):.0f} s; max peak {max(v[1] for v in done.values())} MB; "
          f"warnings {sum(v[3] for v in done.values())}; sorry warnings {sum(v[4] for v in done.values())}")
    if a.record:
        with open(a.record, "w", newline="", encoding="utf-8") as fh:
            w = csv.writer(fh, delimiter="\t", lineterminator="\n")
            w.writerow(["module", "seconds", "peak_rss_MB", "exit", "warnings", "sorry_warnings", "sha256_original",
                        "patched_lines", "errors"])
            for m in sorted(done):
                s, mb, c, nw, ns, e = done[m]
                e = e.replace(str(work) + "/", "")
                w.writerow([m, f"{s:.1f}", mb, c, nw, ns, shas[m], ",".join(map(str, patched.get(m, []))), e])
    return 0 if not failed else 1


running_procs = {}
if __name__ == "__main__":
    sys.exit(main())
