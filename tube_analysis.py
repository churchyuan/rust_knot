import argparse
import csv
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, List, Optional, Sequence, Tuple

import numpy as np

import ropelength_slit


@dataclass
class Frame:
    coords: np.ndarray
    source: str


def _as_array(pts: Sequence[Tuple[float, float, float]]) -> np.ndarray:
    a = np.asarray(pts, dtype=float)
    if a.ndim != 2 or a.shape[1] != 3:
        raise ValueError(f"Expected Nx3 points, got shape {a.shape}")
    return a


def _iter_xyz_frames(path: Path) -> Iterable[np.ndarray]:
    text = path.read_text(encoding="utf-8").splitlines()
    lines = [ln.strip() for ln in text if ln.strip()]
    if not lines:
        raise ValueError(f"Empty file: {path}")

    def _try_int(s: str) -> Optional[int]:
        try:
            return int(s)
        except ValueError:
            return None

    n0 = _try_int(lines[0])
    if n0 is None or n0 <= 0:
        yield _as_array(ropelength_slit.read_xyz_like(path))
        return

    i = 0
    while i < len(lines):
        n = _try_int(lines[i])
        if n is None or n <= 0:
            raise ValueError(f"Unrecognized XYZ frame header at line {i+1} in {path}")
        i += 1
        if i >= len(lines):
            raise ValueError(f"Missing XYZ comment line in {path}")
        i += 1
        if i + n > len(lines):
            raise ValueError(f"Truncated XYZ frame in {path}")
        pts: List[Tuple[float, float, float]] = []
        for j in range(n):
            parts = lines[i + j].split()
            if len(parts) >= 4:
                x, y, z = float(parts[-3]), float(parts[-2]), float(parts[-1])
            elif len(parts) == 3:
                x, y, z = float(parts[0]), float(parts[1]), float(parts[2])
            else:
                raise ValueError(f"Bad XYZ line in {path}: {lines[i+j]!r}")
            pts.append((x, y, z))
        i += n
        yield _as_array(pts)


def _iter_lammps_dump_frames(path: Path) -> Iterable[np.ndarray]:
    frames: List[List[Tuple[float, float, float]]] = []
    current: List[Tuple[float, float, float]] = []
    expect_atoms = False
    x_col = y_col = z_col = None
    with path.open("r", encoding="utf-8") as f:
        for raw in f:
            line = raw.strip()
            if not line:
                continue
            if line.startswith("ITEM:"):
                if expect_atoms:
                    frames.append(current)
                    current = []
                    expect_atoms = False
                if line.startswith("ITEM: ATOMS"):
                    parts = line.split()[2:]
                    col = {name: idx for idx, name in enumerate(parts)}
                    for key in ("x", "y", "z"):
                        if key not in col:
                            raise ValueError(f"Dump {path} missing column {key!r} in: {line!r}")
                    x_col, y_col, z_col = col["x"], col["y"], col["z"]
                    expect_atoms = True
                continue

            if expect_atoms:
                parts = line.split()
                if x_col is None or y_col is None or z_col is None:
                    raise RuntimeError("Internal error: dump columns not set")
                x = float(parts[x_col])
                y = float(parts[y_col])
                z = float(parts[z_col])
                current.append((x, y, z))

    if expect_atoms and current:
        frames.append(current)

    if not frames:
        raise ValueError(f"No frames found in dump: {path}")

    for pts in frames:
        yield _as_array(pts)


def _kabsch_align(mobile: np.ndarray, target: np.ndarray) -> Tuple[np.ndarray, float, np.ndarray]:
    if mobile.shape != target.shape:
        raise ValueError(f"Shape mismatch: mobile {mobile.shape}, target {target.shape}")
    if mobile.shape[1] != 3:
        raise ValueError(f"Expected Nx3, got {mobile.shape}")

    m0 = mobile.mean(axis=0)
    t0 = target.mean(axis=0)
    M = mobile - m0
    T = target - t0
    cov = M.T @ T
    U, S, Vt = np.linalg.svd(cov)
    R = U @ Vt
    if np.linalg.det(R) < 0.0:
        U[:, -1] *= -1.0
        R = U @ Vt
    aligned = (mobile - m0) @ R + t0
    diff = aligned - target
    rmsd = float(math.sqrt((diff * diff).sum(axis=1).mean()))
    return aligned, rmsd, R


@dataclass
class AlignResult:
    aligned: np.ndarray
    rmsd: float
    shift: int
    reversed: bool


def _best_cyclic_align(mobile: np.ndarray, target: np.ndarray, allow_reverse: bool) -> AlignResult:
    n = mobile.shape[0]
    if n != target.shape[0]:
        raise ValueError(f"Point count mismatch: mobile {n}, target {target.shape[0]}")

    best: Optional[AlignResult] = None

    def _try(m: np.ndarray, shift: int, rev: bool) -> None:
        nonlocal best
        aligned, rmsd, _ = _kabsch_align(m, target)
        cand = AlignResult(aligned=aligned, rmsd=rmsd, shift=shift, reversed=rev)
        if best is None or cand.rmsd < best.rmsd:
            best = cand

    for k in range(n):
        _try(np.roll(mobile, shift=k, axis=0), shift=k, rev=False)
    if allow_reverse:
        rev0 = mobile[::-1].copy()
        for k in range(n):
            _try(np.roll(rev0, shift=k, axis=0), shift=k, rev=True)

    if best is None:
        raise RuntimeError("Failed to align")
    return best


def _align(mobile: np.ndarray, target: np.ndarray, cyclic: bool, allow_reverse: bool) -> AlignResult:
    if cyclic:
        return _best_cyclic_align(mobile, target, allow_reverse=allow_reverse)
    aligned, rmsd, _ = _kabsch_align(mobile, target)
    return AlignResult(aligned=aligned, rmsd=rmsd, shift=0, reversed=False)


def _load_frames(
    inputs: Sequence[Path],
    globs: Sequence[str],
    dump: Optional[Path],
) -> List[Frame]:
    paths: List[Path] = []
    for p in inputs:
        paths.append(p)
    for pat in globs:
        for p in sorted(Path().glob(pat)):
            if p.is_file():
                paths.append(p)
    if dump is not None:
        paths.append(dump)

    if not paths:
        raise ValueError("No input frames provided. Use --input/--glob/--dump.")

    frames: List[Frame] = []
    for p in paths:
        if p.suffix.lower() in (".lammpstrj", ".dump", ".traj"):
            for idx, a in enumerate(_iter_lammps_dump_frames(p)):
                frames.append(Frame(coords=a, source=f"{p}:{idx}"))
        else:
            multi = list(_iter_xyz_frames(p))
            if len(multi) == 1:
                frames.append(Frame(coords=multi[0], source=str(p)))
            else:
                for idx, a in enumerate(multi):
                    frames.append(Frame(coords=a, source=f"{p}:{idx}"))

    n = frames[0].coords.shape[0]
    for fr in frames:
        if fr.coords.shape[0] != n:
            raise ValueError(
                "All frames must have the same point count. "
                f"Expected {n}, got {fr.coords.shape[0]} in {fr.source}"
            )
    return frames


def _compute_average(
    frames: Sequence[Frame],
    reference: np.ndarray,
    cyclic: bool,
    allow_reverse: bool,
    max_iter: int,
    tol: float,
) -> Tuple[np.ndarray, List[AlignResult]]:
    avg = reference.copy()
    last_rmsds: List[AlignResult] = []
    for _ in range(max_iter):
        aligned_results = [_align(fr.coords, avg, cyclic=cyclic, allow_reverse=allow_reverse) for fr in frames]
        aligned = np.stack([r.aligned for r in aligned_results], axis=0)
        new_avg = aligned.mean(axis=0)
        delta = float(np.sqrt(((new_avg - avg) ** 2).sum(axis=1).mean()))
        avg = new_avg
        last_rmsds = aligned_results
        if delta < tol:
            break
    final_results = [_align(fr.coords, avg, cyclic=cyclic, allow_reverse=allow_reverse) for fr in frames]
    return avg, final_results


def _write_outputs(
    outdir: Path,
    avg: np.ndarray,
    aligned_results: Sequence[AlignResult],
    frames: Sequence[Frame],
) -> None:
    outdir.mkdir(parents=True, exist_ok=True)

    ropelength_slit.write_xyz(outdir / "average.xyz", [(float(x), float(y), float(z)) for x, y, z in avg])

    rmsd2 = np.array([r.rmsd * r.rmsd for r in aligned_results], dtype=float)
    tube_rms = float(math.sqrt(rmsd2.mean())) if len(rmsd2) else float("nan")

    aligned = np.stack([r.aligned for r in aligned_results], axis=0)
    msd_atom = ((aligned - avg[None, :, :]) ** 2).sum(axis=2).mean(axis=0)
    tube_atom_rms = float(math.sqrt(float(msd_atom.mean())))
    tube_atom_max = float(math.sqrt(float(msd_atom.max())))

    with (outdir / "summary.txt").open("w", encoding="utf-8", newline="\n") as f:
        f.write(f"frames {len(frames)}\n")
        f.write(f"points {avg.shape[0]}\n")
        f.write(f"tube_rmsd_mean {tube_rms:.10g}\n")
        f.write(f"tube_atom_rms {tube_atom_rms:.10g}\n")
        f.write(f"tube_atom_max {tube_atom_max:.10g}\n")

    with (outdir / "rmsd.csv").open("w", encoding="utf-8", newline="") as f:
        w = csv.writer(f)
        w.writerow(["frame", "source", "rmsd", "shift", "reversed"])
        for i, (fr, r) in enumerate(zip(frames, aligned_results)):
            w.writerow([i, fr.source, f"{r.rmsd:.10g}", r.shift, int(r.reversed)])

    with (outdir / "msd_atom.csv").open("w", encoding="utf-8", newline="") as f:
        w = csv.writer(f)
        w.writerow(["atom_index", "msd", "rms"])
        for i, msd in enumerate(msd_atom.tolist()):
            w.writerow([i, f"{msd:.10g}", f"{math.sqrt(msd):.10g}"])


def main(argv: Optional[Sequence[str]] = None) -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", type=Path, action="append", default=[])
    ap.add_argument("--glob", type=str, action="append", default=[])
    ap.add_argument("--dump", type=Path, default=None)
    ap.add_argument("--reference", type=Path, default=None)
    ap.add_argument("--outdir", type=Path, default=Path("tube_out"))
    ap.add_argument("--no-cyclic", action="store_true")
    ap.add_argument("--no-reverse", action="store_true")
    ap.add_argument("--max-iter", type=int, default=10)
    ap.add_argument("--tol", type=float, default=1e-8)
    args = ap.parse_args(argv)

    frames = _load_frames(inputs=args.input, globs=args.glob, dump=args.dump)
    cyclic = not args.no_cyclic
    allow_reverse = not args.no_reverse

    if args.reference is None:
        reference = frames[0].coords
    else:
        reference = next(_iter_xyz_frames(args.reference))
        if reference.shape[0] != frames[0].coords.shape[0]:
            raise ValueError(
                f"Reference point count mismatch: ref {reference.shape[0]}, frames {frames[0].coords.shape[0]}"
            )

    avg, aligned_results = _compute_average(
        frames=frames,
        reference=reference,
        cyclic=cyclic,
        allow_reverse=allow_reverse,
        max_iter=max(1, args.max_iter),
        tol=max(0.0, args.tol),
    )
    _write_outputs(args.outdir, avg=avg, aligned_results=aligned_results, frames=frames)


if __name__ == "__main__":
    main()
