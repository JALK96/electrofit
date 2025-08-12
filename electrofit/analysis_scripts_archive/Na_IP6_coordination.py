#!/usr/bin/env python3
"""
analyze_na_ip6_coordination.py

Loop over all micro-state folders in …/process/, build an
(N_Na, 6, N_frames) Boolean coordination tensor telling whether each
Na⁺ is closer than RCUTOFF to any *peripheral* oxygen of P1, P2, …, P6.
Afterward, aggregate to time-series of counts and plot them.

Author: Arthur Laux
Date  : 2025-07-30
"""

from __future__ import annotations
import os, sys, logging, pathlib, mmap
from typing import Dict, List, Tuple

import numpy as np
import MDAnalysis as mda
from MDAnalysis.lib.distances import distance_array         #  (https://userguide.mdanalysis.org/examples/analysis/distances_and_contacts/distances_between_selections.html)
from MDAnalysis.analysis.rdf import InterRDF

import matplotlib.pyplot as plt
import seaborn as sns
sns.set_context("talk")

# ╔═══════════════════════════════════════════════════════════════╗
# ┃                         CONSTANTS                             ┃
# ╚═══════════════════════════════════════════════════════════════╝
PROJECT_NAME       = "electrofit"
PROCESS_DIR_NAME   = "data/results/final/process_IP6_11ps"
RUN_DIR_NAME       = "run_final_gmx_simulation"
PLOT_DIR_NAME      = "analyze_final_sim/NaP_coordination"
PLOT_NAME_TEMPLATE = "NaP_coordination_counts.png"

RCUTOFF_NM         = 0.32       # coordination threshold  (change if desired)
SAVE_3D_BOOL       = True       # turn off if disk space is an issue
BOOL_FILENAME      = "NaP_coordination_bool.npy"   # stored as numpy.memmap

# Mapping: phosphate → peripheral O atom names (no ester O)
PHOS_OXYGENS: Dict[str, Tuple[str, ...]] = {
    "P" : ("O6",  "O7",  "O8"),
    "P1": ("O9",  "O10", "O11"),
    "P2": ("O12", "O13", "O14"),
    "P3": ("O15", "O16", "O17"),
    "P4": ("O18", "O19", "O20"),
    "P5": ("O21", "O22", "O23"),
}
PHOS_LABELS = tuple(PHOS_OXYGENS.keys())   # ('P', 'P1', …)
PHOS_LABELS_CORRECTED = ('P1', 'P2', 'P3', 'P4', 'P5', 'P6')  # re-order to match the paper

# Set logging level
LOGLEVEL = logging.INFO

# ── Utility: first‑shell boundary from an RDF curve ──────────────────────────
def first_shell_end(r: np.ndarray, g: np.ndarray) -> float:
    """
    Return the radius r (nm) at which the first coordination shell ends,
    i.e. the first *local* minimum of g(r) that follows the first *local*
    maximum (peak).  The algorithm is deliberately conservative: if no
    clear peak‑minimum pair can be detected within 0 < r < 0.8 nm it
    returns NaN so that the calling code can ignore the value.

    Notes
    -----
    *Local* extrema are detected by finite‑difference sign changes:
      peak   : g[i-1] < g[i] >= g[i+1]
      minimum: g[i-1] > g[i] <= g[i+1]

    The search is limited to r ≤ 0.8 nm, which safely brackets the
    Na–O first‑shell region observed in all micro‑states.
    """
    if r.ndim != 1 or g.ndim != 1:
        raise ValueError("r and g must be one‑dimensional arrays")
    if len(r) != len(g):
        raise ValueError("r and g must have the same length")

    N = len(r)
    rlim = 0.8
    # Only consider indices where 0 < r[i] <= 0.8 and 1 <= i < N-1
    valid = [i for i in range(1, N-1) if (r[i] > 0) and (r[i] <= rlim)]

    # Find the first local maximum (peak)
    peak_idx = None
    for i in valid:
        if g[i] > g[i-1] and g[i] >= g[i+1]:
            peak_idx = i
            break
    if peak_idx is None:
        return np.nan

    # Find the first local minimum after the peak
    min_idx = None
    for j in range(peak_idx + 1, N-1):
        if r[j] > rlim:
            break
        if g[j] < g[j-1] and g[j] <= g[j+1]:
            min_idx = j
            break
    if min_idx is None:
        return np.nan

    # Sanity: ensure minimum is after peak and not at the very end
    if min_idx <= peak_idx or min_idx >= N-1:
        return np.nan
    return float(r[min_idx])
# ╔═══════════════════════════════════════════════════════════════╝


def find_project_root(start: pathlib.Path, name: str = PROJECT_NAME) -> pathlib.Path:
    """Walk up to the filesystem root, track all folders named `name`, and return the outermost one."""
    root_candidate = None
    here = start.resolve()

    while True:
        if here.name == name:
            root_candidate = here
        if here == here.parent:  # reached filesystem root
            break
        here = here.parent

    if root_candidate is not None:
        return root_candidate

    raise FileNotFoundError(
        f'Could not locate project root "{name}" starting from {start}'
    )

def setup_logging(logfile: pathlib.Path | None = None, level: int = LOGLEVEL) -> None:
    handlers = [logging.FileHandler(logfile)] if logfile else [logging.StreamHandler()]
    logging.basicConfig(
        handlers=handlers,
        level=level,
        format="%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%H:%M:%S",
    )

# ╔═══════════════════════════════════════════════════════════════╗
# ┃              CORE ANALYSIS ── ONE MICRO-STATE                ┃
# ╚═══════════════════════════════════════════════════════════════╝
def analyse_one_microstate(
        traj_file: pathlib.Path,
        top_file : pathlib.Path,
        dest_dir : pathlib.Path,
        r_cut_nm: float = RCUTOFF_NM,
) -> None:
    """
    For a single trajectory/topology pair:
      1. Build Boolean coordination tensor  A[na, phosph, frame].
      2. Sum over Na dimension → counts_per_frame[phosph, frame].
      3. Plot time-series.
    """

    # ── Load MD trajectory ───────────────────────────────────
    u = mda.Universe(top_file.as_posix(), traj_file.as_posix())
    na_atoms = u.select_atoms("name NA")
    n_na = len(na_atoms)
    if n_na == 0:
        logging.warning("No Na⁺ atoms found – skipped.")
        return
    logging.info("Loaded %d Na⁺ atoms and %d frames", n_na, len(u.trajectory))

    print("Unique atom names in IP6:")
    print(sorted(set(u.select_atoms("resname I*").names)))

    phos_ag = {p: u.select_atoms("resname I* and name " + " ".join(names))
               for p, names in PHOS_OXYGENS.items()}

    # Sanity-check selections
    for p, ag in phos_ag.items():
        if len(ag) != 3:
            raise ValueError(f"Expected 3 peripheral O atoms for {p}, got {len(ag)}")

    # ── Prepare Boolean memmap (optionally) ─────────────────
    n_frames   = len(u.trajectory)
    n_phos     = len(PHOS_LABELS)
    bool_shape = (n_na, n_phos, n_frames)          # (i, j, k)
    coord_mem  = None
    if SAVE_3D_BOOL:
        mmap_path = dest_dir / BOOL_FILENAME
        # open_memmap writes a valid .npy header; np.memmap does NOT
        coord_mem = np.lib.format.open_memmap(
            mmap_path, mode='w+', dtype=np.bool_, shape=bool_shape
        )
        coord_mem[:] = False      # initialise so file is fully allocated
        coord_mem.flush()
        logging.info(
            "Boolean tensor mapped to %s  (≈ %.1f MB).",
            mmap_path, coord_mem.nbytes / 1e6
        )

    # ── Frame-by-frame loop (vectorised in C under the hood) ──
    counts_ts = np.zeros((n_phos, n_frames), dtype=np.int16)  # for final plot

    for k, ts in enumerate(u.trajectory):
        na_pos = na_atoms.positions               # (N_Na, 3)
        min_dists = np.full((n_na, n_phos), np.inf, dtype=np.float32)

        # Compute per-phosphate minimum distances Naᵢ ↔ O_peripheral
        for j, p in enumerate(PHOS_LABELS):
            d = distance_array(na_pos, phos_ag[p].positions, box=ts.dimensions)  # (N_Na, 3)
            min_dists[:, j] = d.min(axis=1)   # (N_Na,)

        # For each Na, find closest phosphate (index of min along axis=1)
        closest_phos_idx = np.argmin(min_dists, axis=1)           # (N_Na,)
        closest_dist     = min_dists[np.arange(n_na), closest_phos_idx]

        # Boolean mask: coordinated if distance < cutoff
        coordin_mask = closest_dist < (r_cut_nm * 10)   # MDAnalysis uses Å
        # Note: *10 converts nm → Å, because Universe positions are in Å.

        # Update Boolean tensor and counts
        if SAVE_3D_BOOL:
            coord_mem[np.arange(n_na), closest_phos_idx, k] = coordin_mask

        # Aggregate counts per phosphate for this frame
        for j in range(n_phos):
            counts_ts[j, k] = np.count_nonzero(coordin_mask & (closest_phos_idx == j))

    # Flush memmap to disk
    if SAVE_3D_BOOL:
        coord_mem.flush()

    # ── Plot counts vs. time ─────────────────────────────────
    plot_counts_subplots(
        counts_ts=counts_ts,
        timestep_ps=u.trajectory.dt,
        out_png  = dest_dir / PLOT_NAME_TEMPLATE,
        title    = f"Na⁺ coordination counts – {dest_dir.parent.name}",
    )

    snapshot_png = dest_dir / f"network_frame_3d.png"
    plot_frame_network_3d_fixed_view(u, frame=521, out_png=snapshot_png, r_cut_nm=r_cut_nm)
    snapshot_png = dest_dir / f"network_frame.png"
    plot_frame_network_plane(
        u, frame=525, out_png=snapshot_png,
        reference_triplet=("C", "C2", "C4"), r_cut_nm=r_cut_nm
    )

    rdf_png = dest_dir / "rdf_Na_periphO.pdf"
    # global y max evaluated once (see outcommentd code below in main function)
    plot_rdf_periphO_Na(u, out_png=rdf_png, r_max=1.2, nbins=240, y_max_global=1800.0) 

    logging.info("Finished %s", dest_dir.parent.name)

# ╔═══════════════════════════════════════════════════════════════╗
# ┃                         PLOTTING                             ┃
# ╚═══════════════════════════════════════════════════════════════╝
def plot_counts_subplots(
        counts_ts : np.ndarray,          # shape (6, n_frames)
        timestep_ps: float,
        out_png   : pathlib.Path,
        title     : str,
) -> None:
    """Re-implementation of your 2×3 layout but plotting counts."""
    sns.set_context("talk")
    n_phos, n_frames = counts_ts.shape
    time_ns = np.arange(n_frames) * timestep_ps / 1000.0

    fig, axes = plt.subplots(2, 3, figsize=(8, 6), sharex=True, sharey=True)
    axes = axes.flatten()

    for idx, (ax, label) in enumerate(zip(axes, PHOS_LABELS)):
        y = counts_ts[idx]
        ax.plot(time_ns, y, lw=0.7, color="black")
        mean_y = y.mean()
        ax.axhline(mean_y, ls="--", lw=1.4, color="red")
        ax.text(0.95, 0.9, f"⟨count⟩ = {mean_y:.1f}",
                ha="right", va="top", transform=ax.transAxes, color="red", fontsize=11)
        ax.set_title(f"{label} – Na⁺")

        if idx % 3 == 0:
            ax.set_ylabel("No. bound Na⁺")
        if idx >= 3:
            ax.set_xlabel("Time / ns")

    fig.suptitle(title, y=1.02, fontsize=16)
    plt.tight_layout()
    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.close(fig)
    logging.info("Plot written to %s", out_png)

# ╔══════════════════════════════════════════════════════════╗
# ┃                 RDF  of Na⁺ vs. peripheral O             ┃
# ╚══════════════════════════════════════════════════════════╝

from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from typing import Optional
# -----------------------------------------------------------------------------
# MAIN PLOTTING FUNCTION -------------------------------------------------------
# -----------------------------------------------------------------------------


def plot_rdf_periphO_Na(
        u: mda.Universe,
        out_png: pathlib.Path,
        r_max: float = 1.2,      # nm
        nbins: int   = 240,      # → dr ≈ 0.005 nm
        *,
        y_max_global: Optional[float] = None,
) -> None:
    """Plot Na⁺–peripheral‑oxygen RDFs for **one** microstate.

    Six rows (P1–P6), shared x‑axis and *global* y‑axis that must be supplied
    via *y_max_global* so **all** microstates use exactly the same scale.

    Parameters
    ----------
    u : MDAnalysis.Universe
        Trajectory of the microstate to analyse.
    out_png : pathlib.Path
        Output filename.
    r_max : float, optional
        Upper bound of the RDF in **nm**.
    nbins : int, optional
        Number of histogram bins.
    y_max_global : float, keyword‑only, required
        *Upper* limit of the y‑axis for **all** figures.  Compute this once
        across the full data set and pass it to every call.  If *None*, an
        exception is raised so the user does not forget.
    """

    if y_max_global is None:
        raise ValueError("plot_rdf_periphO_Na requires *y_max_global* so all figures share the same scale.")

    # ---- compute RDFs -------------------------------------------------------
    na_ag = u.select_atoms("name NA")
    periph_dict = {
        label: u.select_atoms("resname I* and name " + " ".join(PHOS_OXYGENS[label]))
        for label in PHOS_LABELS
    }

    rdf_results = {}
    for label, oxy_ag in periph_dict.items():
        rdf = InterRDF(na_ag, oxy_ag, range=(0, r_max * 10), nbins=nbins)
        rdf.run()
        rdf_results[label] = (rdf.bins / 10.0, rdf.rdf)  # convert Å → nm

    # ---- figure & axes ------------------------------------------------------
    sns.set_context("talk")
    n = len(PHOS_LABELS)
    fig, axes = plt.subplots(n, 1, figsize=(6, 1.2 * n), sharex=True, sharey=True)

    # ---- binary code → colours ---------------------------------------------
    try:
        binary_code = out_png.parents[2].name
        if binary_code.startswith("IP_"):
            binary_code = binary_code[3:]
    except Exception:
        binary_code = "000000"
    if len(binary_code) != 6 or any(c not in "01" for c in binary_code):
        raise ValueError(f"Unexpected binary code '{binary_code}' in output path")
    colours = ["blue" if bit == "1" else "black" for bit in binary_code]

    # ---- constants for inset -----------------------------------------------
    X_LOWER, X_UPPER = 0, 1.2     # nm

    # ---- loop over P1–P6 ----------------------------------------------------
    for ax, (label, col) in zip(axes, zip(PHOS_LABELS, colours)):
        r, g = rdf_results[label]
        g_per_phos = g * 3                    # 3 O atoms → per phosphate

        # --- determine first‑shell boundary ---------------------------------
        shell_end = first_shell_end(r, g_per_phos)
        if not np.isnan(shell_end):
            # visual guide: green dotted line at first minimum
            ax.axvline(shell_end, color="green", ls=":", lw=1.0)
            ax.text(shell_end, 0.95 * y_max_global,
                    f"{shell_end:.2f} nm",
                    rotation=90, va="top", ha="right",
                    color="green", fontsize=8)
            logging.info("First‑shell end for %s: %.3f nm", label, shell_end)

        if col == "blue":
            ax.fill_between(r, g_per_phos, 0.0, color=col, alpha=0.3)
        ax.plot(r, g_per_phos, lw=1.2, color="black")
        ax.axvline(RCUTOFF_NM, ymax=0.8, color="black", ls="--", lw=1.0)

        # global y‑scale ------------------------------------------------------
        ax.set_ylim(0.0, y_max_global)

        # y‑axis label --------------------------------------------------------
        ax.set_ylabel(PHOS_LABELS_CORRECTED[PHOS_LABELS.index(label)],
                      rotation=0, labelpad=40, va="center")
        ax.yaxis.set_major_locator(plt.MaxNLocator(2))
        if ax is not axes[-1]:
            ax.tick_params(labelbottom=False)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

        # coordination number -------------------------------------------------
        mask_cn = r <= RCUTOFF_NM
        vol_nm3 = np.prod(u.dimensions[:3]) / 1000.0  # Å³ → nm³
        rho_Na  = len(na_ag) / vol_nm3
        coord_num = 4.0 * np.pi * rho_Na * np.trapz(g_per_phos[mask_cn] * r[mask_cn]**2, x=r[mask_cn])
        ax.text(0.03, 0.8, f"{coord_num:.2f}", transform=ax.transAxes,
                ha="left", va="top", fontsize=14, color=col, alpha=0.7)

        # --------  OPTIONAL INSET for "blue" (bit == 1), "black" (bit ==0) -----------------
        if col == "blue" or col == "black":
            # restrict data to the desired x-window
            mask_inset = (r >= X_LOWER) & (r <= X_UPPER)
            r_inset    = r[mask_inset]
            g_inset    = (g*3)[mask_inset]        # already per-phosphate

            # create inset: upper-right corner (loc='upper right')
            axins = inset_axes(ax,
                               width = 1.5,   # 35 % of parent
                               height = 0.5,
                               loc = "upper right",
                               borderpad = 1.0)

            # plot and format
            if col == "blue":
                axins.fill_between(r_inset, g_inset, 0,
                               color=col, alpha=0.4)
            axins.plot(r_inset, g_inset, lw=1.0, color="black")

            # set limits so curve fills inset tightly
            axins.set_xlim(X_LOWER, X_UPPER)
            y_min, y_max = g_inset.min(), g_inset.max()
            x_min, x_max = r_inset.min(), r_inset.max()
            axins.set_ylim(y_min*0.98, y_max*1.02)

            # --- show only the top (maximum) y tick ---------------------------
            axins.set_yticks([y_max])                 # one tick at ymax
            axins.set_yticklabels([f"{y_max:.0f}"], fontsize=12)   # optional: format label
            axins.tick_params(axis="y",               # keep the tick mark
                            which="both",
                            direction="out",
                            left=False,
                            right=False,
                            labelright=False)
            if ax is not axes[0]:
                axins.tick_params(axis="x",               # hide x ticks/labels
                                bottom=False,
                                labelbottom=False)
            else:
                axins.set_xticks([x_min, x_max])                 # one tick at ymax
                axins.set_xticklabels([f"{x_min:.2f}", f"{x_max:.2f}"], fontsize=12)   # optional: format label
                axins.tick_params(axis="x",               # show x ticks/labels
                                top=False,
                                labeltop=True,
                                bottom=False,
                                labelbottom=False)
        
        #ax.ticklabel_format(axis='y', style='sci', scilimits=(0,0))  # scientific notation for large y-values

    # ---- common labels & title ---------------------------------------------
    axes[-1].set_xlabel("r / nm")
    fig.suptitle(binary_code, fontsize=20)

    plt.tight_layout()
    fig.subplots_adjust(hspace=0.1)
    fig.savefig(out_png, dpi=300)
    plt.close(fig)
    logging.info("RDF figure written to %s", out_png)

# ╔═══════════════════════════════════════════════════════════════╗
# ┃         2-D PROJECTION  PERPENDICULAR  TO  A  REFERENCE PLANE ┃
# ╚═══════════════════════════════════════════════════════════════╝
def plot_frame_network_plane(
        u: mda.Universe,
        frame: int,
        out_png: pathlib.Path,
        reference_triplet: Tuple[str, str, str] = ("P", "P2", "P3"),
        r_cut_nm: float = RCUTOFF_NM,
) -> None:
    """
    Draw a 2-D network diagram of Na⁺ ions and phosphate centres by
    projecting onto the plane defined by *any* three IP6 atoms.

    Parameters
    ----------
    u : MDAnalysis.Universe
        MD trajectory already loaded.
    frame : int
        Index of the frame to visualise.
    out_png : pathlib.Path
        Output filename (.png).
    reference_triplet : (str, str, str)
        Names of *three* distinct IP6 atoms whose plane you wish to view
        (e.g. ("P","O9","P4") or any others in the same residue).
    r_cut_nm : float
        Only draw a dashed line when Na–peripheral‐O distance < r_cut_nm,
        or `None` to draw all links.
    """
    # 1. frame
    n_frames = len(u.trajectory)
    if not (0 <= frame < n_frames):
        raise IndexError(f"Frame {frame} out of range (0…{n_frames-1})")
    u.trajectory[frame]

    # 2. selections
    na_atoms = u.select_atoms("name NA")
    if len(na_atoms) == 0:
        logging.warning("No Na⁺ ions – snapshot skipped.")
        return

    # phosphate centres for plotting
    phos_atoms = {lab: u.select_atoms(f"resname I* and name {lab}").atoms[0]
                  for lab in PHOS_LABELS}

    # reference atoms for the plane (any three IP6 atom names)
    reference_atoms = {}
    for name in reference_triplet:
        sel = u.select_atoms(f"resname I* and name {name}")
        if len(sel) != 1:
            raise ValueError(f"Expected exactly one IP6 atom named '{name}', found {len(sel)}")
        reference_atoms[name] = sel.atoms[0]
    A, B, C = (reference_atoms[n] for n in reference_triplet)

    # 3. build basis
    rA = A.position
    v1 = B.position - rA
    v2 = C.position - rA
    normal = np.cross(v1, v2)
    if np.linalg.norm(normal) < 1e-6:
        raise ValueError("Reference atoms are collinear – cannot define a unique plane")
    e1 = v1 / np.linalg.norm(v1)
    e2 = np.cross(normal, e1)
    e2 /= np.linalg.norm(e2)

    # 4. project coords
    def project(X: np.ndarray) -> np.ndarray:
        R = X - rA
        return np.vstack((R @ e1, R @ e2)).T

    na_proj   = project(na_atoms.positions)
    phos_proj = project(np.array([phos_atoms[l].position for l in PHOS_LABELS]))

    # compute Na–peripheral-O minimum distances for dashed links
    periph_ag = {
        lab: u.select_atoms("resname I* and name " + " ".join(PHOS_OXYGENS[lab]))
        for lab in PHOS_LABELS
    }
    N = len(na_atoms)
    min_periph = np.full((N, len(PHOS_LABELS)), np.inf, dtype=np.float32)
    for j, lab in enumerate(PHOS_LABELS):
        d = distance_array(na_atoms.positions,
                           periph_ag[lab].positions,
                           box=u.dimensions)
        min_periph[:, j] = d.min(axis=1)
    closest_idx  = min_periph.argmin(axis=1)
    closest_dist = min_periph[np.arange(N), closest_idx]
    r_cut_A = None if r_cut_nm is None else r_cut_nm * 10.0

    # 5. plot
    fig, ax = plt.subplots(figsize=(6, 6))
    ax.set_aspect("equal")

    cmap = plt.cm.get_cmap("tab10", len(PHOS_LABELS))
    # legend for phosphates
    handles = []
    for j, (x, y) in enumerate(phos_proj):
        sc = ax.scatter(x, y, s=90, marker="o", edgecolor="k",
                        facecolor=cmap(j), zorder=3, label=PHOS_LABELS_CORRECTED[j])
        handles.append(sc)
    ax.legend(handles=handles, loc="best", frameon=True, fontsize=10)

    # Na ions
    ax.scatter(na_proj[:, 0], na_proj[:, 1],
               s=20, marker=".", color="black", zorder=2)

    # dashed links
    for i, (x_i, y_i) in enumerate(na_proj):
        j = closest_idx[i]
        if r_cut_A is not None and closest_dist[i] > r_cut_A:
            continue
        x_ph, y_ph = phos_proj[j]
        ax.plot([x_i, x_ph], [y_i, y_ph],
                ls="--", lw=0.6, color="grey", zorder=1)

    # labels & title
    ax.set_xlabel("e₁ projection / Å")
    ax.set_ylabel("e₂ projection / Å")
    tri = "–".join(reference_triplet)
    ax.set_title(f"plane ({tri}) – frame {frame}")
    ax.grid(True, ls=":", lw=0.4)

    plt.tight_layout()
    fig.savefig(out_png, dpi=300)
    plt.close(fig)
    logging.info("Plane-projected snapshot (%s) written to %s", tri, out_png)

# ╔═══════════════════════════════════════════════════════════════╗
# ┃   3-D NETWORK  – ORIGINAL COORDS, CAMERA = normal(P,P2,P4)   ┃
# ╚═══════════════════════════════════════════════════════════════╝
from mpl_toolkits.mplot3d import Axes3D   # noqa: F401 (import is needed)

def plot_frame_network_3d_fixed_view(
        u: mda.Universe,
        frame: int,
        out_png: pathlib.Path,
        anchors: Tuple[str, str, str] = ("P", "P2", "P4"),
        r_cut_nm: float = RCUTOFF_NM,
        look_along: int = +1,          # +1 → camera along  +n,  -1 → along –n
) -> None:
    """
    Draw a 3-D snapshot in laboratory coordinates but orient the Matplotlib
    camera so that it looks *perpendicularly* onto the plane defined by
    the three anchor phosphates (default: P, P2, P4).

    Parameters
    ----------
    u          : MDAnalysis.Universe
    frame      : int          Frame index (0 ≤ frame < n_frames).
    out_png    : pathlib.Path Output image path.
    anchors    : (str,str,str)  Three phosphate labels fixing the plane.
    r_cut_nm   : float | None   Show dashed link only if Na–P distance
                                < r_cut_nm (None → show all).
    look_along : {+1,-1}        Choose whether the camera points parallel
                                (+1) or antiparallel (-1) to the normal.
    """
    # ── 1. go to requested frame ─────────────────────────────────
    if not (0 <= frame < len(u.trajectory)):
        raise IndexError(f"Frame {frame} out of range.")
    u.trajectory[frame]

    # ── 2. selections ────────────────────────────────────────────
    na_atoms = u.select_atoms("name NA")
    if len(na_atoms) == 0:
        logging.warning("No Na⁺ ions — skipping 3-D plot.")
        return

    phos_atoms = {lab: u.select_atoms(f"resname I* and name {lab}").atoms[0]
                  for lab in PHOS_LABELS}

    try:
        r_P  = phos_atoms[anchors[0]].position
        r_P2 = phos_atoms[anchors[1]].position
        r_P4 = phos_atoms[anchors[2]].position
    except KeyError as e:
        raise ValueError(f"Anchor label {e.args[0]} not found among phosphates.")

    # ── 3. normal vector of anchor plane ─────────────────────────
    v1 = r_P2 - r_P
    v2 = r_P4 - r_P
    n = np.cross(v1, v2)
    if np.linalg.norm(n) < 1e-6:
        raise ValueError("Anchor atoms are colinear – plane undefined.")
    n /= np.linalg.norm(n)
    n *= look_along              # flip if the user wants the opposite side

    # convert normal into spherical view angles for Matplotlib
    elev_deg = np.degrees(np.arcsin(n[2]))              # arcsin(z/|n|)
    azim_rad = np.arctan2(n[1], n[0])                  # atan2(y,x)
    azim_deg = np.degrees(azim_rad)

    # ── 4. coordinates (original) & PERIPHERAL-oxygen distances ────
    na_pos   = na_atoms.positions                       # (N,3) Å
    phos_pos = np.array([phos_atoms[l].position for l in PHOS_LABELS])

    # Build a dict AtomGroup → peripheral Oxygens (three each)
    periph_ag = {label: u.select_atoms(
                    "resname I* and name " + " ".join(PHOS_OXYGENS[label]))
                for label in PHOS_LABELS}

    # matrix of minimum Na–O distances  (N_Na, 6)
    min_periph = np.full((len(na_atoms), 6), np.inf, dtype=np.float32)
    for j, label in enumerate(PHOS_LABELS):
        d = distance_array(na_pos, periph_ag[label].positions, box=u.dimensions)
        min_periph[:, j] = d.min(axis=1)          # min over the three O atoms

    closest_idx  = min_periph.argmin(axis=1)                      # phosphate index
    closest_dist = min_periph[np.arange(len(na_atoms)), closest_idx]
    r_cut_A      = r_cut_nm * 10.0                                # nm → Å
    # -------------------------------------------------------------------------

    # ── 5. plotting ─────────────────────────────────────────────
    fig = plt.figure(figsize=(6, 6))
    ax: Axes3D = fig.add_subplot(projection="3d")

    # ---- make all three coordinate panes pure white -----------------
    for axis in (ax.xaxis, ax.yaxis, ax.zaxis):
        axis.set_pane_color((1.0, 1.0, 1.0, 1.0))   # RGBA
        axis._axinfo['tick']['inward_factor']  = 0.0
        axis._axinfo['tick']['outward_factor'] = 0.0
        axis._axinfo['grid']['linewidth'] = 0.3
    # ----------------------------------------------------------------

    cmap = plt.cm.get_cmap("tab10", 6)
    # phosphates
    legend_handles = []                 # collect one handle per phosphate
    for j, (x, y, z) in enumerate(phos_pos):
        sc = ax.scatter(x, y, z, s=90, marker="o",
                        edgecolor="k", facecolor=cmap(j), depthshade=False,
                        label=PHOS_LABELS_CORRECTED[j])      # ← label for legend
        legend_handles.append(sc)                  # keep the handle
    
    ax.legend(handles=legend_handles,
          loc="lower right", frameon=True, fontsize=12) # title="Phosphates"

    # Na⁺ ions
    ax.scatter(na_pos[:, 0], na_pos[:, 1], na_pos[:, 2],
               s=15, c="black", marker=".", depthshade=False)

    # dashed connectors — draw *only* if Na is coordinated via peripheral O
    for i, (x_i, y_i, z_i) in enumerate(na_pos):
        if closest_dist[i] > r_cut_A:
            continue                    # not coordinated → no line
        j        = closest_idx[i]       # phosphate to which it is bound
        x_p, y_p, z_p = phos_pos[j]
        ax.plot([x_i, x_p], [y_i, y_p], [z_i, z_p],
                ls="--", lw=0.6, c="grey")

    # axis limits & labels
    span = np.ptp(np.vstack((na_pos, phos_pos)), axis=0)
    centre = np.mean(np.vstack((na_pos, phos_pos)), axis=0)
    max_half = span.max() / 2
    for axis, c in zip("xyz", centre):
        getattr(ax, f"set_{axis}lim")(c - max_half, c + max_half)

    ax.set_xlabel("x / Å")
    ax.set_ylabel("y / Å")
    ax.set_zlabel("z / Å")

    # --- fixed camera orientation ---
    ax.view_init(elev=elev_deg, azim=azim_deg)



    #ax.set_title(f"3-D Na–phosphate network (frame {frame})")
    ax.set_title(f"frame: {frame}", fontsize=14, loc="left")
    fig.savefig(out_png, dpi=300)
    plt.close(fig)
    logging.info("Fixed-view 3-D snapshot written to %s  (elev %.1f°, azim %.1f°)",
                 out_png, elev_deg, azim_deg)
    
# ╔═══════════════════════════════════════════════════════════════╗
# ┃                 DIRECTORY TRAVERSAL  (unchanged)             ┃
# ╚═══════════════════════════════════════════════════════════════╝
def main() -> None:
    script_dir = pathlib.Path(__file__).resolve()
    project_root = find_project_root(script_dir)
    process_dir  = project_root / PROCESS_DIR_NAME
    print(f"Process Dir: {process_dir}")
    setup_logging()

    ## ---------------------------------- determin global Y ----------------------
#
    ## ---------------------------------------------------------------------------
    ## 0.  parameters for the scan
    ## ---------------------------------------------------------------------------
    #r_max_nm = 1.2          # same value you use in the plotting function
    #nbins    = 240
    #Y_GLOBAL_MAX = 0.0      # will hold the global peak height
#
    ## ---------------------------------------------------------------------------
    ## 1.  loop over all micro-state directories
    ## ---------------------------------------------------------------------------
    #for microstate in sorted(process_dir.iterdir()):
    #    if not microstate.is_dir():
    #        continue                                           # skip files
#
    #    run_dir = microstate / RUN_DIR_NAME
    #    if not run_dir.is_dir():
    #        logging.debug("Skipping %s (no %s)", microstate.name, RUN_DIR_NAME)
    #        continue
#
    #    traj = run_dir / "md_center.xtc"
    #    top  = run_dir / "md.tpr"
    #    if not traj.is_file() or not top.is_file():
    #        logging.warning("Missing files in %s – skipped.", run_dir)
    #        continue                                           # nothing to analyse
#
    #    # -----------------------------------------------------------------------
    #    # 2.  build an MDAnalysis universe (no need to load the full trajectory)
    #    # -----------------------------------------------------------------------
    #    u = mda.Universe(top, traj, topology_format="TPR")
#
    #    # pre-select Na⁺ once
    #    na_ag = u.select_atoms("name NA")
#
    #    # -----------------------------------------------------------------------
    #    # 3.  loop over the six phosphates and update the global max
    #    # -----------------------------------------------------------------------
    #    for label in PHOS_LABELS:
    #        oxy_ag = u.select_atoms("resname I* and name " +
    #                                " ".join(PHOS_OXYGENS[label]))
#
    #        rdf = InterRDF(na_ag, oxy_ag,
    #                    range=(0, r_max_nm * 10),  # Å
    #                    nbins=nbins)
    #        rdf.run()                              # one frame → fast
#
    #        peak_height = (rdf.rdf * 3).max()      # ×3 → per-phosphate
    #        Y_GLOBAL_MAX = max(Y_GLOBAL_MAX, peak_height)
#
    ## add a little head-room so the tallest peak is not cut off
    #Y_GLOBAL_MAX *= 1.05
    #print(f"Global y-limit to use in all figures: {Y_GLOBAL_MAX:.2f}")
#
    ## --------------------------------- end determine global Y -------------------

    for microstate in sorted(process_dir.iterdir()):
        if not microstate.is_dir():
            continue

        run_dir = microstate / RUN_DIR_NAME
        if not run_dir.is_dir():
            logging.debug("Skipping %s (no %s)", microstate.name, RUN_DIR_NAME)
            continue

        dest_dir = microstate / PLOT_DIR_NAME
        dest_dir.mkdir(exist_ok=True)

        traj = run_dir / "md_center.xtc"
        top  = run_dir / "md.tpr"

        if not traj.is_file() or not top.is_file():
            logging.warning("Missing files in %s – skipped.", run_dir)
            continue

        logfile = dest_dir / "NaP_coordination.log"
        setup_logging(logfile)

        analyse_one_microstate(
            traj_file = traj,
            top_file  = top,
            dest_dir  = dest_dir,
            r_cut_nm  = RCUTOFF_NM,
        )


if __name__ == "__main__":
    main()