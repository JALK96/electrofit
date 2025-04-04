#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
For each λ->λ+ transition in base_dir, plots 3 curves:

1) P(ΔU)                 [black, solid]
2) exp(-β ΔU)            [black, dashed]
3) P(ΔU)*exp(-β ΔU)      [red, solid]

We use a 4-column layout. Each subplot is one transition.

Features:
- --skip-time: discard data prior to that time
- --trim: clamp x-range to [p-th, (100-p)-th] percentile
- --nsigma: override x-range to [mean - N*σ, mean + N*σ]
  (this takes priority over --trim if both are given)
- --normalize: scale all 3 curves so their global max = 1
- automatically sets y-limit so that P(ΔU) is visible
  in the non-normalized case (exponential might go off-scale).
"""

import os
import re
import sys
import argparse
import math
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import gaussian_kde
from matplotlib.offsetbox import AnchoredOffsetbox, TextArea, VPacker

sns.set_context("talk", font_scale=1.2)

kB = 0.0083144626  # kJ/(mol*K) - Boltzmann constant approx.

#####################################################################
#  Parsing & Data Gathering
#####################################################################

def parse_dhdl_xvg_current_and_targets(dhdl_file, debug=False, skip_time=0.0):
    """
    Parse a GROMACS dhdl.xvg to find current_lambda, a map of
    target_lambda->column_index, plus time & data arrays.
    
    If skip_time>0, discard rows with time<skip_time.
    """
    current_lambda = None
    targets_map = {}
    legends_map = {}

    time_list = []
    data_list = []

    dh_regex = re.compile(r'@ s(\d+)\s+legend\s+"\\xD\\f\{\}H\s+\\xl\\f\{\}\s+to\s+([^"]+)"')
    subtitle_lambda_regex = re.compile(r'@ subtitle ".*fep-lambda\s*=\s*([\d\.]+)"')
    fep_lambda_line_regex = re.compile(r'fep-lambda\s*=\s*([\d\.]+)')

    with open(dhdl_file, 'r', encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.rstrip("\n")
            if not line or line.startswith('#'):
                continue
            if line.startswith('@'):
                # 1) "fep-lambda" from subtitle
                match_sub = subtitle_lambda_regex.search(line)
                if match_sub:
                    maybe_val = match_sub.group(1)
                    try:
                        current_lambda = float(maybe_val)
                    except ValueError:
                        pass

                # 2) "fep-lambda =" in the line
                match_fep = fep_lambda_line_regex.search(line)
                if match_fep:
                    maybe_val = match_fep.group(1)
                    try:
                        current_lambda = float(maybe_val)
                    except ValueError:
                        pass

                # 3) lines: @ s(\d+) legend "\xD\f{}H \xl\f{} to <float>"
                dh_match = dh_regex.search(line)
                if dh_match:
                    s_index_str = dh_match.group(1)
                    val_str = dh_match.group(2).strip()
                    try:
                        s_index = int(s_index_str)
                        t_val = float(val_str)
                        targets_map[t_val] = s_index
                    except ValueError:
                        pass

                # keep all legend lines
                leg_match = re.search(r'@ s(\d+)\s+legend\s+"([^"]+)"', line)
                if leg_match:
                    s_i = int(leg_match.group(1))
                    legend_txt = leg_match.group(2)
                    legends_map[s_i] = legend_txt

                continue

            # numeric data
            if line.startswith('@'):
                continue

            parts = line.split()
            floats = [float(x) for x in parts]
            time_list.append(floats[0])
            data_list.append(floats[1:])

    time_arr = np.array(time_list)
    data_arr = np.array(data_list)

    if skip_time>0:
        mask = (time_arr>=skip_time)
        time_arr = time_arr[mask]
        data_arr = data_arr[mask,:]

    if debug:
        print(f"[DEBUG] parse_dhdl_xvg_current_and_targets({dhdl_file}) => "
              f"λ={current_lambda}, shape={data_arr.shape}, skip_time={skip_time}, targets={targets_map}")
    return current_lambda, targets_map, time_arr, data_arr, legends_map


def find_next_lambda_and_column(current_lambda, targets_map, data_shape, offset_fix=0, debug=False):
    bigger = [x for x in targets_map.keys() if x>current_lambda]
    if not bigger:
        return None, None
    next_lambda = min(bigger)
    s_index = targets_map[next_lambda]
    col_index = s_index + offset_fix
    if col_index<0 or col_index>=data_shape[1]:
        if debug:
            print(f"[DEBUG] col_index={col_index} out of range")
        return None, None
    return next_lambda, col_index


def parse_single_subdir(dhdl_file, offset_fix=0, debug=False, skip_time=0.0):
    cLam, targets_map, time_arr, data_arr, legends_map = parse_dhdl_xvg_current_and_targets(
        dhdl_file, debug=debug, skip_time=skip_time)
    if cLam is None or data_arr.size==0:
        return None,None,None,None
    nLam,col_idx = find_next_lambda_and_column(cLam, targets_map, data_arr.shape,
                                               offset_fix=offset_fix, debug=debug)
    if nLam is None or col_idx is None:
        return cLam,None,None,None
    deltaU = data_arr[:,col_idx]
    return cLam, nLam, deltaU, col_idx


def gather_dU_per_subdir(base_dir, offset_fix=0, debug=False, skip_time=0.0):
    subdirs = []
    for name in os.listdir(base_dir):
        if name.startswith("lambda_"):
            try:
                idx = float(name.split("_")[-1])
            except ValueError:
                idx = None
            fpath = os.path.join(base_dir, name,"dhdl.xvg")
            if os.path.isfile(fpath):
                subdirs.append((idx,fpath))
    subdirs.sort(key=lambda x: (x[0] if x[0] is not None else 1e9))

    out_list=[]
    for idx,fpath in subdirs:
        cLam,nLam,deltaU,col_idx = parse_single_subdir(
            fpath, offset_fix=offset_fix, debug=debug, skip_time=skip_time)
        if cLam is None:
            print(f"[WARN] No λ found in {fpath}")
            continue
        if nLam is None or deltaU is None:
            print(f"[INFO] No next-lambda or no data for λ={cLam:.3f} in {fpath}, skip.")
            continue
        out_list.append((cLam,nLam,deltaU,col_idx,idx))
    out_list.sort(key=lambda x:x[0])
    return out_list


#####################################################################
#  Plotting
#####################################################################

def compute_kde(deltaU):
    """Return a (kde_func, (data_mean, data_sigma)) so we can evaluate later."""
    if len(deltaU)<2:
        return None,(np.nan, np.nan)
    kde_ = gaussian_kde(deltaU)
    mean_ = np.mean(deltaU)
    sigma_ = np.std(deltaU)
    return kde_, (mean_,sigma_)


def plot_3curves_for_each_transition(base_dir,
                                     offset_fix=0,
                                     debug=False,
                                     normalize=False,
                                     trim=None,
                                     nsigma=None,
                                     skip_time=0.0,
                                     temperature=310.0,
                                     outname=None):
    """
    For each transition:
       1) P(ΔU)  -> black solid
       2) exp(-β ΔU) -> black dashed
       3) product -> red solid

    x-range logic:
      - If nsigma is given => [mean - nsigma*sigma, mean + nsigma*sigma]
      - else if trim is given => [trim%-ile, (100-trim)%-ile]
      - else => [min(ΔU)-5%, max(ΔU)+5%]
    
    if --normalize => scale all 3 so their combined max=1
    else => raw scale but y-limit is 5% above P(ΔU)'s max

    skip_time => discard data prior to that time
    """
    dU_list = gather_dU_per_subdir(base_dir, offset_fix=offset_fix,
                                   debug=debug, skip_time=skip_time)
    n_plots = len(dU_list)
    if n_plots==0:
        print("[WARN] No transitions found or no data.")
        return

    ncols=5
    nrows=math.ceil(n_plots/ncols)
    fig,axs = plt.subplots(nrows=nrows,ncols=ncols,
                           figsize=(4*ncols,3*nrows),
                           sharex=False, sharey=False)
    axs=np.atleast_2d(axs)

    beta = 1./(kB*temperature)

    for i,(cLam,nLam,deltaU,col_idx,subdir_idx) in enumerate(dU_list):
        row = i//ncols
        col = i%ncols
        ax=axs[row,col]

        if len(deltaU)<2:
            ax.set_title(f"λ={cLam:.3f}→{nLam:.3f} [no data]", fontsize=10)
            continue

        # 1) compute KDE
        kde_,(mu,sigma)=compute_kde(deltaU)
        if kde_ is None:
            ax.set_title(f"λ={cLam:.3f}→{nLam:.3f} [no data]", fontsize=10)
            continue

        # define x-range
        if nsigma is not None and nsigma>0 and not math.isnan(mu):
            # x in [mu - nsigma*sigma, mu + nsigma*sigma]
            left = mu - nsigma*sigma
            right= mu + nsigma*sigma
            # fallback if sigma=0
            if sigma==0:
                left-=0.5
                right+=0.5
        elif trim is not None and trim>0 and trim<50:
            # percentile approach
            p_lo = np.percentile(deltaU, trim)
            p_hi = np.percentile(deltaU, 100-trim)
            if p_lo==p_hi:
                p_lo-=0.5
                p_hi+=0.5
            left,right = p_lo,p_hi
        else:
            # fallback to min..max +5%
            m = np.min(deltaU)
            M = np.max(deltaU)
            if m==M:
                m-=0.5
                M+=0.5
            span = M-m
            left = m-0.05*span
            right= M+0.05*span

        xvals = np.linspace(left,right,200)
        pvals = kde_(xvals)  # black solid
        expvals = np.exp(-beta*xvals) # black dashed
        productvals = pvals*expvals   # red

        # do we scale?
        if normalize:
            pvals = pvals / np.max(pvals)
            expvals = expvals / np.max(expvals)
            productvals = productvals / np.max(productvals)
            ax.plot(xvals, pvals, 'k-', label="P(ΔU)")
            ax.plot(xvals, expvals,'k--',label="exp(-βΔU)")
            ax.plot(xvals, productvals,'b-',label="P(ΔU)*exp(-βΔU)")
            ax.set_ylim(bottom=0, top=1.05)
        else:
            # raw scale
            ax.plot(xvals, pvals, 'k-', label="P(ΔU)")
            ax.plot(xvals, expvals,'k--',label="exp(-βΔU)")
            ax.plot(xvals, productvals,'b-',label="P(ΔU)*exp(-βΔU)")
            pmax=np.max(pvals)
            if pmax>0:
                ax.set_ylim(bottom=0, top=1.05*pmax)
            else:
                ax.set_ylim(bottom=0, top=1.0)

        ax.set_title(f"λ={cLam:.3f}→{nLam:.3f}")
        #ax.legend(fontsize='x-small')

    # hide extra subplots
    total_subplots=nrows*ncols
    for j in range(n_plots,total_subplots):
        r = j//ncols
        c = j%ncols
        axs[r,c].set_visible(False)

    # remove per-axis labels
    for ax in axs.flatten():
        ax.set_xlabel("")
        ax.set_ylabel("")

    # add common x-label
    fig.text(0.5,0.08, r"$\Delta U$ (kJ/mol)", ha="center", va="center", fontsize=30)

    # add a y-label
    if normalize:
        # Composite label: black for P and exp, red for product
        ybox1 = TextArea(r"Normalized: $P(\Delta U)$",
                            textprops=dict(color="k", size=30, rotation=90,
                                        ha='center', va='bottom'))
        ybox2 = TextArea(r", $\exp[-\beta \Delta U],$",
                            textprops=dict(color="k", size=30, rotation=90,
                                        ha='center', va='bottom'))
        ybox3 = TextArea(r"$P(\Delta U) \times \exp[-\beta \Delta U]$",
                            textprops=dict(color="blue", size=30, rotation=90,
                                        ha='center', va='bottom'))
        ybox = VPacker(children=[ybox3, ybox2, ybox1],
                        align="center", pad=0, sep=2)

        anchored_ybox = AnchoredOffsetbox(
            loc='center left',
            child=ybox,
            pad=0.,
            frameon=False,
            bbox_to_anchor=(0.08, 0.5),
            bbox_transform=fig.transFigure,
            borderpad=0.
        )
        fig.add_artist(anchored_ybox)
    else:
        # Composite label: black for P and exp, red for product
        ybox1 = TextArea(r"$P(\Delta U)$",
                            textprops=dict(color="k", size=30, rotation=90,
                                        ha='center', va='bottom'))
        ybox2 = TextArea(r", $\exp[-\beta \Delta U],$",
                            textprops=dict(color="k", size=30, rotation=90,
                                        ha='center', va='bottom'))
        ybox3 = TextArea(r"$P(\Delta U) \times \exp[-\beta \Delta U]$",
                            textprops=dict(color="blue", size=30, rotation=90,
                                        ha='center', va='bottom'))
        ybox = VPacker(children=[ybox3, ybox2, ybox1],
                        align="center", pad=0, sep=2)

        anchored_ybox = AnchoredOffsetbox(
            loc='center left',
            child=ybox,
            pad=0.,
            frameon=False,
            bbox_to_anchor=(0.08, 0.5),
            bbox_transform=fig.transFigure,
            borderpad=0.
        )
        fig.add_artist(anchored_ybox)

    plt.tight_layout(rect=[0.1,0.1,1,0.95])
    # Determine default output name if not provided
    if outname is None:
        outname = "plot_3curves_each_transition_nsigma.pdf"
    plt.savefig(outname)
    print(f"[INFO] saved figure to {outname}")


#####################################################################
#  Main
#####################################################################

def main():
    parser = argparse.ArgumentParser(
        description="Plot P(ΔU), exp(-βΔU), & product for each transition in base_dir."
    )
    parser.add_argument("base_dir", help="Directory with subfolders 'lambda_*' each containing dhdl.xvg.")
    parser.add_argument("--offset", type=int, default=0,
                        help="Offset for sX->data column.")
    parser.add_argument("--skip-time", type=float, default=0.0,
                        help="Discard data with time<skip_time (ps).")
    parser.add_argument("--trim", type=float, default=None,
                        help="Clamp x-range to [trim%, (100-trim)%]. e.g. '--trim 2' => [2%,98%].")
    parser.add_argument("--nsigma", type=float, default=None,
                        help="Override x-range => mean ± nsigma*std. If set, overrides --trim.")
    parser.add_argument("--normalize", action="store_true",
                        help="If set, scale all 3 curves so combined max=1.")
    parser.add_argument("--temperature", type=float, default=310,
                        help="Temperature in K for β=1/(kB*T). Default=310.")
    parser.add_argument("--debug", action="store_true",
                        help="Enable debug prints.")
    parser.add_argument("--out_name", type=str, default=None,
                        help=("Output filename for the plot. If not provided, "
                              "the script will default to 'plot_3curves_each_transition_nsigma_<basename>.pdf', "
                              "where <basename> is derived from the base_dir."))
    args=parser.parse_args()
    base = os.path.abspath(args.base_dir)

    # If the user did not specify an output filename, generate one using the base directory's basename.
    if args.out_name:
        outname = args.out_name
    else:
        base_label = os.path.basename(os.path.normpath(base))
        outname = f"plot_3curves_each_transition_nsigma_{base_label}.pdf"

    plot_3curves_for_each_transition(base_dir=base,
                                     offset_fix=args.offset,
                                     debug=args.debug,
                                     normalize=args.normalize,
                                     trim=args.trim,
                                     nsigma=args.nsigma,
                                     skip_time=args.skip_time,
                                     temperature=args.temperature,
                                     outname=outname)

if __name__=="__main__":
    main()