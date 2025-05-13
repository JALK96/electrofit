#!/usr/bin/env python
"""
Script: abfe_workflow_no_class.py
Description: This script reproduces the workflow used in the ABFE class but does so by directly 
             calling the underlying functions from alchemlyb. It performs the following:
             1. Collects file names matching a specified pattern.
             2. Parses simulation output (u_nk and dHdl data) using a chosen parser (here, GROMACS).
             3. Performs preprocessing (skipping equilibration and decorrelating the time series).
             4. Runs free energy estimators (MBAR, BAR, TI) separately on:
                  a) Decorrelated (uncorrelated) data
                  b) Raw (correlated) data 
             5. Generates summary tables and diagnostic plots for both cases.
             6. Logs which state corresponds to which lambda value, so you know, e.g. that state 0--1 is from lambda = 0 to 0.1.
             
The workflow steps mimic the operations inside the ABFE class but are executed here directly.
"""

# =============================================================================
# 1. Import standard libraries and set up logging and paths
# =============================================================================
import os
import warnings
import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import joblib
from loguru import logger
import logging
import seaborn as sns
from matplotlib.font_manager import FontProperties as FP
from sklearn.base import BaseEstimator


try:
    from electrofit.helper.set_logging import setup_logging
except ImportError:
    # If not available, use a basic logging configuration.
    def setup_logging(log_path):
        logging.basicConfig(
            filename=log_path,
            level=logging.DEBUG,
            format="%(asctime)s %(levelname)s: %(message)s"
        )



sns.set_context("talk", font_scale=0.9)

# ─────────────────────────────────────────────────────────────
# MBAR solver settings – tweak here if convergence is difficult
MBAR_OPTS: dict[str, float] = {
    # stop Newton–Raphson after this many iterations (default 10000)
    "maximum_iterations": 2000,
    # fractional change in ΔG vector required for convergence (default 1e‑12)
    "relative_tolerance": 1e-6,
    "verbose": False,
}
# ─────────────────────────────────────────────────────────────

# =============================================================================
# 2. Import functions from alchemlyb used in the ABFE workflow
# =============================================================================
# Parsers: Choose the one matching your simulation software (here GROMACS is assumed)
from alchemlyb.parsing import gmx  # Alternatively: amber, parquet if needed

# Estimators and estimator lists for free energy calculations
from alchemlyb.estimators import MBAR, BAR, TI, FEP_ESTIMATORS, TI_ESTIMATORS

# Convergence analysis functions
from alchemlyb.convergence import forward_backward_convergence

# Unit conversion for free energy outputs
from alchemlyb.postprocessors.units import get_unit_converter

# Decorrelate and subsample the time series
from alchemlyb.preprocessing.subsampling import decorrelate_dhdl, decorrelate_u_nk

# Plotting functions for diagnostics
from alchemlyb.visualisation import plot_mbar_overlap_matrix

# Function to concatenate a list of DataFrames
from alchemlyb import concat

class _EstimatorMixOut:
    """This class creates view for the attributes: states_, delta_f_, d_delta_f_,
    delta_h_, d_delta_h_, delta_sT_, d_delta_sT_ for the estimator class to consume.
    """

    _d_delta_f_ = None
    _delta_f_ = None
    _states_ = None
    _d_delta_h_ = None
    _delta_h_ = None
    _d_delta_sT_ = None
    _delta_sT_ = None

    @property
    def d_delta_f_(self):
        return self._d_delta_f_

    @property
    def delta_f_(self):
        return self._delta_f_

    @property
    def d_delta_h_(self):
        return self._d_delta_h_

    @property
    def delta_h_(self):
        return self._delta_h_

    @property
    def d_delta_sT_(self):
        return self._d_delta_sT_

    @property
    def delta_sT_(self):
        return self._delta_sT_

    @property
    def states_(self):
        return self._states_

class TI(BaseEstimator, _EstimatorMixOut):
    """Thermodynamic integration (TI).

    Parameters
    ----------

    verbose : bool, optional
        Set to True if verbose debug output is desired.

    Attributes
    ----------

    delta_f_ : DataFrame
        The estimated dimensionless free energy difference between each state.

    d_delta_f_ : DataFrame
        The estimated statistical uncertainty (one standard deviation) in
        dimensionless free energy differences.

    states_ : list
        Lambda states for which free energy differences were obtained.

    dhdl : DataFrame
        The estimated dhdl of each state.


    .. versionchanged:: 1.0.0
       `delta_f_`, `d_delta_f_`, `states_` are view of the original object.

    """

    def __init__(self, verbose=False):
        self.verbose = verbose

    def fit(self, dHdl):
        """
        Compute free energy differences between each state by integrating
        dHdl across lambda values.

        Parameters
        ----------
        dHdl : DataFrame
            dHdl[n,k] is the potential energy gradient with respect to lambda
            for each configuration n and lambda k.

        """

        # sort by state so that rows from same state are in contiguous blocks,
        # and adjacent states are next to each other
        dHdl = dHdl.sort_index(level=dHdl.index.names[1:])

        # obtain the mean and variance of the mean for each state
        # variance calculation assumes no correlation between points
        # used to calculate mean
        means = dHdl.groupby(level=dHdl.index.names[1:]).mean()
        variances = np.square(dHdl.groupby(level=dHdl.index.names[1:]).sem())

        # get the lambda names
        l_types = dHdl.index.names[1:]

        # obtain vector of delta lambdas between each state
        # Fix issue #148, where for pandas == 1.3.0
        # dl = means.reset_index()[list(means.index.names[:])].diff().iloc[1:].values
        dl = means.reset_index()[means.index.names[:]].diff().iloc[1:].values

        # apply trapezoid rule to obtain DF between each adjacent state
        deltas = (dl * (means.iloc[:-1].values + means.iloc[1:].values) / 2).sum(axis=1)

        # build matrix of deltas between each state
        adelta = np.zeros((len(deltas) + 1, len(deltas) + 1))
        ad_delta = np.zeros_like(adelta)

        for j in range(len(deltas)):
            out = []
            dout = []
            for i in range(len(deltas) - j):
                out.append(deltas[i] + deltas[i + 1 : i + j + 1].sum())

                # Define additional zero lambda
                a = [0.0] * len(l_types)

                # Define dl series' with additional zero lambda on the left and right
                dll = np.insert(dl[i : i + j + 1], 0, [a], axis=0)
                dlr = np.append(dl[i : i + j + 1], [a], axis=0)

                # Get a series of the form: x1, x1 + x2, ..., x(n-1) + x(n), x(n)
                dllr = dll + dlr

                # Append deviation of free energy difference between state i and i+j+1
                dout.append(
                    (dllr**2 * variances.iloc[i : i + j + 2].values / 4)
                    .sum(axis=1)
                    .sum()
                )
            adelta += np.diagflat(np.array(out), k=j + 1)
            ad_delta += np.diagflat(np.array(dout), k=j + 1)

        # yield standard delta_f_ free energies between each state
        self._delta_f_ = pd.DataFrame(
            adelta - adelta.T, columns=means.index.values, index=means.index.values
        )
        self.dhdl = means

        # yield standard deviation d_delta_f_ between each state
        self._d_delta_f_ = pd.DataFrame(
            np.sqrt(ad_delta + ad_delta.T),
            columns=variances.index.values,
            index=variances.index.values,
        )

        self._states_ = means.index.values.tolist()

        self._delta_f_.attrs = dHdl.attrs
        self._d_delta_f_.attrs = dHdl.attrs
        self.dhdl.attrs = dHdl.attrs

        return self

    def separate_dhdl(self):
        """
        For transitions with multiple lambda, the attr:`dhdl` would return
        a :class:`~pandas.DataFrame` which gives the dHdl for all the lambda
        states, regardless of whether it is perturbed or not. This function
        creates a list of :class:`pandas.Series` for each lambda, where each
        :class:`pandas.Series` describes the potential energy gradient for the
        lambdas state that is perturbed.

        Returns
        ----------
        dHdl_list : list
            A list of :class:`pandas.Series` such that ``dHdl_list[k]`` is the
            potential energy gradient with respect to lambda for each
            configuration that lambda k is perturbed.
        """
        if len(self.dhdl.index.names) == 1:
            name = self.dhdl.columns[0]
            return [
                self.dhdl[name],
            ]
        dhdl_list = []
        # get the lambda names
        l_types = self.dhdl.index.names
        # obtain bool of changed lambdas between each state
        # Fix issue #148, where for pandas == 1.3.0
        # lambdas = self.dhdl.reset_index()[list(l_types)]
        lambdas = self.dhdl.reset_index()[l_types]
        diff = lambdas.diff().to_numpy(dtype="bool")
        # diff will give the first row as NaN so need to fix that
        diff[0, :] = diff[1, :]
        # Make sure that the start point is set to true as well
        diff[:-1, :] = diff[:-1, :] | diff[1:, :]
        for i in range(len(l_types)):
            if any(diff[:, i]):
                new = self.dhdl.iloc[diff[:, i], i]
                # drop all other index
                for l in l_types:
                    if l != l_types[i]:
                        new = new.reset_index(l, drop=True)
                new.attrs = self.dhdl.attrs
                dhdl_list.append(new)
        return dhdl_list

def plot_convergence(dataframe, units=None, final_error=None, ax=None):
    """
    Plot the forward and backward convergence.

    Parameters
    ----------
    dataframe : pandas.DataFrame
         Must have columns 'Forward' and 'Backward'. Optionally,
         'Forward_Error' and 'Backward_Error' for error bars.
         The dataframe.attrs should include information required for unit conversion.
    units : str, optional
         The desired units. If provided, the dataframe is converted with a unit converter.
    final_error : float, optional
         If provided, uses this as the error band around the final backward estimate.
         Otherwise, uses the error in the last value of 'Backward_Error'.
    ax : matplotlib.axes.Axes, optional
         An existing axes on which to plot the data; if not provided, a new figure and axis are created.

    Returns
    -------
    ax : matplotlib.axes.Axes
         Axes with the convergence plot drawn.
    """
    # Apply unit conversion if requested.
    if units is not None:
        dataframe = get_unit_converter(units)(dataframe)
    
    forward = dataframe["Forward"].to_numpy()
    if "Forward_Error" in dataframe:
        forward_error = dataframe["Forward_Error"].to_numpy()
    else:
        forward_error = np.zeros(len(forward))

    backward = dataframe["Backward"].to_numpy()
    if "Backward_Error" in dataframe:
        backward_error = dataframe["Backward_Error"].to_numpy()
    else:
        backward_error = np.zeros(len(backward))
    
    # Create new figure and axis if none is provided.
    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 6))
    

    
    # Generate time fractions for plotting.
    f_ts = np.linspace(0, 1, len(forward) + 1)[1:]
    r_ts = np.linspace(0, 1, len(backward) + 1)[1:]
    
    # Use the final backward error if provided, otherwise take the last backward error.
    if final_error is None:
        final_error = backward_error[-1]
    
    # Plot a horizontal band at the final backward value.
    if np.isfinite(backward[-1]) and np.isfinite(final_error):
        ax.fill_between([0, 1],
                        backward[-1] - final_error,
                        backward[-1] + final_error,
                        color="#D2B9D3", zorder=1)
    
    # Plot forward convergence with errorbars.
    line_forward = ax.errorbar(
        f_ts, forward, yerr=forward_error,
        color="#736AFF", lw=3, marker="o",
        mfc="w", mew=2.5, mec="#736AFF", ms=12, zorder=2)
    
    # Plot backward convergence with errorbars.
    line_backward = ax.errorbar(
        r_ts, backward, yerr=backward_error,
        color="#C11B17", lw=3, marker="o",
        mfc="w", mew=2.5, mec="#C11B17", ms=12, zorder=3)
    
    # Set x-ticks; here we space them out and format as two-decimal strings.
    xticks_spacing = max(1, len(r_ts) // 10)
    xticks = r_ts[::xticks_spacing]
    ax.set_xticks(xticks)
    ax.set_xticklabels([f"{tick:.1f}" for tick in xticks])
    ax.tick_params(axis="y")

    # Set spine properties.
    ax.spines["bottom"].set_color("#D2B9D3")
    ax.spines["bottom"].set_linewidth(2)
    ax.spines["left"].set_color("#D2B9D3")
    ax.spines["left"].set_linewidth(2)
    
    # Remove top and right spines.
    for side in ["top", "right"]:
        ax.spines[side].set_color("none")
    
    ax.xaxis.set_ticks_position("bottom")
    ax.yaxis.set_ticks_position("left")

    # Create legend using custom font properties.
    fp = FP(size=18)
    ax.legend([line_forward[0], line_backward[0]], ["Forward", "Reverse"],
              loc="best", prop=fp, frameon=False)
    
    ax.set_xlabel("Fraction of simulation time", color="k", fontsize=20) #color="#151B54"
    ax.set_ylabel(r"$\Delta G$ ({})".format(units if units else ""), color="k", fontsize=20)

    plt.tick_params(axis="x", color="#D2B9D3")
    plt.tick_params(axis="y", color="#D2B9D3")
    
    plt.tight_layout()
    return ax

def plot_ti_dhdl(dhdl_data, labels=None, colors=None, units=None, ax=None):
    """Plot the dhdl of TI.

    Parameters
    ----------
    dhdl_data : :class:`~alchemlyb.estimators.TI` or list
        One or more :class:`~alchemlyb.estimators.TI` estimator, where the
        dhdl value will be taken from.
    labels : List
        list of labels for labelling all the alchemical transformations.
    colors : List
        list of colors for plotting all the alchemical transformations.
        Default: ['r', 'g', '#7F38EC', '#9F000F', 'b', 'y']
    units : str
        The unit of the estimate. The default is `None`, which is to use the
        unit in the input. Setting this will change the output unit.
    ax : matplotlib.axes.Axes
        Matplotlib axes object where the plot will be drawn on. If ``ax=None``,
        a new axes will be generated.

    Returns
    -------
    matplotlib.axes.Axes
        An axes with the TI dhdl drawn.

    Note
    ----
    The code is taken and modified from
    `Alchemical Analysis <https://github.com/MobleyLab/alchemical-analysis>`_.


    .. versionchanged:: 1.0.0
        If no units is given, the `units` in the input will be used.

    .. versionchanged:: 0.5.0
        The `units` will be used to change the underlying data instead of only
        changing the figure legend.

    .. versionadded:: 0.4.0
    """
    # Make it into a list
    # separate_dhdl method is used so that the input for the actual plotting
    # Function are a uniformed list of series object which only contains one
    # lambda.
    if not isinstance(dhdl_data, list):
        dhdl_list = dhdl_data.separate_dhdl()
    else:
        dhdl_list = []
        for dhdl in dhdl_data:
            dhdl_list.extend(dhdl.separate_dhdl())

    # Convert unit
    if units is None:
        units = dhdl_list[0].attrs["energy_unit"]

    new_unit = []
    convert = get_unit_converter(units)
    for dhdl in dhdl_list:
        new_unit.append(convert(dhdl))

    dhdl_list = new_unit
    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 6))

    ax.spines["bottom"].set_position("zero")
    ax.spines["top"].set_color("none")
    ax.spines["right"].set_color("none")
    ax.xaxis.set_ticks_position("bottom")
    ax.yaxis.set_ticks_position("left")

    for k, spine in ax.spines.items():
        spine.set_zorder(12.2)

    # Make level names
    if labels is None:
        lv_names2 = []
        for dhdl in dhdl_list:
            # Assume that the dhdl has only one columns
            lv_names2.append(dhdl.name.capitalize())
    else:
        if len(labels) == len(dhdl_list):
            lv_names2 = labels
        else:  # pragma: no cover
            raise ValueError(
                "Length of labels ({}) should be the same as the number of data ({})".format(
                    len(labels), len(dhdl_list)
                )
            )

    if colors is None:
        colors = ["#6698FF", "r", "g", "#7F38EC", "#9F000F", "b", "y"]
    else:
        if len(colors) >= len(dhdl_list):
            pass
        else:  # pragma: no cover
            raise ValueError(
                "Number of colors ({}) should be larger than the number of data ({})".format(
                    len(labels), len(dhdl_list)
                )
            )

    # Get the real data out
    xs, ndx, dx = [0], 0, 0.001
    min_y, max_y = 0, 0
    for dhdl in dhdl_list:
        x = dhdl.index.values
        y = dhdl.values.ravel()

        min_y = min(y.min(), min_y)
        max_y = max(y.max(), max_y)

        for i in range(len(x) - 1):
            if i % 2 == 0:
                ax.fill_between(
                    x[i : i + 2] + ndx, 0, y[i : i + 2], color=colors[ndx], alpha=1.0
                )
            else:
                ax.fill_between(
                    x[i : i + 2] + ndx, 0, y[i : i + 2], color=colors[ndx], alpha=0.5
                )

        xlegend = [-100 * wnum for wnum in range(len(lv_names2))]
        ax.plot(
            xlegend,
            [0 * wnum for wnum in xlegend],
            ls="-",
            color=colors[ndx],
            label=lv_names2[ndx],
        )
        xs += (x + ndx).tolist()[1:]
        ndx += 1

    # Make sure the tick labels are not overcrowded.
    xs = np.array(xs)
    dl_mat = np.array([xs - i for i in xs])
    ri = range(len(xs))

    def getInd(r=ri, z=[0]):
        primo = r[0]
        min_dl = ndx * 0.02 * 2 ** (primo > 0)
        if dl_mat[primo].max() < min_dl:
            return z
        for i in r:  # pragma: no cover
            for j in range(len(xs)):
                if dl_mat[i, j] > min_dl:
                    z.append(j)
                    return getInd(ri[j:], z)

    xt = []
    for i in range(len(xs)):
        if i in getInd():
            xt.append(i)
        else:
            xt.append("")

    plt.xticks(xs[1:], xt[1:], fontsize=15)
    ax.yaxis.label.set_size(16)

    # Remove the abscissa ticks and set up the axes limits.
    for tick in ax.get_xticklines():
        tick.set_visible(False)
    ax.set_xlim(0, ndx)
    min_y *= 1.01
    max_y *= 1.01

    # Modified so that the x label won't conflict with the lambda label
    min_y -= (max_y - min_y) * 0.1

    ax.set_ylim(min_y, max_y)

    for i, j in zip(xs[1:], xt[1:]):
        ax.annotate(
            ("{:.2f}".format(i - 1.0 if i > 1.0 else i) if not j == "" else ""),
            xy=(i, 0),
            size=15,
            rotation=90,
            va="bottom",
            ha="center",
            color="#151B54",
        )
    if ndx > 1:
        lenticks = len(ax.get_ymajorticklabels()) - 1
        if min_y < 0:
            lenticks -= 1
        if lenticks < 5:  # pragma: no cover
            from matplotlib.ticker import AutoMinorLocator as AML

            ax.yaxis.set_minor_locator(AML())
    ax.grid(which="both", color="w", lw=0.25, axis="y", zorder=12)
    ax.set_ylabel(
        r"$\langle{\frac{\partial U}{\partial\lambda}}\rangle_{\lambda}$"
        + "({})".format(units),
        fontsize=20,
        color="#151B54",
    )
    ax.annotate(
        r"$\mathit{\lambda}$",
        xy=(0, 0),
        xytext=(0.5, -0.05),
        size=18,
        textcoords="axes fraction",
        va="top",
        ha="center",
        color="#151B54",
    )
    #lege = ax.legend(prop=FP(size=16), frameon=False, loc=1)
    #for l in lege.legend_handles:
    #    l.set_linewidth(10)
    return ax

def plot_dF_state(
    estimators, labels=None, colors=None, units=None, orientation="portrait", nb=10
):
    """Plot the dhdl of TI.

    Parameters
    ----------
    estimators : :class:`~alchemlyb.estimators` or list
        One or more :class:`~alchemlyb.estimators`, where the
        dhdl value will be taken from. For more than one estimators
        with more than one alchemical transformation, a list of list format
        is used.
    labels : List
        list of labels for labelling different estimators.
    colors : List
        list of colors for plotting different estimators.
    units : str
        The unit of the estimate. The default is `None`, which is to use the
        unit in the input. Setting this will change the output unit.
    orientation : string
        The orientation of the figure. Can be `portrait` or `landscape`
    nb : int
        Maximum number of dF states in one row in the `portrait` mode

    Returns
    -------
    matplotlib.figure.Figure
        An Figure with the dF states drawn.

    Note
    ----
    The code is taken and modified from
    `Alchemical Analysis <https://github.com/MobleyLab/alchemical-analysis>`_.


    .. versionchanged:: 1.0.0
        If no units is given, the `units` in the input will be used.

    .. versionchanged:: 0.5.0
        The `units` will be used to change the underlying data instead of only
        changing the figure legend.

    .. versionadded:: 0.4.0
    """
    try:
        len(estimators)
    except TypeError:
        estimators = [
            estimators,
        ]

    formatted_data = []
    for dhdl in estimators:
        try:
            len(dhdl)
            formatted_data.append(dhdl)
        except TypeError:
            formatted_data.append(
                [
                    dhdl,
                ]
            )

    if units is None:
        units = formatted_data[0][0].delta_f_.attrs["energy_unit"]

    estimators = formatted_data

    # Get the dF
    dF_list = []
    error_list = []
    max_length = 0
    convert = get_unit_converter(units)
    for dhdl_list in estimators:
        len_dF = sum([len(dhdl.delta_f_) - 1 for dhdl in dhdl_list])
        if len_dF > max_length:
            max_length = len_dF
        dF = []
        error = []
        for dhdl in dhdl_list:
            for i in range(len(dhdl.delta_f_) - 1):
                dF.append(convert(dhdl.delta_f_).iloc[i, i + 1])
                error.append(convert(dhdl.d_delta_f_).iloc[i, i + 1])

        dF_list.append(dF)
        error_list.append(error)

    # Get the determine orientation
    if orientation == "landscape":
        if max_length < 8:
            fig, ax = plt.subplots(figsize=(8, 6))
        else:
            fig, ax = plt.subplots(figsize=(max_length, 6))
        axs = [
            ax,
        ]
        xs = [
            np.arange(max_length),
        ]
    elif orientation == "portrait":
        if max_length < nb:
            xs = [
                np.arange(max_length),
            ]
            fig, ax = plt.subplots(figsize=(8, 6))
            axs = [
                ax,
            ]
        else:
            xs = np.array_split(np.arange(max_length), max_length / nb + 1)
            fig, axs = plt.subplots(nrows=len(xs), figsize=(8, 6))
        mnb = max([len(i) for i in xs])
    else:
        raise ValueError(
            "Not recognising {}, only supports 'landscape' or 'portrait'.".format(
                orientation
            )
        )

    # Sort out the colors
    if colors is None:
        colors_dict = {
            "TI": "#C45AEC",
            "TI-CUBIC": "#33CC33",
            "DEXP": "#F87431",
            "IEXP": "#FF3030",
            "GINS": "#EAC117",
            "GDEL": "#347235",
            "BAR": "#6698FF",
            "UBAR": "#817339",
            "RBAR": "#C11B17",
            "MBAR": "#F9B7FF",
        }
        colors = []
        for dhdl in estimators:
            dhdl = dhdl[0]
            if isinstance(dhdl, TI):
                colors.append(colors_dict["TI"])
            elif isinstance(dhdl, BAR):
                colors.append(colors_dict["BAR"])
            elif isinstance(dhdl, MBAR):
                colors.append(colors_dict["MBAR"])
    else:
        if len(colors) >= len(estimators):
            pass
        else:
            raise ValueError(
                "Number of colors ({}) should be larger than the number of data ({})".format(
                    len(colors), len(estimators)
                )
            )

    # Sort out the labels
    if labels is None:
        labels = []
        for dhdl in estimators:
            dhdl = dhdl[0]
            if isinstance(dhdl, TI):
                labels.append("TI")
            elif isinstance(dhdl, BAR):
                labels.append("BAR")
            elif isinstance(dhdl, MBAR):
                labels.append("MBAR")
    else:
        if len(labels) == len(estimators):
            pass
        else:
            raise ValueError(
                "Length of labels ({}) should be the same as the number of data ({})".format(
                    len(labels), len(estimators)
                )
            )

    # Plot the figure
    width = 1.0 / (len(estimators) + 1)
    elw = 30 * width
    ndx = 1
    for x, ax in zip(xs, axs):
        lines = []
        for i, (dF, error) in enumerate(zip(dF_list, error_list)):
            y = [dF[j] for j in x]
            ye = [error[j] for j in x]
            if orientation == "landscape":
                lw = 0.1 * elw
            elif orientation == "portrait":
                lw = 0.05 * elw
            line = ax.bar(
                x + len(lines) * width,
                y,
                width,
                color=colors[i],
                yerr=ye,
                lw=lw,
                error_kw=dict(elinewidth=elw, ecolor="black", capsize=0.5 * elw),
            )
            lines += (line[0],)
        for dir in ["left", "right", "top", "bottom"]:
            if dir == "left":
                ax.yaxis.set_ticks_position(dir)
            else:
                ax.spines[dir].set_color("none")

        if orientation == "landscape":
            plt.yticks(fontsize=12)
            ax.set_xlim(x[0] - width, x[-1] + len(lines) * width)
            plt.xticks(
                x + 0.5 * width * len(estimators),
                tuple([f"{i}--{i+1}" for i in x]),
                fontsize=12,
            )
        elif orientation == "portrait":
            #plt.yticks(fontsize=10)
            ax.xaxis.set_ticks([])
            for i in x + 0.5 * width * len(estimators):
                ax.annotate(
                    r"$\mathrm{%d-%d}$" % (i, i + 1),
                    xy=(i, 0),
                    xycoords=("data", "axes fraction"),
                    xytext=(0, -2),
                    size=10,
                    textcoords="offset points",
                    va="top",
                    ha="center",
                    rotation=15
                )
            ax.set_xlim(x[0] - width, x[-1] + len(lines) * width + (mnb - len(x)))
        ndx += 1
    x = np.arange(max_length)

    ax = plt.gca()

    for tick in ax.get_xticklines():
        tick.set_visible(False)
    if orientation == "landscape":
        leg = plt.legend(lines, labels, loc=3, ncol=2, prop=FP(size=18), fancybox=True)
        plt.title("The free energy change breakdown", fontsize=18)
        plt.xlabel("States", fontsize=18, color="#151B54")
        plt.ylabel(r"$\Delta G$ ({})".format(units), fontsize=18, color="#151B54")
    elif orientation == "portrait":
        leg = ax.legend(
            lines,
            labels,
            loc=0,
            ncol=2,
            #title=r"$\Delta G$ ({})".format(units),
            prop=FP(size=12),
            fancybox=True,
        )
        #plt.ylabel(r"$\Delta G$ ({})".format(units), fontsize=18, color="#151B54")
        fig.supylabel(r"$\Delta G$ ({})".format(units), fontsize=18, color="#151B54")

    leg.get_frame().set_alpha(0.5)
    return fig

# =============================================================================
# Command-Line Argument Parsing
# =============================================================================
def parse_args():
    parser = argparse.ArgumentParser(
        description="FEP analysis script that runs free energy estimations and generates diagnostic plots."
    )
    parser.add_argument(
        "-i", "--inputdir", required=True,
        help="Input directory containing simulation output files."
    )
    parser.add_argument(
        "-o", "--outputdir", required=True,
        help="Output directory for saving images and summary results."
    )
    parser.add_argument(
        "--temp", type=float, default=310,
        help="Simulation temperature in Kelvin (default: 310)."
    )
    parser.add_argument(
        "--units", default="kcal/mol",
        help="Energy units for reporting free energy results (default: kcal/mol)."
    )
    parser.add_argument(
        "--skiptime", type=float, default=300,
        help="Skiptime in ps."
    )
    parser.add_argument(
        "--raw-only", action="store_true",
        help="If set, skip all decorrelation and use only raw data."
    )

    return parser.parse_args()

# =============================================================================
# Main analysis function
# =============================================================================
def main():
    args = parse_args()
    raw_only = args.raw_only
    
    # Use command-line arguments for the base parameters.
    dir_input = args.inputdir
    dir_output = args.outputdir
    temperature = args.temp
    energy_units = args.units

    # Default file patterns and analysis options.
    prefix_pattern = 'dhdl'      # File prefix to match
    suffix_pattern = 'xvg'       # File suffix (extension)
    skiptime = args.skiptime     # Time to skip as equilibration, defaults to 300
    uncorr = 'dhdl'              # Observable for autocorrelation analysis
    threshold = 50               # Minimal number of decorrelated samples required
    estimators_to_use = ("MBAR", "BAR", "TI")  # Which estimators to run
    convergence_num_points = 10  # Number of points for convergence analysis

    os.makedirs(dir_output, exist_ok=True)
    log_file_path = os.path.join(dir_output, "analyse.log")
    setup_logging(log_file_path)
    logging.getLogger().setLevel(logging.INFO)
    logging.info(f"Logging initialized. Log file: {log_file_path}")

    # =============================================================================
    # 4. Discover data files using pathlib globbing
    # =============================================================================
    logging.info("Scanning for input files.")
    pattern = "**/" + prefix_pattern + "*" + suffix_pattern
    file_list = list(map(str, Path(dir_input).glob(pattern)))
    if not file_list:
        raise ValueError(f"No files matched the pattern '{pattern}' under directory {dir_input}")
    logging.info(f"Found {len(file_list)} files:")
    for f in file_list:
        logging.info(f" - {f}")

    # =============================================================================
    # 5. Select parser functions (assuming GROMACS)
    # =============================================================================
    _extract_u_nk = gmx.extract_u_nk
    _extract_dHdl = gmx.extract_dHdl

    # =============================================================================
    # 6. Define helper functions for reading data from files
    # =============================================================================
    def read_u_nk(file, T):
        try:
            u_nk = _extract_u_nk(file, T)
            logging.info(f"Read {len(u_nk)} lines of u_nk from {file}")
            return u_nk
        except Exception as exc:
            msg = f"Error reading u_nk from {file}"
            logging.error(msg)
            raise OSError(msg) from exc

    def read_dHdl(file, T):
        try:
            dhdl = _extract_dHdl(file, T)
            logging.info(f"Read {len(dhdl)} lines of dHdl from {file}")
            return dhdl
        except Exception as exc:
            msg = f"Error reading dHdl from {file}"
            logging.error(msg)
            raise OSError(msg) from exc

    # Read the data files in parallel.
    n_jobs = 1  # Set to -1 to use all available cores.
    u_nk_list = joblib.Parallel(n_jobs=n_jobs)(
        joblib.delayed(read_u_nk)(file, temperature) for file in file_list
    )
    dHdl_list = joblib.Parallel(n_jobs=n_jobs)(
        joblib.delayed(read_dHdl)(file, temperature) for file in file_list
    )

    # =============================================================================
    # 7. Sort the file lists (and corresponding data lists) by state order.
    # =============================================================================
    if u_nk_list:
        sort_keys = [u_nk.reset_index("time").index.values[0] for u_nk in u_nk_list]
        sorted_indices = sorted(range(len(file_list)), key=lambda i: sort_keys[i])
    elif dHdl_list:
        sort_keys = [d.reset_index("time").index.values[0] for d in dHdl_list]
        sorted_indices = sorted(range(len(file_list)), key=lambda i: sort_keys[i])
    else:
        raise ValueError("No data files were successfully read.")
    file_list = [file_list[i] for i in sorted_indices]
    u_nk_list = [u_nk_list[i] for i in sorted_indices]
    dHdl_list = [dHdl_list[i] for i in sorted_indices]
    logging.info("Files reordered according to state order.")

    # =============================================================================
    # 8. Preprocessing: Remove equilibration and decorrelate the data.
    # =============================================================================
    logging.info(f"Preprocessing data (skiptime={skiptime}, uncorr={uncorr}, threshold={threshold}, raw_only={raw_only}).")
    u_nk_sample_list = []
    if u_nk_list:
        for index, df in enumerate(u_nk_list):
            df_proc = df[df.index.get_level_values("time") >= skiptime]
            if raw_only:
                logging.info(f"State {index}: raw‑only flag set, skipping u_nk decorrelation.")
                u_nk_sample_list.append(df_proc)
            else:
                subsample = decorrelate_u_nk(df_proc, uncorr, remove_burnin=True)
                if len(subsample) < threshold:
                    logging.warning(f"u_nk for state {index} has only {len(subsample)} uncorrelated samples (< {threshold}). Using full data.")
                    subsample = df_proc
                else:
                    logging.info(f"State {index}: using {len(subsample)} uncorrelated u_nk samples.")
                u_nk_sample_list.append(subsample)
    else:
        logging.info("No u_nk data available for decorrelation.")

    dHdl_sample_list = []
    if dHdl_list:
        for index, df in enumerate(dHdl_list):
            df_proc = df[df.index.get_level_values("time") >= skiptime]
            if raw_only:
                logging.info(f"State {index}: raw‑only flag set, skipping dHdl decorrelation.")
                dHdl_sample_list.append(df_proc)
            else:
                subsample = decorrelate_dhdl(df_proc, remove_burnin=True)
                if len(subsample) < threshold:
                    logging.warning(f"dHdl for state {index} has only {len(subsample)} uncorrelated samples (< {threshold}). Using full data.")
                    subsample = df_proc
                else:
                    logging.info(f"State {index}: using {len(subsample)} uncorrelated dHdl samples.")
                dHdl_sample_list.append(subsample)
    else:
        logging.info("No dHdl data available for decorrelation.")

    # =============================================================================
    # 9. Concatenate data for estimation.
    # =============================================================================
    if u_nk_sample_list:
        u_nk_concat = concat(u_nk_sample_list)
        logging.info(f"Concatenated decorrelated u_nk data: {len(u_nk_concat)} lines.")
    else:
        u_nk_concat = concat(u_nk_list)
        logging.warning("Using original u_nk data as decorrelated data is not available.")

    if dHdl_sample_list:
        dHdl_concat = concat(dHdl_sample_list)
        logging.info(f"Concatenated decorrelated dHdl data: {len(dHdl_concat)} lines.")
    else:
        dHdl_concat = concat(dHdl_list)
        logging.warning("Using original dHdl data as decorrelated data is not available.")

    u_nk_raw_concat = concat(u_nk_list)
    logging.info(f"Concatenated raw u_nk data: {len(u_nk_raw_concat)} lines.")
    dHdl_raw_concat = concat(dHdl_list)
    logging.info(f"Concatenated raw dHdl data: {len(dHdl_raw_concat)} lines.")

    # =============================================================================
    # 10. Estimation: Compute free energy differences using selected estimators.
    # =============================================================================
    # Estimation for decorrelated data.
    estimators_results = {}
    for estimator in estimators_to_use:
        if estimator in FEP_ESTIMATORS:
            if estimator == "MBAR":
                logging.info("Running MBAR estimator on decorrelated data.")
                estimators_results["MBAR"] = MBAR(**MBAR_OPTS).fit(u_nk_concat)
            elif estimator == "BAR":
                logging.info("Running BAR estimator on decorrelated data.")
                estimators_results["BAR"] = BAR().fit(u_nk_concat)
            else:
                raise ValueError(f"FEP estimator {estimator} is not recognized.")
        elif estimator in TI_ESTIMATORS:
            if estimator == "TI":
                logging.info("Running TI estimator on decorrelated data.")
                estimators_results["TI"] = TI().fit(dHdl_concat)
            else:
                raise ValueError(f"TI estimator {estimator} is not recognized.")
        else:
            msg = f"Estimator {estimator} is not supported."
            logging.error(msg)
            raise ValueError(msg)

    # Log lambda mapping for decorrelated data.
    ref_estimator = estimators_results.get("MBAR", next(iter(estimators_results.values())))
    num_states = len(ref_estimator.states_)
    logging.info("Lambda values for states (decorrelated):")
    for i, state in enumerate(ref_estimator.states_):
        logging.info(f"  State {i}: Lambda = {state}")
    for i in range(num_states - 1):
        logging.info(f"  Transition {i}--{i+1}: Lambda from {ref_estimator.states_[i]} to {ref_estimator.states_[i+1]}")

    # Estimation for raw data.
    estimators_results_raw = {}
    for estimator in estimators_to_use:
        if estimator in FEP_ESTIMATORS:
            if estimator == "MBAR":
                logging.info("Running MBAR estimator on raw data.")
                estimators_results_raw["MBAR"] = MBAR(**MBAR_OPTS).fit(u_nk_raw_concat)
            elif estimator == "BAR":
                logging.info("Running BAR estimator on raw data.")
                estimators_results_raw["BAR"] = BAR().fit(u_nk_raw_concat)
            else:
                raise ValueError(f"FEP estimator {estimator} is not recognized.")
        elif estimator in TI_ESTIMATORS:
            if estimator == "TI":
                logging.info("Running TI estimator on raw data.")
                estimators_results_raw["TI"] = TI().fit(dHdl_raw_concat)
            else:
                raise ValueError(f"TI estimator {estimator} is not recognized.")
        else:
            msg = f"Estimator {estimator} is not supported."
            logging.error(msg)
            raise ValueError(msg)

    ref_estimator_raw = estimators_results_raw.get("MBAR", next(iter(estimators_results_raw.values())))
    num_states_raw = len(ref_estimator_raw.states_)
    logging.info("Lambda values for states (raw data):")
    for i, state in enumerate(ref_estimator_raw.states_):
        logging.info(f"  State {i}: Lambda = {state}")
    for i in range(num_states_raw - 1):
        logging.info(f"  Transition {i}--{i+1}: Lambda from {ref_estimator_raw.states_[i]} to {ref_estimator_raw.states_[i+1]}")

    # =============================================================================
    # 11. Generate summary DataFrame for decorrelated data.
    # =============================================================================
    summary_dict = {"name": [], "state": []}
    for i in range(num_states - 1):
        summary_dict["name"].append(f"{i} -- {i+1}")
        summary_dict["state"].append("States")
    try:
        stages = u_nk_list[0].reset_index("time").index.names
        logging.info("Stage names taken from u_nk data.")
    except Exception:
        stages = dHdl_list[0].reset_index("time").index.names
        logging.info("Stage names taken from dHdl data.")
    for stage in stages:
        summary_dict["name"].append(stage.split("-")[0])
        summary_dict["state"].append("Stages")
    summary_dict["name"].append("TOTAL")
    summary_dict["state"].append("Stages")

    col_order = []
    for estimator_name, estimator_obj in estimators_results.items():
        logging.info(f"Assembling results from estimator: {estimator_name} (decorrelated)")
        summary_dict[estimator_name] = []
        summary_dict[f"{estimator_name}_Error"] = []
        col_order.extend([estimator_name, f"{estimator_name}_Error"])
        delta_f_ = estimator_obj.delta_f_
        d_delta_f_ = estimator_obj.d_delta_f_
        for j in range(1, num_states):
            summary_dict[estimator_name].append(delta_f_.iloc[j - 1, j])
            summary_dict[f"{estimator_name}_Error"].append(d_delta_f_.iloc[j - 1, j])
        for idx, stage in enumerate(stages):
            if len(stages) == 1:
                start, end = 0, num_states - 1
            else:
                state_values = [state[idx] for state in estimator_obj.states_]
                lambda_min, lambda_max = min(state_values), max(state_values)
                if lambda_min == lambda_max:
                    start, end = 0, 0
                else:
                    start = num_states - list(reversed(state_values)).index(lambda_min) - 1
                    end = state_values.index(lambda_max)
            result = delta_f_.iloc[start, end]
            if estimator_name != "BAR":
                error = d_delta_f_.iloc[start, end]
            else:
                error = np.sqrt(sum(d_delta_f_.iloc[i, i + 1] ** 2 for i in range(start, end)))
            summary_dict[estimator_name].append(result)
            summary_dict[f"{estimator_name}_Error"].append(error)
        result_total = delta_f_.iloc[0, -1]
        if estimator_name != "BAR":
            error_total = d_delta_f_.iloc[0, -1]
        else:
            error_total = np.sqrt(sum(d_delta_f_.iloc[i, i + 1] ** 2 for i in range(num_states - 1)))
        summary_dict[estimator_name].append(result_total)
        summary_dict[f"{estimator_name}_Error"].append(error_total)

    summary_df = pd.DataFrame.from_dict(summary_dict)
    summary_df = summary_df.set_index(["state", "name"])
    summary_df = summary_df[col_order]
    summary_df.attrs["temperature"] = temperature
    summary_df.attrs["energy_unit"] = 'kT'
    converter = get_unit_converter(energy_units)
    summary_df = converter(summary_df)
    logging.info("Free energy summary (decorrelated):")
    logging.info(summary_df.to_string())
    summary_csv_file = os.path.join(dir_output, "free_energy_summary.csv")
    summary_df.to_csv(summary_csv_file)
    logging.info(f"Summary (decorrelated) saved to {summary_csv_file}")

    # =============================================================================
    # 12. Generate diagnostic plots for decorrelated data.
    # =============================================================================
    if not raw_only:
        # 12a. MBAR Overlap Matrix.
        if "MBAR" in estimators_results:
            logging.info("Plotting MBAR overlap matrix (decorrelated).")
            ax_overlap = plot_mbar_overlap_matrix(estimators_results["MBAR"].overlap_matrix)
            overlap_file = os.path.join(dir_output, "O_MBAR.pdf")
            ax_overlap.figure.savefig(overlap_file)
            plt.close(ax_overlap.figure)
            logging.info(f"Overlap matrix (decorrelated) saved to {overlap_file}")
        else:
            logging.warning("MBAR estimator not available for decorrelated data; skipping overlap plot.")
        
        # 12b. TI dH/dλ Plot.
        if "TI" in estimators_results:
            logging.info("Plotting TI dH/dλ curve (decorrelated).")
            # Here, we pass the TI estimator object; its separate_dhdl() method will be used.
            ax_ti = plot_ti_dhdl(estimators_results["TI"], units=energy_units)
            ti_file = os.path.join(dir_output, "dhdl_TI.pdf")
            ax_ti.figure.tight_layout()
            ax_ti.figure.savefig(ti_file)
            plt.close(ax_ti.figure)
            logging.info(f"TI dH/dλ plot (decorrelated) saved to {ti_file}")
        else:
            logging.warning("TI estimator not available for decorrelated data; skipping TI dH/dλ plot.")
        
        # 12c. Free Energy State Differences.
        logging.info("Plotting free energy state differences (decorrelated).")
        print(estimators_results.values())
        print(estimators_results)
        fig_df = plot_dF_state(estimators_results.values(), units=energy_units, orientation="portrait", nb=10)
        df_state_file = os.path.join(dir_output, "dF_state.pdf")
        fig_df.tight_layout()
        fig_df.savefig(df_state_file)
        plt.close(fig_df)
        logging.info(f"Free energy state differences (portrait, decorrelated) saved to {df_state_file}")
        fig_df_land = plot_dF_state(estimators_results.values(), units=energy_units, orientation="landscape", nb=10)
        df_state_long_file = os.path.join(dir_output, "dF_state_long.pdf")
        fig_df_land.tight_layout()
        fig_df_land.savefig(df_state_long_file)
        plt.close(fig_df_land)
        logging.info(f"Free energy state differences (landscape, decorrelated) saved to {df_state_long_file}")
        
        # 12d. Convergence Analysis.

        if "MBAR" in estimators_results:
            data_for_conv = u_nk_sample_list if u_nk_sample_list else u_nk_list
            conv_estimator = "MBAR"
        elif "TI" in estimators_results:
            data_for_conv = dHdl_sample_list if dHdl_sample_list else dHdl_list
            conv_estimator = "TI"
        else:
            raise ValueError("No suitable estimator available for decorrelated convergence analysis.")
        
        #------------override for testing purposes
        #data_for_conv = dHdl_sample_list if dHdl_sample_list else dHdl_list
        #conv_estimator = "TI"
        #------------override for testing purposes

        logging.info(f"Performing convergence analysis (decorrelated) with {conv_estimator}.")
        if conv_estimator == "MBAR":
            convergence_df = forward_backward_convergence(data_for_conv, estimator=conv_estimator, num=convergence_num_points, **MBAR_OPTS)
        elif conv_estimator == "TI":
            convergence_df = forward_backward_convergence(data_for_conv, estimator=conv_estimator, num=convergence_num_points)
        else:
            raise ValueError(f"Estimator {conv_estimator} is not recognized for convergence analysis.")
        logging.info("Convergence analysis results:")
        logging.info(convergence_df.to_string())
        converted_conv = get_unit_converter(energy_units)(convergence_df)
        converted_conv["data_fraction"] = convergence_df["data_fraction"]
        ax_conv = plot_convergence(converted_conv, units=energy_units)
        conv_file = os.path.join(dir_output, "dF_t.pdf")
        ax_conv.figure.tight_layout()
        ax_conv.figure.savefig(conv_file)
        plt.close(ax_conv.figure)
        logging.info(f"Convergence plot (decorrelated) saved to {conv_file}")

    # =============================================================================
    # 13. Raw (correlated) Data: Repeat summary and plots for raw data.
    # =============================================================================
    summary_dict_raw = {"name": [], "state": []}
    for i in range(num_states_raw - 1):
        summary_dict_raw["name"].append(f"{i} -- {i+1}")
        summary_dict_raw["state"].append("States")
    try:
        stages = u_nk_list[0].reset_index("time").index.names
        logging.info("Stage names taken from u_nk data for raw analysis.")
    except Exception:
        stages = dHdl_list[0].reset_index("time").index.names
        logging.info("Stage names taken from dHdl data for raw analysis.")
    for stage in stages:
        summary_dict_raw["name"].append(stage.split("-")[0])
        summary_dict_raw["state"].append("Stages")
    summary_dict_raw["name"].append("TOTAL")
    summary_dict_raw["state"].append("Stages")
    
    col_order_raw = []
    for estimator_name, estimator_obj in estimators_results_raw.items():
        logging.info(f"Assembling results from estimator: {estimator_name} (raw data)")
        summary_dict_raw[estimator_name] = []
        summary_dict_raw[f"{estimator_name}_Error"] = []
        col_order_raw.extend([estimator_name, f"{estimator_name}_Error"])
        delta_f_ = estimator_obj.delta_f_
        d_delta_f_ = estimator_obj.d_delta_f_
        for j in range(1, num_states_raw):
            summary_dict_raw[estimator_name].append(delta_f_.iloc[j - 1, j])
            summary_dict_raw[f"{estimator_name}_Error"].append(d_delta_f_.iloc[j - 1, j])
        for idx, stage in enumerate(stages):
            if len(stages) == 1:
                start, end = 0, num_states_raw - 1
            else:
                state_values = [state[idx] for state in estimator_obj.states_]
                lambda_min, lambda_max = min(state_values), max(state_values)
                if lambda_min == lambda_max:
                    start, end = 0, 0
                else:
                    start = num_states_raw - list(reversed(state_values)).index(lambda_min) - 1
                    end = state_values.index(lambda_max)
            result = delta_f_.iloc[start, end]
            if estimator_name != "BAR":
                error = d_delta_f_.iloc[start, end]
            else:
                error = np.sqrt(sum(d_delta_f_.iloc[i, i + 1] ** 2 for i in range(start, end)))
            summary_dict_raw[estimator_name].append(result)
            summary_dict_raw[f"{estimator_name}_Error"].append(error)
        result_total = delta_f_.iloc[0, -1]
        if estimator_name != "BAR":
            error_total = d_delta_f_.iloc[0, -1]
        else:
            error_total = np.sqrt(sum(d_delta_f_.iloc[i, i + 1] ** 2 for i in range(num_states_raw - 1)))
        summary_dict_raw[estimator_name].append(result_total)
        summary_dict_raw[f"{estimator_name}_Error"].append(error_total)
    
    summary_df_raw = pd.DataFrame.from_dict(summary_dict_raw)
    summary_df_raw = summary_df_raw.set_index(["state", "name"])
    summary_df_raw = summary_df_raw[col_order_raw]
    summary_df_raw.attrs["temperature"] = temperature
    summary_df_raw.attrs["energy_unit"] = 'kT'
    converter = get_unit_converter(energy_units)
    summary_df_raw = converter(summary_df_raw)
    logging.info("Free energy summary (raw data):")
    logging.info(summary_df_raw.to_string())
    summary_csv_file_raw = os.path.join(dir_output, "free_energy_summary_raw.csv")
    summary_df_raw.to_csv(summary_csv_file_raw)
    logging.info(f"Summary (raw data) saved to {summary_csv_file_raw}")
    
    # Raw data diagnostic plots.
    if "MBAR" in estimators_results_raw:
        logging.info("Plotting MBAR overlap matrix (raw data).")
        ax_overlap_raw = plot_mbar_overlap_matrix(estimators_results_raw["MBAR"].overlap_matrix)
        overlap_file_raw = os.path.join(dir_output, "O_MBAR_raw.pdf")
        ax_overlap_raw.figure.savefig(overlap_file_raw)
        plt.close(ax_overlap_raw.figure)
        logging.info(f"Overlap matrix (raw data) saved to {overlap_file_raw}")
    else:
        logging.warning("MBAR estimator not available for raw data; skipping overlap plot.")
    
    if "TI" in estimators_results_raw:
        logging.info("Plotting TI dH/dλ curve (raw data).")
        ax_ti_raw = plot_ti_dhdl(estimators_results_raw["TI"], units=energy_units)
        ax_ti_raw.figure.tight_layout()
        ti_file_raw = os.path.join(dir_output, "dhdl_TI_raw.pdf")
        ax_ti_raw.figure.tight_layout()
        ax_ti_raw.figure.savefig(ti_file_raw)
        plt.close(ax_ti_raw.figure)
        logging.info(f"TI dH/dλ plot (raw data) saved to {ti_file_raw}")
    else:
        logging.warning("TI estimator not available for raw data; skipping TI dH/dλ plot.")
    
    logging.info("Plotting free energy state differences (raw data).")
    fig_df_raw = plot_dF_state(estimators_results_raw.values(), units=energy_units, orientation="portrait", nb=10)
    df_state_file_raw = os.path.join(dir_output, "dF_state_raw.pdf")
    fig_df_raw.tight_layout()
    fig_df_raw.savefig(df_state_file_raw)
    plt.close(fig_df_raw)
    logging.info(f"Free energy state differences (portrait, raw data) saved to {df_state_file_raw}")
    fig_df_raw_land = plot_dF_state(estimators_results_raw.values(), units=energy_units, orientation="landscape", nb=10)
    df_state_long_file_raw = os.path.join(dir_output, "dF_state_long_raw.pdf")
    fig_df_raw_land.tight_layout()
    fig_df_raw_land.savefig(df_state_long_file_raw)
    plt.close(fig_df_raw_land)
    logging.info(f"Free energy state differences (landscape, raw data) saved to {df_state_long_file_raw}")
    
    if "MBAR" in estimators_results_raw:
        data_for_conv_raw = u_nk_list
        conv_estimator_raw = "MBAR"
    elif "TI" in estimators_results_raw:
        data_for_conv_raw = dHdl_list
        conv_estimator_raw = "TI"
    else:
        raise ValueError("No suitable estimator available for raw convergence analysis.")
    
    #------------override for testing purposes
    #data_for_conv_raw = dHdl_list
    #conv_estimator_raw = "TI"
    #------------override for testing purposes

    logging.info(f"Performing convergence analysis (raw data) with {conv_estimator_raw}.")
    if conv_estimator_raw == "MBAR":
        convergence_df_raw = forward_backward_convergence(data_for_conv_raw, estimator=conv_estimator_raw, num=convergence_num_points, **MBAR_OPTS)
    elif conv_estimator_raw == "TI":
        convergence_df_raw = forward_backward_convergence(data_for_conv_raw, estimator=conv_estimator_raw, num=convergence_num_points)
    else:
        raise ValueError(f"Estimator {conv_estimator_raw} is not recognized for convergence analysis.")
    logging.info("Convergence analysis results:")
    logging.info(convergence_df_raw.to_string())

    converted_conv_raw = get_unit_converter(energy_units)(convergence_df_raw)
    converted_conv_raw["data_fraction"] = convergence_df_raw["data_fraction"]
    ax_conv_raw = plot_convergence(converted_conv_raw, units=energy_units)
    conv_file_raw = os.path.join(dir_output, "dF_t_raw.pdf")
    ax_conv_raw.figure.tight_layout()
    ax_conv_raw.figure.savefig(conv_file_raw)
    plt.close(ax_conv_raw.figure)
    logging.info(f"Convergence plot (raw data) saved to {conv_file_raw}")
    
    # =============================================================================
    # 14. End of Workflow
    # =============================================================================
    logging.info("Free energy analysis workflow complete.")

if __name__ == '__main__':
    main()