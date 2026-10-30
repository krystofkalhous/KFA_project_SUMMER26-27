"""
Maps of the GAM teleconnection results, in the style of Skala (2025).

Adapted from Petr Skala's grid_plot_discrete2.py (discrete_dataset_plot) using Claude Opus 5.5: Robinson projection,
discrete colour classes, hatching where there is not enough data.

All statistics are computed in R (section 7 of KFA_project_final_defin): the BH correction, the
significance classes on Skala's sigma scale and the agreement with his classes. R writes each field
as a headerless 90 x 180 CSV in Skala's format (rows = GISTEMP latitudes, columns = longitudes,
NaN = no result) to maps/; this script only draws them, so every number on a map is the same as in
the R summary. Skala's own classes are read from the digitised files in Data/.

Maps (one PDF each, for NAO and NINO3.4), all in maps/:
  1. nonlinearity on Skala's scale (classes 0-1 ... >4 sigma), three maps drawn identically:
     Skala's own result (his figs. 5.4 and 5.2, bottom panels, redrawn from the digitised file),
     our pooled model (his set-up: all months pooled) and our main model (after the seasonality
     correction). Uncorrected per-cell values on all three, as in his maps.
  2. agreement with Skala (>= 2 sigma on both sides): neither / GAM only / Skala only / both,
     for the pooled and for the main model.
  3. nonlinear edf (0 = straight line), pooled and main model; grey where BH q >= ALPHA.
  4. bend of the main model's curve (K); grey where BH q >= ALPHA.
  5. strength in K per 4 SD (the scale of Miksovsky et al. 2016): the straight line through all
     months (grey where BH q >= ALPHA; comparable with his Fig. 5 - for ENSO he used the SOI, so
     his signs are flipped) and one straight line per season, 2 x 2 panels (grey where the joint
     test of "no dependence in any season" is not rejected).
  6. the AR(1) check of the main model: lag-1 autocorrelation left after whitening.

AR(1) marker (switch AR1_MARK below, or the environment variable KFA_AR1_MARK=0 to turn it off): on
the main model's nonlinearity maps (1-4), a black dot marks every cell that the map shows as a
result (>= 2 sigma, flagged by the GAM, or BH q < ALPHA) but where the AR(1) check of the main model
fails (|lag-1 autocorrelation left| > AR1_LIMIT, default 0.1; KFA_AR1_LIMIT=0.05 marks every failure). There the autocorrelation is not fully removed, so
the p-value is probably too small: the reader sees directly which results rest on it.

Requires: numpy matplotlib xarray netCDF4 cartopy
    conda install -c conda-forge cartopy xarray netcdf4
The first run downloads Natural Earth coastlines, so it needs internet once. Without internet,
set the environment variable KFA_NO_COAST=1 to draw the maps without coastlines.
"""

import os
import numpy as np
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.colors import BoundaryNorm, ListedColormap
from matplotlib.lines import Line2D
import xarray as xr
import cartopy.crs as ccrs

# ---- Paths and settings -----------------------------------------------------------
BASE    = os.environ.get("KFA_BASE", "/Users/krystofkalhous/Desktop/Projekt KFA")
DATA    = os.path.join(BASE, "Data")
GISTEMP = os.path.join(DATA, "gistemp1200_GHCNv4_ERSSTv5.nc")
MAPDIR  = os.path.join(BASE, "maps")       # the fields from R (section 7); the PDFs are saved here too
SKALA   = {"nao":  os.path.join(DATA, "skala_nao_nonlin_class.csv"),
           "nino": os.path.join(DATA, "skala_nino34_nonlin_class.csv")}

ALPHA      = 0.05                          # as ALPHA in the R script
SKALA_SIG  = 2                             # as SKALA_SIG in the R script (>= 2 sigma = significant)
COASTLINES = os.environ.get("KFA_NO_COAST") is None
AR1_MARK   = os.environ.get("KFA_AR1_MARK", "1") != "0"   # True: dot the results where the AR(1) check fails
AR1_LIMIT  = float(os.environ.get("KFA_AR1_LIMIT", "0.1"))  # dot where |lag-1 autocorrelation left| > this:
                                          # 0.1 = clear failures only; 0.05 (= 2/sqrt(1728)) = every failure
NAMES      = {"nao": "NAO", "nino": "NINO3.4"}
SEASONS    = ["DJF", "MAM", "JJA", "SON"]


# ---- Grid axes straight from GISTEMP, as Skala does ----------------------------------
with xr.open_dataset(GISTEMP) as ds:
    lats = ds["lat"].values
    lons = ds["lon"].values
lon_grid, lat_grid = np.meshgrid(lons, lats)          # both 90 x 180


def load(name):
    """A field written by R: maps/gam_<name>.csv -> 90 x 180 array (NaN = no result)."""
    path = os.path.join(MAPDIR, f"gam_{name}.csv")
    if not os.path.exists(path):
        raise FileNotFoundError(f"{path} not found - run section 7 of the R script first")
    F = np.loadtxt(path, delimiter=",")
    assert F.shape == (len(lats), len(lons)), (path, F.shape)
    return F


def not_sig(F, Q):
    """Cells with a result where BH q >= ALPHA (drawn light grey)."""
    return ~np.isnan(F) & ~(np.nan_to_num(Q, nan=1.0) < ALPHA)


RHO_LEFT = load("rho_left_sw") if AR1_MARK else None


def ar1_mark(shown):
    """The cells among `shown` (the results a map highlights) where the AR(1) check of the main
    model fails; None if the marker is switched off."""
    if not AR1_MARK:
        return None
    return shown & (np.abs(np.nan_to_num(RHO_LEFT, nan=0.0)) > AR1_LIMIT)


# ---- Drawing (Skala's discrete_dataset_plot, generalised) ------------------------------
def colormap(cmap, n):
    """n discrete classes from a matplotlib colour map name, or a list of n colours."""
    cm = mpl.colormaps[cmap].resampled(n) if isinstance(cmap, str) else ListedColormap(cmap)
    cm.set_bad(color="white", alpha=1.0)
    return cm


def draw(ax, F, bins, cmap, grey_mask=None, mark=None):
    """One field on a Robinson map; mark: cells to dot (AR(1) check failed). Returns the colour
    mesh (for the colour bar)."""
    num_bins = len(bins) - 1
    nan_mask = np.isnan(F)
    hide = nan_mask if grey_mask is None else (nan_mask | grey_mask)
    c = ax.pcolormesh(lon_grid, lat_grid, np.ma.masked_where(hide, F),
                      cmap=colormap(cmap, num_bins), norm=BoundaryNorm(bins, num_bins, clip=True),
                      shading="auto", transform=ccrs.PlateCarree())
    if grey_mask is not None:
        grey = np.ma.masked_where(~grey_mask, np.ones_like(F))
        ax.pcolormesh(lon_grid, lat_grid, grey, cmap=ListedColormap(["0.85"]),
                      shading="auto", transform=ccrs.PlateCarree())
    # Hatching where there is no result. Skala uses alpha=0.0 with edgecolor="none", which stopped
    # drawing the hatch in matplotlib 3.11; facecolor="none" with an explicit edgecolor renders on both.
    ax.pcolor(lon_grid, lat_grid, np.ma.masked_where(~nan_mask, nan_mask), hatch="///",
              facecolor="none", edgecolor="black", linewidth=0.0, transform=ccrs.PlateCarree())
    if mark is not None and mark.any():
        ax.scatter(lon_grid[mark], lat_grid[mark], s=2.5, c="black", marker="o", linewidths=0,
                   transform=ccrs.PlateCarree(), zorder=3)
    if COASTLINES:
        ax.coastlines(linewidth=0.6)
    ax.set_global()
    return c


def colorbar(c, cax, bins, label, tick_labels=None, extend_low=False):
    """Colour bar. Categories (tick_labels given, e.g. Skala's sigma classes): one label per class,
    at its centre. Numeric classes: labels at the class edges, so each colour reads as the interval
    between its two labels; the outermost classes are open-ended (the top one always, the bottom one
    if extend_low), so their outer edge gets no label."""
    if tick_labels is not None:
        centers = [(bins[k] + bins[k + 1]) / 2 for k in range(len(bins) - 1)]
        cb = plt.colorbar(c, cax=cax, orientation="vertical", ticks=centers)
        cb.ax.set_yticklabels(tick_labels)
    else:
        edges = list(bins[1:-1]) if extend_low else list(bins[:-1])
        cb = plt.colorbar(c, cax=cax, orientation="vertical", ticks=edges)
        cb.ax.set_yticklabels([f"{b:g}" for b in edges])
    cb.set_label(label)


def legend_handles(grey_label=None, mark=False):
    h = [plt.Rectangle((0, 0), 40, 40, facecolor="white", hatch="///", edgecolor="black",
                       linewidth=1.0, label="Not enough data")]
    if grey_label:
        h.append(plt.Rectangle((0, 0), 40, 40, facecolor="0.85", edgecolor="black",
                               linewidth=1.0, label=grey_label))
    if mark:
        h.append(Line2D([], [], marker="o", markersize=2.5, color="black", linestyle="none",
                        label=f"AR(1) check failed (> {AR1_LIMIT:g})"))
    return h


def save(fig, name):
    path = os.path.join(MAPDIR, name)
    fig.savefig(path, format="pdf", bbox_inches="tight")
    plt.close(fig)
    print("saved", path)


def map_plot(F, bins, label, save_name, cmap="YlGnBu", grey_mask=None, grey_label=None,
             tick_labels=None, extend_low=False, title=None, figsize=(12, 6), mark=None):
    """One map with its colour bar and legend, as in Skala's figures. mark: cells to dot."""
    fig = plt.figure(figsize=figsize)
    gs = gridspec.GridSpec(1, 2, width_ratios=[20, 1], wspace=0.05)
    ax = plt.subplot(gs[0], projection=ccrs.Robinson())
    c = draw(ax, F, bins, cmap, grey_mask, mark)
    if title:
        ax.set_title(title, loc="left", fontsize="medium")
    colorbar(c, plt.subplot(gs[1]), bins, label, tick_labels, extend_low)
    handles = legend_handles(grey_label if grey_mask is not None else None, mark is not None)
    ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.88, 1.12 + 0.05 * (len(handles) - 2)),
              frameon=False, fontsize="medium",
              handlelength=3.2, handletextpad=0.5, borderpad=0.8, labelspacing=0.6)
    save(fig, save_name)


def season_panels(Fs, bins, label, save_name, cmap, grey_masks, grey_label, title,
                  extend_low=False):
    """The four seasons (DJF, MAM, JJA, SON) as 2 x 2 maps with one colour bar."""
    fig = plt.figure(figsize=(16, 9))
    gs = gridspec.GridSpec(2, 3, width_ratios=[20, 20, 1], wspace=0.05, hspace=0.15)
    for k, (F, gm, sn) in enumerate(zip(Fs, grey_masks, SEASONS)):
        ax = plt.subplot(gs[k // 2, k % 2], projection=ccrs.Robinson())
        c = draw(ax, F, bins, cmap, gm)
        ax.set_title(sn, loc="left", fontsize="medium")
    colorbar(c, plt.subplot(gs[:, 2]), bins, label, extend_low=extend_low)
    fig.suptitle(title, x=0.12, ha="left", fontsize="large")
    fig.legend(handles=legend_handles(grey_label), loc="upper right", ncol=2, frameon=False)
    save(fig, save_name)


# ---- 1. Nonlinearity on Skala's scale -------------------------------------------------
SIGMA_BINS   = [-0.5, 0.5, 1.5, 2.5, 3.5, 4.5]
SIGMA_LABELS = ["0-1", "1-2", "2-3", "3-4", ">4"]
SIGMA_LABEL  = "Significance level ($\\sigma$)"

for tag, name in NAMES.items():
    if os.path.exists(SKALA[tag]):
        map_plot(np.loadtxt(SKALA[tag], delimiter=","), SIGMA_BINS, SIGMA_LABEL,
                 f"skala_{tag}_nonlin.pdf", tick_labels=SIGMA_LABELS,
                 title=f"{name}: Skala (2025), nonlinear contribution (MI, IAAFT surrogates)")
    map_plot(load(f"sigma_{tag}"), SIGMA_BINS, SIGMA_LABEL, f"gam_{tag}_sigma_pooled.pdf",
             tick_labels=SIGMA_LABELS,
             title=f"{name}: GAM, all months pooled (Skala (2025) set-up), nonlinear part")
    S = load(f"sigma_{tag}_sw")
    map_plot(S, SIGMA_BINS, SIGMA_LABEL, f"gam_{tag}_sigma_main.pdf",
             tick_labels=SIGMA_LABELS, mark=ar1_mark(np.nan_to_num(S, nan=-1) >= SKALA_SIG),
             title=f"{name}: GAM, main model (seasonality corrected), nonlinear part")

# ---- 2. Agreement with Skala ----------------------------------------------------------
AGREE_BINS   = [-0.5, 0.5, 1.5, 2.5, 3.5]
AGREE_LABELS = ["neither", "GAM only", "Skala only", "both"]
AGREE_COLORS = ["#f4f4f2", "#2a78d6", "#eda100", "#d55181"]   # checked: distinguishable with colour blindness

for tag, name in NAMES.items():
    for suf, model in (("", "pooled"), ("_sw", "main")):
        A = load(f"agree_{tag}{suf}")
        map_plot(A, AGREE_BINS, "Nonlinearity $\\geq$ 2 $\\sigma$",
                 f"gam_{tag}_agree_{model}.pdf", cmap=AGREE_COLORS, tick_labels=AGREE_LABELS,
                 mark=ar1_mark(np.isin(A, [1, 3])) if suf == "_sw" else None,   # GAM only / both
                 title=f"{name}: GAM ({'all months pooled' if suf == '' else 'main model'}) vs Skala (2025)")

# ---- 3. Nonlinear edf ------------------------------------------------------------------
EDF_BINS = [0, 0.5, 1, 1.5, 2, 3, 4, 5, 8]

for tag, name in NAMES.items():
    for suf, model in (("", "pooled"), ("_sw", "main")):
        E = load(f"edf_{tag}{suf}")
        ns = not_sig(E, load(f"q_{tag}{suf}"))
        map_plot(E, EDF_BINS, "Nonlinear edf (0 = straight line)", f"gam_{tag}_edf_{model}.pdf",
                 cmap="YlOrRd", grey_mask=ns, grey_label=f"Straight line (BH q $\\geq$ {ALPHA})",
                 mark=ar1_mark(~np.isnan(E) & ~ns) if suf == "_sw" else None,
                 title=f"{name}: GAM, {'all months pooled' if suf == '' else 'main model'}, nonlinear edf")

# ---- 4. Bend of the main model's curve ---------------------------------------------------
BEND_BINS = [-1, -0.5, -0.25, -0.1, 0.1, 0.25, 0.5, 1]

for tag, name in NAMES.items():
    B = load(f"bend_{tag}")
    ns = not_sig(B, load(f"q_{tag}_sw"))
    map_plot(B, BEND_BINS, "Bend f(+2SD) - 2 f(mean) + f(-2SD) (K)", f"gam_{tag}_bend.pdf",
             cmap="RdBu_r", extend_low=True, grey_mask=ns,
             grey_label=f"Straight line (BH q $\\geq$ {ALPHA})", mark=ar1_mark(~np.isnan(B) & ~ns),
             title=f"{name}: GAM, main model, size of the nonlinearity")

# ---- 5. Strength, K per 4 SD ---------------------------------------------------------------
LIN4_BINS = [-6, -3, -2, -1, -0.5, -0.25, 0.25, 0.5, 1, 2, 3, 6]
LIN4_LABEL = "Response to +4 SD of the index (K)"

for tag, name in NAMES.items():
    L = load(f"lin4_{tag}")
    map_plot(L, LIN4_BINS, LIN4_LABEL, f"gam_{tag}_strength_all_months.pdf", cmap="RdBu_r",
             extend_low=True, grey_mask=not_sig(L, load(f"q_lin_{tag}")),
             grey_label=f"BH q $\\geq$ {ALPHA}",
             title=f"{name}: straight line through all months"
                   + ("\ntheir ENSO index, the SOI, has the opposite sign" if tag == "nino" else ""))
    Qd = load(f"q_dep_{tag}")
    Ls = [load(f"lin4_{sn}_{tag}") for sn in SEASONS]
    season_panels(Ls, LIN4_BINS, LIN4_LABEL, f"gam_{tag}_strength_seasons.pdf", "RdBu_r",
                  [not_sig(F, Qd) for F in Ls], f"No dependence in any season (BH q $\\geq$ {ALPHA})",
                  f"{name}: one straight line per season, main model", extend_low=True)

# ---- 6. AR(1) check of the main model --------------------------------------------------------
RHO_BINS   = [-0.2, -0.05, 0.05, 0.1, 0.2, 0.3, 0.6]
RHO_COLORS = ["#4575b4", "#f7f7f7", "#fee090", "#fdae61", "#f46d43", "#d73027"]

map_plot(load("rho_left_sw"), RHO_BINS, "Lag-1 autocorrelation left after whitening",
         "gam_ar1_check.pdf", cmap=RHO_COLORS, extend_low=True,
         title="AR(1) check of the main model (adequate within $\\pm$0.05)")
