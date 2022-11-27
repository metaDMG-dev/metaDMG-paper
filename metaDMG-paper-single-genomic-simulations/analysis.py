#%%

from importlib import reload
from multiprocessing.spawn import get_executable
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages
from tqdm import tqdm

import utils

#%%

plt.style.use("plotstyle.mplstyle")


#%%

# save_plots = False
# save_plots = True


make_plots = False
# make_plots = True


#%%

# D_DAMAGE = utils.D_DAMAGE


#%%

all_species = [
    "homo",
    "betula",
    "GC-low",
    "GC-mid",
    "GC-high",
    "contig1k",
    "contig10k",
    "contig100k",
]

#%%

reload(utils)
df_all = utils.load_multiple_species(all_species)

#%%

df = df_all.query(
    f"sim_length == 60 and sim_seed < 100 and sim_species in {all_species[:-3]}"
)

#%%


# reload(utils)
# df_damaged_reads = utils.load_multiple_damaged_reads(all_species)
df_damaged_reads = None

#%%

# reload(utils)
df_known_damage = utils.get_df_known_damage()


#%%

# x = x

#%%

reload(utils)
if make_plots:
    print("Plotting individual damage")
    utils.plot_individual_damage_results(
        df=df,
        df_known_damage=df_known_damage,
        df_damaged_reads=df_damaged_reads,
    )


#%%

reload(utils)
df_aggregated = utils.get_df_aggregated(
    df_in=df,
    df_known_damage=df_known_damage,
    df_damaged_reads=df_damaged_reads,
)


# %%

# reload(utils)
if make_plots:
    print("Plotting combined damage")
    utils.plot_combined_damage_results(
        df=df,
        df_aggregated=df_aggregated,
        df_known_damage=df_known_damage,
        df_damaged_reads=df_damaged_reads,
    )


#%%


df_aggregated_homo = df_aggregated.query("sim_species == 'homo'")

reload(utils)
if make_plots:
    fig = utils.plot_combined_MAEs(
        df_aggregated_homo=df_aggregated_homo,
        df_known_damage=df_known_damage,
        method="Bayesian",
        # ylim=(0, 2),
    )


# %%

reload(utils)
if make_plots:
    print("Plotting contour lines")
    utils.plot_contour_lines(df)


#%%


# reload(utils)
if make_plots:

    reload(utils)
    fig, ax = plt.subplots(figsize=(3.5, 3))

    utils.plot_contour_lines_on_ax(
        df,
        ax=ax,
        sim_species="homo",
        method="Bayesian",
        title="",
        frac_legend_pos=(1, 0.60),
        fill_contour=True,
        fill_cut_value=4,
        # fill_cut_value=2,
        fill_level_value=0.95,
        # fill_level_value=0.5,
    )

    ax.set(xlim=(20, 1.2 * 10**5))

    # fig.tight_layout()

    fig.savefig("figures/contour_homo.pdf", bbox_inches="tight")


#%%

# reload(utils)

df_zero_damage = df_all.query(
    "sim_species == 'homo' and sim_damage == 0 and sim_length == 60"
)

if make_plots:
    utils.plot_zero_damage_groups(df_zero_damage)


#%%

df_homo_99 = df_all.query("sim_species == 'homo' and sim_seed < 100")

reload(utils)
if make_plots:

    print("Plotting individual damage lengths")
    utils.plot_individual_damage_results_lengths(
        df_homo_99=df_homo_99,
        df_known_damage=df_known_damage,
        df_damaged_reads=df_damaged_reads,
    )


#%%

reload(utils)

df_aggregated_lengths = utils.get_df_aggregated(
    df_in=df_homo_99,
    df_known_damage=df_known_damage,
    df_damaged_reads=df_damaged_reads,
)

if make_plots:

    print("Plotting combined damage lengths")
    utils.plot_combined_damage_results_lengths(
        df_homo_99=df_homo_99,
        df_aggregated_lengths=df_aggregated_lengths,
        df_known_damage=df_known_damage,
        df_damaged_reads=df_damaged_reads,
    )


#%%


df_contigs = df_all.query(f"sim_species in {all_species[-3:]}")

reload(utils)
if make_plots:

    print("Plotting individual damage: contigs")
    utils.plot_individual_damage_results(
        df=df_contigs,
        df_known_damage=df_known_damage,
        df_damaged_reads=df_damaged_reads,
        suffix="_contigs",
    )

# %%

reload(utils)

df_aggregated_contigs = utils.get_df_aggregated(
    df_in=df_contigs,
    df_known_damage=df_known_damage,
    df_damaged_reads=df_damaged_reads,
)

# reload(utils)
if make_plots:
    print("Plotting combined damage: contigs")
    utils.plot_combined_damage_results(
        df=df_contigs,
        df_aggregated=df_aggregated_contigs,
        df_known_damage=df_known_damage,
        df_damaged_reads=df_damaged_reads,
        suffix="_contigs",
    )


#%%

sim_species = "homo"
sim_damage = 0.31
sim_N_reads = 100
sim_length = 60
min_seed = 60
max_seed = 80

group = df.query(
    f"sim_species == '{sim_species}'"
    f" and sim_damage == {sim_damage}"
    f" and sim_length == {sim_length}"
    f" and sim_N_reads == {sim_N_reads}"
    f" and {min_seed} <= sim_seed < {max_seed}"
)


group_agg = df_aggregated.query(
    f"sim_species == '{sim_species}'"
    f" and sim_damage == {sim_damage}"
    f" and sim_length == {sim_length}"
    # f" and sim_N_reads == {sim_N_reads}"
    # f" and sim_seed < {max_seed}"
)


#%%

import arviz as az
from scipy.stats import beta as sp_beta
from scipy.stats import betabinom as sp_betabinom
from scipy.stats import norm as sp_norm

posterior = az.from_netcdf("sim-homo-0.31-100-60-69.nc").posterior


A = posterior["A"].values[0]
phi = posterior["phi"].values[0]
N = 16
mu = np.mean(A)
stds = np.sqrt(A * (1 - A) * (phi + N) / ((phi + 1) * N))
std = np.mean(stds)


#%%

fig, ax = plt.subplots()
xmin, xmax = 0, 0.18
Nbins = 50
ax.hist(A, Nbins, range=(xmin, xmax), density=True, histtype="step", label=r"$\bar{D}$")
ax.hist(
    stds, Nbins, range=(xmin, xmax), density=True, histtype="step", label=r"$\sigma_D$"
)
ax.set(xlim=(xmin, xmax))
ax.legend()

#%%
if make_plots:

    fig, (ax1, ax2) = plt.subplots(figsize=(10, 3.5), ncols=2)

    ymin, ymax = -0.0001, 0.25

    reload(utils)
    utils.plot_individual_damage_result(
        df_in=df,
        group_all_keys=group,
        df_known_damage=df_known_damage,
        df_damaged_reads=df_damaged_reads,
        method="Bayesian",
        splitby="species",
        keys=["homo"],
        xlim=(min_seed - 0.5, max_seed - 0.1),
        ylim=(ymin, ymax),
        # fig_title=f"Simulation, {sim_N_reads} reads",
        fig_title=f"",
        ax_titles=False,
        ax_in=ax1,
        # loc="upper left",
        bbox_to_anchor=None,
        ncols=1,
        markerscale=0.7,
    )

    reload(utils)
    utils.plot_combined_damage_result(
        df_in=df,
        group_agg_all_keys=group_agg,
        df_known_damage=df_known_damage,
        df_damaged_reads=None,
        method="Bayesian",
        splitby="species",
        keys=["homo"],
        fig_title=f"Simulation",
        # xlim=(0.7 * 10**2, 1.3 * 10**5),
        xlim=(0.7 * 10, 1.3 * 10**5),
        ylim=(ymin, ymax),
        ax_titles=False,
        delta=0.1,
        ax_in=ax2,
        loc="upper right",
        markerscale=0.7,
    )

    # ax1.annotate("A)", (0.02, 0.9), xycoords="axes fraction", fontsize=14)
    ax1.set_title(
        r"\textbf{A}) Homo, $\delta_\mathrm{ss}$ = "
        f"{sim_damage:.2f}, L = {sim_length}, {sim_N_reads} reads",
        pad=15,
    )
    ax2.set_title(
        r"\textbf{B}) Homo, $\delta_\mathrm{ss}$ = "
        f"{sim_damage:.2f}, L = {sim_length}",
        pad=15,
    )

    fig.tight_layout()

    fig.savefig("figures/ngsngs_overview.pdf", bbox_inches="tight")


#%%


group_zero_damage = df_zero_damage.query("sim_N_reads == 1000")

if make_plots:

    reload(utils)

    g = utils.plot_zero_damage_group(
        group_zero_damage,
        method="Bayesian",
        title="",
        xlim=(0.1, 0.9),
        ylim=(0.0, 0.008),
    )

    g.savefig("figures/zero_damage_1000_reads.pdf")
# %%


#%%


def parse_pydamage_name(name):
    _, species, damage, N_reads, L, seed = name.split("-")
    return species, float(damage), int(N_reads), int(L), int(seed.split(".")[0])


def load_pydamage_results():

    if Path("pydamage/pydamage.parquet").exists():
        df = pd.read_parquet("pydamage/pydamage.parquet")

    else:

        paths = Path("pydamage") / "homo"

        dfs = []
        for path in tqdm(list(paths.glob("*.csv"))):
            df = pd.read_csv(path)

            sim_data = parse_pydamage_name(path.stem)
            for i, col in enumerate(utils.simulation_columns):
                df[col] = sim_data[i]

            dfs.append(df)

        df = (
            pd.concat(dfs)
            .drop(columns=["reference"])
            .sort_values(["sim_damage", "sim_N_reads", "sim_seed"])
            .reset_index(drop=True)
        )
        df["sim_damage_percent_approx"] = df["sim_damage"].map(utils.D_DAMAGE_APPROX)
        df["D_max"] = df["damage_model_pmax"]
        df["D_max_std"] = df["damage_model_pmax_stdev"]
        df["significance"] = df["D_max"] / df["D_max_std"]

        df.to_parquet("pydamage/pydamage.parquet")

    return df


df_pydamage = load_pydamage_results().query("sim_length == 60")


#%%


# %%


sim_species = "homo"
sim_length = 60

df_pydamage_100 = df_pydamage.query("sim_seed < 100").copy()
df_metaDMG_100 = df_all.query(
    f"sim_species == '{sim_species}' & sim_length == {sim_length} & sim_seed < 100"
).copy()


df_pydamage_100["method"] = "pydamage"
df_metaDMG_100["method"] = "metaDMG"

cols = utils.simulation_columns + [
    "sim_damage_percent_approx",
    "D_max",
    "D_max_std",
    "significance",
    "method",
]


df_combined = pd.concat(
    [
        df_pydamage_100[cols + ["predicted_accuracy", "qvalue"]],
        df_metaDMG_100[cols],
    ]
)


df_combined_wide = pd.concat(
    [
        df_metaDMG_100.reset_index(drop=True)[cols[:-1]],
        df_pydamage_100.reset_index(drop=True)[
            ["predicted_accuracy", "damage_model_pmax", "qvalue"]
        ],
    ],
    axis=1,
)

#%%


x = x

#%%

filename = Path(f"figures/pydamage_comparison_special_axis.pdf")
filename.parent.mkdir(parents=True, exist_ok=True)

with PdfPages(filename) as pdf:

    it = tqdm(df_combined.groupby(["sim_damage", "sim_N_reads"]))

    for (sim_damage, sim_N_reads), group in it:
        # break

        if sim_N_reads == 100 and sim_damage == 0.472:
            break

        reload(utils)
        fig = utils.plot_pydamage_comparison(
            group=group,
            df_known_damage=df_known_damage,
            sim_damage=sim_damage,
            sim_species=sim_species,
            sim_N_reads=sim_N_reads,
            sim_length=sim_length,
            use_special_axis=True,
        )

        pdf.savefig(fig, bbox_inches="tight")
        plt.close()


#%%

df_pydamage_zero_damage = df_pydamage.query(
    "sim_species == 'homo' and sim_damage == 0 and sim_length == 60"
)

df_metaDMG_zero_damage = df_all.query(
    "sim_species == 'homo' and sim_damage == 0 and sim_length == 60"
)


df_pydamage_zero_damage.sort_values("significance")
df_pydamage_zero_damage.sort_values("pvalue")


df_pydamage_zero_damage_cut = df_pydamage_zero_damage.query(
    "pvalue < 0.05 and predicted_accuracy > 0.67"
)

df_pydamage_zero_damage_cut = df_pydamage_zero_damage.query(
    "pvalue < 0.05 and predicted_accuracy > 0.50 and D_max > 0.01"
    # "pvalue < 0.05 and predicted_accuracy > 0.50"
)


df_pydamage_zero_damage_cut100 = df_pydamage_zero_damage.query("sim_N_reads == 100")
df_pydamage_zero_damage_cut500 = df_pydamage_zero_damage.query("sim_N_reads == 500")
df_pydamage_zero_damage_cut500.sort_values("pvalue")


# sns.scatterplot(data=df_pydamage_zero_damage_cut, x="predicted_accuracy", y="damage_model_pmax")


#%%


# %%

filename = Path(f"figures/pydamage_comparison_zero_damage.pdf")
filename.parent.mkdir(parents=True, exist_ok=True)

with PdfPages(filename) as pdf:

    it = tqdm(df_combined.query("sim_damage == 0").groupby("sim_N_reads"))

    for sim_N_reads, group_zero_damage in it:
        fig = utils.plot_pydamage_comparison_zero_damage(sim_N_reads, group_zero_damage)

        pdf.savefig(fig, bbox_inches="tight")
        plt.close()


#%%

x = x

reload(utils)


filename = Path(f"figures/pydamage_comparison2.pdf")
filename.parent.mkdir(parents=True, exist_ok=True)

with PdfPages(filename) as pdf:

    it = df_combined_wide.groupby(["sim_damage", "sim_N_reads"])

    for (sim_damage, sim_N_reads), group in tqdm(it):

        group = utils.add_colors(group)

        fig = utils.make_parallel_plots(
            group,
            sim_damage,
            sim_N_reads,
            smooth_lines=True,
        )

        pdf.savefig(fig, bbox_inches="tight", transparent=True)
        plt.close()


filename = Path(f"figures/pydamage_comparison2-without-10-reads.pdf")
filename.parent.mkdir(parents=True, exist_ok=True)

with PdfPages(filename) as pdf:

    it = df_combined_wide.query("sim_N_reads > 11").groupby(
        ["sim_damage", "sim_N_reads"]
    )

    for (sim_damage, sim_N_reads), group in tqdm(it):

        group = utils.add_colors(group)

        fig = utils.make_parallel_plots(
            group,
            sim_damage,
            sim_N_reads,
            smooth_lines=True,
        )

        pdf.savefig(fig, bbox_inches="tight", transparent=True)
        plt.close()

#%%

reload(utils)
utils.plot_all_aggregate_groups(
    df_aggregated,
    df_aggregated_lengths,
    df_aggregated_contigs,
    df_known_damage=df_known_damage,
)

# %%
