#%%

from importlib import reload
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import parse
import seaborn as sns
from tqdm import tqdm

import utils

#%%

plt.style.use("plotstyle.mplstyle")


#%%


parser_template = "{sample}__{simulation_method:SimName}__{simulated_N_reads:Int}"
parser = parse.compile(
    parser_template,
    dict(
        Int=int,
        SimName=lambda s: utils.fix_sim_name(s, utils.D_SIM_NAME_TRANSLATE),
    ),
)


sim_columns = ["sample", "simulation_method", "simulated_N_reads"]


def split_simulation_name(simulation_name):
    result = parser.parse(simulation_name).named
    return result["sample"], result["simulation_method"], result["simulated_N_reads"]


def split_name_pd(name):
    return pd.Series(split_simulation_name(name), index=sim_columns)


def load_results():

    df = (
        pd.read_parquet("data/results/")
        .reset_index(drop=True)
        .rename(columns={"sample": "simulation_name"})
        .astype({"tax_id": int})
    ).copy()

    df["Bayesian_significance"] = df["Bayesian_D_max"] / df["Bayesian_D_max_std"]
    df["Bayesian_prob_not_zero_damage"] = 1 - df["Bayesian_prob_zero_damage"]
    df["Bayesian_prob_gt_1p_damage"] = 1 - df["Bayesian_prob_lt_1p_damage"]
    df[sim_columns] = df["simulation_name"].apply(split_name_pd)

    columns = [
        "simulation_name",
        "sample",
        "simulation_method",
        "simulated_N_reads",
        "tax_id",
        "tax_name",
        "tax_rank",
        "N_reads",
        "N_alignments",
        "D_max",
        "D_max_std",
        "Bayesian_D_max",
        "Bayesian_D_max_std",
        "significance",
        "Bayesian_significance",
        "Bayesian_prob_not_zero_damage",
        "Bayesian_prob_gt_1p_damage",
        "mean_L",
        "mean_GC",
        "q",
        "A",
        "c",
        "phi",
        "rho_Ac",
        "valid",
        "asymmetry",
        "std_L",
        "std_GC",
        "lambda_LR",
        "q_std",
        "phi_std",
        "A_std",
        "c_std",
        "N_x=1_forward",
        # "N_x=1_reverse",
        "N_sum_total",
        "N_sum_forward",
        # "N_sum_reverse",
        "N_min",
        "k_sum_total",
        "k_sum_forward",
        # "k_sum_reverse",
        "Bayesian_D_max_median",
        "Bayesian_D_max_confidence_interval_1_sigma_low",
        "Bayesian_D_max_confidence_interval_1_sigma_high",
        "Bayesian_D_max_confidence_interval_2_sigma_low",
        "Bayesian_D_max_confidence_interval_2_sigma_high",
        "Bayesian_D_max_confidence_interval_3_sigma_low",
        "Bayesian_D_max_confidence_interval_3_sigma_high",
        "Bayesian_D_max_confidence_interval_95_low",
        "Bayesian_D_max_confidence_interval_95_high",
        "Bayesian_A",
        "Bayesian_A_std",
        "Bayesian_A_median",
        "Bayesian_A_confidence_interval_1_sigma_low",
        "Bayesian_A_confidence_interval_1_sigma_high",
        "Bayesian_q",
        "Bayesian_q_std",
        "Bayesian_q_median",
        "Bayesian_q_confidence_interval_1_sigma_low",
        "Bayesian_q_confidence_interval_1_sigma_high",
        "Bayesian_c",
        "Bayesian_c_std",
        "Bayesian_c_median",
        "Bayesian_c_confidence_interval_1_sigma_low",
        "Bayesian_c_confidence_interval_1_sigma_high",
        "Bayesian_phi",
        "Bayesian_phi_std",
        "Bayesian_phi_median",
        "Bayesian_phi_confidence_interval_1_sigma_low",
        "Bayesian_phi_confidence_interval_1_sigma_high",
        "Bayesian_rho_Ac",
        # "Bayesian_significance",
        "var_L",
        "var_GC",
        # "f+1",
        # "f+15",
        # "f-1",
        # "f-15",
    ]

    df = df.loc[:, columns]

    return df


#%%

df_results = load_results()
df_results_species = df_results.query("tax_rank == 'species'")


#%%

reload(utils)

directory = Path("input-data") / "data-pre-mapping"
df_simulation = utils.get_simulation_details(directory)


#%%

# x = x


#%%

good_samples = [
    "Lake-9",
    # "Lake-7",
    "Lake-7-forward",
    "Cave-22",
    # "Cave-100",
    "Cave-100-forward",
    "Cave-102",
    "Pitch-6",
]


df_good = df_results_species.query("sample in @good_samples")


dfs_ancient = []
dfs_modern = []
dfs_non_simulation = []

for (sample, simulated_N_reads), group in df_good.groupby(
    ["sample", "simulated_N_reads"]
):
    # break

    if "forward" in sample:
        sample = sample.replace("-forward", "")

    query = f"sample == '{sample}' and simulated_N_reads == {simulated_N_reads}"
    df_simulation_group = df_simulation.query(query)
    drop_cols = ["sample", "simulated_N_reads"]

    # ancient

    query = "simulated_only_ancient == 'True'"
    df_simulation_group_ancient = df_simulation_group.query(query).copy()
    df_simulation_group_ancient["type"] = "Ancient"

    dfs_ancient.append(
        group.merge(
            df_simulation_group_ancient.drop(columns=drop_cols),
            on="tax_id",
        )
    )

    # modern

    query = "simulated_only_ancient != 'True'"
    df_simulation_group_modern = df_simulation_group.query(query).copy()
    df_simulation_group_modern["type"] = "Non-ancient"

    dfs_modern.append(
        group.merge(
            df_simulation_group_modern.drop(columns=drop_cols),
            on="tax_id",
        ),
    )

    # non-simulation

    tax_ids_non_simulation = set(group.tax_id) - set(df_simulation_group.tax_id)
    df_non_simulation = group.query(f"tax_id in @tax_ids_non_simulation").copy()
    df_non_simulation["type"] = "Non-simulated"
    dfs_non_simulation.append(df_non_simulation)


df_ancient = pd.concat(dfs_ancient).reset_index(drop=True)
df_modern = pd.concat(dfs_modern).reset_index(drop=True)
df_non_simulation = pd.concat(dfs_non_simulation).reset_index(drop=True)

#%%

df_all = pd.concat([df_ancient, df_modern, df_non_simulation], axis=0)

df_all.to_parquet("df_all.parquet")


#%%

import matplotlib.ticker as mtick
from matplotlib.backends.backend_pdf import PdfPages


def plot_overview(
    df_in,
    title="",
    xlims=None,
    ylims=None,
    figsize=(10, 4),
    types=None,
    loc="upper left",
    legend_i=1,
):

    if types is None:
        types = df_all["type"].unique()

    samples = [
        "Cave-100-forward",
        "Cave-102",
        "Pitch-6",
        "Cave-22",
        "Lake-9",
        "Lake-7-forward",
    ]
    d_colors = {sample: f"C{i}" for i, sample in enumerate(df_all["sample"].unique())}

    symbols = [
        "o",
        "s",
        "D",
        "v",
        "^",
        ">",
        "<",
        "*",
        "x",
    ]

    fig, axes = plt.subplots(
        figsize=figsize,
        ncols=len(types),
        sharey=True,
    )

    for i_type, (type_, ax) in enumerate(zip(types, axes)):
        # break

        df_damage_type = df_in.query("type == @type_")

        for i_sample, sample in enumerate(samples):

            group = df_damage_type.query("sample == @sample")

            ax.scatter(
                group["Bayesian_significance"],
                group["Bayesian_D_max"],
                s=2 + np.sqrt(group["N_reads"]) / 5,
                alpha=0.8,
                color=d_colors[sample],
                label=sample,
                marker=symbols[i_sample],
            )

        ax.set(
            title=type_,
            xlabel="Significance",
            ylabel="Damage",
            xlim=xlims[i_type] if xlims is not None else None,
            ylim=ylims[i_type] if ylims is not None else None,
        )

        ax.yaxis.set_major_formatter(mtick.PercentFormatter(1.0))

        if i_type == legend_i:
            ax.yaxis.set_tick_params(labelbottom=True)
            leg = ax.legend(markerscale=5, loc=loc)
            for handle in leg.legendHandles:
                handle.set_sizes([30.0])
                handle.set_alpha(1)

    if title != "":
        fig.suptitle(title, fontsize=16)
        fig.subplots_adjust(top=0.85)

    return fig


#%%


x = x

#%%


simulation_methods = ["art", "deam", "frag"]

filename = Path("figures/overview_bayesian_all.pdf")
filename.parent.mkdir(exist_ok=True, parents=True)
with PdfPages(filename) as pdf:

    for simulation_method in simulation_methods:
        # break

        fig = plot_overview(
            df_all.query(f"simulation_method == '{simulation_method}'"),
            title="Simulation method: " + simulation_method,
            figsize=(10, 4),
            xlims=[(-0.1, 26), (0, 2), (0, 4.1)],
            ylims=[(-0.01, 0.7), (-0.01, 0.7), (-0.01, 0.7)],
        )

        fig.tight_layout()
        fig

        pdf.savefig(fig)
        plt.close()

#%%

fig = plot_overview(
    df_all.query(f"simulation_method == 'art'"),
    figsize=(7, 3),
    xlims=[(0, 26), (0, 2)],
    ylims=[(0, 0.7), (0, 0.7)],
    types=["Ancient", "Non-ancient"],
)

fig.tight_layout()
fig.savefig("figures/overview_bayesian_art.pdf")


#%%


fig = plot_overview(
    df_all.query(f"simulation_method == 'frag'"),
    figsize=(7, 3),
    xlims=[(0, 2), (0, 2)],
    ylims=[(0.0, 0.15), (0.0, 0.15)],
    types=["Ancient", "Non-ancient"],
    # legend_i=0,
    loc="upper right",
)

fig.tight_layout()
fig.savefig("figures/overview_bayesian_frag.pdf")


#%%

# df_modern.sort_values("Bayesian_significance", ascending=False).head(10)
# df_non_simulation.sort_values("Bayesian_significance", ascending=False).head(10)
# df_results.query("tax_id == 134927")

# %%


# %%

from adjustText import adjust_text
from engineering_notation import EngNumber


def map_simulated_N_reads_to_x_axis(xs, delta=0):
    d_translate = {
        1000000: 1,
        5000000: 2,
        10000000: 3,
        50000000: 4,
        100000000: 5,
    }
    return [d_translate[x] + delta for x in xs]


def plot_damage_for_df_tax_id(
    df_tax_id,
    tax_id,
    use_bayesian=True,
    title=None,
):

    if use_bayesian:
        D_max_col = "Bayesian_D_max"
        D_max_std_col = "Bayesian_D_max_std"
        # D_max_std_col = [
        #     "Bayesian_D_max_confidence_interval_1_sigma_low",
        #     "Bayesian_D_max_confidence_interval_1_sigma_high",
        # ]
    else:
        D_max_col = "D_max"
        D_max_std_col = "D_max_std"

    sample = df_tax_id["sample"].iloc[0]
    tax_name = df_tax_id["tax_name"].iloc[0]
    type_ = df_tax_id["type"].iloc[0]

    d_colors = {"frag": "C2", "deam": "C0", "art": "C1"}
    d_y_diffs = {"frag": 0.3 / 100, "deam": 0.3 / 100, "art": 0.5 / 100}

    DELTA = 0.20

    d_method_names = {
        "frag": "Fragmentation",
        "deam": "Deamination",
        "art": "Sequencing Noise",
    }

    methods = ["frag", "deam", "art"]
    deltas = [-DELTA, 0, DELTA]

    fig, ax = plt.subplots()

    for method, delta in zip(methods, deltas):
        # break

        df_method = df_tax_id.query(f"simulation_method == '{method}'").set_index(
            "simulated_N_reads",
            drop=False,
        )

        if len(df_method) == 0:
            continue

        xx = map_simulated_N_reads_to_x_axis(df_method.simulated_N_reads, delta)
        yy = df_method[D_max_col]

        sy = df_method[D_max_std_col]

        ax.errorbar(
            xx,
            yy,
            sy,
            fmt=".",
            label=d_method_names[method],
            color=d_colors[method],
        )

    if type_ == "Ancient":

        simulated_D = df_tax_id["simulated_D_max"].iloc[0]
        ax.axhline(
            simulated_D,
            label="Simulated Damage",
            color="k",
            alpha=0.5,
            ls="--",
            zorder=0,
        )
        ax.text(
            5.55,
            simulated_D,
            f"{simulated_D*100:.1f}" + r"\%",
            horizontalalignment="left",
            verticalalignment="center",
            fontsize=10,
        )

    ax.set_xticks(
        np.arange(5) + 1,
        labels=["1M", "5M", "10M", "50M", "100M"],
    )

    ax.set(
        xlim=(0.5, 5.5),
        ylim=(0, ax.get_ylim()[1] * 1.25),
        ylabel="Damage" if use_bayesian else "Damage (MAP)",
        xlabel="Simulation size",
    )

    for method, delta in zip(methods, deltas):

        df_method = df_tax_id.query(f"simulation_method == '{method}'").set_index(
            "simulated_N_reads",
            drop=False,
        )

        if len(df_method) == 0:
            continue

        xx = map_simulated_N_reads_to_x_axis(df_method.simulated_N_reads, delta)
        yy = df_method[D_max_col]
        sy = df_method[D_max_std_col]

        for xi, yi, si, N_i in zip(xx, yy, sy, df_method["N_reads"]):

            # ypos = yi + si + 0.55 / 100
            ypos = yi + si + 3 * ax.get_ylim()[1] / 100
            if method == "art" and yi - si - 10 * ax.get_ylim()[1] / 100 > 0:
                ypos = yi - si - 5 * ax.get_ylim()[1] / 100

            ax.text(
                xi,
                ypos,
                f"{EngNumber(N_i, precision=0)}",
                horizontalalignment="center",
                verticalalignment="center",
                fontsize=10,
                color=d_colors[method],
            )

    ax.yaxis.set_major_formatter(mtick.PercentFormatter(1.0))

    if title == "":
        pass
    else:
        if title is None:
            title = f"{sample}\nTax ID: {tax_id}, {tax_name}\n{type_}\n"
        fig.suptitle(
            title,
            fontsize=12,
            ha="left",
            x=0.15,
            y=0.975,
        )
        fig.subplots_adjust(top=0.80)

    ax.legend(
        # title="Legend",
        ncol=1,
        loc="upper right",
        markerscale=0.7,
        bbox_to_anchor=(1.07, 1.15),
        # fontsize=10,
    )

    return fig


# %%


def plot_damages(df_all, use_bayesian=True):

    for sample, df_sample in df_all.groupby("sample", sort=False):

        tax_ids_in_all = (
            df_sample.groupby("tax_id")
            .sum(numeric_only=True)["N_reads"]
            .sort_values(ascending=False)
            .index
        )

        if len(tax_ids_in_all) == 0:
            continue

        suffix = "bayesian" if use_bayesian else "MAP"

        fig_name = (
            Path("figures") / suffix / "damages" / f"{suffix}.{sample}.damage.pdf"
        )
        fig_name.parent.mkdir(exist_ok=True, parents=True)

        with PdfPages(fig_name) as pdf:

            for tax_id in tqdm(tax_ids_in_all):
                # break

                df_tax_id = df_sample.query(f"tax_id == {tax_id}")

                if df_tax_id.N_reads.max() < 10:
                    continue

                fig = plot_damage_for_df_tax_id(
                    df_tax_id,
                    tax_id,
                    use_bayesian=use_bayesian,
                )
                pdf.savefig(fig)
                plt.close()


plot_damages(df_all, use_bayesian=True)
plot_damages(df_all, use_bayesian=False)
# %%


df_ancient_art = df_ancient.query("simulation_method == 'art'")

print("\n\nLoose cut:")
for sample, df_ancient_sample in df_ancient_art.groupby("sample", sort=False):
    # break

    if sample == "Lake-7-forward":
        break

    mask = (df_ancient_sample.Bayesian_D_max > 0.01) & (
        df_ancient_sample.Bayesian_significance > 2
    )

    df_ancient_sample[mask]

    mask_N_reads_100 = df_ancient_sample.N_reads > 100

    mask_conditional = (df_ancient_sample[mask_N_reads_100].Bayesian_D_max > 0.01) & (
        df_ancient_sample[mask_N_reads_100].Bayesian_significance > 2
    )

    s = (
        f"{sample+':':20s}"
        f" Total: {len(df_ancient_sample):4d}."
        f" Total pass: {np.sum(mask):3d}, {np.mean(mask):5.1%}"
        f" >100 reads: {np.sum(mask_N_reads_100):3d}, {np.mean(mask_N_reads_100):5.1%}"
        f" Conditional: {np.sum(mask_conditional):3d}, {np.mean(mask_conditional):5.1%}"
    )
    print(s)

# %%


tax_id = 129488
df_tax_id = df_all.query(f"tax_id == {tax_id} and sample == 'Pitch-6'")

fig = plot_damage_for_df_tax_id(
    df_tax_id,
    tax_id,
    use_bayesian=True,
    title="",
)

fig.savefig("figures/damage.Pitch-6.129488.pdf", bbox_inches="tight")

# %%
