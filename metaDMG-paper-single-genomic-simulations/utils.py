#%%
from pathlib import Path

import matplotlib as mpl
import matplotlib.patches as patches
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import LogNorm
from matplotlib.path import Path as mpl_Path
from scipy.ndimage import gaussian_filter
from scipy.stats import norm as sp_norm
from tqdm import tqdm

#%%


#%%


# D_DAMAGE = {
#     0.0: 0.0,  # np.mean([0.00676 - 0.00526, 0.00413 - 0.00137]),
#     0.014: 0.01,  # np.mean([0.01127 - 0.00526, 0.00841 - 0.00137]),
#     0.047: 0.02,  # np.mean([0.02163 - 0.00526, 0.01881 - 0.00137]),
#     0.138: 0.05,  # np.mean([0.05111 - 0.00524, 0.04824 - 0.00137]),
#     0.303: 0.10,  # np.mean([0.10149 - 0.00523, 0.09855 - 0.00137]),
#     0.466: 0.15,  # np.mean([0.15183 - 0.00518, 0.14900 - 0.00137]),
#     0.626: 0.20,
#     0.96: 0.30,  # np.mean([0.30046 - 0.00518, 0.29910 - 0.00141]),
# }


D_DAMAGE_APPROX = {
    0.0: 0.0,
    0.035: 0.01,
    0.065: 0.02,
    0.162: 0.05,
    0.31: 0.10,
    0.472: 0.15,
    0.633: 0.20,
    0.96: 0.30,
}


#%%


def sdom(x):
    return x.std() / np.sqrt(len(x))


#%%


simulation_columns = [
    "sim_species",
    "sim_damage",
    "sim_N_reads",
    "sim_length",
    "sim_seed",
]


def split_name(name):
    splitted = name.split("-")
    damage, N_reads, length, seed = splitted[-4:]
    specie = "-".join(splitted[1:-4])
    return specie, float(damage), int(N_reads), int(length), int(seed)


def split_name_pd(name):
    return pd.Series(split_name(name), index=simulation_columns)


#%%


def get_df_known_damage_single_path(path):

    specie, damage, N_reads, length, seed = split_name(path.stem)

    d = {
        "sim_species": specie,
        "sim_damage": damage,
        "sim_N_reads": N_reads,
        "sim_length": length,
        "sim_seed": seed,
    }

    df_forward_reverse = pd.read_csv(path, sep="\t")
    df_forward = df_forward_reverse.iloc[:50].astype(float)
    df_reverse = df_forward_reverse.iloc[51:].astype(float).reset_index(drop=True)

    d["C>T(1)"] = df_forward["C>T"][1 - 1]
    d["C>T(15)"] = df_forward["C>T"][15 - 1]
    d["C>T(50)"] = df_forward["C>T"][50 - 1]
    d["C>T_diff"] = d["C>T(1)"] - d["C>T(15)"]

    d["G>A(1)"] = df_reverse["G>A"][1 - 1]
    d["G>A(15)"] = df_reverse["G>A"][15 - 1]
    d["G>A(50)"] = df_reverse["G>A"][50 - 1]
    d["G>A_diff"] = d["G>A(1)"] - d["G>A(15)"]

    d["known_damage"] = np.mean([d["C>T_diff"], d["G>A_diff"]])

    df_known_damage = pd.DataFrame.from_dict(d, orient="index").T

    return df_known_damage


def get_df_known_damage():

    paths = Path("true-damage").glob("truedamage-*.txt")

    dfs = []
    for path in paths:
        dfs.append(get_df_known_damage_single_path(path))

    df_known_damage = (
        pd.concat(dfs, axis=0, ignore_index=True)
        .sort_values(simulation_columns)
        .reset_index(drop=True)
    )

    return df_known_damage


def get_known_damage(
    df_known_damage,
    sim_damage,
    sim_species,
    sim_length,
):

    query = (
        "sim_damage == @sim_damage"
        " & sim_species == @sim_species"
        " & sim_length == @sim_length"
    )
    return df_known_damage.query(query).iloc[0]["known_damage"]


# df_known_damage = get_df_known_damage()


#%%


def load_df(data_dir, path_parquet):

    if isinstance(data_dir, str):
        data_dir = Path(data_dir)

    results_dir = data_dir / "results"

    all_paths = list(results_dir.glob("*.parquet"))

    if path_parquet.exists():
        df_previous = pd.read_parquet(path_parquet)

        previous_names = set(df_previous["sample"])

        paths_to_load = []
        for path in all_paths:
            name = path.stem.removesuffix(".results")
            if name not in previous_names:
                paths_to_load.append(path)
    else:
        df_previous = None
        paths_to_load = all_paths

    if len(paths_to_load) == 0:
        return df_previous

    out = []
    for path in tqdm(paths_to_load):
        out.append(pd.read_parquet(path))

    if df_previous is None:
        df = pd.concat(out, ignore_index=True)
    else:
        df = pd.concat([df_previous, pd.concat(out)], ignore_index=True)

    return df


ALL_SPECIES = [
    "homo",
    "betula",
    "GC-low",
    "GC-mid",
    "GC-high",
    "contig1k",
    "contig10k",
    "contig100k",
]


def get_data_dir(specie):
    if specie in ALL_SPECIES:
        return Path("data") / specie
    raise AssertionError(f"Unknown specie: {specie}")


def get_damaged_reads_path(specie):
    if specie in ALL_SPECIES:
        return f"damaged_reads_{specie}.txt"
    raise AssertionError(f"Unknown specie: {specie}")


def load_results(specie=ALL_SPECIES, use_columns_subset=True):

    data_dir = get_data_dir(specie)
    path_parquet = data_dir / "df.parquet"
    df = load_df(data_dir, path_parquet)
    if df is None:
        return None
    df.to_parquet(path_parquet)

    df["Bayesian_significance"] = df["Bayesian_D_max"] / df["Bayesian_D_max_std"]

    columns = [
        "sample",
        "tax_id",
        "N_reads",
        "D_max",
        "D_max_std",
        "Bayesian_D_max",
        "Bayesian_D_max_std",
        "significance",
        "Bayesian_significance",
        "Bayesian_prob_lt_5p_damage",
        "Bayesian_prob_lt_2p_damage",
        "Bayesian_prob_lt_1p_damage",
        "Bayesian_prob_lt_0.1p_damage",
        "Bayesian_prob_zero_damage",
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
        "N_x=1_reverse",
        "N_sum_total",
        "N_sum_forward",
        "N_sum_reverse",
        "N_min",
        "k_sum_total",
        "k_sum_forward",
        "k_sum_reverse",
        "non_CT_GA_damage_frequency_mean",
        "non_CT_GA_damage_frequency_std",
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
        "Bayesian_z",
        "var_L",
        "var_GC",
        "f+1",
        "f+15",
        "f-1",
        "f-15",
    ]

    if use_columns_subset:
        df = df.loc[:, columns]

    df[simulation_columns] = (
        df["sample"]
        .apply(split_name_pd)
        .astype(
            {
                "sim_N_reads": "int",
                "sim_length": "int",
                "sim_seed": "int",
            }
        )
    )

    df = df.sort_values(simulation_columns).reset_index(drop=True)

    # for col in simulation_columns:
    #     df[col] = df[col].astype("category")

    df["sim_species"] = df["sim_species"].astype("category")

    df["sim_damage_percent_approx"] = df["sim_damage"].map(D_DAMAGE_APPROX)

    df["Bayesian_prob_not_zero_damage"] = 1 - df["Bayesian_prob_zero_damage"]
    df["Bayesian_prob_gt_1p_damage"] = 1 - df["Bayesian_prob_lt_1p_damage"]

    return df


#%%


#%%


def load_multiple_species(species=ALL_SPECIES):

    if not (isinstance(species, list) or isinstance(species, tuple)):
        species = [species]

    dfs = []
    for specie in species:
        df = load_results(specie)
        if df is not None:
            dfs.append(df)
    return pd.concat(dfs, axis=0, ignore_index=True)


#%%


def get_df_damaged_reads(path):

    if isinstance(path, str):
        path = Path(path)

    if not path.exists():
        return None

    with open(path) as f:
        d = {}
        for line in f:
            line = line.strip()
            if line.startswith("HEADER"):
                continue

            if line.startswith("sim"):
                filename = Path(line)

                species, damage, N_reads, length, seed = split_name(filename.stem)

                d[str(filename)] = {
                    "mod0000": 0,
                    "mod1000": 0,
                    "sim_species": species,
                    "sim_damage": damage,
                    "sim_N_reads": N_reads,
                    "sim_length": length,
                    "sim_seed": seed,
                }
            else:
                counts, key, _ = line.split(" ")
                d[str(filename)][key] = int(counts)

    df_damaged_reads = (
        pd.DataFrame(d).T.sort_values(simulation_columns).reset_index(drop=True)
    )

    df_damaged_reads["frac_damaged"] = df_damaged_reads["mod1000"] / (
        df_damaged_reads["mod1000"] + df_damaged_reads["mod0000"]
    )

    return df_damaged_reads


def load_multiple_damaged_reads(species):
    if not (isinstance(species, list) or isinstance(species, tuple)):
        species = [species]

    dfs = []
    for specie in species:
        damaged_reads_path = get_damaged_reads_path(specie)
        dfs.append(get_df_damaged_reads(damaged_reads_path))

    return pd.concat(dfs, axis=0, ignore_index=True)


#%%


#%%


def mean_of_CI_halfrange(group):
    s_low = "Bayesian_D_max_confidence_interval_1_sigma_low"
    s_high = "Bayesian_D_max_confidence_interval_1_sigma_high"
    return np.mean((group[s_high] - group[s_low]) / 2)


def mean_of_CI_range_low(group):
    return group["Bayesian_D_max_confidence_interval_1_sigma_low"].mean()


def mean_of_CI_range_high(group):
    return group["Bayesian_D_max_confidence_interval_1_sigma_"].mean()


def MAE(x, y):
    return np.mean(np.abs(x - y))


def MAPE(actual, predicted):
    return np.mean(np.abs((actual - predicted) / actual))


def get_df_aggregated(
    df_in,
    df_known_damage,
    df_damaged_reads=None,
):

    dfg = df_in.groupby(["sim_species", "sim_length", "sim_damage", "sim_N_reads"])

    out = []
    for (sim_species, sim_length, sim_damage, sim_N_reads), group in dfg:
        # break

        prefix = "Bayesian_D_max"
        prefix_MAP = "D_max"

        known_damage = get_known_damage(
            df_known_damage=df_known_damage,
            sim_damage=sim_damage,
            sim_species=sim_species,
            sim_length=sim_length,
        )

        d = {
            "sim_damage": sim_damage,
            "sim_damage_percent_approx": D_DAMAGE_APPROX[sim_damage],
            "sim_species": sim_species,
            "sim_length": sim_length,
            "sim_N_reads": sim_N_reads,
            "known_damage": known_damage,
            "N_simulations": len(group),
            # Bayesian
            f"{prefix}_mean_of_mean": group[f"{prefix}"].mean(),
            f"{prefix}_mean_of_median": group[f"{prefix}_median"].mean(),
            f"{prefix}_median_of_median": group[f"{prefix}_median"].median(),
            f"{prefix}_std_of_mean": group[f"{prefix}"].std(),
            f"{prefix}_mean_of_std": group[f"{prefix}_std"].mean(),
            f"{prefix}_mean_of_CI_halfrange": mean_of_CI_halfrange(group),
            f"{prefix}_median_of_CI_range_low": group[
                f"{prefix}_confidence_interval_1_sigma_low"
            ].median(),
            f"{prefix}_median_of_CI_range_high": group[
                f"{prefix}_confidence_interval_1_sigma_high"
            ].median(),
            f"{prefix}_MAE_mean": MAE(known_damage, group[f"{prefix}"]),
            f"{prefix}_MAE_median": MAE(known_damage, group[f"{prefix}_median"]),
            f"{prefix}_MAPE_mean": MAPE(known_damage, group[f"{prefix}"]),
            f"{prefix}_MAPE_median": MAPE(known_damage, group[f"{prefix}_median"]),
            # MAP
            f"{prefix_MAP}_mean_of_mean": group[f"{prefix_MAP}"].mean(),
            f"{prefix_MAP}_std_of_mean": group[f"{prefix_MAP}"].std(),
            f"{prefix_MAP}_mean_of_std": group[f"{prefix_MAP}_std"].mean(),
            f"{prefix_MAP}_MAE_mean": MAE(group[f"{prefix}"], known_damage),
            # Fit quality, Bayesian
            f"Bayesian_z_mean": group[f"Bayesian_z"].mean(),
            f"Bayesian_z_std": group[f"Bayesian_z"].std(),
            f"Bayesian_z_sdom": sdom(group[f"Bayesian_z"]),
            # Fit quality, MAP
            f"lambda_LR_mean": group[f"lambda_LR"].mean(),
            f"lambda_LR_std": group[f"lambda_LR"].std(),
            f"lambda_LR_sdom": sdom(group[f"lambda_LR"]),
            # damaged reads
            # f"reads_damaged": df_damaged_reads.,
        }

        if df_damaged_reads is not None:

            series_damaged_reads = get_damaged_reads(
                df_damaged_reads,
                sim_species,
                sim_damage,
                sim_N_reads,
                sim_length,
            )
            series = series_damaged_reads

            if len(series) > 1:
                d["reads_damaged"] = series["mod1000"].median()
                d["reads_non_damaged"] = series["mod0000"].median()
                d["reads_damaged_fraction"] = series["frac_damaged"].mean()
            else:
                d["reads_damaged"] = np.nan
                d["reads_non_damaged"] = np.nan
                d["reads_damaged_fraction"] = np.nan
        out.append(d)

    df_aggregated = pd.DataFrame(out)

    return df_aggregated


#%%


def _from_low_high_to_errors(x_low, x_high):
    yerr = np.vstack([x_low, x_high])
    yerr_mean = yerr.mean(axis=0)
    yerr2 = yerr[1, :] - yerr_mean
    return yerr_mean, yerr2


def from_low_high_to_errors_corrected(x_low, x_center, x_high):
    # xerr = np.vstack([x_low, x_high])

    x_low = np.array(x_low)
    x_center = np.array(x_center)
    x_high = np.array(x_high)

    xerr = np.zeros((2, len(x_center)))
    xerr[0, :] = x_center - x_low
    xerr[1, :] = x_high - x_center
    return x_center, xerr


# y, sy = from_low_high_to_errors_corrected([0, 0, 0], [0.5, 0.25, 0], [1, 1, 1])
# plt.errorbar([0, 1, 2], y, yerr=sy, fmt="o")


#%%


def get_damaged_reads(
    df_damaged_reads,
    sim_species=None,
    sim_damage=None,
    sim_N_reads=None,
    sim_length=None,
):

    query = ""

    if sim_species is not None:
        query += f"and sim_species == '{sim_species}' "

    if sim_damage is not None:
        query += f"and sim_damage == {sim_damage} "

    if sim_N_reads is not None:
        query += f"and sim_N_reads == {sim_N_reads} "

    if sim_length is not None:
        query += f"and sim_length == {sim_length} "

    return df_damaged_reads.query(query[4:])


#%%


# y_limits_individual_damage = {
#     0.0: (0, 0.10),
#     0.014: (0, 0.15),
#     0.047: (0, 0.15),
#     0.138: (0, 0.15),
#     0.303: (0, 0.20),
#     0.466: (0, 0.25),
#     0.626: (0, 0.35),
#     0.96: (0, 0.60),
# }


y_limits_individual_damage = {
    0.0: (0, 0.10),
    0.035: (0, 0.15),
    0.065: (0, 0.15),
    0.162: (0, 0.15),
    0.31: (0, 0.20),
    0.472: (0, 0.25),
    0.633: (0, 0.35),
    0.96: (0, 0.60),
}


def plot_individual_damage_result(
    df_in,
    group_all_keys,
    df_known_damage,
    df_damaged_reads=None,
    xlim=None,
    ylim=None,
    figsize=(15, 15),
    method="Bayesian",
    keys=None,
    fig_title=None,
    ax_titles=True,
    ax_in=None,
    splitby="species",
    loc="upper right",
    bbox_to_anchor=(1, 1.16),
    ncols=3,
    markerscale=0.75,
):

    sim_damage = group_all_keys["sim_damage"].iloc[0]
    sim_damage_percent_approx = D_DAMAGE_APPROX[sim_damage]
    sim_N_reads = group_all_keys["sim_N_reads"].iloc[0]

    if method.lower() == "bayesian":
        prefix = "Bayesian_"
        ylabel = r"Damage"
    else:
        prefix = ""
        ylabel = r"Damage (MAP)"

    # delta = 0.1
    delta = 0.0

    if splitby.lower() == "species":
        sim_length = 60
        if keys is None:
            keys = df_in["sim_species"].unique()

    elif splitby.lower() == "lengths":
        sim_species = "homo"
        if keys is None:
            keys = df_in["sim_length"].unique()
    else:
        raise ValueError(f"splitby must be 'species' or 'lengths', not {splitby}")

    N_keys = len(keys)

    if ax_in is None:
        fig, axes = plt.subplots(figsize=figsize, nrows=N_keys, sharex=True)
    else:
        axes = ax_in

    if N_keys == 1:
        axes = [axes]

    for i, (key, ax) in enumerate(zip(keys, axes)):
        # break

        if splitby.lower() == "species":
            sim_species = key
        elif splitby.lower() == "lengths":
            sim_length = key

        query = f"sim_species == '{sim_species}' and sim_length == {sim_length}"
        group = group_all_keys.query(query)

        known_damage = get_known_damage(
            df_known_damage=df_known_damage,
            sim_damage=sim_damage,
            sim_species=sim_species,
            sim_length=sim_length,
        )

        str_mean_damaged_reads = ""
        if df_damaged_reads is not None:
            series_damaged_reads = get_damaged_reads(
                df_damaged_reads,
                sim_species=sim_species,
                sim_damage=sim_damage,
                sim_N_reads=sim_N_reads,
                sim_length=sim_length,
            )
            if len(series_damaged_reads) > 0:
                mean_damaged_reads = series_damaged_reads["frac_damaged"].mean()
                str_mean_damaged_reads = (
                    f", {mean_damaged_reads*100:.1f}"
                    r"\% damaged reads (mean) in fasta file"
                )

        ax.set(
            ylabel=ylabel,
            title=f"{splitby.capitalize()} = {key}{str_mean_damaged_reads}"
            if ax_titles
            else "",
            ylim=y_limits_individual_damage[sim_damage] if ylim is None else ylim,
        )

        ax.yaxis.set_major_formatter(mtick.PercentFormatter(1.0))

        ax.axhline(
            known_damage,
            color="k",
            linestyle="--",
            label=r"$D_\mathrm{known} = " f"{known_damage*100:.1f}" r"\%$",
        )

        if len(group) == 0:
            continue

        x = group["sim_seed"]

        ax.errorbar(
            x - delta,
            group[f"{prefix}D_max"],
            group[f"{prefix}D_max_std"],
            fmt="o",
            # color="C0",
            color="grey",
            # alpha=0.5,
            label=r"Mean $\pm$ std.",
            capsize=0,
        )

        # if method.lower() == "bayesian":

        #     y, sy = from_low_high_to_errors_corrected(
        #         group["Bayesian_D_max_confidence_interval_1_sigma_low"],
        #         group["Bayesian_D_max_median"],
        #         group["Bayesian_D_max_confidence_interval_1_sigma_high"],
        #     )

        #     mask = sy[1, :] >= 0
        #     not_mask = np.logical_not(mask)

        #     ax.errorbar(
        #         x[mask] + delta,
        #         y[mask],
        #         sy[:, mask],
        #         fmt="s",
        #         label=r"Median $\pm$ 68\% C.I.",
        #         color="C1",
        #     )

        #     y2, sy2 = _from_low_high_to_errors(
        #         group["Bayesian_D_max_confidence_interval_1_sigma_low"],
        #         group["Bayesian_D_max_confidence_interval_1_sigma_high"],
        #     )

        #     ax.plot(
        #         x[not_mask] + delta,
        #         group["Bayesian_D_max_median"].values[not_mask],
        #         "s",
        #         color="C1",
        #     )

        #     ax.errorbar(
        #         x[not_mask] + delta,
        #         y2[not_mask],
        #         sy2[not_mask],
        #         fmt="None",
        #         color="C1",
        #     )

    ax.set(
        xlabel="Iteration",
        xlim=(-0.9, 100 - 0.1) if xlim is None else xlim,
    )

    if ax_in is None:
        if fig_title is None:
            fig.suptitle(
                f"Individual damages: \n"
                f"{sim_N_reads} reads\n"
                f"Briggs damage = {sim_damage}\n"
                f"Damage percent (approx) = {sim_damage_percent_approx*100:.0f}"
                r"\% "
            )
        else:
            fig.suptitle(fig_title)
        fig.subplots_adjust(top=0.85)

    # try:
    leg_kws = dict(
        markerscale=markerscale,
        bbox_to_anchor=bbox_to_anchor,
        loc=loc,
        ncols=ncols,
    )
    handles, labels = axes[0].get_legend_handles_labels()
    order = [1, 0]
    # order = [0, 1]
    axes[0].legend(
        [handles[idx] for idx in order],
        [labels[idx] for idx in order],
        **leg_kws,
    )

    # if method.lower() == "bayesian":
    #     handles, labels = axes[0].get_legend_handles_labels()
    #     order = [1, 2, 0]
    #     axes[0].legend(
    #         [handles[idx] for idx in order],
    #         [labels[idx] for idx in order],
    #         **leg_kws,
    #     )
    # else:
    #     axes[0].legend(**leg_kws)
    # except IndexError:
    #     pass

    if ax_in is None:
        return fig


#%%


def plot_individual_damage_results(
    df,
    df_known_damage,
    df_damaged_reads=None,
    suffix="",
):

    filename = Path(f"figures/bayesian/bayesian_individual_damage_results{suffix}.pdf")
    filename.parent.mkdir(parents=True, exist_ok=True)

    with PdfPages(filename) as pdf:

        for _, group_all_species in tqdm(df.groupby(["sim_damage", "sim_N_reads"])):
            # break

            fig = plot_individual_damage_result(
                df_in=df,
                group_all_keys=group_all_species,
                df_known_damage=df_known_damage,
                df_damaged_reads=df_damaged_reads,
                splitby="species",
                method="Bayesian",
            )

            pdf.savefig(fig)
            plt.close()

    filename = Path(f"figures/MAP/MAP_individual_damage_results{suffix}.pdf")
    filename.parent.mkdir(parents=True, exist_ok=True)

    with PdfPages(filename) as pdf:

        for _, group_all_species in tqdm(df.groupby(["sim_damage", "sim_N_reads"])):

            fig = plot_individual_damage_result(
                df_in=df,
                group_all_keys=group_all_species,
                df_known_damage=df_known_damage,
                df_damaged_reads=df_damaged_reads,
                splitby="species",
                method="MAP",
            )

            pdf.savefig(fig)
            plt.close()


#%%


def plot_combined_damage_result(
    df_in,
    group_agg_all_keys,
    df_known_damage,
    df_damaged_reads=None,
    sim_length=60,
    method="Bayesian",
    all_species=None,
    figsize=(15, 15),
    fig_title=None,
    keys=None,
    xlim=None,
    ylim=None,
    ax_titles=True,
    ax_in=None,
    # delta=0.07,
    delta=0.0,
    splitby="species",
    loc="upper right",
    markerscale=0.75,
):

    sim_damage = group_agg_all_keys["sim_damage"].iloc[0]
    sim_damage_percent_approx = D_DAMAGE_APPROX[sim_damage]

    if method.lower() == "bayesian":
        prefix = "Bayesian_"
        ylabel = r"Damage"
    else:
        prefix = ""
        ylabel = r"Damage (MAP)"

    if splitby.lower() == "species":
        sim_length = 60
        if keys is None:
            keys = df_in["sim_species"].unique()

    elif splitby.lower() == "lengths":
        sim_species = "homo"
        if keys is None:
            keys = df_in["sim_length"].unique()
    else:
        raise ValueError(f"splitby must be 'species' or 'lengths', not {splitby}")

    N_keys = len(keys)

    # if all_species is None:
    # all_species = df["sim_species"].unique()
    # N_species = len(all_species)

    if ax_in is None:
        fig, axes = plt.subplots(figsize=figsize, nrows=N_keys, sharex=True)
    else:
        axes = ax_in

    if N_keys == 1:
        axes = [axes]

    for i, (key, ax) in enumerate(zip(keys, axes)):
        # break

        if splitby.lower() == "species":
            sim_species = key
            title = sim_species
        elif splitby.lower() == "lengths":
            sim_length = key
            title = sim_length

        query = f"sim_species == '{sim_species}' and sim_length == {sim_length}"
        group_agg = group_agg_all_keys.query(query)

        known_damage = get_known_damage(
            df_known_damage=df_known_damage,
            sim_damage=sim_damage,
            sim_species=sim_species,
            sim_length=sim_length,
        )

        ax.axhline(
            known_damage,
            color="k",
            linestyle="--",
            label=r"$D_\mathrm{known} = " f"{known_damage*100:.1f}" r"\%$",
        )
        ax.set_xscale("log")

        ax.set(
            title=f"{splitby.capitalize()} = {title}" if ax_titles else None,
            ylabel=ylabel,
            xlim=(0.8 * 10**1, 1.2 * 10**5) if xlim is None else xlim,
            ylim=(0, 0.48) if ylim is None else ylim,
        )

        ax.yaxis.set_major_formatter(mtick.PercentFormatter(1.0))

        if len(group_agg) == 0:
            continue

        x = group_agg["sim_N_reads"]

        ax.errorbar(
            x * (1 - delta),
            group_agg[f"{prefix}D_max_mean_of_mean"],
            group_agg[f"{prefix}D_max_mean_of_std"],
            fmt="o",
            label="Mean of mean ± mean of std",
            # color="C0",
            color="grey",
            capsize=0,
        )

        # if method.lower() == "bayesian":

        #     y, sy = from_low_high_to_errors_corrected(
        #         group_agg["Bayesian_D_max_median_of_CI_range_low"],
        #         group_agg["Bayesian_D_max_median_of_median"],
        #         group_agg["Bayesian_D_max_median_of_CI_range_high"],
        #     )

        #     mask = sy[1, :] >= 0
        #     not_mask = np.logical_not(mask)

        #     ax.errorbar(
        #         x.values[mask] * (1 + delta),
        #         y[mask],
        #         sy[:, mask],
        #         fmt="s",
        #         label=r"Median of median ± median of CI (68\%)",
        #         color="C1",
        #     )

        #     y2, sy2 = _from_low_high_to_errors(
        #         group_agg["Bayesian_D_max_median_of_CI_range_low"],
        #         group_agg["Bayesian_D_max_median_of_CI_range_high"],
        #     )

        #     ax.plot(
        #         x[not_mask] * (1 + delta),
        #         y[not_mask],
        #         "s",
        #         color="C1",
        #         # label="Median of median",
        #     )

        #     ax.errorbar(
        #         x[not_mask] * (1 + delta),
        #         y2[not_mask],
        #         sy2[not_mask],
        #         fmt="None",
        #         color="C1",
        #         # label="Median of CI (16%-84%)",
        #     )

        if df_damaged_reads is not None:

            for sim_N_reads, group_agg_N_reads in group_agg.groupby("sim_N_reads"):
                # break

                series_damaged_reads = get_damaged_reads(
                    df_damaged_reads,
                    sim_species=sim_species,
                    sim_damage=sim_damage,
                    sim_length=sim_length,
                    sim_N_reads=sim_N_reads,
                )

                str_mean_damaged_reads = ""
                if len(series_damaged_reads) > 0:
                    mean_damaged_reads = series_damaged_reads["frac_damaged"].mean()
                    str_mean_damaged_reads = f"{mean_damaged_reads*100:.1f}" + r"\%"

                top1 = (
                    group_agg_N_reads[f"{prefix}D_max_mean_of_mean"]
                    + group_agg_N_reads[f"{prefix}D_max_mean_of_std"]
                )
                if method.lower() == "bayesian":
                    top2 = group_agg_N_reads["Bayesian_D_max_median_of_CI_range_high"]
                    y = max([top1.iloc[0], top2.iloc[0]])
                else:
                    y = top1.iloc[0]

                ax.text(
                    sim_N_reads,
                    y * 1.02,
                    str_mean_damaged_reads,
                    ha="center",
                    va="bottom",
                    fontsize=6,
                )

    ax.set(xlabel="Number of reads")
    if ax_in is None:
        if fig_title is None:
            fig.suptitle(
                f"{ylabel}\n"
                f"Briggs damage = {sim_damage}\n"
                f"Damage percent (approx) = {sim_damage_percent_approx*100:.0f}"
                + r"\%",
            )
        else:
            fig.suptitle(fig_title)

    leg_kws = dict(
        loc=loc,
        markerscale=markerscale,
    )
    handles, labels = axes[0].get_legend_handles_labels()
    order = [1, 0]
    axes[0].legend(
        [handles[idx] for idx in order],
        [labels[idx] for idx in order],
        **leg_kws,
    )
    # axes[0].legend(**leg_kws)
    # try:
    #     if method.lower() == "bayesian":
    #         handles, labels = axes[0].get_legend_handles_labels()
    #         order = [1, 2, 0]
    #         axes[0].legend(
    #             [handles[idx] for idx in order],
    #             [labels[idx] for idx in order],
    #             **leg_kws,
    #         )
    #     else:
    #         handles, labels = axes[0].get_legend_handles_labels()
    #         order = [1, 0]
    #         axes[0].legend(
    #             [handles[idx] for idx in order],
    #             [labels[idx] for idx in order],
    #             **leg_kws,
    #         )
    # except IndexError:
    #     pass

    if ax_in is None:
        return fig


#%%


def plot_combined_damage_results(
    df,
    df_aggregated,
    df_known_damage,
    df_damaged_reads=None,
    suffix="",
):

    filename = Path(f"figures/bayesian/bayesian_combined_damage_results{suffix}.pdf")
    filename.parent.mkdir(parents=True, exist_ok=True)

    with PdfPages(filename) as pdf:

        for sim_damage, group_agg_all_species in tqdm(
            df_aggregated.groupby("sim_damage")
        ):
            # break

            fig = plot_combined_damage_result(
                df_in=df,
                group_agg_all_keys=group_agg_all_species,
                df_known_damage=df_known_damage,
                df_damaged_reads=df_damaged_reads,
                sim_length=60,
                method="Bayesian",
            )
            pdf.savefig(fig)
            plt.close()

    filename = Path(f"figures/MAP/MAP_combined_damage_results{suffix}.pdf")
    filename.parent.mkdir(parents=True, exist_ok=True)

    with PdfPages(filename) as pdf:

        for sim_damage, group_agg_all_species in tqdm(
            df_aggregated.groupby("sim_damage")
        ):
            # break

            fig = plot_combined_damage_result(
                df_in=df,
                group_agg_all_keys=group_agg_all_species,
                df_known_damage=df_known_damage,
                df_damaged_reads=df_damaged_reads,
                sim_length=60,
                method="MAP",
            )
            pdf.savefig(fig)
            plt.close()


#%%


def plot_combined_MAEs(
    df_aggregated_homo,
    df_known_damage,
    method="Bayesian",
    sim_length=60,
    figsize=(6, 6),
    fig_title=None,
    xlim=None,
    ylim=None,
    ax_in=None,
    # delta=0.07,
    delta=0.0,
    relative=False,
):

    if not relative:
        title = f"Mean Absolute Error"
        ylabel = "MAE"
    else:
        title = f"Mean Absolute Percentage Error"
        ylabel = "MAPE"

    if method.lower() == "bayesian":
        prefix = "Bayesian_"
        ylabel = ylabel
    else:
        prefix = ""
        ylabel = ylabel + " (MAP)"

    sim_species = "homo"
    sim_damages = df_aggregated_homo["sim_damage"].unique()

    # N_keys = len(sim_damages)

    if ax_in is None:
        fig, ax = plt.subplots(
            # figsize=figsize,
        )
    else:
        ax = ax_in

    for (i, sim_damage) in enumerate(sim_damages):
        # break

        query = f"sim_damage == {sim_damage} and sim_length == {sim_length}"
        group_agg = df_aggregated_homo.query(query)

        sim_damage_percent_approx = D_DAMAGE_APPROX[sim_damage]

        known_damage = get_known_damage(
            df_known_damage=df_known_damage,
            sim_damage=sim_damage,
            sim_species=sim_species,
            sim_length=sim_length,
        )

        ax.set_xscale("log")

        x = group_agg["sim_N_reads"]

        ax.plot(
            x * (1 - delta),
            group_agg[f"{prefix}D_max_{ylabel}_mean"],
            "-",
            color=f"C{i}",
            label=f"{sim_damage} ({sim_damage_percent_approx*100:.0f}" + r"\%)",
        )

    ax.set(
        title=title,
        ylabel=ylabel,
        xlim=(0.8 * 10**1, 1.2 * 10**5) if xlim is None else xlim,
        ylim=(0, None) if ylim is None else ylim,
    )

    if relative:
        ax.yaxis.set_major_formatter(mtick.PercentFormatter(1.0))

    ax.set(xlabel="Number of reads")

    leg_kws = dict(loc="upper right", markerscale=0.75)
    ax.legend(**leg_kws)

    fig.tight_layout()

    return fig


#%%


def get_df_frac(df_in, column, cut):

    out = []
    for (sim_species, sim_damage, sim_N_reads), group in df_in.groupby(
        ["sim_species", "sim_damage", "sim_N_reads"]
    ):

        numerator = (group[column] > cut).sum()
        denominator = len(group)
        frac = numerator / denominator

        out.append(
            {
                "sim_species": sim_species,
                "sim_damage": sim_damage,
                "sim_damage_percent_approx": D_DAMAGE_APPROX[sim_damage],
                "sim_N_reads": sim_N_reads,
                "log10_sim_N_reads": np.log10(sim_N_reads),
                "frac": frac,
            }
        )

    df_fracs = pd.DataFrame(out)
    return df_fracs


#%%


def get_n_sigma_probability(n_sigma):
    return sp_norm.cdf(n_sigma) - sp_norm.cdf(-n_sigma)


CUTS = [2, 3, 4]


def get_contour_settings(cut_type, cuts=None, method="bayesian"):

    if method.lower() == "bayesian":
        column = "Bayesian_significance"
        title = "Bayesian Significance"
    else:
        column = "significance"
        title = "Significance (MAP)"

    if cut_type == "significance":
        contour_settings = {
            "column": column,
            "label_title": "Significance cut:",
            "figure_title": title,
            "cuts": CUTS if cuts is None else cuts,
            "cut_transform": lambda x: x,
            "label_template": lambda cut: f"{cut} " r"$\sigma$",
        }

    elif cut_type == "prob_not_zero_damage":
        contour_settings = {
            "column": "Bayesian_prob_not_zero_damage",
            "label_title": r"$P(D > 0\%)$ cut:",
            "figure_title": r"$P(D > 0\%)$",
            "cuts": CUTS if cuts is None else cuts,
            "cut_transform": get_n_sigma_probability,
            "label_template": lambda cut: f"{cut} " r"$\sigma$",
            # "cuts": [0.9, 0.95, 0.99, 0.999],
            # "label_template": lambda cut: f"{cut:.1%}",
        }

    elif cut_type == "prob_gt_1p_damage":
        contour_settings = {
            "column": "Bayesian_prob_gt_1p_damage",
            "label_title": r"$P(D > 1\%)$ cut:",
            "figure_title": r"$P(D > 1\%)$",
            "cuts": CUTS if cuts is None else cuts,
            "cut_transform": get_n_sigma_probability,
            "label_template": lambda cut: f"{cut} " r"$\sigma$",
            # "cuts": [0.9, 0.95, 0.99, 0.999],
            # "cuts": [get_n_sigma_probability(cut) for cut in [2, 3, 4]],
            # "label_template": lambda cut: f"{cut:.1%}",
        }

    else:
        raise ValueError(f"Unknown cut_type: {cut_type}")

    contour_settings["colors"] = ["C0", "C1", "C2", "C3", "C4"]
    contour_settings["levels"] = [0.5, 0.95]
    contour_settings["alphas"] = [0.5, 1.0]
    contour_settings["linestyles"] = ["dashed", "solid"]
    return contour_settings


#%%


# fig, ax = plt.subplots()

# contour_settings = get_contour_settings(
#     cut_type,
#     cuts=cuts,
#     method=method,
# )


def plot_contour_lines_on_ax(
    df,
    ax,
    contour_settings=None,
    sim_species="homo",
    gaussian_noise=None,
    method="Bayesian",
    title=None,
    frac_legend_pos=(1, 0.74),
    fill_contour=False,
    fill_cut_value=4,
    fill_level_value=0.95,
):

    if contour_settings is None:
        contour_settings = get_contour_settings(
            "significance",
            method="Bayesian",
        )

    for cut, color in zip(contour_settings["cuts"], contour_settings["colors"]):
        # break

        cut_transformed = contour_settings["cut_transform"](cut)

        df_fracs = get_df_frac(
            df.query(f"sim_species == '{sim_species}'"),
            column=contour_settings["column"],
            cut=cut_transformed,
        )

        df_wide = pd.pivot(
            df_fracs,
            index="sim_damage_percent_approx",
            columns="sim_N_reads",
            values="frac",
        )

        df_wide

        if gaussian_noise is None:
            data = df_wide.values

        else:
            data = gaussian_filter(df_wide.values, gaussian_noise)

        for level, alpha, linestyle in zip(
            contour_settings["levels"],
            contour_settings["alphas"],
            contour_settings["linestyles"],
        ):
            CS_ = ax.contour(
                df_wide.columns,
                df_wide.index,
                data,
                levels=[level],
                alpha=alpha,
                colors=color,
                linestyles=linestyle,
            )
            if cut == fill_cut_value and level == fill_level_value:
                CS = CS_
                fill_color = color

        ax.plot(
            [np.nan, np.nan],
            [np.nan, np.nan],
            color=color,
            label=contour_settings["label_template"](cut),
        )

    if fill_contour:
        xs, ys = CS.allsegs[0][0].T
        d = 0.3 * np.ones(len(ys))
        ax.fill_between(xs, ys, d, interpolate=True, color=fill_color, alpha=0.1)

    ax.set(
        ylabel="Damage",
        xlabel="Number of reads",
        xscale="log",
        title=f"Species = {sim_species}" if title is None else title,
    )
    ax.yaxis.set_major_formatter(mtick.PercentFormatter(1.0))

    # if i == N_species - 1:
    # if i == 0:
    # if i >= 0:

    ax.legend(
        loc="upper right",
        bbox_to_anchor=(1, 0.99),
        frameon=False,
        title=contour_settings["label_title"],
        alignment="right",
    )

    ax2 = ax.twinx()
    for level, alpha, linestyle in zip(
        contour_settings["levels"],
        contour_settings["alphas"],
        contour_settings["linestyles"],
    ):
        ax2.plot(
            np.nan,
            np.nan,
            ls=linestyle,
            label=f"{level*100:.0f}" + r"\%",
            c="black",
            alpha=alpha,
        )

    ax2.get_yaxis().set_visible(False)

    ax2.legend(
        loc="upper right",
        bbox_to_anchor=frac_legend_pos,
        frameon=False,
        title="Sim. fraction:",
        alignment="right",
    )


# %%


def plot_contour_line(
    df,
    cut_type,
    gaussian_noise=None,
    cuts=None,
    method="Bayesian",
):

    contour_settings = get_contour_settings(
        cut_type,
        cuts=cuts,
        method=method,
    )

    all_species = df["sim_species"].unique()
    N_species = len(all_species)

    fig, axes = plt.subplots(figsize=(20, 5), ncols=N_species, sharey=True)
    for i, (sim_species, ax) in enumerate(zip(all_species, axes)):
        # break

        try:
            plot_contour_lines_on_ax(
                df,
                ax=ax,
                contour_settings=contour_settings,
                sim_species=sim_species,
                gaussian_noise=gaussian_noise,
                method=method,
            )
        except TypeError:
            pass

    fig.suptitle(contour_settings["figure_title"], fontsize=16)
    fig.subplots_adjust(top=0.85)
    return fig


#%%


def plot_contour_lines(df):

    filename = Path(f"figures/bayesian/bayesian_contours.pdf")
    filename.parent.mkdir(parents=True, exist_ok=True)

    with PdfPages(filename) as pdf:

        cut_types = [
            "significance",
            "prob_not_zero_damage",
            "prob_gt_1p_damage",
        ]

        for cut_type in cut_types:
            # break

            fig = plot_contour_line(
                df,
                cut_type,
                method="Bayesian",
                # gaussian_noise=0.5,
                # cuts=[1, 2, 3, 4],
            )

            pdf.savefig(fig)
            plt.close()

    filename = Path(f"figures/MAP/MAP_contours.pdf")
    filename.parent.mkdir(parents=True, exist_ok=True)

    with PdfPages(filename) as pdf:

        cut_types = [
            "significance",
            # "prob_not_zero_damage",
            # "prob_gt_1p_damage",
        ]

        for cut_type in cut_types:
            # break

            fig = plot_contour_line(
                df,
                cut_type,
                method="MAP",
                # gaussian_noise=0.5,
                # cuts=[1, 2, 3, 4],
            )

            pdf.savefig(fig)
            plt.close()


#%%


def plot_zero_damage_group(
    group_zero_damage,
    method="Bayesian",
    title=None,
    xlim=None,
    ylim=None,
    ax=None,
):

    if method.lower() == "bayesian":
        prefix = "Bayesian_"
        xlabel = "Significance"
        ylabel = "Damage"
        ncols = 4
    else:
        prefix = ""
        ncols = 3
        xlabel = "Significance (MAP)"
        ylabel = "Damage (MAP)"

    g = sns.jointplot(
        data=group_zero_damage,
        x=f"{prefix}significance",
        y=f"{prefix}D_max",
        height=4,
        ratio=5,
        color="k",
        alpha=0.5,
        marker="+",
        s=25,
        marginal_kws=dict(bins=30, fill=False),
        xlim=xlim,
        ylim=ylim,
        # marginal_ticks=True,
        # space=0.5,
    )
    g.plot_joint(
        sns.kdeplot,
        color="C1",
        zorder=0,
        levels=5,
        alpha=0.5,
    )
    # g.plot_marginals(
    #     sns.rugplot,
    #     color="k",
    #     height=-0.1,
    #     clip_on=False,
    #     alpha=0.92,
    # )

    g.set_axis_labels(xlabel=xlabel, ylabel=ylabel)
    g.ax_joint.yaxis.set_major_formatter(mtick.PercentFormatter(1.0))

    # g.ax_joint.grid()
    # g.ax_marg_x.set_yscale("log")
    # g.ax_marg_y.set_xscale("log")

    if title is None or title == "":
        pass
    else:
        g.fig.suptitle(title, fontsize=16)
        g.fig.subplots_adjust(top=0.85)

    g.fig.tight_layout()

    return g


#%%


def plot_zero_damage_groups(df_zero_damage):

    filename = Path(f"figures/bayesian/bayesian_zero_damage_plots.pdf")
    filename.parent.mkdir(parents=True, exist_ok=True)

    with PdfPages(filename) as pdf:

        for sim_N_reads, group_zero_damage in tqdm(
            df_zero_damage.groupby("sim_N_reads")
        ):
            # break

            fig = plot_zero_damage_group(
                group_zero_damage,
                method="Bayesian",
                title=rf"sim_N_reads = {sim_N_reads}, \# = {len(group_zero_damage)}",
            )
            pdf.savefig(fig.fig, transparent=True)
            plt.close()

    filename = Path(f"figures/MAP/MAP_zero_damage_plots.pdf")
    filename.parent.mkdir(parents=True, exist_ok=True)

    with PdfPages(filename) as pdf:

        for sim_N_reads, group_zero_damage in tqdm(
            df_zero_damage.groupby("sim_N_reads")
        ):

            # if sim_N_reads == 25000:
            #     break

            fig = plot_zero_damage_group(
                group_zero_damage,
                method="MAP",
                title=rf"sim_N_reads = {sim_N_reads}, \# = {len(group_zero_damage)}",
            )
            pdf.savefig(fig.fig, transparent=True)
            plt.close()


#%%


def plot_individual_damage_results_lengths(
    df_homo_99,
    df_known_damage,
    df_damaged_reads=None,
):

    filename = Path(f"figures/bayesian/bayesian_individual_damage_results_lengths.pdf")
    filename.parent.mkdir(parents=True, exist_ok=True)

    with PdfPages(filename) as pdf:

        for (sim_damage, sim_N_reads), group_all_lengths in tqdm(
            df_homo_99.groupby(["sim_damage", "sim_N_reads"])
        ):

            # break

            fig = plot_individual_damage_result(
                df_in=df_homo_99,
                group_all_keys=group_all_lengths,
                df_known_damage=df_known_damage,
                df_damaged_reads=df_damaged_reads,
                splitby="lengths",
                method="Bayesian",
            )

            pdf.savefig(fig)
            plt.close()

    filename = Path(f"figures/MAP/MAP_individual_damage_results_lengths.pdf")
    filename.parent.mkdir(parents=True, exist_ok=True)

    with PdfPages(filename) as pdf:

        for (sim_damage, sim_N_reads), group_all_lengths in tqdm(
            df_homo_99.groupby(["sim_damage", "sim_N_reads"])
        ):

            fig = plot_individual_damage_result(
                df_in=df_homo_99,
                group_all_keys=group_all_lengths,
                df_known_damage=df_known_damage,
                df_damaged_reads=df_damaged_reads,
                splitby="lengths",
                method="MAP",
            )

            pdf.savefig(fig)
            plt.close()


#%%


def plot_combined_damage_results_lengths(
    df_homo_99,
    df_aggregated_lengths,
    df_known_damage,
    df_damaged_reads=None,
):

    filename = Path(f"figures/bayesian/bayesian_combined_damage_results_lengths.pdf")
    filename.parent.mkdir(parents=True, exist_ok=True)

    with PdfPages(filename) as pdf:

        for sim_damage, group_agg_all_lengths in tqdm(
            df_aggregated_lengths.groupby("sim_damage")
        ):
            # break

            fig = plot_combined_damage_result(
                df_in=df_homo_99,
                group_agg_all_keys=group_agg_all_lengths,
                df_known_damage=df_known_damage,
                df_damaged_reads=df_damaged_reads,
                splitby="lengths",
                method="Bayesian",
            )
            pdf.savefig(fig)

            plt.close()

    filename = Path(f"figures/MAP/MAP_combined_damage_results_lengths.pdf")
    filename.parent.mkdir(parents=True, exist_ok=True)

    with PdfPages(filename) as pdf:

        for sim_damage, group_agg_all_lengths in tqdm(
            df_aggregated_lengths.groupby("sim_damage")
        ):

            fig = plot_combined_damage_result(
                df_in=df_homo_99,
                group_agg_all_keys=group_agg_all_lengths,
                df_known_damage=df_known_damage,
                df_damaged_reads=df_damaged_reads,
                splitby="lengths",
                method="MAP",
            )
            pdf.savefig(fig)

            plt.close()


# %%


#%%

# # np.set_printoptions(suppress=True)
# from scipy.interpolate import griddata


# def get_CS_locations(points, values, y_axis_positions, levels):

#     # N_reads = np.linspace(points[:, 1].min(), points[:, 1].max(), 1000)
#     # N_reads = np.linspace(0, 500, 10 + 1)
#     N_reads = np.logspace(
#         np.log10(points[:, 1].min()), np.log10(points[:, 1].max()), 1000
#     )

#     grid_x, grid_y = np.meshgrid(
#         y_axis_positions,
#         N_reads,
#         # np.log10(N_reads),
#     )

#     grid_z1 = griddata(points, values, (grid_x, grid_y), method="cubic")

#     manual_locations = []
#     for i, level in enumerate(levels):
#         # break
#         N_read_position = N_reads[np.abs(grid_z1[:, i] - level).argmin()]
#         manual_locations.append((N_read_position, y_axis_positions[i]))
#     manual_locations

#     return manual_locations


#%%

# # significance_cut = 3

# df["log10_sim_N_reads"] = np.log10(df["sim_N_reads"])
# df["log10_Bayesian_D_max_significance"] = np.log10(df["Bayesian_significance"])
# df["log10_Bayesian_prob_zero_damage"] = np.log10(df["Bayesian_prob_zero_damage"])
# df["log10_Bayesian_prob_lt_1p_damage"] = np.log10(df["Bayesian_prob_lt_1p_damage"])

# #%%


# xys = [
#     ("Bayesian_significance", "Bayesian_D_max"),
#     ("Bayesian_prob_lt_1p_damage", "Bayesian_D_max"),
#     ("Bayesian_prob_zero_damage", "Bayesian_D_max"),
#     ("Bayesian_prob_lt_1p_damage", "Bayesian_significance"),
#     ("Bayesian_prob_zero_damage", "Bayesian_significance"),
#     ("Bayesian_prob_lt_1p_damage", "Bayesian_prob_zero_damage"),
# ]


# xys = [
#     ("Bayesian_significance", "Bayesian_D_max"),
#     ("log10_Bayesian_prob_lt_1p_damage", "Bayesian_D_max"),
#     ("log10_Bayesian_prob_zero_damage", "Bayesian_D_max"),
#     ("log10_Bayesian_prob_lt_1p_damage", "Bayesian_significance"),
#     ("log10_Bayesian_prob_zero_damage", "Bayesian_significance"),
#     ("log10_Bayesian_prob_lt_1p_damage", "log10_Bayesian_prob_zero_damage"),
# ]


# for xy in tqdm(xys):

#     fig, ax = plt.subplots(figsize=(10, 6))
#     sns.scatterplot(
#         data=df,
#         x=xy[0],
#         y=xy[1],
#         hue="sim_damage_percent_approx",
#         palette="deep",
#         size="sim_N_reads",
#         legend=False,
#         sizes=(2, 100),
#         alpha=0.5,
#         ax=ax,
#     )

#     x_str = xy[0].replace("Bayesian_", "")
#     y_str = xy[1].replace("Bayesian_", "")

#     ax.set(title=f"Bayesian, {x_str} vs {y_str}", xlabel=x_str, ylabel=y_str)

#     fig.savefig(f"figures/comparison_{species}_{xy[0]}_vs_{xy[1]}.pdf")

#     # plt.close("all")


# #%%

# columns = [
#     "Bayesian_significance",
#     "log10_Bayesian_prob_lt_1p_damage",
#     "log10_Bayesian_prob_zero_damage",
#     "Bayesian_D_max",
# ]

# g = sns.PairGrid(
#     df,
#     vars=columns,
#     hue="sim_damage_percent_approx",
#     palette="deep",
#     diag_sharey=False,
#     corner=True,
# )

# g.map_diag(
#     sns.histplot,
#     log_scale=(False, True),
#     element="step",
#     fill=False,
# )
# # g.map_diag(sns.kdeplot, log_scale=(False, True))
# g.map_lower(
#     sns.scatterplot,
#     size=df["sim_N_reads"],
#     sizes=(2, 100),
#     alpha=0.5,
# )

# # g.add_legend()
# g.add_legend(
#     title="Legend:",
#     adjust_subtitles=True,
# )

# # g.tight_layout()

# g.figure.savefig(f"figures/comparison_{species}_pairgrid.pdf")


#%%


d_xlim = {
    0.0: (0, 40),
    0.01: (0, 40),
    0.02: (0, 40),
    0.05: (0, 50),
    0.1: (0, 50),
    0.15: (0, 60),
    0.20: (0, 100),
    0.30: (0, 200),
}
d_ylim = {
    0.0: (0, 0.1),
    0.01: (0, 0.2),
    0.02: (0, 0.2),
    0.05: (0, 0.2),
    0.1: (0.0, 0.4),
    0.15: (0.0, 0.6),
    0.20: (0.0, 0.6),
    0.30: (0.0, 1.0),
}


def plot_pydamage_comparison(
    group,
    df_known_damage,
    sim_damage,
    sim_species,
    sim_N_reads,
    sim_length,
    use_special_axis=False,
):

    sim_damage_percent_approx = D_DAMAGE_APPROX[sim_damage]

    xlim = d_xlim[sim_damage_percent_approx]
    ylim = d_ylim[sim_damage_percent_approx]

    known_damage = get_known_damage(
        df_known_damage=df_known_damage,
        sim_damage=sim_damage,
        sim_species=sim_species,
        sim_length=sim_length,
    )

    fig, ax = plt.subplots()

    data = group.query(f"method == 'metaDMG'")
    mask = (data.D_max >= 0.01) & (data.significance >= 2)

    ax.scatter(
        data[mask][f"significance"],
        data[mask][f"D_max"],
        s=10,
        marker="*",
        label="metaDMG (damaged)",
        color="C0",
    )
    ax.scatter(
        data[~mask][f"significance"],
        data[~mask][f"D_max"],
        s=10,
        marker="x",
        label="metaDMG (non-damaged)",
        color="C0",
    )

    data = group.query(f"method == 'pydamage'")
    mask = (data.predicted_accuracy > 0.5) & (data.qvalue < 0.05)

    ax.scatter(
        data[mask][f"significance"],
        data[mask][f"D_max"],
        s=10,
        marker="*",
        label="PyDamage (damaged)",
        color="C2",
    )

    ax.scatter(
        data[~mask][f"significance"],
        data[~mask][f"D_max"],
        s=10,
        marker="x",
        label="PyDamage (non-damaged)",
        color="C1",
    )

    ax.axhline(
        known_damage,
        color="k",
        linestyle="--",
        label=r"$D_\mathrm{known} = " f"{known_damage*100:.1f}" r"\%$",
    )

    if use_special_axis:
        ax.set_xscale("function", functions=functions_significance)
        ax.set_xticks(cuts_significance[:-1])

        ax.set_yscale("function", functions=functions_damage)
        ax.set_yticks(cuts_damage[:-1])

        xlim = (0, 1000)
        ylim = (0, 1.0)

    ax.yaxis.set_major_formatter(mtick.PercentFormatter(1.0))

    ax.set(
        xlabel="Significance (MAP)",
        ylabel="Damage (MAP)",
        xlim=xlim,
        ylim=ylim,
    )

    title = f"{sim_N_reads} reads\n" f"Briggs damage = {sim_damage}\n"
    ax.set_title(title, pad=30, fontsize=12, loc="left")

    ax.spines.right.set_visible(True)
    ax.spines.top.set_visible(True)

    leg_kws = dict(
        markerscale=1.5,
        bbox_to_anchor=(1, 1.35),
        loc="upper right",
        ncols=1,
    )

    ax.legend(**leg_kws)

    if use_special_axis:
        kwargs = dict(
            color="C2",
            linestyle="--",
            alpha=0.3,
            linewidth=1,
        )
        xlim = ax.get_xlim()
        ylim = ax.get_ylim()
        ax.plot([2, xlim[1]], [0.01, 0.01], **kwargs)
        ax.plot([2, 2], [0.01, ylim[1]], **kwargs)
        ax.set(xlim=xlim)

        ax.fill_between(
            [2, xlim[1]],
            [0.01, 0.01],
            [ylim[1], ylim[1]],
            color="C2",
            alpha=0.1,
        )

        ax.fill_between(
            [0, 2, 2, xlim[1]],
            [ylim[1], ylim[1], 0.01, 0.01],
            color="C1",
            alpha=0.1,
        )

        d = 0.5  # proportion of vertical to horizontal extent of the slanted line
        kwargs = dict(
            markersize=5,
            linestyle="none",
            color="k",
            mec="k",
            mew=1,
            clip_on=False,
        )
        ax.plot(
            [0, 0],
            [0.060, 0.055],
            marker=[(-1, d), (1, -d)],
            **kwargs,
        )

        ax.plot(
            [5, 5.5],
            [0.0, 0.0],
            marker=[(-0.5, d), (0.5, -d)],
            **kwargs,
        )

    # handles, labels = ax.get_legend_handles_labels()
    # order = [1, 0, 2]
    # ax.legend(
    #     [handles[idx] for idx in order],
    #     [labels[idx] for idx in order],
    #     **leg_kws,
    # )

    return fig


#%%


def f_forward(xs, cuts):
    out = []

    for x in xs:
        y = x
        for i_cut in range(len(cuts) - 1):
            if cuts[i_cut] <= x < cuts[i_cut + 1]:
                y = i_cut + (x - cuts[i_cut]) / (cuts[i_cut + 1] - cuts[i_cut])
                break

        out.append(y)
    return np.array(out)


def f_reverse(xs, cuts):

    cuts_x_forward = f_forward(cuts[:-1], cuts)

    out = []
    for x in xs:
        y = x
        for i_cut in range(len(cuts_x_forward) - 1):
            if cuts_x_forward[i_cut] <= x < cuts_x_forward[i_cut + 1]:
                y = (x - i_cut) * (cuts[i_cut + 1] - cuts[i_cut]) + cuts[i_cut]
                break
        out.append(y)
    return out


cuts_significance = np.array([0, 1, 2, 3, 4, 5, 10, 100, 1000, np.inf])

functions_significance = [
    lambda xs: f_forward(xs, cuts_significance),
    lambda xs: f_reverse(xs, cuts_significance),
]


cuts_damage = np.array([0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.10, 0.20, 0.5, 1, np.inf])

functions_damage = [
    lambda xs: f_forward(xs, cuts_damage),
    lambda xs: f_reverse(xs, cuts_damage),
]


#%%


def plot_pydamage_comparison_zero_damage(sim_N_reads, group_zero_damage):

    fig, axes = plt.subplots(ncols=2, sharey=False, figsize=(8, 3))
    ax1, ax2 = axes

    for method, ax in zip(["metaDMG", "pydamage"], axes):

        data = group_zero_damage.query(f"method == '{method}'")
        mask = (data[f"significance"] > 2) & (data[f"D_max"] > 0.01)

        ax.scatter(data[f"significance"][mask], data[f"D_max"][mask], s=10, color="C1")
        ax.scatter(
            data[f"significance"][~mask], data[f"D_max"][~mask], s=10, color="C2"
        )

        title = f"{method}: {sim_N_reads} reads"
        ax.set_title(title, pad=10, fontsize=12, loc="left")

        ax.set_xscale("function", functions=functions_significance)
        ax.set_xticks(cuts_significance[:-1])

        ax.set_yscale("function", functions=functions_damage)
        ax.set_yticks(cuts_damage[:-1])
        ax.yaxis.set_major_formatter(mtick.PercentFormatter(1.0))

        # ax.grid()

    xlim = (0, 1000)
    ylim = (0, 0.5)
    ax1.set(
        xlabel="Significance (MAP)",
        ylabel="Damage (MAP)",
        xlim=xlim,
        ylim=ylim,
    )
    ax2.set(
        xlabel="Significance (MAP)",
        xlim=xlim,
        ylim=ylim,
    )

    for ax in axes:

        kwargs = dict(
            color="C1",
            linestyle="--",
            alpha=0.3,
            linewidth=1,
        )
        xlim = ax.get_xlim()
        ylim = ax.get_ylim()
        ax.plot([2, xlim[1]], [0.01, 0.01], **kwargs)
        ax.plot([2, 2], [0.01, ylim[1]], **kwargs)
        ax.set(xlim=xlim)

        ax.fill_between(
            [2, xlim[1]],
            [0.01, 0.01],
            [ylim[1], ylim[1]],
            color="C1",
            alpha=0.1,
        )

        ax.fill_between(
            [0, 2, 2, xlim[1]],
            [ylim[1], ylim[1], 0.01, 0.01],
            color="C2",
            alpha=0.1,
        )

        d = 0.5  # proportion of vertical to horizontal extent of the slanted line
        kwargs = dict(
            markersize=5,
            linestyle="none",
            color="k",
            mec="k",
            mew=1,
            clip_on=False,
        )
        ax.plot(
            [0, 0],
            [0.060, 0.055],
            marker=[(-1, d), (1, -d)],
            **kwargs,
        )

        ax.plot(
            [5, 5.5],
            [0.0, 0.0],
            marker=[(-0.5, d), (0.5, -d)],
            **kwargs,
        )

    return fig


#%%


def add_colors(group):
    bad_color = -1
    mask = (group.D_max < 0.01) | (group.significance < 2)
    colors = mask * bad_color + (1 - mask) * (group.significance)
    group["colors"] = colors
    return group.sort_values("colors", ascending=True, inplace=False)


#%%


def make_parallel_plots(
    group,
    sim_damage,
    sim_N_reads,
    percentage_columns=None,
    smooth_lines=True,
    d_names=None,
    cmap=None,
    norm=None,
):

    sim_damage_percent = D_DAMAGE_APPROX[sim_damage]
    title = rf"sim_N_reads: {sim_N_reads}, Damage: {sim_damage_percent*100:.1f}\%"

    if d_names is None:
        d_names = {
            "D_max": r"$D$",
            "damage_model_pmax": r"$p_\mathrm{max}$",
            "significance": r"$Z$",
            "predicted_accuracy": r"$\mathrm{Acc}_\mathrm{pred}$",
            "qvalue": r"$q$",
        }

    if cmap is None:
        cmap = plt.get_cmap("cividis")
        cmap.set_bad(color="k", alpha=0.2)
    if norm is None:
        norm = LogNorm(vmin=2, vmax=50)

    ys = group[d_names.keys()].values
    N = len(ys)
    ys.shape

    # organize the data
    ymins = ys.min(axis=0)
    ymaxs = ys.max(axis=0)
    dys = ymaxs - ymins
    ymins -= dys * 0.05  # add 5% padding below and above
    ymaxs += dys * 0.05

    if percentage_columns is None:
        percentage_columns = {
            "D_max",
            "damage_model_pmax",
            "predicted_accuracy",
            "qvalue",
        }
    percentage_columns_indices = []
    for i, column in enumerate(d_names.keys()):
        if column in percentage_columns:
            ymins[i] = 0
            ymaxs[i] = 1
            percentage_columns_indices.append(i)

    dys = ymaxs - ymins

    # transform all data to be compatible with the main axis
    zs = np.zeros_like(ys)
    zs[:, 0] = ys[:, 0]
    zs[:, 1:] = (ys[:, 1:] - ymins[1:]) / dys[1:] * dys[0] + ymins[0]

    fig, (host, ax_c) = plt.subplots(
        figsize=(8, 4),
        nrows=2,
        gridspec_kw={"height_ratios": [20, 1]},
    )

    axes = [host] + [host.twinx() for i in range(ys.shape[1] - 1)]
    for i, ax in enumerate(axes):
        ax.set_ylim(ymins[i], ymaxs[i])
        ax.spines["top"].set_visible(False)
        ax.spines["bottom"].set_visible(False)
        if ax != host:
            ax.spines["left"].set_visible(False)
            ax.yaxis.set_ticks_position("right")
            ax.spines["right"].set_position(("axes", i / (ys.shape[1] - 1)))

    host.set_xlim(0, ys.shape[1] - 1)
    host.set_xticks(range(ys.shape[1]))
    host.set_xticklabels(d_names.values(), fontsize=14)
    host.tick_params(axis="x", which="major", pad=7)
    host.spines["right"].set_visible(False)
    host.xaxis.tick_top()
    host.set_title(title, fontsize=18, pad=30)

    for i in percentage_columns_indices:
        axes[i].yaxis.set_major_formatter(mtick.PercentFormatter(1.0))

    colors_masked = np.ma.masked_where(group["colors"] < 0, group["colors"])

    for j in range(N):
        color = cmap(norm(colors_masked[j]))

        if not smooth_lines:

            # to just draw straight lines between the axes:
            host.plot(
                range(ys.shape[1]),
                zs[j, :],
                color=color,
            )
        else:
            # create bezier curves
            # for each axis, there will a control vertex at the point itself, one at 1/3rd towards the previous and one
            #   at one third towards the next axis; the first and last axis have one less control vertex
            # x-coordinate of the control vertices: at each integer (for the axes) and two inbetween
            # y-coordinate: repeat every point three times, except the first and last only twice
            xs = np.linspace(0, len(ys) - 1, len(ys) * 3 - 2, endpoint=True)
            verts = list(zip([x for x in xs], np.repeat(zs[j, :], 3)[1:-1]))

            codes = [mpl_Path.MOVETO] + [mpl_Path.CURVE4 for _ in range(len(verts) - 1)]
            path = mpl_Path(verts, codes)
            patch = patches.PathPatch(
                path,
                facecolor="none",
                lw=1,
                edgecolor=color,
            )
            host.add_patch(patch)

            # to show the control points of the beziers
            # for x, y in verts:
            #     host.plot(x, y, "go")

    colorbar_ticks = [2, 5, 10, 20, 50]
    cb1 = mpl.colorbar.ColorbarBase(
        ax_c,
        cmap=cmap,
        norm=norm,
        orientation="horizontal",
        ticks=colorbar_ticks,
    )
    cb1.set_label(r"Significance, $Z$")
    cb1.ax.set_xticklabels(colorbar_ticks)

    fig.tight_layout()

    return fig


#%%


def plot_single_aggregate_group(
    group_agg_specific,
    sim_species,
    sim_length,
    df_known_damage,
    title="",
    method="Bayesian",
    loc="upper right",
    markerscale=0.75,
    xlim=None,
    ylim=None,
):

    if method.lower() == "bayesian":
        prefix = "Bayesian_"
        ylabel = r"Damage"
    else:
        prefix = ""
        ylabel = r"Damage (MAP)"

    d_ylim = {
        0: (0, 0.2),
        0.035: (0, 0.2),
        0.065: (0, 0.2),
        0.162: (0, 0.3),
        0.310: (0, 0.3),
        0.472: (0, 0.4),
        0.633: (0, 0.4),
        0.960: (0, 0.5),
    }

    figsize = (9, 9 * np.sqrt(2))

    fig, axes = plt.subplots(figsize=figsize, nrows=4, ncols=2)
    axes = axes.flatten()

    groups = group_agg_specific.groupby("sim_damage")
    for ax, (sim_damage, group_agg) in zip(axes, groups):
        # break

        sim_damage_percent_approx = D_DAMAGE_APPROX[sim_damage]

        known_damage = get_known_damage(
            df_known_damage=df_known_damage,
            sim_damage=sim_damage,
            sim_species=sim_species,
            sim_length=sim_length,
        )

        ax.axhline(
            known_damage,
            color="k",
            linestyle="--",
            label=r"$D_\mathrm{known} = " f"{known_damage*100:.1f}" r"\%$",
        )
        ax.set_xscale("log")

        ax.set(
            # title=f"{splitby.capitalize()} = {title}" if ax_titles else None,
            title=r"$\delta_\mathrm{SS}=" f"{sim_damage:.3f}$",
            ylabel=ylabel,
            xlim=(0.8 * 10**1, 1.2 * 10**5) if xlim is None else xlim,
            ylim=d_ylim[sim_damage] if ylim is None else ylim,
        )

        ax.yaxis.set_major_formatter(mtick.PercentFormatter(1.0))

        if len(group_agg) == 0:
            continue

        x = group_agg["sim_N_reads"]

        ax.errorbar(
            x,
            group_agg[f"{prefix}D_max_mean_of_mean"],
            group_agg[f"{prefix}D_max_mean_of_std"],
            fmt="o",
            label="Mean of mean ± mean of std",
            # color="C0",
            color="grey",
            capsize=0,
        )

        ax.set_ylim((0, ax.get_ylim()[1]))

        leg_kws = dict(
            loc=loc,
            markerscale=markerscale,
        )
        handles, labels = ax.get_legend_handles_labels()
        order = [1, 0]
        ax.legend(
            [handles[idx] for idx in order],
            [labels[idx] for idx in order],
            **leg_kws,
        )

    axes[6].set(xlabel="Number of reads")
    axes[7].set(xlabel="Number of reads")

    fig.suptitle(title, y=1.01, fontsize=20)

    fig.tight_layout()

    return fig


#%%


def _plot_all_aggregate_groups_iterator(
    df_aggregated,
    df_aggregated_lengths,
    df_aggregated_contigs,
    df_known_damage,
    method="Bayesian",
):

    suffix = "" if method.lower() == "bayesian" else " (MAP)"

    for sim_species in ["homo", "betula", "GC-low", "GC-mid", "GC-high"]:

        group_agg_specific = df_aggregated.query(f"sim_species == '{sim_species}'")

        title = (
            sim_species.capitalize()
            if sim_species in ["homo", "betula"]
            else sim_species
        )

        fig = plot_single_aggregate_group(
            group_agg_specific,
            sim_species=sim_species,
            sim_length=60,
            df_known_damage=df_known_damage,
            title=title + suffix,
            method=method,
        )
        yield fig

    for sim_length in [35, 60, 90]:

        group_agg_specific = df_aggregated_lengths.query(f"sim_length == {sim_length}")

        fig = plot_single_aggregate_group(
            group_agg_specific,
            sim_species="homo",
            sim_length=sim_length,
            df_known_damage=df_known_damage,
            title=f"Fragment Length Average: {sim_length}" + suffix,
            method=method,
        )
        yield fig

    d_contig_translations = {
        "contig1k": "Contig length: 1 000",
        "contig10k": "Contig length: 10 000",
        "contig100k": "Contig length: 100 000",
    }

    for sim_species in ["contig1k", "contig10k", "contig100k"]:

        group_agg_specific = df_aggregated_contigs.query(
            f"sim_species == '{sim_species}'"
        )

        fig = plot_single_aggregate_group(
            group_agg_specific,
            sim_species=sim_species,
            sim_length=60,
            df_known_damage=df_known_damage,
            title=d_contig_translations[sim_species] + suffix,
            method=method,
        )
        yield fig


def plot_all_aggregate_groups(
    df_aggregated,
    df_aggregated_lengths,
    df_aggregated_contigs,
    df_known_damage,
):

    filename = Path(f"figures/bayesian/bayesian_combined_damage_ALL.pdf")
    filename.parent.mkdir(parents=True, exist_ok=True)

    with PdfPages(filename) as pdf:
        figs = _plot_all_aggregate_groups_iterator(
            df_aggregated,
            df_aggregated_lengths,
            df_aggregated_contigs,
            df_known_damage=df_known_damage,
            method="Bayesian",
        )

        for fig in tqdm(figs, desc="Saving Bayesian figures", total=11):
            pdf.savefig(fig, bbox_inches="tight", transparent=True)
            plt.close()

    filename = Path(f"figures/MAP/MAP_combined_damage_ALL.pdf")
    filename.parent.mkdir(parents=True, exist_ok=True)

    with PdfPages(filename) as pdf:
        figs = _plot_all_aggregate_groups_iterator(
            df_aggregated,
            df_aggregated_lengths,
            df_aggregated_contigs,
            df_known_damage=df_known_damage,
            method="MAP",
        )

        for fig in tqdm(figs, desc="Saving MAP figures", total=11):
            pdf.savefig(fig, bbox_inches="tight", transparent=True)
            plt.close()
