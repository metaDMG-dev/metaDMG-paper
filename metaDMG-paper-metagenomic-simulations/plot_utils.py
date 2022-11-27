#%%

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import parse
import plotly.express as px
import plotly.graph_objects as go
from engineering_notation import EngNumber
from matplotlib.backends.backend_pdf import PdfPages
from tqdm import tqdm

#%%

plt.rcParams["font.size"] = "16"


#%%


def plot_single_group(
    group,
    tax_id,
    sample,
    N_reads_simulated,
    simulation_method,
):

    fig, (ax1, ax2) = plt.subplots(nrows=2, figsize=(10, 6))

    delta = 0.15

    markersize = 12

    ax1.plot(
        1 - delta,
        group["A) f_CT (x=1)"],
        color="#C00000",
        marker="o",
        markersize=5 * np.log10(1 + group["A) N_C (x=1)"]),
        # label="A) f_CT",
        linestyle="None",
    )
    ax1.plot(
        1 + delta,
        group["A) f_GA (x=-1)"],
        color="#C00000",
        marker="^",
        markersize=5 * np.log10(1 + group["A) N_G (x=-1)"]),
        # label="A) f_GA",
        linestyle="None",
    )

    ax1.plot(
        2 - delta,
        group["A\B) f_CT (x=1)"],
        color="#802540",
        marker="o",
        markersize=5 * np.log10(1 + group["A\B) N_C (x=1)"]),
        # label="A\B) f_CT",
        linestyle="None",
    )
    ax1.plot(
        2 + delta,
        group["A\B) f_GA (x=-1)"],
        color="#802540",
        marker="^",
        markersize=5 * np.log10(1 + group["A\B) N_G (x=-1)"]),
        # label="A\B) f_GA",
        linestyle="None",
    ),

    ax1.plot(
        3 - delta,
        group["B) f_CT (x=1)"],
        color="#0070C0",
        marker="o",
        markersize=5 * np.log10(1 + group["B) N_C (x=1)"]),
        # label="B) f_CT",
        linestyle="None",
    )
    ax1.plot(
        3 + delta,
        group["B) f_GA (x=-1)"],
        color="#0070C0",
        marker="^",
        markersize=5 * np.log10(1 + group["B) N_G (x=-1)"]),
        # label="B) f_GA",
        linestyle="None",
    )

    ax1.plot(
        4 - delta,
        group["B\C) f_CT (x=1)"],
        color="#00859B",
        marker="o",
        markersize=5 * np.log10(1 + group["B\C) N_C (x=1)"]),
        # label="B\C) f_CT",
        linestyle="None",
    )
    ax1.plot(
        4 + delta,
        group["B\C) f_GA (x=-1)"],
        color="#00859B",
        marker="^",
        markersize=5 * np.log10(1 + group["B\C) N_G (x=-1)"]),
        # label="B\C) f_GA",
        linestyle="None",
    )

    ax1.plot(
        5 - delta,
        group["C) f_CT (x=1)"],
        color="#00B050",
        marker="o",
        markersize=5 * np.log10(1 + group["C) N_C (x=1)"]),
        # label="C) f_CT",
        linestyle="None",
    )
    ax1.plot(
        5 + delta,
        group["C) f_GA (x=-1)"],
        color="#00B050",
        marker="^",
        markersize=5 * np.log10(1 + group["C) N_G (x=-1)"]),
        # label="C) f_GA",
        linestyle="None",
    )

    ax1.plot(
        6 - delta,
        group["D) f_CT (x=1)"],
        color="#FFC000",
        marker="o",
        markersize=5 * np.log10(1 + group["D) N_C (x=1)"]),
        # label="D) f_CT",
        linestyle="None",
    )
    ax1.plot(
        6 + delta,
        group["D) f_GA (x=-1)"],
        color="#FFC000",
        marker="^",
        markersize=5 * np.log10(1 + group["D) N_G (x=-1)"]),
        # label="D) f_GA",
        linestyle="None",
    )

    ax1.plot(
        0,
        0,
        color="k",
        marker="o",
        markersize=10,
        linestyle="None",
        label=r"$C \rightarrow T$",
    )

    ax1.plot(
        0,
        0,
        color="k",
        marker="^",
        markersize=10,
        linestyle="None",
        label=r"$G \rightarrow A$",
    )

    ax1.set(
        ylabel="Fraction",
        ylim=(0, ax1.get_ylim()[1] * 1.1),
        xlim=(0.5, 6.5),
    )
    ax1.set_xticks(
        [1, 2, 3, 4, 5, 6],
        labels=["A", "A\B", "B", "B\C", "C", "D"],
    )

    ax1.axhline(
        group.simulated_D_max,
        label="Simulated D-max",
        color="k",
        ls="-",
    )
    ax1.axhline(
        group.Bayesian_D_max,
        label="Bayesian D-max",
        color="grey",
        ls="--",
    )
    ax1.axhline(
        group.D_max,
        label="MAP D-max",
        color="grey",
        ls="-.",
    )

    width = 0.15
    ax2.bar(
        1 - delta,
        height=group["A) N_C (x=1)"],
        color="#C00000",
        width=width,
    )
    ax2.bar(
        1 + delta,
        height=group["A) N_G (x=-1)"],
        color="#C00000",
        width=width,
    )

    ax2.bar(
        2 - delta,
        height=group["A\B) N_C (x=1)"],
        color="#802540",
        width=width,
    )
    ax2.bar(
        2 + delta,
        height=group["A\B) N_G (x=-1)"],
        color="#802540",
        width=width,
    )

    ax2.bar(
        3 - delta,
        height=group["B) N_C (x=1)"],
        color="#0070C0",
        width=width,
    )
    ax2.bar(
        3 + delta,
        height=group["B) N_G (x=-1)"],
        color="#0070C0",
        width=width,
    )

    ax2.bar(
        4 - delta,
        height=group["B\C) N_C (x=1)"],
        color="#00859B",
        width=width,
    )
    ax2.bar(
        4 + delta,
        height=group["B\C) N_G (x=-1)"],
        color="#00859B",
        width=width,
    )

    ax2.bar(
        5 - delta,
        height=group["C) N_C (x=1)"],
        color="#00B050",
        width=width,
    )
    ax2.bar(
        5 + delta,
        height=group["C) N_G (x=-1)"],
        color="#00B050",
        width=width,
    )

    ax2.bar(
        6 - delta,
        height=group["D) N_C (x=1)"],
        color="#FFC000",
        width=width,
    )
    ax2.bar(
        6 + delta,
        height=group["D) N_G (x=-1)"],
        color="#FFC000",
        width=width,
    )

    ax2.set(
        ylabel="Counts (N)",
        ylim=(0, ax2.get_ylim()[1] * 1.1),
        xlim=(0.5, 6.5),
    )

    ax2.set_xticks(
        [1, 2, 3, 4, 5, 6],
        labels=["A", "A\B", "B", "B\C", "C", "D"],
    )

    ax1.legend(
        ncol=1,
        loc="upper right",
        bbox_to_anchor=(1.1, 1.65),
        fontsize=10,
    )

    N_ancient = EngNumber(int(group.simulated_seq_depth_ancient))
    N_modern = EngNumber(int(group.simulated_seq_depth_modern))

    s = (
        f"Only ancient: {group.simulated_only_ancient}"
        "\n"
        f"Ancient:Modern = {N_ancient} : {N_modern}"
    )

    ax1.text(
        -0.146,
        1.45,
        s,
        horizontalalignment="left",
        transform=ax1.transAxes,
        fontsize=10,
    )

    N_reads_simulated_str = EngNumber(int(N_reads_simulated))

    title = (
        f"{sample}, {N_reads_simulated_str}, {simulation_method}"
        "\n"
        f"tax_id={tax_id}"
        "\n"
        f"tax_name={group['tax_name']}"
    )
    fig.suptitle(
        title,
        fontsize=16,
    )
    fig.subplots_adjust(
        top=0.80,
    )

    return fig


#%%


def plot_df_comparison_plt(
    df_comparison,
    sample,
    N_reads_simulated,
    simulation_method,
    use_tqdm=False,
):

    fig_name = (
        Path("figures")
        / "individual-comparison"
        / f"{sample}.{N_reads_simulated}.{simulation_method}.individual-comparison.pdf"
    )
    fig_name.parent.mkdir(exist_ok=True)

    groupby = df_comparison.groupby("tax_id", sort=False)
    if use_tqdm:
        groupby = tqdm(groupby)

    with mpl.rc_context({"text.usetex": False}):

        with PdfPages(fig_name) as pdf:

            for tax_id, group in groupby:
                # break

                if len(group) != 1:
                    raise ValueError(f"{tax_id} has {len(group)} rows")

                group = group.iloc[0]

                fig = plot_single_group(
                    group,
                    tax_id,
                    sample,
                    N_reads_simulated,
                    simulation_method,
                )

                pdf.savefig(fig)
                plt.close()


#%%


def map_N_reads_simulated_to_x_axis(xs, delta=0):
    d_translate = {
        1000000: 1,
        5000000: 2,
        10000000: 3,
        50000000: 4,
        100000000: 5,
    }
    return [d_translate[x] + delta for x in xs]


def plot_single_comparison_across_N_reads_simulated_and_sim_method(
    df_tmp,
    tax_id,
    use_bayesian=True,
):

    if use_bayesian:
        D_max_col = "Bayesian_D_max"
        # D_max_std_col = "Bayesian_D_max_std"
        D_max_std_col = [
            "Bayesian_D_max_confidence_interval_1_sigma_low",
            "Bayesian_D_max_confidence_interval_1_sigma_high",
        ]
    else:
        D_max_col = "D_max"
        D_max_std_col = "D_max_std"

    d_colors = {"frag": "C2", "deam": "C0", "art": "C1"}

    fig, ax = plt.subplots(figsize=(10, 6))

    DELTA = 0.20

    methods = ["frag", "deam", "art"]
    deltas = [-DELTA, 0, DELTA]

    for method, delta in zip(methods, deltas):
        # break

        df_method = df_tmp.query(f"simulation_method == '{method}'").set_index(
            "N_reads_simulated",
            drop=False,
        )

        if len(df_method) == 0:
            continue

        xx = map_N_reads_simulated_to_x_axis(df_method.N_reads_simulated, delta)
        yy = df_method[D_max_col]

        if use_bayesian:
            sy = (df_method[D_max_std_col].T - yy).abs()
        else:
            sy = df_method[D_max_std_col]

        ax.errorbar(
            xx,
            yy,
            sy,
            fmt="o",
            capsize=10,
            label=method,
            color=d_colors[method],
        )

        if use_bayesian:
            sy = sy.values[1, :]

        for xi, yi, si, N_i in zip(xx, yy, sy, df_method["|D|"]):
            # break
            ax.text(
                xi,
                (yi + si) * 1.001 + 0.001,
                f"{EngNumber(N_i)}",
                horizontalalignment="center",
                fontsize=8,
                color=d_colors[method],
            )

    simulated_D_max = df_tmp["simulated_D_max"].iloc[0]
    ax.axhline(
        simulated_D_max,
        label="simulated D-max",
        color="k",
        alpha=0.5,
        ls="--",
        zorder=0,
    )
    ax.text(
        5.55,
        simulated_D_max,
        f"{simulated_D_max:.3f}",
        horizontalalignment="left",
        verticalalignment="center",
        fontsize=10,
    )

    s = f"Only ancient: {df_tmp.simulated_only_ancient.iloc[0]}"
    ax.text(
        -0.1,
        1.1,
        s,
        horizontalalignment="left",
        verticalalignment="center",
        transform=ax.transAxes,
        fontsize=10,
    )

    ax.set_xticks(
        np.arange(5) + 1,
        labels=["1M", "5M", "10M", "50M", "100M"],
    )

    sample = df_tmp["sample"].iloc[0]
    tax_name = df_tmp["tax_name"].iloc[0]
    title = f"\n{sample}, Tax ID: {tax_id} \n{tax_name}\n"

    ax.set(
        xlim=(0.5, 5.5),
        ylim=(0, ax.get_ylim()[1] * 1.1),
        ylabel="Bayesian D-max" if use_bayesian else "D-max (MAP)",
        xlabel="Number of simulated reads",
    )

    # ax.legend(loc="upper right")

    ax.legend(
        ncol=1,
        loc="upper right",
        bbox_to_anchor=(1.1, 1.24),
        fontsize=10,
    )

    fig.suptitle(
        title,
        fontsize=16,
    )
    fig.subplots_adjust(
        top=0.83,
    )

    return fig


def plot_comparison_across_N_reads_simulated_and_sim_method(
    df_comparisons,
    use_bayesian=True,
):

    for sample, df_sample in tqdm(df_comparisons.groupby("sample", sort=False)):
        # break

        # if sample == "Cave-102":
        #     break

        # value_counts = df_sample["tax_id"].value_counts()
        # tax_ids_in_all = value_counts[value_counts > 10].index
        # tax_ids_in_all = df_sample.tax_id.unique()

        tax_ids_in_all = (
            df_sample.groupby("tax_id")
            .sum(numeric_only=True)["|D|"]
            .sort_values(ascending=False)
            .index
        )

        if len(tax_ids_in_all) > 0:

            suffix = "bayesian" if use_bayesian else "MAP"

            fig_name = (
                Path("figures")
                / "comparison-across-N_reads_simulated-and-sim_method"
                / f"{sample}.comparison-across-N_reads_simulated-and-sim_method.{suffix}.pdf"
            )

            fig_name.parent.mkdir(exist_ok=True, parents=True)

            groupby = df_sample.groupby("tax_id", sort=False)

            with PdfPages(fig_name) as pdf:

                for tax_id in tax_ids_in_all:
                    # break

                    df_tmp = df_sample.query(f"tax_id == {tax_id}")

                    fig = (
                        plot_single_comparison_across_N_reads_simulated_and_sim_method(
                            df_tmp,
                            tax_id,
                            use_bayesian=use_bayesian,
                        )
                    )
                    pdf.savefig(fig)
                    plt.close()


# %%
