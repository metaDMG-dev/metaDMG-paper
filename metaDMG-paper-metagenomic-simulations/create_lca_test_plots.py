#%%
from importlib import reload
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import parse
import plotly.express as px
import plotly.graph_objects as go
from tqdm import tqdm

import plot_utils
import utils

# reload(plot_utils)


#%%

plt.style.use("plotstyle.mplstyle")


#%%

parser_comparison = parse.compile(
    "{sample}.{N_reads:Int}.{sim_name}.comparison",
    dict(Int=int),
)


def get_sample_N_reads_simulation_method(path):
    d = parser_comparison.parse(path.stem).named
    return d["sample"], d["N_reads"], d["sim_name"]


def load_df_comparisons():

    paths = sorted(Path("data/analysis_comparison").glob("*.csv"))

    dfs = []
    for path in tqdm(paths):
        # break

        sample, N_reads, simulation_method = get_sample_N_reads_simulation_method(path)

        df_comparison = pd.read_csv(path).rename(
            columns={"N_reads": "N_reads_simulated"}
        )

        try:

            df_metaDMG_results = (
                utils.load_df_metaDMG_results_all(
                    sample,
                    N_reads,
                    simulation_methods=[simulation_method],
                )[simulation_method]
                .rename(columns={"sample": "sample_name"})
                .astype({"tax_id": int})
            )

        except FileNotFoundError:
            continue

        if not (
            "Bayesian_D_max_confidence_interval_1_sigma_low"
            in df_metaDMG_results.columns
        ):
            raise AssertionError("Expected this file to be found")

        drop_cols = [
            "tax_name",
            "D_max",
            "Bayesian_D_max",
            # "lambda_LR",
            # "Bayesian_z",
            "Bayesian_significance",
            "Bayesian_prob_not_zero_damage",
            "Bayesian_prob_gt_1p_damage",
            "Bayesian_D_max_CI_low",
            "Bayesian_D_max_CI_high",
        ]

        if "f-15" in df_metaDMG_results.columns:
            drop_cols.extend(df_metaDMG_results.loc[:, "k+1":"f-15"].columns)
        else:
            drop_cols.extend(df_metaDMG_results.loc[:, "k+1":"f+15"].columns)

        tmp = df_metaDMG_results.drop(columns=drop_cols)
        df_comparison = pd.merge(df_comparison, tmp, on="tax_id")

        dfs.append(df_comparison)

    dfs = pd.concat(dfs).sort_values(
        by=["sample", "N_reads_simulated", "simulation_method", "|A|"],
        ascending=[True, True, False, False],
    )

    return dfs


#%%


df_comparisons = load_df_comparisons()


#%%

df_all = pd.read_parquet("df_all.parquet")


#%%

x = x


#%%


df_comparisons.query("tax_id == 127401")
df_all.query("tax_id == 127401")


#%%


for sample, df_sample in tqdm(df_comparisons.groupby("sample", sort=False)):
    if sample == "Pitch-6":
        break


df_sample

#%%


for sample, df_sample in tqdm(df_all.groupby("sample", sort=False)):
    if sample == "Pitch-6":
        break


df_sample

#%%

reload(plot_utils)

plot_utils.plot_comparison_across_N_reads_simulated_and_sim_method(
    df_comparisons,
    use_bayesian=True,
)
plot_utils.plot_comparison_across_N_reads_simulated_and_sim_method(
    df_comparisons,
    use_bayesian=False,
)


#%%

if False:

    reload(plot_utils)
    groups = df_comparisons.groupby(
        by=["sample", "N_reads_simulated", "simulation_method"]
    )
    for (sample, N_reads_simulated, simulation_method), df_comparison in tqdm(groups):
        plot_utils.plot_df_comparison_plt(
            df_comparison, sample, N_reads_simulated, simulation_method
        )

# %%
