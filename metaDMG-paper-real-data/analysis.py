#%%
from importlib import reload
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages
from tqdm import tqdm

#%%

save_plots = False

#%%

plt.style.use("plotstyle.mplstyle")


#%%


def load_df():

    df = (
        pd.read_parquet("data/results/")
        .reset_index(drop=True)
        .astype({"tax_id": int})
        .copy()
    )

    df["Bayesian_significance"] = df["Bayesian_D_max"] / df["Bayesian_D_max_std"]
    df["Bayesian_prob_not_zero_damage"] = 1 - df["Bayesian_prob_zero_damage"]
    df["Bayesian_prob_gt_1p_damage"] = 1 - df["Bayesian_prob_lt_1p_damage"]
    # df[sim_columns] = df["simulation_name"].apply(split_name_pd)

    return df


df_all = load_df()

good_samples = [
    "Cave-100-forward",
    # "Cave-100",
    "Cave-102",
    "Cave-22",
    "Lake-1",
    "Lake-7-forward",
    # "Lake-7",
    "Lake-9",
    "Library-0",
    "Pitch-6",
    "Shelter-39",
]

df = df_all.query("sample in @good_samples")

df_species = df.query("tax_rank == 'species'").copy()

#%%


def plot_overview(
    df_species,
    title="",
    xlims=None,
    ylims=None,
    size_normalization=5,
    figsize=(10, 6),
    samples=None,
    color_order=None,
):

    fig, ax = plt.subplots(figsize=figsize)

    # colors = [f"C{i}" for i in range(10)]

    if samples is None:
        samples = [
            "Cave-100-forward",
            "Cave-102",
            "Pitch-6",
            "Cave-22",
            "Lake-9",
            "Lake-7-forward",
            "Lake-1",
            "Shelter-39",
            "Library-0",
        ]

    if color_order is None:
        color_order = [
            "Cave-100-forward",
            "Cave-102",
            "Cave-22",
            "Lake-7-forward",
            "Lake-9",
            "Pitch-6",
            "Lake-1",
            "Shelter-39",
            "Library-0",
        ]

    d_colors = {sample: f"C{i}" for i, sample in enumerate(color_order)}

    symbols = ["o", "s", "D", "v", "^", ">", "<", "*", "x"]
    # samples = df_species["sample"].unique()

    for i_sample, sample in enumerate(samples):

        group = df_species.query("sample == @sample")

        ax.scatter(
            group["Bayesian_significance"],
            group["Bayesian_D_max"],
            s=2 + np.sqrt(group["N_reads"]) / size_normalization,
            marker=symbols[i_sample],
            # alpha=0.5,
            color=d_colors[sample],
            label=sample,
        )

    ax.set(
        # title=type_,
        xlabel="Significance",
        ylabel="Damage",
        xlim=xlims,
        ylim=ylims,
    )

    ax.yaxis.set_major_formatter(mtick.PercentFormatter(1.0))

    ax.yaxis.set_tick_params(labelbottom=True)

    handles, labels = ax.get_legend_handles_labels()
    leg = fig.legend(
        handles,
        labels,
        loc="upper center",
        ncol=3,
        bbox_to_anchor=(0.5, 1.2),
        fontsize=11,
    )

    for handle in leg.legendHandles:
        handle.set_sizes([30.0])
        handle.set_alpha(1)

    fig.suptitle(title, fontsize=16)
    fig.subplots_adjust(top=0.85)

    fig.tight_layout()

    return fig


#%%


fig = plot_overview(
    df_species,
    title="",
    xlims=(0, 30),
    ylims=(0, 0.7),
    figsize=(6, 3),
)
if save_plots:
    fig.savefig("overview-real-data.pdf", bbox_inches="tight")


fig_zoom = plot_overview(
    df_species, title="", xlims=(1, 4), ylims=(0, 0.1), size_normalization=1
)
if save_plots:
    fig_zoom.savefig("overview-real-data-zoom.pdf")


#%%


# fig_shelter = plot_overview(
#     df_species.query("sample == 'Shelter-39'"),
#     title="",
#     xlims=(0, 5),
#     ylims=(0, 0.4),
#     figsize=(6, 3),
# )
# if save_plots:
#     fig_shelter.savefig("overview-real-data-shelter.pdf", bbox_inches="tight")


#%%

df_species_library0 = df_species.query("sample == 'Library-0'")


mask = (df_species_library0.Bayesian_D_max > 0.01) & (
    df_species_library0.Bayesian_significance > 2
)
mask.sum()


mask2 = (df_species_library0.Bayesian_D_max > 0.02) & (
    df_species_library0.Bayesian_significance > 3
)
mask2.sum()

df_species_library0[mask]

mask.mean()


#%%


def compare_Bayesian_MAP(
    df_species_in,
    title="",
    damage_limits=None,
    significance_limits=None,
    size_normalization=5,
    figsize=(9, 4),
):

    samples = [
        "Cave-100-forward",
        "Cave-102",
        "Pitch-6",
        "Cave-22",
        "Lake-9",
        "Lake-7-forward",
        "Lake-1",
        "Shelter-39",
        "Library-0",
    ]

    color_order = [
        "Cave-100-forward",
        "Cave-102",
        "Cave-22",
        "Lake-7-forward",
        "Lake-9",
        "Pitch-6",
        "Lake-1",
        "Shelter-39",
        "Library-0",
    ]

    d_colors = {sample: f"C{i}" for i, sample in enumerate(color_order)}

    symbols = ["o", "s", "D", "v", "^", ">", "<", "*", "x"]

    fig, axes = plt.subplots(figsize=figsize, ncols=2)
    for i_sample, sample in enumerate(samples):

        group = df_species_in.query("sample == @sample")

        axes[0].scatter(
            group["Bayesian_D_max"],
            group["D_max"],
            s=2 + np.sqrt(group["N_reads"]) / size_normalization,
            marker=symbols[i_sample],
            # alpha=0.5,
            color=d_colors[sample],
            label=sample,
        )

        axes[1].scatter(
            group["Bayesian_significance"],
            group["significance"],
            s=2 + np.sqrt(group["N_reads"]) / size_normalization,
            marker=symbols[i_sample],
            # alpha=0.5,
            color=d_colors[sample],
            label=sample,
        )

    corr_matrix = np.corrcoef(df_species_in["Bayesian_D_max"], df_species_in["D_max"])
    rho = corr_matrix[0, 1]
    axes[0].annotate(
        r"$\rho = " f"{rho*100:.2f}" + r"\%$",
        xy=(0.03, 0.99),
        xycoords="axes fraction",
        ha="left",
        va="top",
        size=14,
    )

    corr_matrix = np.corrcoef(
        df_species_in["Bayesian_significance"], df_species_in["significance"]
    )
    rho = corr_matrix[0, 1]
    axes[1].annotate(
        r"$\rho = " f"{rho*100:.2f}" + r"\%$",
        xy=(0.03, 0.99),
        xycoords="axes fraction",
        ha="left",
        va="top",
        size=14,
    )

    axes[0].set(
        xlabel="Damage",
        ylabel="Damage (MAP)",
        xlim=damage_limits,
        ylim=damage_limits,
    )
    axes[0].plot(
        (damage_limits[0], damage_limits[1]),
        (damage_limits[0], damage_limits[1]),
        "--",
        color="grey",
        zorder=0,
    )

    axes[1].set(
        xlabel="Significance",
        ylabel="Significance (MAP)",
        xlim=significance_limits,
        ylim=significance_limits,
    )
    axes[1].plot(
        (significance_limits[0], significance_limits[1]),
        (significance_limits[0], significance_limits[1]),
        "--",
        color="grey",
        zorder=0,
    )

    axes[0].xaxis.set_major_formatter(mtick.PercentFormatter(1.0))
    axes[0].yaxis.set_major_formatter(mtick.PercentFormatter(1.0))

    handles, labels = axes[0].get_legend_handles_labels()
    leg = fig.legend(
        handles,
        labels,
        loc="upper center",
        ncol=5,
        bbox_to_anchor=(0.5, 1.2),
        fontsize=11,
    )

    for handle in leg.legendHandles:
        handle.set_sizes([30.0])
        handle.set_alpha(1)

    fig.tight_layout()

    return fig


fig_cut = compare_Bayesian_MAP(
    df_species.query(
        "Bayesian_D_max > 0.01 & Bayesian_significance > 2 & N_reads > 100"
    ),
    title="",
    damage_limits=(0, 0.8),
    significance_limits=(0, 30),
    size_normalization=10,
    figsize=(8, 3),
)

if save_plots:
    fig_cut.savefig("comparison-real-data-cutted.pdf", bbox_inches="tight")


fig_all = compare_Bayesian_MAP(
    df_species,
    title="",
    damage_limits=(0, 0.8),
    significance_limits=(0, 30),
    size_normalization=10,
    figsize=(8, 3),
)

if save_plots:
    fig_all.savefig("comparison-real-data-all.pdf", bbox_inches="tight")


#%%


d_name_to_age = {
    "Library-0": 10 / 1000,
    "Cave-100-forward": 100,
    "Cave-102": 102,
    "Cave-22": 22,
    "Lake-1": 1,
    "Lake-7-forward": 7,
    "Lake-9": 9,
}


df_species.loc[:, "age"] = df_species["sample"].map(d_name_to_age)


samples = [
    "Library-0",
    # "Cave-100-forward",
    "Cave-102",
    "Cave-22",
    # "Lake-1",
    # "Lake-7-forward",
    # "Lake-9",
]

df_samples = pd.concat(
    [
        df_species.query(
            f"sample in {samples[1:]}"
            " and Bayesian_D_max > 0.01"
            " and Bayesian_significance > 3"
            " and N_reads > 200"
        ),
        df_species.query("sample == 'Library-0' and N_reads > 1000"),
    ],
)


if False:
    for path in df_samples.tax_path:
        print(path)


fig_samples = plot_overview(
    df_samples,
    title="",
    xlims=(
        df_samples.Bayesian_significance.min() * 0.9,
        df_samples.Bayesian_significance.max() * 1.1,
    ),
    ylims=(df_samples.Bayesian_D_max.min(), df_samples.Bayesian_D_max.max() * 1.1),
    figsize=(6, 3),
    samples=samples,
    color_order=samples,
)


#%%

fig, ax = plt.subplots()
ages = []

xmax = df_samples.Bayesian_D_max.max() * 1.1
ymax = max(df_samples.age) * 1.1

for sample, group in df_samples.groupby("sample"):
    damage = group.Bayesian_D_max
    age = group.age
    ax.plot(damage, age, "o", label=sample)
ax.xaxis.set_major_formatter(mtick.PercentFormatter(1.0))
ax.set(xlabel="Damage", ylabel="Age (kyr)", ylim=(0, ymax), xlim=(0, xmax))

from sklearn.linear_model import LinearRegression

w = 1 / df_samples.Bayesian_D_max_std**2

reg = LinearRegression(fit_intercept=True).fit(
    df_samples.Bayesian_D_max.array.reshape(-1, 1),
    df_samples.age,
    # sample_weight=w.values,
)
reg.coef_
reg.intercept_

ax.plot(
    [0, xmax],
    [reg.intercept_, reg.intercept_ + xmax * reg.coef_[0]],
    "--",
    color="gray",
    label=f"$A = {reg.coef_[0]:.1f} " r"\cdot D" f" + {reg.intercept_:.1f}$",
)

ax.legend()

#%%

fig, ax = plt.subplots(figsize=(5, 3))
ax.plot(df_samples.Bayesian_D_max, df_samples.age, "o")

beta = 0.1
alpha = 50

xmin = df_samples.Bayesian_D_max.min()
xmax = df_samples.Bayesian_D_max.max()
ymin = df_samples.age.min()
ymax = df_samples.age.max()

damages = np.linspace(xmin, xmax, 1000)
y = beta * np.exp(alpha * damages)
y = beta * np.exp(alpha * damages)


ax.plot(damages, y, "-")

ax.set(
    xlim=(0, xmax * 1.1),
    ylim=(0, ymax * 1.1),
)

#%%


fig, ax = plt.subplots(figsize=(5, 3))
ax.plot(df_samples.Bayesian_D_max, df_samples.age, "o")

beta = 0.01
alpha = 150

xmin = df_samples.Bayesian_D_max.min()
xmax = df_samples.Bayesian_D_max.max()
ymin = df_samples.age.min()
ymax = df_samples.age.max()

damages = np.linspace(xmin, xmax, 1000)
y = beta * np.exp(alpha * damages)
y = beta * np.exp(alpha * damages)


ax.plot(damages, y, "-")

ax.set(
    xlim=(0, xmax * 1.1),
    ylim=(0, ymax * 1.1),
)


# %%


import jax.numpy as jnp
from jax import jit, random, vmap
from jax.random import PRNGKey as Key
from numba import njit
from numpyro import distributions as dist
from numpyro.infer import MCMC, NUTS, Predictive, log_likelihood
from scipy.special import logsumexp
from scipy.stats import beta as sp_beta
from scipy.stats import betabinom as sp_betabinom
from scipy.stats import norm as sp_norm

numpyro.enable_x64()


#%%


def numpyro_model_linear(damage, age=None):
    alpha = numpyro.sample("alpha", dist.Normal(1000, 1000))
    beta = numpyro.sample("beta", dist.Normal(0, 100))
    y = alpha * damage + beta
    numpyro.sample("obs", dist.Normal(y, y / 10), obs=age)


def numpyro_model(damage, age=None):
    # alpha = numpyro.sample("alpha", dist.Normal(200, 10))
    alpha = numpyro.sample("alpha", dist.Normal(100, 10))
    # beta = numpyro.sample("beta", dist.Normal(0, 0.01))
    beta = numpyro.sample("beta", dist.Normal(0, 0.1))
    # c = numpyro.sample("c", dist.Normal(0, 10))
    # y = jnp.exp(alpha * damage)
    y = beta * jnp.exp(alpha * damage)
    # y = c + beta * jnp.exp(alpha * damage)
    # y = alpha * damage + beta
    numpyro.sample("obs", dist.Normal(y, y / 10 + 1), obs=age)


def _init_mcmc(model, **kwargs):

    mcmc_kwargs = dict(
        progress_bar=True,
        num_warmup=1000,
        num_samples=1000,
        num_chains=4,  # problem when setting to 2
        chain_method="vectorized",
        # http://num.pyro.ai/en/stable/_modules/numpyro/infer/mcmc.html#MCMC
    )

    return MCMC(NUTS(model), jit_model_args=True, **mcmc_kwargs, **kwargs)


def init_mcmc(**kwargs):
    mcmc = _init_mcmc(numpyro_model, **kwargs)
    return mcmc


def fit_mcmc(mcmc, data, seed=0):
    mcmc.run(Key(seed), **data)


mcmc = init_mcmc()

data = {"damage": df_samples.Bayesian_D_max.values, "age": df_samples.age.values}
fit_mcmc(mcmc, data)

mcmc.print_summary(prob=0.68)

samples_1 = mcmc.get_samples()

#%%

from numpyro.diagnostics import hpdi


def plot_regression(x, y_mean, y_hpdi, title=""):
    # Sort values for plotting by x axis
    idx = jnp.argsort(x)
    damage = x[idx]
    mean = y_mean[idx]
    hpdi = y_hpdi[:, idx]
    age = df_samples.age.values[idx]

    # Plot
    fig, ax = plt.subplots(figsize=(5, 3))
    ax.plot(damage, mean)
    ax.plot(damage, age, "o")
    ax.fill_between(damage, hpdi[0], hpdi[1], alpha=0.3, interpolate=True)
    ax.set(
        xlabel="Damage",
        ylabel="Age",
        title=title,
        xlim=(0, None),
        ylim=(0, age.max() * 1.1),
    )
    return ax


rng_key = random.PRNGKey(0)
rng_key, rng_key_ = random.split(rng_key)
predictive = Predictive(numpyro_model, samples_1)
predictions = predictive(rng_key_, damage=df_samples.Bayesian_D_max.values)["obs"]

mean_pred = jnp.mean(predictions, axis=0)
hpdi_pred = hpdi(predictions, 0.9)

ax = plot_regression(
    df_samples.Bayesian_D_max.values,
    mean_pred,
    hpdi_pred,
    title="Posterior with 90\% CI",
)


#%%

# Compute empirical posterior distribution over mu
posterior_mu = (
    jnp.expand_dims(samples_1["beta"], -1)
    + jnp.expand_dims(samples_1["alpha"], -1) * df_samples.Bayesian_D_max.values
)

mean_mu = jnp.mean(posterior_mu, axis=0)
hpdi_mu = hpdi(posterior_mu, 0.9)

ax = plot_regression(
    df_samples.Bayesian_D_max.values,
    mean_mu,
    hpdi_mu,
    title="Mu with 90\% CI",
)

#%%

rng_key = random.PRNGKey(0)
rng_key, rng_key_ = random.split(rng_key)
prior_predictive = Predictive(numpyro_model, num_samples=100)
prior_predictions = prior_predictive(
    rng_key_,
    damage=df_samples.Bayesian_D_max.values,
)["obs"]

mean_prior_pred = jnp.mean(prior_predictions, axis=0)
hpdi_prior_pred = hpdi(prior_predictions, 0.9)

ax = plot_regression(
    df_samples.Bayesian_D_max.values,
    mean_prior_pred,
    hpdi_prior_pred,
    title="Prior with 90\% CI",
)

#%%


#%%

import arviz as az

dims = {"damage": ["damage"], "damage": ["damage"]}
idata_kwargs = {
    "dims": dims,
    "constant_data": {"damage": df_samples.Bayesian_D_max.values},
}

idata = az.from_numpyro(
    mcmc,
    prior=prior_predictive(
        rng_key_,
        damage=df_samples.Bayesian_D_max.values,
    ),
    posterior_predictive=predictive(rng_key_, damage=df_samples.Bayesian_D_max.values),
    **idata_kwargs,
)

az.plot_trace(idata, compact=False)

# %%
az.plot_pair(
    idata,
    # var_names=["mu", "theta"],
    kind=["scatter", "kde"],
    kde_kwargs={"fill_last": False},
    marginals=True,
    # coords=coords,
    point_estimate="median",
    figsize=(11.5, 5),
)

# %%
