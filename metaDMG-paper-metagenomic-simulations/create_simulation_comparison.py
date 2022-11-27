#%%

import gzip
import multiprocessing as mp
import platform
from importlib import reload
from pathlib import Path

import numpy as np
import pandas as pd
from tqdm import tqdm

import utils

#%%

# reload(utils)

#%%

MAX_WORKERS = 20
MAX_WORKERS = 1

is_mac = platform.system() == "Darwin"
if is_mac:
    MAX_WORKERS = 1


path_data = utils.PATH_DATA

path_alignment_files = (
    Path("/willerslev")
    / "edna"
    / "antonio"
    / "projects"
    / "metaDMG-sims"
    / "20221010"
    / "synthetic-data"
    / "results"
)
if is_mac:
    path_alignment_files = Path("input-data") / "data-pre-mapping"

simulation_methods = ["frag", "deam", "art"]

path_alignment_parquet = path_data / "analysis_alignment"
path_alignment_parquet.mkdir(exist_ok=True)

path_analysis_lca = path_data / "analysis_lca"
path_analysis_lca.mkdir(exist_ok=True)

path_comparison = path_data / "analysis_comparison"
path_comparison.mkdir(exist_ok=True)

path_genome_fasta = (
    Path("/willerslev")
    / "edna"
    / "antonio"
    / "geogenetics"
    / "DBs"
    / "gtdb"
    / "r202"
    / "flavours"
    / "vanilla-organelles-virus"
    / "pkg"
    / "genomes"
    / "fasta"
)

if is_mac:
    path_genome_fasta = Path("genomes")


#%%

# reload(utils)

sample = "Pitch-6"
N_reads = 1_000_000
simulation_method = "deam"

# sample = "Lake-9"
# N_reads = 1_000_000
# simulation_method = "frag"

# pbar = tqdm(it)


#%%


_inputs = list(utils.get_sample_N_reads_simulation_method(path_data))

if is_mac:
    _inputs = [
        ("Pitch-6", 1_000_000, "deam"),
        # ("Lake-9", 1000000, "frag"),
        # ("Cave-100-forward", 1000000, "art"),
    ]


inputs = []
for (sample, N_reads, simulation_method) in _inputs:
    inputs.append(
        (
            sample,
            N_reads,
            simulation_method,
            path_alignment_files,
            path_alignment_parquet,
            path_analysis_lca,
            path_genome_fasta,
            path_comparison,
            simulation_methods,
        )
    )


if __name__ == "__main__":
    if MAX_WORKERS == 1:
        for p in tqdm(inputs):
            if is_mac:
                x = x
            utils.main(p)

    else:
        with mp.Pool(processes=MAX_WORKERS) as pool:
            results = tqdm(
                pool.imap_unordered(utils.main, inputs),
                total=len(inputs),
            )  # 'total' is redundant here but can be useful
            # when the size of the iterable is unobvious

            for result in results:
                pass


#%%

## Pitch-6.communities_read-abundances.tsv
# comm	    taxon	                                    frag_type	seq_depth	seq_depth_rounded
# Pitch-6	Homo_sapiens----NC_012920.1	                ancient	    20844	    108872.0
# Pitch-6	Homo_sapiens----NC_012920.1	                modern	    0	        0.0
# Pitch-6	Geofilum_rhodophaeum----GCF_002210225.1	    ancient	    16	        84.0
# Pitch-6	Geofilum_rhodophaeum----GCF_002210225.1	    modern	    0	        0.0


## Pitch-6.genome-compositions.tsv
# Taxon	                                    Community	Coverage	        Read_length	    Read_length_std	    Read_length_min	    Read_length_max	onlyAncient	    D_max
# Homo_sapiens----NC_012920.1	            Pitch-6	    72.96559633027523	81	            16.263479571	    30	                81	            True	        0.14912868
# Geofilum_rhodophaeum----GCF_002210225.1	Pitch-6	    0.0072481496380986	81	            17.183398659	    30	                81	            True	        0.11171601

## Pitch-6.communities.tsv
# Community	    Taxon	                                    Rank	Read_length	    Read_length_std	    Read_length_min	    Read_length_max	    Perc_rel_abund
# Pitch-6	    Homo_sapiens----NC_012920.1	                1	    81	            16.2634795716	    30	                81	                23.169172
# Pitch-6	    Geofilum_rhodophaeum----GCF_002210225.1	    83	    81	            17.1833986596	    30	                81	                0.00234860


## Pitch-6.filepaths.tsv
# Taxon	                                    TaxId	Accession	    Fasta
# Homo_sapiens----NC_012920.1	            134313	NC_012920.1	    /maps/projects/lundbeck/scratch/eDNA/DBs/gtdb/r202/flavours/vanilla-organelles-virus/pkg/genomes/fasta/NC_012920.1_genomic.fna.gz
# Geofilum_rhodophaeum----GCF_002210225.1	121400	GCF_002210225.1	/maps/projects/lundbeck/scratch/eDNA/DBs/gtdb/r202/flavours/vanilla-organelles-virus/pkg/genomes/fasta/GCF_002210225.1_genomic.fna.gz


## Pitch-6.communities.json
# "comm": "Pitch-6",
# "taxon": "Homo_sapiens----NC_012920.1",
# "rel_abund": 23.169172561638344,
# "genome_size": 16569,
# "accession": "NC_012920.1",
# "onlyAncient": true,
# "fragments_ancient": {
#     "fragments": {
#         "length": [
#             30,
#             31,
#             76,
#             81
#         ],
#         "freq": [
#             0.005098572399728076,
#             0.007672137515781295,
#             0.008594736330970186,
#             0.24084684859667865
#         ]
#     },
#     "dist_params": {
#         "scale": null,
#         "sigma": null,
#         "rnd_seed": null
#     },
#     "avg_len": 58.695882295814314,
#     "seq_depth": 20844,
#     "seq_depth_original": 20844,
#     "fold": 72.96469310157524,
#     "fold_original": 72.96559633027523
# },
# "fragments_modern": null,
# "coverage_enforced": true,
# "seq_depth": 108872

#%%

# %%


# tax_id_homo = 134_313
# # taxon_accession = "Homo_sapiens----NC_012920.1"

# df_simulation_alignment_frag.query(f"tax_id == {tax_id_homo}")
# df_simulation_alignment_deam.query(f"tax_id == {tax_id_homo}")

# df_simulation_mismatch_deam.query(f"tax_id == {tax_id_homo} & abs(position) == 1")[
#     ["position", "k", "N", "f"]
# ]

# df_metaDMG_mismatch_deam.query(f"tax_id == '{tax_id_homo}' & abs(position) == 1")["k"]
# df_metaDMG_mismatch_deam.query(f"tax_id == '{tax_id_homo}' & abs(position) == 1")["N"]

# # N_reads_alignment_deam =
# len(df_simulation_alignment_deam.query(f"tax_id == {tax_id_homo}"))

# # N_reads_lca_deam =
# df_metaDMG_results_deam.query(f"tax_id == '{tax_id_homo}'")["N_reads"]

# df_metaDMG_results_deam.query(f"tax_id == '{tax_id_homo}'")["N_alignments"]


#%%

# reload(utils)


#%%

# reload(utils)


#%%

# if False:

#     row = df_simulation_alignment_deam.query(f"tax_id == {121_400}").iloc[2]

#     seq = Seq(row.seq)
#     ref = d_reference_sequences[row.accession][row.contig_num][
#         row.reference_start : row.reference_end
#     ]
#     ref, seq, ref.reverse_complement()


# LCA filtering analysis

if False:

    bam_file = (
        Path("input-data") / "data" / f"{sample}__deamSim__{N_reads}.dedup.filtered.bam"
    )

    lca_csv_file = Path("data") / "analysis_lca" / f"{sample}.{N_reads}.deam.lca.csv"
    lca_csv_file.parent.mkdir(exist_ok=True, parents=True)

    if lca_csv_file.exists():
        df_lca_stats_min_similarity_score = pd.read_csv(lca_csv_file)

    else:
        df_lca_stats_min_similarity_score = (
            utils.compute_lca_stats_min_similarity_score(
                bam_file=bam_file,
                df_simulated_alignment=df_simulation_alignment_deam,
            )
        )
        df_lca_stats_min_similarity_score.to_csv(lca_csv_file, index=False)

    fig, ax = plt.subplots(figsize=(16, 10))
    for tax_id, df_ in df_lca_stats_min_similarity_score.groupby("tax_id"):
        ax.plot(df_["fraction_lca_stat"].values)
    ax.set_xticks(
        np.arange(len(utils.MIN_SIMILARITY_SCORES)),
        labels=utils.MIN_SIMILARITY_SCORES,
    )
    ax.set(
        xlabel="min_similarity_score",
        ylabel="fraction of reads after LCA (stat)",
        title=f"{sample} {N_reads} deam (stat)",
    )
    fig.savefig(f"figures/{sample}.{N_reads}.deam.fraction_reads_after_lca_stats.pdf")

    fig, ax = plt.subplots(figsize=(16, 10))
    for tax_id, df_ in df_lca_stats_min_similarity_score.groupby("tax_id"):
        ax.plot(df_["fraction_lca_full"].values)
    ax.set_xticks(
        np.arange(len(utils.MIN_SIMILARITY_SCORES)),
        labels=utils.MIN_SIMILARITY_SCORES,
    )
    ax.set(
        xlabel="min_similarity_score",
        ylabel="fraction of reads after LCA (full)",
        title=f"{sample} {N_reads} deam (full)",
    )
    fig.savefig(f"figures/{sample}.{N_reads}.deam.fraction_reads_after_lca_full.pdf")
