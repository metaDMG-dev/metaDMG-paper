#%%
import gzip
import shlex
import subprocess
from collections import Counter, defaultdict
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
import parse
from Bio import SeqIO
from Bio.Seq import Seq
from tqdm import tqdm

#%%

PATH_DATA = Path("data")


D_SIM_NAME_TRANSLATE = {
    "frag": "fragSim",
    "deam": "deamSim",
    "art": "art",
}

#%%


def nth_repl(s, old, new, n):
    """helper function to find the nth occurrence of a substring `old` in the
    original string `s` and replace it with `new`."""
    find = s.find(old)
    # If find is not -1 we have found at least one match for the substring
    i = find != -1
    # loop until we find the nth or we find no match
    while find != -1 and i != n:
        # find + 1 means we start searching from after the last match
        find = s.find(old, find + 1)
        i += 1
    # If i is equal to n we found nth match so replace
    if i == n and i <= len(s.split(old)) - 1:
        return s[:find] + new + s[find + len(old) :]
    return s


# "Abiotrophia_sp001815865----GCF_001815865.1"

parser_template = (
    "{sample}___"
    "{taxon}____"
    "{accession:Accession}__"
    "{contig_num}---"
    "{read_num:Int}:"
    "{ancient_modern}:"
    "{strand}:"
    "{reference_start:Int}:"
    "{reference_end:Int}:"
    "{fragment_length:Int}:"
    "{damaged_positions_in_fragment:Fragment}"
)


def fix_accession(accession):
    if accession.startswith("_"):
        accession = accession[1:]
    accession_fixed = nth_repl(accession, "_", ".", 2)
    return accession_fixed


def fragment_parser(damaged_positions):
    if damaged_positions == "None":
        return []
    elif "," not in damaged_positions:
        # return str([int(damaged_positions)])
        return [int(damaged_positions)]
    else:
        # return str([int(string_pos) for string_pos in damaged_positions.split(",")])
        return [int(string_pos) for string_pos in damaged_positions.split(",")]


parser = parse.compile(
    parser_template,
    dict(Accession=fix_accession, Int=int, Fragment=fragment_parser),
)


#%%


def strip(str_: str) -> str:
    """
    :param str_: a string
    """
    return str_.strip()


# def load_names(names_file: str | Path) -> pd.DataFrame:
def load_names(names_file):
    """
    load names.dmp and convert it into a pandas.DataFrame.
    Taken from https://github.com/zyxue/ncbitax2lin/blob/master/ncbitax2lin/data_io.py
    """

    df_data = pd.read_csv(
        names_file,
        sep="|",
        header=None,
        index_col=False,
        names=["tax_id", "name_txt", "unique_name", "name_class"],
    )

    return (
        df_data.assign(
            name_txt=lambda df: df["name_txt"].apply(strip),
            unique_name=lambda df: df["unique_name"].apply(strip),
            name_class=lambda df: df["name_class"].apply(strip),
        )
        .loc[lambda df: df["name_class"] == "scientific name"]
        .reset_index(drop=True)
    )


names_file = Path("names.dmp")
df_names = load_names(names_file)
df_names


# def load_nodes(nodes_file: str | Path) -> pd.DataFrame:
def load_nodes(nodes_file):
    """
    load nodes.dmp and convert it into a pandas.DataFrame
    """

    df_data = pd.read_csv(
        nodes_file,
        sep="|",
        header=None,
        index_col=False,
        names=["tax_id", "parent_tax_id", "rank"],
    )

    return df_data.assign(rank=lambda df: df["rank"].apply(strip))


nodes_file = Path("nodes.dmp")
df_nodes = load_nodes(nodes_file)


# def load_acc2taxid(acc2taxid_file: str | Path) -> pd.DataFrame:
def load_acc2taxid(acc2taxid_file):
    """
    load acc2taxid.map and convert it into a pandas.DataFrame
    """
    df_data = pd.read_csv(acc2taxid_file, sep="\t")
    return df_data


acc2taxid_file = Path("acc2taxid.map.gz")
df_acc2taxid = load_acc2taxid(acc2taxid_file)


def get_key2val_dict(df_acc2taxid, key_col, val_col):
    return df_acc2taxid[[key_col, val_col]].set_index(key_col).to_dict()[val_col]


def propername(name):
    return "_".join(name.split(" "))


d_acc2taxid = get_key2val_dict(df_acc2taxid, "accession", "taxid")
d_taxid2acc = get_key2val_dict(df_acc2taxid, "taxid", "accession")
d_taxid2name = get_key2val_dict(df_names, "tax_id", "name_txt")


def acc2name(accesion):
    taxid = d_acc2taxid[accesion]
    name = d_taxid2name[taxid]
    return propername(name)


#%%


def get_simulation_alignment_paths(path_alignment_files, name, N_reads):
    path = path_alignment_files / name / "single" / str(N_reads) / "reads"

    paths = {
        "frag": path / f"{name}_fragSim.fa.gz",
        "deam": path / f"{name}_deamSim.fa.gz",
        "art": path / f"{name}_art.fq.gz",
    }

    return paths


#%%


def fix_sim_name(
    sim_name,
    d_sim_name_translate=D_SIM_NAME_TRANSLATE,
):
    d_translate = {value: key for key, value in d_sim_name_translate.items()}
    return d_translate[sim_name]


def get_sample_N_reads_simulation_method(
    path_data=PATH_DATA,
    d_sim_name_translate=D_SIM_NAME_TRANSLATE,
):

    parser_results_path = parse.compile(
        "{sample}__{sim_name:SimName}__{N_reads:Int}.results",
        dict(Int=int, SimName=lambda s: fix_sim_name(s, d_sim_name_translate)),
    )

    for path in (path_data / "results").glob("*.parquet"):
        d = parser_results_path.parse(path.stem).named
        yield d["sample"], d["N_reads"], d["sim_name"]


#%%


def extract_unique_read_name(result):
    return (
        f"{result['sample']}___"
        f"{result['taxon']}____"
        f"{result['accession'].replace('.1', '_1')}__"
        f"{result['contig_num']}---"
        f"{result['read_num']}"
    )


# path_alignment = path_alignment_deam


def load_simulation_alignment(
    path_alignment,
    # d_reference_sequences,
    feather_file,
    max_position=15,
    use_tqdm=False,
    fix_forward=True,
):

    # feather_file_mismatch = str(feather_file).replace(".feather", ".mismatch.feather")

    try:
        df = pd.read_feather(feather_file)
        # df_mismatch = pd.read_feather(feather_file_mismatch)
        # return df, df_mismatch
        return df

    except FileNotFoundError:
        pass

    if ".fq" in path_alignment.name:
        alignment_type = "fastq"
    elif ".fa" in path_alignment.name:
        alignment_type = "fasta"
    else:
        raise AssertionError(f"Unknown alignment type: {path_alignment.suffix}")

    if fix_forward and "-forward" in path_alignment.name:
        path_alignment = Path(str(path_alignment).replace("-forward", ""))

    d_counter = defaultdict(lambda: Counter())
    results = []
    with gzip.open(path_alignment, "rt") as handle:

        records = SeqIO.parse(handle, alignment_type)
        if use_tqdm:
            records = tqdm(records)

        for record in records:
            # break

            name = record.name
            seq = record.seq
            result = parser.parse(name).named
            result["seq"] = str(seq)

            result["tax_id"] = d_acc2taxid[result["accession"]]
            result["taxon_accession"] = f"{result['taxon']}----{result['accession']}"
            result["read_name"] = extract_unique_read_name(result)
            results.append(result)

            # ref = d_reference_sequences[result["accession"]][result["contig_num"]][
            #     result["reference_start"] : result["reference_end"]
            # ]

            # update_d_counter_bang(
            #     d_counter=d_counter[result["tax_id"]],
            #     seq=seq,
            #     ref=ref,
            #     strand=result["strand"],
            #     max_position=max_position,
            # )

    df = pd.DataFrame(results)

    categories = [
        "sample",
        "taxon",
        "accession",
        "contig_num",
        "ancient_modern",
        "strand",
        "tax_id",
        "taxon_accession",
    ]
    for category in categories:
        df[category] = df[category].astype("category")

    df.loc[:, "index_original"] = df.index

    df.to_feather(feather_file)

    # mismatches = []
    # for tax_id, d in d_counter.items():
    #     mismatches.append(_d_counter_to_dataframe(d, tax_id))
    # df_mismatch = pd.concat(mismatches).reset_index(drop=True)
    # df_mismatch = fix_mismatch_df(df_mismatch)
    # df_mismatch.to_feather(feather_file_mismatch)

    # return df, df_mismatch
    return df


#%%


def load_simulation_alignment_all(
    path_simulation_alignment_all,
    path_alignment_parquet,
    simulation_methods,
    sample,
    N_reads,
):

    df_simulation_alignment_all = {}

    for sim_method in simulation_methods:

        df_simulation_alignment = load_simulation_alignment(
            path_alignment=path_simulation_alignment_all[sim_method],
            feather_file=path_alignment_parquet
            / f"{sample}.{N_reads}.{sim_method}.feather",
        )

        df_simulation_alignment_all[sim_method] = df_simulation_alignment

    return df_simulation_alignment_all


#%%


def load_df_metaDMG_mismatch_all(
    sample,
    N_reads,
    simulation_methods,
    path_data=PATH_DATA,
    d_sim_name_translate=D_SIM_NAME_TRANSLATE,
):

    df_metaDMG_mismatch_all = {}

    for sim_method in simulation_methods:

        sim_name = d_sim_name_translate[sim_method]

        df_metaDMG_mismatch = pd.read_parquet(
            path_data
            / "mismatches"
            / f"{sample}__{sim_name}__{N_reads}.mismatches.parquet"
        )

        df_metaDMG_mismatch_all[sim_method] = df_metaDMG_mismatch

    return df_metaDMG_mismatch_all


#%%


def load_df_metaDMG_results_all(
    sample,
    N_reads,
    simulation_methods,
    path_data=PATH_DATA,
    d_sim_name_translate=D_SIM_NAME_TRANSLATE,
):

    df_metaDMG_results_all = {}

    for sim_method in simulation_methods:

        sim_name = d_sim_name_translate[sim_method]

        df_metaDMG_results = pd.read_parquet(
            path_data / "results" / f"{sample}__{sim_name}__{N_reads}.results.parquet"
        )

        df_metaDMG_results["significance"] = (
            df_metaDMG_results["D_max"] / df_metaDMG_results["D_max_std"]
        )
        df_metaDMG_results["Bayesian_significance"] = (
            df_metaDMG_results["Bayesian_D_max"]
            / df_metaDMG_results["Bayesian_D_max_std"]
        )

        df_metaDMG_results["Bayesian_prob_not_zero_damage"] = (
            1 - df_metaDMG_results["Bayesian_prob_zero_damage"]
        )
        df_metaDMG_results["Bayesian_prob_gt_1p_damage"] = (
            1 - df_metaDMG_results["Bayesian_prob_lt_1p_damage"]
        )

        df_metaDMG_results["Bayesian_D_max_CI_low"] = df_metaDMG_results[
            "Bayesian_D_max_confidence_interval_1_sigma_low"
        ]
        df_metaDMG_results["Bayesian_D_max_CI_high"] = df_metaDMG_results[
            "Bayesian_D_max_confidence_interval_1_sigma_high"
        ]

        df_metaDMG_results_all[sim_method] = df_metaDMG_results

    return df_metaDMG_results_all


#%%


parser_template_lca = (
    parser_template + ":{seq}:{L:Int}:{N_alignments:Int}:{GC:Float}\t{lca}"
)

parser_lca = parse.compile(
    parser_template_lca,
    dict(Accession=fix_accession, Int=int, Fragment=fragment_parser, Float=float),
)


def _read_metaDMG_lca_file(lca_file):

    lca_parsed_results = []
    with gzip.open(lca_file, "rt") as a_file:
        for i, line in enumerate(a_file):
            if i == 0:
                continue

            # break

            stripped_line = line.strip()
            parsed_result = parser_lca.parse(stripped_line).named
            parsed_result["read_name"] = extract_unique_read_name(parsed_result)
            parsed_result["tax_id_simulated"] = d_acc2taxid[parsed_result["accession"]]

            # tax_ids = [s.split(":")[0] for s in parsed_result["lca"].split("\t")]
            # parsed_result["lca_tax_ids"] = "|".join(tax_ids)
            tax_id_lca = int(parsed_result["lca"].split("\t")[0].split(":")[0])
            parsed_result["tax_id_lca"] = tax_id_lca
            del parsed_result["lca"]

            lca_parsed_results.append(parsed_result)

    df_lca = pd.DataFrame(lca_parsed_results)
    df_lca.loc[:, "index_original"] = df_lca.index

    categories = [
        "sample",
        "taxon",
        "accession",
        "contig_num",
        "ancient_modern",
        "strand",
        "N_alignments",
        "tax_id_simulated",
        "tax_id_lca",
    ]
    for category in categories:
        df_lca[category] = df_lca[category].astype("category")

    return df_lca


def get_lca_file_path(
    sample,
    simulation_method,
    N_reads,
    path_data=PATH_DATA,
    d_sim_name_translate=D_SIM_NAME_TRANSLATE,
):
    sim_name = d_sim_name_translate[simulation_method]
    lca_file = path_data / "lca" / f"{sample}__{sim_name}__{N_reads}.lca.txt.gz"
    return lca_file


def read_metaDMG_lca_file(
    sample,
    simulation_method,
    N_reads,
    path_data=PATH_DATA,
    d_sim_name_translate=D_SIM_NAME_TRANSLATE,
):
    lca_file = get_lca_file_path(
        sample,
        simulation_method,
        N_reads,
        path_data,
        d_sim_name_translate,
    )
    df_metaDMG_mapped = _read_metaDMG_lca_file(lca_file)
    return df_metaDMG_mapped


def load_metaDMG_lca_file(
    sample,
    simulation_method,
    N_reads,
    path_analysis_lca,
    path_data=PATH_DATA,
    d_sim_name_translate=D_SIM_NAME_TRANSLATE,
):

    filename = (
        path_analysis_lca / f"{sample}.{N_reads}.{simulation_method}.lca_df.feather"
    )

    try:
        df_metaDMG_mapped = pd.read_feather(filename)
        return df_metaDMG_mapped
    except FileNotFoundError:
        pass

    df_metaDMG_mapped = read_metaDMG_lca_file(
        sample,
        simulation_method,
        N_reads,
        path_data,
        d_sim_name_translate,
    )

    df_metaDMG_mapped.to_feather(filename)

    return df_metaDMG_mapped


#%%


#%%


def _d_counter_to_dataframe(d_counter, tax_id, max_position=15):
    bases = ["A", "C", "G", "T"]
    out = []
    for pos in list(range(1, max_position + 1)) + list(
        range(-1, -max_position - 1, -1)
    ):

        if pos > 0:
            direction = "5'"
        else:
            direction = "3'"

        d_tmp = {"tax_id": tax_id, "direction": direction, "position": pos}

        for ref in bases:
            for obs in bases:
                d_tmp[ref + obs] = d_counter[(pos, ref + obs)]

        out.append(d_tmp)
    return pd.DataFrame(out)


# def update_d_counter_bang_readname_based_old(
#     d_counter,
#     seq,
#     damaged_positions_in_fragment,
#     max_position=15,
# ):

#     seq_forward = Seq(seq)
#     seq_reverse = Seq(seq[::-1])

#     for i, base in enumerate(seq_forward):
#         pos = i + 1
#         ref = base
#         if pos in damaged_positions_in_fragment:
#             ref = "C"
#         d_counter[(pos, ref + base)] += 1

#         if pos >= max_position:
#             break

#     for i, base in enumerate(seq_reverse):
#         pos = -(i + 1)
#         ref = base
#         if pos in damaged_positions_in_fragment:
#             ref = "G"
#         d_counter[(pos, ref + base)] += 1

#         if abs(pos) >= max_position:
#             break


def update_d_counter_bang(
    d_counter,
    seq,
    ref,
    strand,
    max_position=15,
):

    if strand == "-":
        ref = ref.reverse_complement()

    for i, (seq_i, ref_i) in enumerate(zip(seq, ref)):
        pos = i + 1
        d_counter[pos, ref_i + seq_i] += 1
        if pos >= max_position:
            break

    for i, (seq_i, ref_i) in enumerate(zip(reversed(seq), reversed(ref))):
        pos = -(i + 1)
        d_counter[pos, ref_i + seq_i] += 1
        if abs(pos) >= max_position:
            break


def _dataframe_single_group_to_dataframe(
    d_reference_sequences,
    df_tax_id,
    tax_id,
    max_position=15,
):

    d_counter = Counter()

    for row in df_tax_id.itertuples():

        seq = Seq(row.seq)
        ref = d_reference_sequences[row.accession][row.contig_num][
            row.reference_start : row.reference_end
        ]

        update_d_counter_bang(
            d_counter=d_counter,
            seq=seq,
            ref=ref,
            strand=row.strand,
            # damaged_positions_in_fragment=row.damaged_positions_in_fragment,
            max_position=max_position,
        )

    df_mismatch = _d_counter_to_dataframe(d_counter, tax_id)
    return df_mismatch


ACTG = ["A", "C", "G", "T"]


def get_base_columns(df):
    base_columns = []
    for column in df.columns:
        if len(column) == 2 and column[0] in ACTG and column[1] in ACTG:
            base_columns.append(column)
    return base_columns


def get_reference_columns(df, ref):
    ref_columns = []
    for column in get_base_columns(df):
        if column[0] == ref:
            ref_columns.append(column)
    return ref_columns


def add_reference_count(df, ref):
    reference_columns = get_reference_columns(df, ref)
    df[ref] = df[reference_columns].sum(axis=1)
    return df


def compute_fraction_and_uncertainty(x, N, set_zero_to_nan=False):
    f = x / N
    if set_zero_to_nan:
        f = f.mask(x == 0, np.nan)
    sf = np.sqrt(f * (1 - f) / N)
    return f, sf


def compute_error_rates(df, ref, obs):
    s_ref_obs = ref + obs
    x = df[s_ref_obs]
    N_ref = df[ref]
    f, sf = compute_fraction_and_uncertainty(x, N_ref)
    return f, sf


def add_error_rate(df, ref, obs, include_uncertainties=False):
    f, sf = compute_error_rates(df, ref, obs)
    df[f"f_{ref}{obs}"] = f
    if include_uncertainties:
        df[f"sf_{ref}{obs}"] = sf
    return df


def add_k_N_x_names(df):
    mask_forward = df.position > 0
    df["k"] = np.where(mask_forward, df["CT"], df["GA"])
    df["N"] = np.where(mask_forward, df["C"], df["G"])
    df["f"] = df["k"] / df["N"]
    df["|x|"] = np.abs(df["position"])
    return df


def add_N_reads(df_mismatch):
    df_mismatch.loc[:, "N_reads"] = df_mismatch.loc[:, "AA":"TT"].sum(axis=1)
    return df_mismatch


def fix_mismatch_df(df_mismatch):
    return (
        df_mismatch.pipe(add_reference_count, ref="C")
        .pipe(add_reference_count, ref="G")
        .pipe(add_error_rate, ref="C", obs="T")
        .pipe(add_error_rate, ref="G", obs="A")
        .pipe(add_k_N_x_names)
        .pipe(add_N_reads)
    )


def dataframe_to_mismatch(
    d_reference_sequences,
    df,
    max_position=15,
    tax_id_col="tax_id",
):

    mismatches = []
    for tax_id, df_tax_id in df.groupby(tax_id_col, observed=True):
        # break
        df_mismatch = _dataframe_single_group_to_dataframe(
            d_reference_sequences=d_reference_sequences,
            df_tax_id=df_tax_id,
            tax_id=tax_id,
            max_position=max_position,
        )
        mismatches.append(df_mismatch)

    df_mismatch = fix_mismatch_df(pd.concat(mismatches))

    return df_mismatch


#%%


def delete_temp_lca_files():
    for path in Path(".").glob("outnames*"):
        path.unlink()

    for path in Path(".").glob("*.bam.bin"):
        path.unlink()


MIN_SIMILARITY_SCORES = [
    0.5,
    0.7,
    0.8,
    0.9,
    0.93,
    0.94,
    0.95,
    0.955,
    0.96,
    0.965,
    0.97,
    0.975,
    0.98,
    0.985,
    0.99,
    0.9925,
    0.995,
    0.9975,
    1.0,
]


# min_similarity_score = 0.9
def compute_lca_stats_min_similarity_score(
    bam_file,
    df_alignment,
    min_similarity_scores=None,
):

    if min_similarity_scores is None:
        min_similarity_scores = MIN_SIMILARITY_SCORES

    d_out = []
    # min_similarity_score = 0.95
    for min_similarity_score in tqdm(min_similarity_scores):

        delete_temp_lca_files()

        command = (
            f"./metaDMG-cpp lca "
            f"-bam {bam_file} "
            f"-names names.dmp "
            f"-nodes nodes.dmp "
            f"-acc2tax acc2taxid.map.gz "
            f"-simscorelow {min_similarity_score} "
            f"-simscorehigh 1.0 "
            f"-minmapq 0 "
            f"-howmany 15 "
            f"-weighttype 1 "
            f"-fix_ncbi 0 "
            f"-lca_rank species"
        )

        run_lca = subprocess.run(
            shlex.split(command),
            stdout=subprocess.DEVNULL,
            check=True,
        )

        lca_stats_cols = [
            "tax_id",
            "N_reads",
            "L_mean",
            "L_var",
            "GC_mean",
            "GC_var",
            "name",
            "rank",
        ]
        df_lca_stats = pd.read_csv("outnames.stat", sep="\t", names=lca_stats_cols)

        lca_file = "outnames.lca.gz"
        d_lca_full = read_lca_file(lca_file)

        for tax_id, df_species in df_alignment.groupby("tax_id"):
            # break
            N_reads_seq_depth = len(df_species)

            df_lca_stats_tax_id = df_lca_stats.query(f"tax_id == {tax_id}")
            assert len(df_lca_stats_tax_id) == 1
            N_reads_lca_stat = df_lca_stats_tax_id.iloc[0]["N_reads"]
            name = df_lca_stats_tax_id.iloc[0]["name"]

            N_reads_lca_full = len(d_lca_full[tax_id])

            d_out.append(
                {
                    "min_similarity_score": min_similarity_score,
                    "tax_id": tax_id,
                    "name": name,
                    "N_reads_seq_depth": N_reads_seq_depth,
                    "N_reads_lca_stat": N_reads_lca_stat,
                    "N_reads_lca_full": N_reads_lca_full,
                }
            )

    delete_temp_lca_files()

    df = pd.DataFrame(d_out)
    df["fraction_lca_stat"] = df["N_reads_lca_stat"] / df["N_reads_seq_depth"]
    df["fraction_lca_full"] = df["N_reads_lca_full"] / df["N_reads_seq_depth"]
    return df


#%%


# def load_communities_read_abundances(path_alignment_files, name, N_reads):

#     path_communities_read_abundances = (
#         path_alignment_files
#         / name
#         / "single"
#         / str(N_reads)
#         / f"{name}.communities_read-abundances.tsv"
#     )
#     return pd.read_csv(path_communities_read_abundances, sep="\t")


# #%%


# def load_genome_composition(path_alignment_files, name, N_reads):

#     path_genome_composition = (
#         path_alignment_files
#         / name
#         / "single"
#         / str(N_reads)
#         / f"{name}.genome-compositions.tsv"
#     )
#     return pd.read_csv(path_genome_composition, sep="\t")


#%%


# def extract_simulation_parameters(name, N_reads):

#     df_communities_read_abundances = load_communities_read_abundances(name, N_reads)
#     df_genome_composition = load_genome_composition(name, N_reads)

#     out = []
#     for (taxon_accession, _df) in df_communities_read_abundances.groupby(
#         "taxon", sort=False
#     ):
#         # break
#         assert len(_df) == 2

#         taxon, accession = taxon_accession.split("----")
#         tax_id = d_acc2taxid[accession]

#         seq_depth_ancient = _df.query("frag_type == 'ancient'").iloc[0]["seq_depth"]
#         seq_depth_modern = _df.query("frag_type == 'modern'").iloc[0]["seq_depth"]

#         s_genome_composition = df_genome_composition.query(
#             f"Taxon == '{taxon_accession}'"
#         )
#         assert len(s_genome_composition) == 1
#         s_genome_composition = s_genome_composition.iloc[0]
#         only_ancient = s_genome_composition["onlyAncient"]
#         D_max_simulation = s_genome_composition["D_max"]

#         out.append(
#             {
#                 "tax_id": tax_id,
#                 "taxon": taxon,
#                 "accession": accession,
#                 "seq_depth_ancient": seq_depth_ancient,
#                 "seq_depth_modern": seq_depth_modern,
#                 "fraction_ancient": seq_depth_ancient
#                 / (seq_depth_ancient + seq_depth_modern),
#                 "only_ancient": only_ancient,
#                 "D_max_simulation": D_max_simulation,
#             }
#         )
#     return pd.DataFrame(out)


#%%


def load_single_genome(reference_fasta_path):
    reference_sequence = {}
    with gzip.open(reference_fasta_path, "rt") as handle:
        for i, record in enumerate(SeqIO.parse(handle, "fasta")):
            reference_sequence[f"seq-{i}"] = record.seq
    return reference_sequence


# def _load_reference_genomes():
#     d_reference_sequences = {}
#     for reference_fasta_path in Path("genomes").glob("*_genomic.fna.gz"):
#         accession = reference_fasta_path.stem.split("_genomic.fna")[0]
#         d_reference_sequences[accession] = load_single_genome(reference_fasta_path)
#     return d_reference_sequences


def _load_reference_genomes(
    path_genome_fasta,
    df_simulation_alignment_all,
):
    accessions = list(df_simulation_alignment_all["frag"]["accession"].unique())
    d_reference_sequences = {}
    for accession in accessions:
        reference_fasta_path = path_genome_fasta / f"{accession}_genomic.fna.gz"
        d_reference_sequences[accession] = load_single_genome(reference_fasta_path)
    return d_reference_sequences


def load_reference_genomes(
    path_genome_fasta,
    df_simulation_alignment_all,
    sample,
    N_reads,
):

    outdir = Path("genomes")
    outdir.mkdir(exist_ok=True)

    filename = outdir / f"{sample}.{N_reads}.reference_genomes.pkl"
    try:
        d_reference_sequences = joblib.load(filename)
        return d_reference_sequences
    except FileNotFoundError:
        pass

    d_reference_sequences = _load_reference_genomes(
        path_genome_fasta,
        df_simulation_alignment_all,
    )
    joblib.dump(d_reference_sequences, filename)

    return d_reference_sequences


#%%


def _add_counts_to_dict_bang(d, df, name):
    df_forward = df.query("position == 1")
    d[f"{name}) k_CT (x=1)"] = df_forward["k"].iloc[0]
    d[f"{name}) N_C (x=1)"] = df_forward.query("position == 1")["N"].iloc[0]
    d[f"{name}) f_CT (x=1)"] = df_forward.query("position == 1")["f"].iloc[0]

    df_reverse = df.query("position == -1")
    if len(df_reverse) > 0:
        d[f"{name}) k_GA (x=-1)"] = df_reverse["k"].iloc[0]
        d[f"{name}) N_G (x=-1)"] = df_reverse.query("position == -1")["N"].iloc[0]
        d[f"{name}) f_GA (x=-1)"] = df_reverse.query("position == -1")["f"].iloc[0]


def _add_counts_to_dict_bang_empty(d, name):
    d[f"{name}) k_CT (x=1)"] = 0
    d[f"{name}) N_C (x=1)"] = 0
    d[f"{name}) f_CT (x=1)"] = 0
    d[f"{name}) k_GA (x=-1)"] = 0
    d[f"{name}) N_G (x=-1)"] = 0
    d[f"{name}) f_GA (x=-1)"] = 0


def compute_mismatch_and_add_to_dict(
    d,
    df_in,
    name,
    d_reference_sequences,
    tax_id_col,
):

    d[f"|{name}|"] = len(df_in)

    if len(df_in) > 0:

        df_mismatch = dataframe_to_mismatch(
            d_reference_sequences,
            df=df_in,
            tax_id_col=tax_id_col,
        )
        _add_counts_to_dict_bang(d, df_mismatch.query("abs(position) == 1"), name)
        return df_mismatch

    else:
        _add_counts_to_dict_bang_empty(d, name)
        return None


def compute_comparison(
    df_simulation_alignment,
    df_metaDMG_mapped,
    df_metaDMG_mismatch,
    df_metaDMG_results,
    d_reference_sequences,
    use_tqdm=False,
):

    colnames_int = [
        "read_num",
        "reference_start",
        "reference_end",
        "fragment_length",
        "index_original",
    ]
    dtype = {col: "int" for col in colnames_int}

    out = []

    tax_ids = df_simulation_alignment.tax_id.unique()
    if use_tqdm:
        tax_ids = tqdm(tax_ids)

    for tax_id in tax_ids:
        # break

        # tax_id = 134313  # homo

        d = {
            "tax_id": tax_id,
            "tax_name": d_taxid2name[tax_id],
        }

        df_A = df_simulation_alignment.query(f"tax_id == {tax_id}")
        df_A_mismatch = compute_mismatch_and_add_to_dict(
            d=d,
            df_in=df_A,
            name="A",
            d_reference_sequences=d_reference_sequences,
            tax_id_col="tax_id",
        )
        df_A
        len(df_A)

        df_B = df_metaDMG_mapped.query(f"tax_id_simulated == {tax_id}")
        df_B_mismatch = compute_mismatch_and_add_to_dict(
            d=d,
            df_in=df_B,
            name="B",
            d_reference_sequences=d_reference_sequences,
            tax_id_col="tax_id_simulated",
        )
        df_B
        len(df_B)

        df_C = df_metaDMG_mapped.query(f"tax_id_lca == {tax_id}")
        df_C_mismatch = compute_mismatch_and_add_to_dict(
            d=d,
            df_in=df_C,
            name="C",
            d_reference_sequences=d_reference_sequences,
            tax_id_col="tax_id_lca",
        )
        df_C
        len(df_C)

        # in df_A, but not in df_B
        df_A_slash_B = (
            df_A.merge(df_B["read_name"], how="outer", indicator=True)
            .loc[lambda x: x["_merge"] == "left_only"]
            .drop("_merge", axis=1)
        ).astype(dtype)
        df_A_slash_B_mismatch = compute_mismatch_and_add_to_dict(
            d=d,
            df_in=df_A_slash_B,
            name="A\B",
            d_reference_sequences=d_reference_sequences,
            tax_id_col="tax_id",
        )
        df_A_slash_B

        # in df_A, but not in df_C
        df_A_slash_C = (
            df_A.merge(df_C["read_name"], how="outer", indicator=True)
            .loc[lambda x: x["_merge"] == "left_only"]
            .drop("_merge", axis=1)
        ).astype(dtype)
        df_A_slash_C_mismatch = compute_mismatch_and_add_to_dict(
            d=d,
            df_in=df_A_slash_C,
            name="A\C",
            d_reference_sequences=d_reference_sequences,
            tax_id_col="tax_id",
        )
        df_A_slash_C

        # in df_B, but not in df_C
        df_B_slash_C = (
            df_B.merge(df_C["read_name"], how="outer", indicator=True)
            .loc[lambda x: x["_merge"] == "left_only"]
            .drop("_merge", axis=1)
        ).astype(dtype)
        df_B_slash_C_mismatch = compute_mismatch_and_add_to_dict(
            d=d,
            df_in=df_B_slash_C,
            name="B\C",
            d_reference_sequences=d_reference_sequences,
            tax_id_col="tax_id_simulated",
        )
        df_B_slash_C

        # in df_C, but not in df_A
        df_C_slash_A = (
            df_C.merge(df_A["read_name"], how="outer", indicator=True)
            .loc[lambda x: x["_merge"] == "left_only"]
            .drop("_merge", axis=1)
        ).astype(dtype)
        df_C_slash_A_mismatch = compute_mismatch_and_add_to_dict(
            d=d,
            df_in=df_C_slash_A,
            name="C\A",
            d_reference_sequences=d_reference_sequences,
            tax_id_col="tax_id_lca",
        )
        df_C_slash_A

        # mismatch information from metaDMG:

        columns_to_keep = [
            "D_max",
            "Bayesian_D_max",
            "significance",
            "Bayesian_significance",
            "Bayesian_prob_not_zero_damage",
            "Bayesian_prob_gt_1p_damage",
            "Bayesian_D_max_CI_low",
            "Bayesian_D_max_CI_high",
        ]

        series = df_metaDMG_results.query(f"tax_id == '{tax_id}'")
        if len(series) > 1:
            raise AssertionError("series should be of length 1")

        elif len(series) == 1:
            series = series.iloc[0]
            d[f"|D|"] = series["N_reads"]
            _add_counts_to_dict_bang(
                d,
                df_metaDMG_mismatch.query(f"tax_id == '{tax_id}' & abs(position) == 1"),
                "D",
            )

            for column in columns_to_keep:
                d[column] = series[column]

        # if tax id not in metaDMG results
        else:
            d[f"|D|"] = 0
            _add_counts_to_dict_bang_empty(d, "D")

            for column in columns_to_keep:
                d[column] = np.nan

        out.append(d)

    df_comparison = pd.DataFrame(out).sort_values("|A|", ascending=False)
    df_comparison

    df_comparison.loc[:, "|B|/|A|"] = df_comparison["|B|"] / df_comparison["|A|"]
    df_comparison.loc[:, "|C|/|B|"] = df_comparison["|C|"] / df_comparison["|B|"]
    df_comparison.loc[:, "|C|/|A|"] = df_comparison["|C|"] / df_comparison["|A|"]
    df_comparison.loc[:, "|C\A|/|C|"] = df_comparison["|C\A|"] / df_comparison["|C|"]

    return df_comparison.reset_index(drop=True)


#%%


def add_simulation_information_to_df_comparison(
    path_alignment_files,
    df_comparison,
    sample,
    N_reads,
    simulation_method,
    N_reads_col="N_reads",
):

    df = df_comparison.copy()

    df.loc[:, "sample"] = sample
    df.loc[:, N_reads_col] = N_reads
    df.loc[:, "simulation_method"] = simulation_method

    if "-forward" in sample:
        sample = sample.replace("-forward", "")
    prefix = path_alignment_files / sample / "single" / str(N_reads)

    df_communities_read_abundances = pd.read_csv(
        prefix / f"{sample}.communities_read-abundances.tsv",
        sep="\t",
    )
    df_genome_compositions = pd.read_csv(
        prefix / f"{sample}.genome-compositions.tsv",
        sep="\t",
    )

    df_file_paths = pd.read_csv(
        prefix / f"{sample}.filepaths.tsv",
        sep="\t",
    )

    d_tax_id_to_taxon = get_key2val_dict(df_file_paths, "TaxId", "Taxon")

    tmp = []

    for tax_id in df["tax_id"]:
        # df.query(f"tax_id == {tax_id}")

        try:
            taxon = d_tax_id_to_taxon[tax_id]

            d_tmp = {"tax_id": tax_id}

            series1 = df_communities_read_abundances.query(f"taxon == '{taxon}'")
            assert len(series1) == 2

            d_tmp["simulated_seq_depth_ancient"] = series1.query(
                "frag_type == 'ancient'"
            ).iloc[0]["seq_depth"]
            d_tmp["simulated_seq_depth_modern"] = series1.query(
                "frag_type == 'modern'"
            ).iloc[0]["seq_depth"]

            series2 = df_genome_compositions.query(f"Taxon == '{taxon}'")
            assert len(series2) == 1
            series2 = series2.iloc[0]

            d_tmp["simulated_only_ancient"] = series2.onlyAncient
            d_tmp["simulated_D_max"] = series2.D_max

            tmp.append(d_tmp)

        except KeyError:
            pass

    df_tmp = pd.DataFrame(tmp)

    return pd.merge(df, df_tmp, on="tax_id")


#%%


def get_df_comparison_path(path_comparison, sample, N_reads, simulation_method):

    filename = (
        path_comparison / f"{sample}.{N_reads}.{simulation_method}.comparison.csv"
    )

    return filename


def load_df_comparison(
    df_simulation_alignment,
    df_metaDMG_mapped,
    df_metaDMG_mismatch,
    df_metaDMG_results,
    d_reference_sequences,
    path_comparison,
    path_alignment_files,
    sample,
    N_reads,
    simulation_method,
):

    filename = get_df_comparison_path(
        path_comparison,
        sample,
        N_reads,
        simulation_method,
    )

    try:
        df_comparison = pd.read_csv(filename)
        return df_comparison

    except FileNotFoundError:
        pass

    df_comparison = compute_comparison(
        df_simulation_alignment,
        df_metaDMG_mapped,
        df_metaDMG_mismatch,
        df_metaDMG_results,
        d_reference_sequences,
    )

    df_comparison = add_simulation_information_to_df_comparison(
        path_alignment_files,
        df_comparison,
        sample,
        N_reads,
        simulation_method,
    )

    try:
        filename.parent.mkdir(exist_ok=True)
        # print(f"Saving: {filename}")
        df_comparison.to_csv(filename, index=False)
        return df_comparison

    except OSError:
        print(f"Error saving {filename}")


#%%


def main(p):

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
    ) = p

    # pbar.set_description(
    #     f"Processing sample {sample}, {N_reads} reads, {simulation_method} method"
    # )

    filename_comparison = get_df_comparison_path(
        path_comparison,
        sample,
        N_reads,
        simulation_method,
    )
    if filename_comparison.exists():
        return None

    print(f"Processing sample {sample}, {N_reads} reads, {simulation_method} method")
    path_simulation_alignment_all = get_simulation_alignment_paths(
        path_alignment_files=path_alignment_files,
        name=sample,
        N_reads=str(N_reads),
    )

    df_simulation_alignment_all = load_simulation_alignment_all(
        path_simulation_alignment_all,
        path_alignment_parquet,
        simulation_methods,
        sample,
        N_reads,
    )

    try:
        df_metaDMG_mismatch_all = load_df_metaDMG_mismatch_all(
            sample,
            N_reads,
            simulation_methods,
        )
    except FileNotFoundError:
        return None

    df_metaDMG_results_all = load_df_metaDMG_results_all(
        sample,
        N_reads,
        simulation_methods,
    )

    df_metaDMG_mapped = load_metaDMG_lca_file(
        sample,
        simulation_method,
        N_reads,
        path_analysis_lca,
    )

    d_reference_sequences = load_reference_genomes(
        path_genome_fasta,
        df_simulation_alignment_all,
        sample,
        N_reads,
    )

    df_comparison = load_df_comparison(
        df_simulation_alignment=df_simulation_alignment_all[simulation_method],
        df_metaDMG_mapped=df_metaDMG_mapped,
        df_metaDMG_mismatch=df_metaDMG_mismatch_all[simulation_method],
        df_metaDMG_results=df_metaDMG_results_all[simulation_method],
        d_reference_sequences=d_reference_sequences,
        path_comparison=path_comparison,
        path_alignment_files=path_alignment_files,
        sample=sample,
        N_reads=N_reads,
        simulation_method=simulation_method,
    )



#%%


def list_subdirs(directory):
    return [path for path in directory.glob("*/") if path.is_dir()]


def get_simulation_details(directory):

    out = []
    for subdir_sample in tqdm(list_subdirs(directory)):
        # break

        sample = subdir_sample.name

        for subdir_N_reads in list_subdirs(subdir_sample / "single"):
            # break

            N_reads = int(subdir_N_reads.name)

            df_communities_read_abundances = pd.read_csv(
                subdir_N_reads / f"{sample}.communities_read-abundances.tsv",
                sep="\t",
            )  # .rename(columns={"taxon": "Taxon"})

            df_genome_compositions = pd.read_csv(
                subdir_N_reads / f"{sample}.genome-compositions.tsv",
                sep="\t",
            )

            df_file_paths = pd.read_csv(
                subdir_N_reads / f"{sample}.filepaths.tsv",
                sep="\t",
            )

            for row in df_file_paths.itertuples():

                taxon = row.Taxon
                tax_id = row.TaxId

                series1 = df_genome_compositions[df_genome_compositions.Taxon == taxon]
                assert len(series1) == 1
                series1 = series1.iloc[0]

                series2 = df_communities_read_abundances[
                    df_communities_read_abundances.taxon == taxon
                ]
                assert len(series2) == 2

                d_tmp = {
                    "sample": sample,
                    "tax_id": tax_id,
                    "taxon": taxon,
                    "simulated_N_reads": N_reads,
                    "simulated_only_ancient": series1.onlyAncient,
                    "simulated_D_max": series1.D_max,
                    "simulated_seq_depth_ancient": series2.query(
                        "frag_type == 'ancient'"
                    ).iloc[0]["seq_depth"],
                    "simulated_seq_depth_modern": series2.query(
                        "frag_type == 'modern'"
                    ).iloc[0]["seq_depth"],
                }
                out.append(d_tmp)

    return pd.DataFrame(out)
