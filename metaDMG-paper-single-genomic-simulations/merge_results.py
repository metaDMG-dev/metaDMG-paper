#%%

from importlib import reload
from pathlib import Path

import numpy as np
import pandas as pd
from tqdm import tqdm

import utils

#%%

reload(utils)

species = "homo"
df = utils.load_results(species, use_columns_subset=False)


#%%

reload(utils)
data_dir = utils.get_data_dir(species)

# reload(utils)
for p, dfg in tqdm(df.groupby(utils.sim_columns[:-1])):

    sim_species, sim_damage, sim_N_reads, sim_length = p

    name = f"{sim_species}-{sim_damage}-{sim_N_reads}-{sim_length}"
    path = data_dir / f"results-{sim_species}" / f"{name}.parquet"

    dfg["tax_id"] = dfg["sample"]
    dfg["sample"] = name

    path.parent.mkdir(exist_ok=True, parents=True)
    dfg.to_parquet(path)

# %%
