#%%

import numpy as np
import pysam

#%%

contig_sizes = [1_000, 10_000, 100_000]
linelength = 80


#%%

np.random.seed(42)

for contig_size in contig_sizes:

    contig = "".join(
        [np.random.choice(["A", "C", "G", "T"]) for _ in range(contig_size)]
    )

    filename = f"genome/contig{contig_size}.fna"

    with open(filename, "w") as ofile:
        ofile.write(f">contig{contig_size}\n")
        for i in range(0, contig_size, linelength):
            ofile.write(contig[i : i + linelength])
            ofile.write("\n")

# %%
