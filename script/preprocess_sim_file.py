import numpy as np
import pandas as pd

from common import normalize

SIM_FILE = snakemake.input.sim
LANG_LIST_FILE = snakemake.input.lang_list

OUTPUT_FILE = snakemake.output.processed_sim


simil = pd.read_pickle(SIM_FILE)
lang_list = np.load(LANG_LIST_FILE)

simil = (np.sqrt(2) - (1 - simil)) / np.sqrt(2)
for x in simil.columns:
    simil[x][x] = 0
sim_mtx = simil.reindex(columns=lang_list, index=lang_list).to_numpy()

np.save(OUTPUT_FILE, sim_mtx)
