import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr
from tqdm import tqdm

from common import normalize, get_temp_weights, simulate


INIT_VEC_FILE = snakemake.input.init_vec
SIM_FILE = snakemake.input.sim

lr = snakemake.params.lr
repetition = snakemake.params.repetition

OUTPUT_FILE = snakemake.output.simul_result


init_vec = np.load(INIT_VEC_FILE)
sim = np.load(SIM_FILE)
n_lang = sim.shape[0]
rand_weight = rand_weights = normalize(np.random.random((n_lang, n_lang)))



results = []

for i in tqdm(range(repetition)):
    sim_simul, arg_min_index = simulate(init_vec, sim, rand_weight, lr)
    flatten_sim_simul = np.extract(1 - np.eye(len(sim_simul)), sim_simul)
    flatten_sim = np.extract(1 - np.eye(len(sim)), sim)
    p_corr = pearsonr(flatten_sim_simul, flatten_sim)
    r_corr = spearmanr(flatten_sim_simul, flatten_sim)

    results.append([sim_simul, arg_min_index, p_corr, r_corr])

np.save(OUTPUT_FILE, results)
