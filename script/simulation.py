import numpy as np
import pandas as pd
from numpy.linalg import norm
from scipy.stats import pearsonr, spearmanr
from sklearn.metrics.pairwise import euclidean_distances
from tqdm import tqdm

from common import normalize

INIT_VEC_FILE = snakemake.input.init_vec
SIM_FILE = snakemake.input.sim
WEIGHT_FILE = snakemake.input.weight

lr = snakemake.params.lr
repetition = snakemake.params.repetition

OUTPUT_FILE = snakemake.output.simul_result


init_vec = np.load(INIT_VEC_FILE)
sim = np.load(SIM_FILE)
emperical_weight = np.load(WEIGHT_FILE)


def get_temp_weights(weight_mtx):
    temp_weights = np.zeros_like(weight_mtx)
    for i, row in enumerate(weight_mtx):
        temp_weights[i][np.random.choice(len(weight_mtx), p=row)] = 1
    return temp_weights


def simulate(init_vec, sim, weight_mtx, max_iter=5000):
    simul_sims = []
    norms = []
    vec = init_vec.copy()
    for i in range(max_iter):
        # update vectors sequentially
        temp_weight = get_temp_weights(weight_mtx)
        vec += lr * (temp_weight @ vec - vec)
        vec = normalize(vec)

        # calculate sim
        simul_sim = euclidean_distances(vec)
        simul_sim = (np.sqrt(2) - simul_sim) / np.sqrt(2)
        for i in range(len(weight_mtx)):
            simul_sim[i][i] = 0
        simul_sims.append(simul_sim)
        norms.append(norm(sim - simul_sim))
    best_result = simul_sims[np.argmin(norms)]

    return best_result, np.argmin(norms)


results = []

for i in tqdm(range(repetition)):
    sim_simul, arg_min_index = simulate(init_vec, sim, emperical_weight)
    flatten_sim_simul = np.extract(1 - np.eye(len(sim_simul)), sim_simul)
    flatten_sim = np.extract(1 - np.eye(len(sim)), sim)
    p_corr = pearsonr(flatten_sim_simul, flatten_sim)
    r_corr = spearmanr(flatten_sim_simul, flatten_sim)

    results.append([sim_simul, arg_min_index, p_corr, r_corr])

np.save(OUTPUT_FILE, results)
