import numpy as np
from numpy.linalg import norm
from sklearn.metrics.pairwise import euclidean_distances


def normalize(mtx):
    row_sums = mtx.sum(axis=1)
    new_mtx = mtx / row_sums[:, np.newaxis]
    return new_mtx


def get_temp_weights(weight_mtx):
    temp_weights = np.zeros_like(weight_mtx)
    for i, row in enumerate(weight_mtx):
        temp_weights[i][np.random.choice(len(weight_mtx), p=row)] = 1
    return temp_weights


def simulate(init_vec, sim, weight_mtx, lr, max_iter=5000):
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
