import numpy as np

from common import normalize

n_lang = snakemake.params.n_lang
length = snakemake.params.length

OUTPUT_DIR = snakemake.output.init_vecotr_file


def initialize_vector(col, row):
    vec = np.zeros((row, col))
    ratio = int(col / row)
    for rep in range(ratio):
        for i in range(row):
            vec[i][rep * row + i] = 1
    vec = vec[:, np.random.permutation(vec.shape[1])]
    return vec


init_vector = normalize(initialize_vector(length, n_lang))
np.save(OUTPUT_DIR, init_vector)
