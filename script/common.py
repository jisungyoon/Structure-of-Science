import numpy as np


def normalize(mtx):
    row_sums = mtx.sum(axis=1)
    new_mtx = mtx / row_sums[:, np.newaxis]
    return new_mtx
