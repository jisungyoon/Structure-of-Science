import numpy as np
import pandas as pd

from common import normalize

INPUT_FILE = snakemake.input.external_val_file
OUTPUT_FILE = snakemake.output.weight_mtx_file
measure_type = snakemake.wildcards.measure


input_data = pd.read_pickle(INPUT_FILE)

if measure_type == "distance":
    weight_mtx = normalize(1 / np.log10(input_data))

else:
    weights = input_data[1]
    # inpute with one if assigned value is lower than 1
    weights[weights < 1] = 1
    for i in range(len(weights)):
        weights[i][i] = 1
    weight_mtx = normalize(np.log10(weights))

np.save(OUTPUT_FILE, weight_mtx)
