from os.path import join as j
import numpy as np


configfile: "config.yaml"


DATA_DIR = config["data_dir"]
SIMUL_RESULT_DIR = config["simul_result_dir"]


###############################################################################
# RAW_DATA
###############################################################################

SIM_FILE = j(DATA_DIR, "raw_similarity.pkl")
LANG_LIST_FILE = j(DATA_DIR, "lang_list_for_simulation.npy")
EXTERNAL_DATA_FILES = j(
    DATA_DIR, "external_data/country_to_lang_figure/data/{measure}.pkl"
)


###############################################################################
# SIMUL_RESULTS
###############################################################################

INIT_VECTORS_DIR = j(SIMUL_RESULT_DIR, "init_vectors")
INIT_VECTORS_FILES = j(INIT_VECTORS_DIR, "{number}.npy")

WEIGHT_MATRIX_DIR = j(SIMUL_RESULT_DIR, "weight_mtx")
WEIGHT_MATRIX_FILES = j(WEIGHT_MATRIX_DIR, "{measure}.npy")
PROCESSED_SIM_FILE = j(WEIGHT_MATRIX_DIR, "similarity.npy")

SIMPLE_SIMULATION_RESULT_DIR = j(SIMUL_RESULT_DIR, "result", "simple")
SIMPLE_SIMULATION_RESULT_FILE = j(
    SIMPLE_SIMULATION_RESULT_DIR, "{number},{measure}.npy"
)
SIMPLE_SIMULATION_RESULT_FILE_RAND = j(
    SIMPLE_SIMULATION_RESULT_DIR, "{number},random.npy"
)


###############################################################################
# HYPER PARAMETERS
###############################################################################
ratio = 10
n_lang = 52
length = 52 * ratio
lr = 0.001

external_vals = [
    "export",
    "good_un_data",
    "no_citation_pair",
    "no_collab_pair",
    "foriegn_stuedent_all",
    "facebook",
    "distance",
]
number_of_sample = 100
repetition = 10

rule all:
    input:
        expand(INIT_VECTORS_FILES, number=range(number_of_sample)),
        expand(WEIGHT_MATRIX_FILES, measure=external_vals),
        PROCESSED_SIM_FILE,
        expand(SIMPLE_SIMULATION_RESULT_FILE, number=range(number_of_sample), measure=external_vals),
        expand(SIMPLE_SIMULATION_RESULT_FILE_RAND, number=range(number_of_sample)),


rule make_init_vectors:
    params:
        n_lang=n_lang,
        length=length,
    output:
        init_vecotr_file=INIT_VECTORS_FILES,
    script:
        "script/generate_init_vectors.py"


rule make_weight_matrix:
    input:
        external_val_file=EXTERNAL_DATA_FILES,
    output:
        weight_mtx_file=WEIGHT_MATRIX_FILES,
    script:
        "script/generate_weight_mtx.py"


rule preprocess_sim_file:
    input:
        sim=SIM_FILE,
        lang_list=LANG_LIST_FILE,
    output:
        processed_sim=PROCESSED_SIM_FILE,
    script:
        "script/preprocess_sim_file.py"


rule simulation:
    input:
        init_vec=INIT_VECTORS_FILES,
        sim=PROCESSED_SIM_FILE,
        weight=WEIGHT_MATRIX_FILES,
    params:
        lr=lr,
        repetition=repetition,
    output:
        simul_result=SIMPLE_SIMULATION_RESULT_FILE,
    script:
        "script/simulation.py"



rule simulation_random:
    input:
        init_vec=INIT_VECTORS_FILES,
        sim=PROCESSED_SIM_FILE,
    params:
        lr=lr,
        repetition=repetition,
    output:
        simul_result=SIMPLE_SIMULATION_RESULT_FILE_RAND,
    script:
        "script/simulation_random.py"
