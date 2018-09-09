# utilities
import functools
# multiprocessing tools
import multiprocessing as mp
import os
import pickle
import random
from time import time

import numpy as np
# data analysis tools
import pandas as pd
from joblib import Parallel, delayed

# network analysis tools

os.chdir("{}/Science/eQTL_analysis/".format(os.environ["HOME"]))

# visualization tools

import qtls

random.seed(int(time()))


''' Where possible, gene names were converted from systematic to standard notation '''
# mRNA expression and genotypes of strains the data is available for
expression_df = pd.read_table("./data/eQTLs/2011/expression.csv")

# qtls_df = pd.read_table("./data/eQTLs/qtls_2018.csv")
qtls_df = pd.read_table("./data/eQTLs/2011/qtls_albert&bloom_script_output.csv")  # .query("q_value <= 0.05")
qtls_type = "eQTLs_2011_A&B"

interactions_type = "all"
interactome_graphs_dict = {}
with open("./data/interactions/{}_interactions_graph.pkl".format(interactions_type), "rb") as pickle_file:
    interactome_graph = pickle.load(pickle_file)

modules_type = "thebiogrid"
with open("./results/{}/modules_dict.pkl".format(modules_type), "rb") as pickle_file:
    modules_dict = pickle.load(pickle_file)
modules_dict.pop("all", None)
modules_dict = {key: value for key, value in modules_dict.items() if key == "physical"}

print(modules_dict.keys())

gene_pool = expression_df["gene"].values
t0 = time()
with Parallel(n_jobs=min(len(modules_dict.keys()), mp.cpu_count())) as parallel:
    full_graph_test_handler = functools.partial(
        qtls.full_graph_test,
        gene_pool=gene_pool, RANDITER_COUNT=200, qtl_df=qtls_df
    )
    results_df = pd.DataFrame(
        np.vstack(parallel(
            delayed(full_graph_test_handler)(module_genes)
            for module_genes in modules_dict.values())
        ),
        columns=["real", "random", "significant", "lo_ci", "hi_ci"]
    )
    results_df.insert(0, "module_name", modules_dict.keys())
    results_df["significant"] = results_df["significant"].astype("bool")
    results_df.to_csv("./formatted_results/{}/{}/full_graph_test.csv".format(modules_type, qtls_type),
                      sep='\t', index=False)

    ppin_test_handler = functools.partial(
        qtls.ppin_test,
        interactions_type="physical", interactome_graph=interactome_graph,
        RANDITER_COUNT=1024, qtl_df=qtls_df
    )
    results_df = pd.DataFrame(
        np.vstack(parallel(
            delayed(ppin_test_handler)(module_genes)
            for module_genes in modules_dict.values())
        ),
        columns=["real", "random", "significant"]
    )
    results_df.insert(0, "module_name", modules_dict.keys())
    results_df["significant"] = results_df["significant"].astype("bool")
    results_df.to_csv("./formatted_results/{}/{}/ppin_test.csv".format(modules_type, qtls_type),
                      sep='\t', index=False)
    print("{} elapsed".format(time() - t0))
