# utilities
import functools
# multiprocessing tools
import multiprocessing as mp
import os
import pickle
import random
import time

# data analysis tools
import joblib
import numpy as np
import pandas as pd

# network analysis tools

os.chdir("{}/Science/eQTL_analysis/".format(os.environ["HOME"]))

# visualization tools

import qtls

random.seed(int(time.time()))

# 112 segregants genotyped by inherited marker variants
full_genotypes_df = pd.read_table("./data/genotypes/processed_genotypes.csv")

''' Where possible, gene names were converted from systematic to standard notation '''
# mRNA expression and genotypes of strains the data is available for
expression_df = pd.read_table("./data/eQTLs/expression.csv")
genotypes_df = pd.read_table("./data/eQTLs/genotypes_2018.csv")

qtls_df = pd.read_table("./data/eQTLs/qtls_2018.csv")
# qtls_df = pd.read_table("./data/eQTLs/qtls.csv").query("q_value <= 0.05")

interactions_type = "all"
interactome_graphs_dict = {}
with open("./data/interactions/{}_interactions_graph.pkl".format(interactions_type), "rb") as pickle_file:
    interactome_graph = pickle.load(pickle_file)

modules_type = "physical"
with open("./results/{}/modules_dict.pkl".format(modules_type), "rb") as pickle_file:
    modules_dict = pickle.load(pickle_file)
modules_dict.pop("all", None)

gene_pool = expression_df["gene"].values

results_df = pd.DataFrame(
    [qtls.community_graph_test(modules_dict=modules_dict, gene_pool=gene_pool, RANDITER_COUNT=200, qtl_df=qtls_df)],
    columns=["real", "random", "significant", "lo_ci", "hi_ci"]
)
results_df.insert(0, "module_name", [modules_type])
results_df["significant"] = results_df["significant"].astype("bool")
results_df.to_csv("./results/{}/community_graph_test_new_eQTLs.csv".format(modules_type),
                  sep='\t', index=False)

full_graph_test_handler = functools.partial(
    qtls.full_graph_test,
    gene_pool=gene_pool, RANDITER_COUNT=200, qtl_df=qtls_df
)
results_df = pd.DataFrame(
    np.vstack(joblib.Parallel(n_jobs=mp.cpu_count())(
        joblib.delayed(full_graph_test_handler)(module_genes)
        for module_genes in modules_dict.values())
    ),
    columns=["real", "random", "significant", "lo_ci", "hi_ci"]
)
results_df.insert(0, "module_name", modules_dict.keys())
results_df["significant"] = results_df["significant"].astype("bool")
results_df.to_csv("./results/{}/full_graph_test_new_eQTLs.csv".format(modules_type),
                  sep='\t', index=False)

ppin_test_handler = functools.partial(
    qtls.ppin_test,
    interactions_type="physical", interactome_graph=interactome_graph,
    RANDITER_COUNT=1024, qtl_df=qtls_df
)
results_df = pd.DataFrame(
    np.vstack(joblib.Parallel(n_jobs=mp.cpu_count())(
        joblib.delayed(ppin_test_handler)(module_genes)
        for module_genes in modules_dict.values())
    ),
    columns=["real", "random", "significant"]
)
results_df.insert(0, "module_name", modules_dict.keys())
results_df["significant"] = results_df["significant"].astype("bool")
results_df.to_csv("./results/{}/ppin_test_new_eQTLs.csv".format(modules_type),
                  sep='\t', index=False)
