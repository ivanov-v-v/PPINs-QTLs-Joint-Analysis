import os
import pickle
import random
import time

import numpy as np
import pandas as pd

import qtls

os.chdir("{}/Science/eQTL_analysis/".format(os.environ["HOME"]))

random.seed(time.time())

full_genotypes_df = pd.read_table("./data/genotypes/processed_genotypes.csv")
eQTLs_expression_df = pd.read_table("./data/eQTLs/expression.csv")
eQTLs_genotypes_df = pd.read_table("./data/eQTLs/genotypes_2018.csv")
# Protein expression and genotypes of strains the data is available for
pQTLs_expression_df = pd.read_table("./data/pQTLs/expression.csv")
pQTLs_genotypes_df = pd.read_table("./data/pQTLs/genotypes_2018.csv")

# QTLs estimated with MatrixEQTL package for R
eQTLs_df = pd.read_table("./data/eQTLs/qtls.csv").query("q_value < 0.05")
pQTLs_df = pd.read_table("./data/pQTLs/qtls.csv").query("q_value < 0.05")

interactions_type = "all"
with open("./data/interactions/{}_interactions_graph.pkl".format(interactions_type), "rb") as pickle_file:
    interactome_graph = pickle.load(pickle_file)

modules_type = "thebiogrid"
with open("./results/" + modules_type + "/modules_dict.pkl", "rb") as pickle_file:
    modules_dict = pickle.load(pickle_file)
# modules_to_test.pop("all", None)

# for qtl_type, qtl_df, expression_df in [("eQTLs", eQTLs_df, eQTLs_expression_df),
#                                         ("pQTLs", pQTLs_df, pQTLs_expression_df)]:
#     module_analyzer = qtls.LinkageSharingAnalyzer(
#         expression_df=expression_df,
#         interactions_type=interactions_type,
#         interactome_graph=interactome_graph,
#         modules_type=modules_type,
#         modules_dict=modules_dict,
#         qtl_type=qtl_type,
#         qtl_df=qtl_df,
#         q_thresholds=np.logspace(-5, -2, 10),
#         pairwise_test_iter=200,
#         ppi_test_iter=1024
#     )
#     module_analyzer.pairwise_test()
#     module_analyzer.ppi_test()
#

predictions_dict = {}
for module_name, module_genes in modules_dict.items():
    module_graph = interactome_graph.subgraph(
        np.intersect1d(interactome_graph.vs["name"], module_genes)
    )
    if len(qtls.linked_markers(pQTLs_df, module_genes)) == 0 \
            or len(qtls.linked_markers(eQTLs_df, module_genes)) == 0:
        predictions_dict[module_name] = "No linkages"
    else:
        predictions_dict[module_name] = qtls.PqtlPredictor(
            eQTLs_df[eQTLs_df["gene"].isin(module_graph.vs["name"])],
            pQTLs_df[pQTLs_df["gene"].isin(module_graph.vs["name"])],
            eQTLs_expression_df, eQTLs_genotypes_df,
            pQTLs_expression_df, pQTLs_genotypes_df,
            full_genotypes_df,
            module_name, module_graph
        ).predict()

pd.DataFrame(
    ((key, *values) for key, values in predictions_dict.items() if values != "No linkages"),
    columns=["module_name", "common_pQTLs", "overlap_ratio", "old_pQTLs", "predicted_pQTLs", "delta"]
).sort_values(by="overlap_ratio", ascending=False) \
    .to_csv("./results/{}/pQTLs_from_eQTLs.csv".format(modules_type), sep='\t', index=False)
