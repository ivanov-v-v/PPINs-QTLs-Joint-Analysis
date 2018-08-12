import multiprocessing as mp
import os
import pickle
import random
import time

import pandas as pd
from joblib import Parallel, delayed

import qtls
import util

os.chdir("{}/Science/eQTL_analysis/".format(os.environ["HOME"]))

random.seed(time.time())

# 112 segregants genotyped by inherited marker variants

''' Where possible, gene names were converted from systematic to standard notation '''

expression_df, genotypes_df, qtls_df = {}, {}, {}

# expression_df["eQTLs_old"] = pd.read_table("./data/eQTLs/expression.csv")
# expression_df["eQTLs_new"] = pd.read_table("./data/eQTLs/expression_2018.csv")
expression_df["pQTLs"] = pd.read_table("./data/pQTLs/expression.csv")

# genotypes_df["eQTLs_old"] = pd.read_table("./data/eQTLs/genotypes.csv")
# genotypes_df["eQTLs_new"] = pd.read_table("./data/eQTLs/genotypes_2018.csv")
# genotypes_df["eQTLs_interpolated"] = pd.read_table("./data/eQTLs/genotypes_interpolated.csv")
genotypes_df["pQTLs"] = pd.read_table("./data/pQTLs/genotypes.csv")

# QTLs estimated with MatrixEQTL package for R
qtls_df["eQTLs_old"] = pd.read_table("./data/eQTLs/qtls.csv").query("q_value <= 0.05")
# qtls_df["eQTLs_new"] = pd.read_table("./data/eQTLs/qtls_2018.csv")#.query("q_value <= 0.05")
# qtls_df["eQTLs_interpolated"] = pd.read_table("./data/eQTLs/qtls_interpolated.csv")
qtls_df["pQTLs"] = pd.read_table("./data/pQTLs/qtls.csv").query("q_value <= 0.05")

interactions_type = "physical"
with open("./data/interactions/{}_interactions_graph.pkl".format(interactions_type), "rb") as pickle_file:
    interactome_graph = pickle.load(pickle_file)

with Parallel(n_jobs=mp.cpu_count()) as parallel:
    for modules_type in ["thecellmap"]:
        # modules_type = "thebiogrid"
        with open("./results/" + modules_type + "/modules_dict.pkl", "rb") as pickle_file:
            modules_dict = pickle.load(pickle_file)
        # modules_to_test.pop("all", None)
        util.ensure_dir("./data/pQTLs/temp")
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

        t0 = time.time()


        def predict_pQTLs(module_name, module_genes):
            return qtls.PqtlPredictor(
                eqtls_df=qtls_df["eQTLs_old"], pqtls_df=qtls_df["pQTLs"],
                pqtls_expression_df=expression_df["pQTLs"], pqtls_genotypes_df=genotypes_df["pQTLs"],
                module_name=module_name, module_genes=module_genes,
                interactome_graph=interactome_graph
            ).predict()


        modules_for_pQTL_prediction = [module_name
                                       for module_name, module_genes in modules_dict.items()
                                       if len(qtls.linked_markers(qtls_df["pQTLs"], module_genes)) > 0]

        pQTL_prediction_results = \
            parallel(
                delayed(predict_pQTLs)(module_name, modules_dict[module_name])
                for module_name in modules_for_pQTL_prediction
            )

        pd.DataFrame(
            ((key, *values) for key, values in zip(modules_for_pQTL_prediction, pQTL_prediction_results)),
            columns=["module_name", "common_pQTLs", "overlap_ratio", "old_pQTLs", "predicted_pQTLs", "delta"]
        ).sort_values(by="overlap_ratio", ascending=False) \
            .to_csv("./results/{}/pQTLs_from_eQTLs_old.csv".format(modules_type), sep='\t', index=False)
        print("{}: {} elapsed".format(modules_type, time.time() - t0))
