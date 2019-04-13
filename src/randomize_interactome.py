import multiprocessing as mp
import os
import pickle

import igraph as ig
from joblib import Parallel, delayed

os.chdir("..")

RANDITER = 1024

interactions_type = "all"
with open("./data/raw/interactions/{}_interactions_graph.pkl".format(interactions_type), "rb") as pickle_file:
    interactome_graph = pickle.load(pickle_file)

def _produce_randomized_interactome_copy(i):
    randomized_interactome = interactome_graph.Degree_Sequence(
        interactome_graph.degree(),
        method='vl'
    )
    randomized_interactome.vs["name"] = interactome_graph.vs["name"]
    ig.write(
        graph=randomized_interactome,
        filename="./data/raw/randomized_interactome_copies/{0}/{1}.pkl".format(interactions_type, i),
        format="pickle"
    )


Parallel(n_jobs=mp.cpu_count())(
    delayed(_produce_randomized_interactome_copy)(i)
    for i in range(RANDITER)
)
