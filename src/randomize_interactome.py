import multiprocessing as mp
import os
import pickle

import igraph as ig

# os.chdir("/home/vvi/Science/eQTL_analysis/")
os.chdir("/home/ivanov_vv/eQTL_analysis/")

interactions_type = "genetic"
with open("./data/interactions/{}_interactions_graph.pkl".format(interactions_type), "rb") as pickle_file:
    interactome_graph = pickle.load(pickle_file)


def _produce_randomized_interactome_copy(i):
    randomized_interactome = interactome_graph.Degree_Sequence(
        interactome_graph.degree(),
        method='vl'
    )
    randomized_interactome.vs["name"] = interactome_graph.vs["name"]
    ig.write(
        graph=randomized_interactome,
        filename="./data/randomized_interactome_copies/{0}/{1}.pkl".format(interactions_type, i),
        format="pickle"
    )


pool = mp.Pool(mp.cpu_count())
pool.map(_produce_randomized_interactome_copy, range(1024))
pool.close()
pool.join()
