"""
General-purpose graph routines for biological networks analysis.
At the moment supports:
- Conversion between different naming nomenclatures.
- Graph construction from edges list.
- Display of basic information about graphs.
"""

import pickle

import igraph as ig
import numpy as np
import pandas as pd


# Returns a list of genes connected to gene_name in interactions_graph
# where "connected" is interpreted as "reachable in bfs_depth steps".
def convert_to_std(gene_name):
    with open("./data/nomenclature/sys2std.pkl", "rb") as pickle_file:
        sys2std = pickle.load(pickle_file)
    try:
        result = sys2std[gene_name]
    except KeyError:
        result = gene_name
    return result


def convert_to_sys(gene_name):
    with open("./data/nomenclature/std2sys.pkl", "rb") as pickle_file:
        std2sys = pickle.load(pickle_file)
    try:
        result = std2sys[gene_name]
    except KeyError:
        result = gene_name
    return result


# Given a list of gene-gene interactions of some kind (pairs of genes treated as edges),
# construct a graph and, if specified, randomize it, preserving some invariant
def graph_from_edges(edges, directed=False, randomize=False):
    interaction_graph = ig.Graph(directed=directed)
    vertex_names = np.unique(np.ravel(edges))
    interaction_graph.add_vertices(vertex_names)
    interaction_graph.add_edges(edges)

    if randomize:
        if directed:
            interaction_graph = interaction_graph.Degree_Sequence(
                interaction_graph.outdegree(),
                interaction_graph.indegree(),
                method='vl'
            )
        else:
            interaction_graph = interaction_graph.Degree_Sequence(
                interaction_graph.degree(),
                method='vl'
            )

    interaction_graph.vs["name"] = vertex_names

    # Due to domain-related specifics,
    # directed graphs produced by this function
    # will represent marker-gene interactions
    # and must be bipartite therefore
    if directed:
        interaction_graph.vs["part"] = (np.array(interaction_graph.outdegree()) == 0)

    return interaction_graph


def basic_module_info(interactome_graph, modules_dict, modules_type, make_plots=True, simplify=False):
    module_graphs_dict = {}
    for module_name in modules_dict.keys():
        module_graph = interactome_graph.subgraph(
            set(modules_dict[module_name])
            & set(interactome_graph.vs["name"])
        )
        if simplify:
            module_graph = module_graph.simplify()
            module_graph.vs.select(_degree=0).delete()
        if make_plots:
            destdir = '/'.join(["./img/functional_module_graphs",
                                modules_type,
                                "simplified" if simplify else "raw"])
            ig.plot(
                module_graph,
                layout=module_graph.layout_kamada_kawai(),
                target="/".join([destdir, module_name + ".png"]),
                vertex_label=module_graph.vs["name"]
            );
        module_graphs_dict[module_name] = module_graph
    return pd.DataFrame(
        [(mname, mgraph.vcount(), mgraph.ecount()) for mname, mgraph in module_graphs_dict.items()],
        columns=["module_name", "genes count", "interactions count"]
    )
