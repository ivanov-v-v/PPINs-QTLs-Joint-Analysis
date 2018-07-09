''' General-purpose graph routines for biological networks analysis.
    At the moment supports:

    - Conversion between different naming nomenclatures.
    - Graph construction from edges list.

'''


import igraph as ig
import numpy as np


# Returns a list of genes connected to gene_name in interactions_graph
# where "connected" is interpreted as "reachable in bfs_depth steps".
def convert_to_std(gene_name, sys2std):
    try:
        result = sys2std[gene_name]
    except KeyError:
        result = gene_name
    return result


def convert_to_sys(gene_name, std2sys):
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
