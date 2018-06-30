from __future__ import print_function

import sys

import numpy as np

from .networks import assemble_graph_of_interactions


def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

# Returns a list of genes connected to gene_name in interactions_graph
# where "connected" is interpreted as "reachable in bfs_depth steps".
def convert_to_std(gene_name, sys_to_std):
    try:
        result = sys_to_std[gene_name]
    except KeyError:
        result = gene_name
    return result


def convert_to_sys(gene_name, std_to_sys):
    try:
        result = std_to_sys[gene_name]
    except KeyError:
        result = gene_name
    return result


# Given a dataframe of estimated QTL linkages and
# full strain genotype lookup table, plot the number
# of linkages against marker location in genome
def map_linkages_to_genome_location(QTL_df, full_genotypes_df):
    QTL_graph = assemble_graph_of_interactions(
        edges=QTL_df[["SNP", "gene"]].values,
        directed=True
    )
    # Select only the marker vertices
    left_part = QTL_graph.vs.select(part=False)
    # And map them to their absolute genome location
    marker_to_rownum = \
        dict(
            zip(
                full_genotypes_df.iloc[:, 0],
                np.arange(full_genotypes_df.shape[0])
            )
        )
    # In order to combine plots for different QTLs
    # it's wise to plot them against the total set of markers,
    # disregarding whether linkages for the given marker
    # are calculated or not
    QTL_marker_to_linkages = \
        dict(
            zip(
                full_genotypes_df.iloc[:, 0],
                np.zeros(full_genotypes_df.shape[0])
            )
        )
    # Add the available linkage data
    QTL_marker_to_linkages.update(dict(
        zip(
            left_part["name"],
            left_part.outdegree()
        )
    ))
    # Sort the dictionary accordingly to the marker position on the chromosome
    # and unzip it to extract the ordered linkage list

    QTL_x, QTL_y = \
        map(list,
            zip(
                *sorted(QTL_marker_to_linkages.items(),
                        key=lambda p: marker_to_rownum[p[0]])
            )
        )
    return QTL_x, QTL_y