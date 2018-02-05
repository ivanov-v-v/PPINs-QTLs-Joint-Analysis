import igraph as ig


# Given a list of gene-gene interactions of some kind (pairs of genes treated as edges),
# construct a graph and, if specified, randomize it, preserving some invariant
def assemble_graph_of_interactions(edges, directed=False, randomize=False):
    interaction_graph = ig.Graph(directed=directed)
    # For some weird reason, igraph can add multiple copies of the same
    # vertex without even signalling about it, therefore the duplicates
    # require manual removal from the dataset
    vertices = set()
    for source, target in edges:
        vertices |= {source, target}
    vertex_names = list(vertices)
    interaction_graph.add_vertices(vertex_names)
    interaction_graph.add_edges(edges)
    ''' TODO:   There must be hundreds to thousands of 
                randomization iterations over which the
                results will be averaged and used,
                this must be implemented efficiently
    '''
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
    # Due to domain-related specifics,
    # directed graphs produced by this function
    # will represent marker-gene interactions
    # and must be bipartite therefore
    if directed:
        interaction_graph.vs["type"] = \
            [False if deg > 0 else True for deg in interaction_graph.outdegree()]
    interaction_graph.vs["name"] = vertex_names
    return interaction_graph


''' TODO:   Стоит попробовать более совершенные метрики
            подобия, а также добиться лучшей скорости работы 
'''

# Given a pair of graphs, representing gene-gene interactions
# and estimated QTL-linkages, calculate for each pair of interacting genes
# a Jaccard similarity coefficient, and then average it over all edges
def mean_linkage_similarity(interaction_graph, QTL_graph):
    linked_genes = set([vertex["name"] for vertex in QTL_graph.vs])
    interacting_genes = [vertex["name"] for vertex in interaction_graph.vs]

    mean_coeff = 0.
    # Перебрать все рёбра и сопоставить каждому
    # пару множеств: eQTL, которые линкуются с инцидентными вершинами,
    # а затем рассмотреть меру пересечения их объединения с мерой пересечения
    if interaction_graph.ecount():
        for edge in interaction_graph.es:
            s_id, t_id = edge.source, edge.target
            s_name = interacting_genes[s_id]
            t_name = interacting_genes[t_id]
            if s_name in linked_genes and t_name in linked_genes:
                s_neigh = set(QTL_graph.neighbors(s_name, mode="IN"))
                t_neigh = set(QTL_graph.neighbors(t_name, mode="IN"))
                mean_coeff += len(s_neigh & t_neigh) / len(s_neigh | t_neigh)
        mean_coeff /= interaction_graph.ecount()
    return mean_coeff