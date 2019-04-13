from __future__ import print_function

import collections
import gc
import itertools
import multiprocessing as mp
import os
import pickle
import subprocess

import igraph as ig
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.stats as sps
from tqdm import *

import networks
import util


########################################################################################################################
#                                                     UTILITIES                                                        #
########################################################################################################################


def linked_markers(qtl_df, gene_list):
    return np.unique(qtl_df[qtl_df["gene"].isin(np.atleast_1d(gene_list))]["SNP"].values)


def linked_genes(qtl_df, marker_list):
    return np.unique(qtl_df[qtl_df["SNP"].isin(np.atleast_1d(marker_list))]["gene"].values)


########################################################################################################################
#                                                  VISUALIZATIONS                                                      #
########################################################################################################################

def linkages2gencoords(qtl_df):
    '''
    Given a dataframe of estimated QTL linkages and
    full strain genotype lookup table, plot the number
    of linkages against marker location in genome

    :param qtl_df: must contain columns "SNP" and "gene"
    :param full_genotypes_df: a list o genes sorted by genetic coordinate in ascending order
    :return: a list of linkages sorted in the same order
    '''
    with open("./data/raw/genomic_features/order_of_markers.pkl", "rb") as pickle_file:
        marker2location = pickle.load(pickle_file)

    qtl_graph = networks.graph_from_edges(
        edges=qtl_df[["SNP", "gene"]].values,
        directed=True
    )
    # Select only the marker vertices
    marker_vertices = qtl_graph.vs.select(part=False)
    # And map them to their absolute genome location
    # Add the available linkage data
    qtl_marker_to_linkages = dict(
        zip(
            marker_vertices["name"],
            marker_vertices.outdegree()
        )
    )
    # Sort the dictionary accordingly to the marker position on the chromosome
    # and unzip it to extract the ordered linkage list
    qtl_x, qtl_y = \
        map(list, zip(*sorted(qtl_marker_to_linkages.items(),
                              key=lambda p: marker2location[p[0]])))
    return np.array(qtl_x), np.array(qtl_y)


def plot_module_graph(destdir, module_graph, qtl_df, filename="module_graph", format="png"):
    '''
    Plot functional module graph in .svg format with vertices colored in accordance to
    the number of linked QTLs. Vertices with no linkages are not shown.
    All other vertices have number of linked QTLs written in their labels.

    :param module_graph: interactome subgraph
    :param qtl_df: s/e
    :param destdir: s/e
    :param module_name: s/e
    :return: None
    '''
    num_linkages = [len(linked_markers(qtl_df, [gene_name])) for gene_name in module_graph.vs["name"]]
    if np.max(num_linkages) != np.min(num_linkages):
        normalized_num_linkages = (num_linkages - np.min(num_linkages)) / (np.max(num_linkages)) - np.min(num_linkages)
        vertex_colors = plt.cm.coolwarm(normalized_num_linkages)
        vertex_colors[:, -1] = 0.4
        module_graph.vs["color"] = [util.rgba2hex(rgba) for rgba in vertex_colors]

    edge_colors = [util.rgba2hex(rgba) for rgba in [[0, 0, 0, 0.1]] * module_graph.ecount()]
    for i, e in enumerate(module_graph.es):
        if num_linkages[e.source] != 0 and num_linkages[e.target] != 0:
            edge_colors[i] = util.rgba2hex([1, 0, 0, 0.5])

    module_graph.es["color"] = edge_colors
    n = module_graph.vcount()
    module_graph.vs["size"] = [5 if num_linkages[i] == 0 else 10 for i in range(n)]
    module_graph.vs["label_size"] = [0 if num_linkages[i] == 0 else 15 for i in range(n)]
    module_graph.vs["label_dist"] = [0 if num_linkages[i] == 0 else 2 for i in range(n)]

    ig.plot(
        module_graph,
        target=destdir + '/' + filename + "." + format,
        layout=module_graph.layout_fruchterman_reingold(maxiter=1000, area=n ** 3, repulserad=n ** 3),
        bbox=(1440, 720),
        vertex_label=[gene_name + ", " + str(num_linkages[i]) if num_linkages[i] != 0 else ''
                      for i, gene_name in enumerate(module_graph.vs["name"])],
        margin=50,
    )
    plt.close("all")
    gc.collect()


def plot_q_hist(destdir, module_graph, module_name, module_type, qtl_df, qtl_type):
    """ ??? """
    fig, ax = plt.subplots(figsize=(20, 10))
    plt.xscale("log")
    counts, bins, patches = plt.hist(
        qtl_df["q_value"],
        bins=10. ** np.arange(-8, -1),
        alpha=0.5,
        color="xkcd:turquoise",
        edgecolor="black",
        linewidth=1.2,
        label="{0} linkages in total\n"
              "mean q-value: {1:.4f}\n"
              "median q-value: {2:.4f}\n"
              "|V| = {3}, |E| = {4}"
            .format(qtl_df.shape[0],
                    qtl_df["q_value"].mean(),
                    qtl_df["q_value"].median(),
                    module_graph.vcount(),
                    module_graph.ecount())
    )
    ax.grid(linestyle='dotted', alpha=0.8)
    ax.tick_params(labelsize=15)

    ax.set_xticks(bins)
    ax.tick_params(axis='x', labelrotation=90)

    ax.set_title("{} linkage q-value distribution for {}; {}".format(
        qtl_type, module_name, module_type
    ), fontsize=20)
    ax.set_xlabel("q-value", fontsize=20)
    ax.set_ylabel("number of linkages", fontsize=20)
    ax.legend(fontsize=15)

    plt.savefig(destdir + "/q_hist.png")
    plt.close('all')
    gc.collect()


def plot_test_results(module_name, module_type, test_results_df, test_type, q_thresholds, qtl_type):
    """ ??? """
    # for the plot to be informative, don't show the range
    # where there are no linkage with such q-values at all
    fig, ax = plt.subplots(figsize=(20, 10))
    ax.set_title("{} {} test for {}; {}".format(
        qtl_type, test_type, module_name, module_type
    ), fontsize=20)
    plt.xscale('log')

    significant_point_ids = list(np.where(test_results_df["significant"]))
    ax.set_xlabel("linkage q-value threshold, lg-scale", fontsize=20)
    ax.set_ylabel("average linkage similarity", fontsize=20)
    ax.plot(q_thresholds, test_results_df["real"], '-b*',
            label="real data", markevery=significant_point_ids)
    ax.plot(q_thresholds, test_results_df["simulated"], '-r*',
            label="synthesized data", markevery=significant_point_ids)

    ax.fill_between(q_thresholds,
                    test_results_df["lo_ci"], test_results_df["hi_ci"],
                    color="#FFD700", alpha=0.4)

    ax.grid(linestyle='dotted')
    ax.legend(loc=2, fontsize=15)
    plt.setp(ax.get_xticklabels(), fontsize=15)
    plt.setp(ax.get_yticklabels(), fontsize=12)

    destdir = '/'.join(["./results", module_type, module_name, qtl_type])
    util.ensure_dir(destdir)
    plt.savefig(destdir + "/" + test_type + "_test.png")
    plt.close('all')
    gc.collect()


def qtl_overlap_hist(overlap_data, modules_type):
    fig = plt.figure(figsize=(30, 10))
    hist_ax = fig.add_subplot(1, 2, 1)

    hist_ax.hist(
        overlap_data,
        bins=util.bincount_scott(overlap_data),
        alpha=0.8,
        label="mean: {0:.4f}\nmedian: {1:.4f}".format(
            np.mean(overlap_data),
            np.median(overlap_data)),
        density=True
    );
    hist_ax.tick_params(axis='x', labelsize=20)
    hist_ax.tick_params(axis='y', labelsize=20)
    hist_ax.set_xlabel("jaccard coefficient", fontsize=20)
    hist_ax.legend(loc=1, fontsize=20)
    hist_ax.set_title("eQTLs/pQTLs linkage overlap for {}, in %".format(modules_type), fontsize=25);

    sdata_ax = fig.add_subplot(1, 2, 2)
    sdata_ax.plot(sorted(overlap_data))
    # sdata_ax.axvline(x=np.median(range(len(overlap_data))), color='r', linestyle='--')
    sdata_ax.axhline(y=np.median(overlap_data), color='r', alpha=0.8, linestyle='--', label="median jaccard")
    sdata_ax.axhline(y=np.mean(overlap_data), color='g', alpha=0.8, linestyle='--', label="mean jaccard")
    sdata_ax.set_xlabel("order statistics", fontsize=20)
    sdata_ax.set_ylabel("jaccard coefficient", fontsize=20)
    sdata_ax.tick_params(axis='x', labelsize=20)
    sdata_ax.tick_params(axis='y', labelsize=20)
    sdata_ax.set_title("Sorted Jaccard coefficient values", fontsize=25)
    sdata_ax.legend(loc=2, fontsize=20)

    plt.savefig("img/eqtl_pqtl_intersection" + modules_type + "_linkage_overlap.png", dpi=300)
    plt.show()
    plt.close("all")
    gc.collect()


def plot_vignettes(module_type, module_names):
    for module_name in tqdm(module_names, desc="vignettes plotted"):
        plots = []
        for qtl_type in ["eQTLs", "pQTLs"]:
            path = "/".join(["./results", module_type, module_name, qtl_type])
            q_hist = util.try_imread(path + "/q_hist.png", mode="RGBA")
            pairwise_test = util.try_imread(path + "/pairwise_test.png", mode="RGBA")
            interactome_test = util.try_imread(path + "/ppi_test.png", mode="RGBA")
            network = util.try_imread(path + "/module_graph.png", mode="RGBA")

            vignette_side_parts = [img for img in [q_hist, pairwise_test, interactome_test, network]
                                   if img is not None]
            if len(vignette_side_parts) != 0:
                plots.append(np.vstack(vignette_side_parts))

        combined_plot = np.column_stack([img for img in plots if img is not None])
        plt.imsave("/".join(["./results", module_type, module_name + ".png"]), combined_plot)
        del q_hist, pairwise_test, interactome_test, network, vignette_side_parts, combined_plot
        plt.close("all")
        gc.collect()
        
def jaccard(s1, s2):
    """
    Notice: due to interpretation, similarity between empty sets is redefined to be 0 instead of 1.
    :return: Jaccard coefficient for sets s1 and s2
    """
    if len(s1) == 0 and len(s2) == 0:
        return 0
    return len(s1 & s2) / len(s1 | s2)


def linkage_similarity(module_graph, qtl_graph, mode="mean"):
    """
    :param module_graph: simple graph of some functional module
    :param qtl_graph: bipartite graph of linkages
    :param mode: 'full' — return vector of Jaccard coefficients representing all edges
                 'mean' — return only average linkage similarity
    :return: some statistics about linkage similarity specified by "mode"
    """
    if module_graph.ecount() == 0:
        return 0. if mode == "mean" else np.array([])

    results = np.zeros(shape=module_graph.ecount())
    for i, edge in enumerate(module_graph.es):
        source = module_graph.vs[edge.source]
        target = module_graph.vs[edge.target]

        try:
            s_neigh = set(qtl_graph.neighbors(source["name"], mode="IN"))
            t_neigh = set(qtl_graph.neighbors(target["name"], mode="IN"))
            results[i] = jaccard(s_neigh, t_neigh)
        except ValueError:
            '''Sometimes there is no such vertex in qtl graph'''
            pass

    return results.mean() if mode == "mean" else results

########################################################################################################################
#                                            LINKAGE SHARING ESTIMATION                                                #
########################################################################################################################
#----------------------------------------------------------------------------------------------------------------------#
#                                                    SIMPLE TESTS                                                      #     
#----------------------------------------------------------------------------------------------------------------------#

def community_graph_test(modules_dict, gene_pool, RANDITER_COUNT, qtl_df):
    """
    Каждый модуль преобразуется в полный граф, потом аналогичное "созвездие"
    симулируется из случайных наборов генов тех же размеров, измеряется
    среднее подобие и возвращается в качестве ответа.
    Вроде мы хотели этим тестом проверять на осмысленность наши наборы модулей,
    причём это был самый первый, предварительный тест в соотв. разделе.
    При этом я не уверен, что мы в итоге не отказались от этого теста вообще.
    """
    qtl_graph = networks.graph_from_edges(
        qtl_df[["SNP", "gene"]].values,
        directed=True
    )

    community_graph = ig.Graph()
    vertex_names = np.unique(np.hstack(modules_dict.values()))
    community_graph.add_vertices(vertex_names)
    community_graph.vs["name"] = vertex_names
    for module_name, module_genes in modules_dict.items(): # dafuq
        community_graph.add_edges(itertools.product(module_genes, repeat=2))
    community_graph.simplify()
    community_graph.vs.select(_degree=0).delete()
    
    real_avg_link_sim = linkage_similarity(community_graph, qtl_graph, mode='mean')

    rand_link_sim = []
    random_state = np.random.RandomState()
    for _ in itertools.repeat(None, RANDITER_COUNT):
        gene_samples = []
        for module_name, module_genes in modules_dict.items():
            gene_samples.append(random_state.choice(
                gene_pool, size=len(module_genes), replace=False
            ))
            
        community_graph = ig.Graph()
        vertex_names = np.unique(np.hstack(gene_samples))
        community_graph.add_vertices(vertex_names)
        community_graph.vs["name"] = vertex_names
        for sample in gene_samples:
            community_graph.add_edges(itertools.product(sample, repeat=2))
        community_graph.simplify()
        community_graph.vs.select(_degree=0).delete()
        
        rand_link_sim.append(linkage_similarity(community_graph, qtl_graph, mode='mean'))
        
    rand_link_sim = np.array(rand_link_sim)
    rv = sps.rv_discrete(values=(rand_link_sim, np.full(RANDITER_COUNT, 1. / RANDITER_COUNT)))
    ci = rv.interval(0.95)

    empirical_p_value = np.mean(rand_link_sim >= real_avg_link_sim)
    significant = empirical_p_value < 0.05
    random_avg_link_sim = np.mean(rand_link_sim)

    return real_avg_link_sim, random_avg_link_sim, significant, ci[0], ci[1]


def full_graph_test(module_genes, gene_pool, RANDITER_COUNT, qtl_df):
    """
    Compare linkage similarity averaged over all possible edges in a functional
    module and in random gene sets of the same size taken from gene pool.
    ---------------------------------------------------------------------------
    :param module_genes — iterable —
        names of genes the make up the functional module
    :param gene_pool — iterable — 
        set of genes used for sampling (usually the set of expressed genes)
    :param RANDITER_COUNT — uint —
        number of randomization iterations to carry out
    :param qtl_df — pd.DataFrame —
        pandas DataFrame that must contain such columns as "SNP" and "gene"
        "SNP" is a silly name choice, it actually means "marker"
    ---------------------------------------------------------------------------
    :returns (float, float, bool, float, float) —
        -   linkage similarity averaged over all possible pairs of nodes from
            the real functional module
        -   same, but for random samples of the same size, taken from the 
            gene_pool, and then averaged across randomization iterations
        -   whether the empirical p-value is smaller than 0.05 or not
        -   confidence interval
    ---------------------------------------------------------------------------
    :note   this test makes no sense for broad interaction categories
            don't run it on large chunks of the interactome graph
    """
    qtl_graph = networks.graph_from_edges(
        qtl_df[["SNP", "gene"]].values,
        directed=True
    )

    full_graph = ig.Graph.Full(len(module_genes))
    full_graph.vs["name"] = module_genes
    real_avg_link_sim = linkage_similarity(full_graph, qtl_graph, mode='mean')

    rand_link_sim = []
    random_state = np.random.RandomState()
    for _ in itertools.repeat(None, RANDITER_COUNT):
        full_graph.vs["name"] = random_state.choice(gene_pool, len(module_genes), replace=False)
        rand_link_sim.append(linkage_similarity(full_graph, qtl_graph, mode='mean'))

    rand_link_sim = np.array(rand_link_sim)
    rv = sps.rv_discrete(values=(rand_link_sim, np.full(RANDITER_COUNT, 1. / RANDITER_COUNT)))
    ci = rv.interval(0.95)

    empirical_p_value = np.mean(rand_link_sim >= real_avg_link_sim)
    significant = empirical_p_value < 0.05
    random_avg_link_sim = np.mean(rand_link_sim)

    return real_avg_link_sim, random_avg_link_sim, significant, ci[0], ci[1]

from tqdm import tqdm, tqdm_notebook
def ppin_test(module_genes, interactions_type, interactome_graph, qtl_df, 
              n_iters, path_to_randomized, return_pvalues=False):
    """
    Compute average linkage similiary for each edge of the function module.
    Repeat the same computation, but substitute the real interactome graph
    with its simulation obtained through structure-preserving randomization
    of the original graph.
    ---------------------------------------------------------------------------
    :param module_genes — iterable —
        names of genes the make up the functional module
    :param interactions_type — str —
        graphs of all, genetic and physical interactions were randomized
        in advance; this parameter is used to conclude which versions to use
    :param interactome_graph — ig.Graph — 
        original real interactome graph (note: edge type mustn't be the 
        one specified by interactions_type, but must align with the latter: 
        the passed graph must be a subgraph of interactions_type-specified
        interactome graph
    :param qtl_df — pd.DataFrame —
        pandas DataFrame that must contain such columns as "SNP" and "gene"
        "SNP" is a silly name choice, it actually means "marker"
    :param n_iters — uint —
        number of randomized graph to process
    :param path_to_randomized — str —
        path to the directory with the randomized graphs
    :param return_pvalues — bool —
        when set to True forces the function to return the vector of p-values
        instead of the aggregate
    ---------------------------------------------------------------------------
    :returns 
        -   linkage similarity averaged over edges of the real functional module
        -   same, but for randomized versions of the interactome graph
        -   aggregated
    ---------------------------------------------------------------------------
    :note   this test can be used both with the broad interaction categories
            like "all", "physical", "genetic" and with the functional modules
            
    """
    qtl_graph = networks.graph_from_edges(qtl_df[["SNP", "gene"]].values, directed=True)
    real_linksimvec = linkage_similarity(
        interactome_graph.subgraph(
            np.intersect1d(interactome_graph.vs["name"],
                           module_genes)
        ),
        qtl_graph,
        mode='full'
    )
    randomized_j_means = []
    p_values = []
    for iter_num in tqdm(range(n_iters), "n_iters"):
        randomized_interactome_graph = ig.Graph().Read_Pickle(
            os.path.join(path_to_randomized, 
                         "{0}/{1}.pkl".format(interactions_type, iter_num))
        ).simplify()
        randomized_interactome_graph.vs.select(_degree=0).delete()
        randomized_linksimvec = linkage_similarity(
            randomized_interactome_graph.subgraph(
                np.intersect1d(randomized_interactome_graph.vs["name"],
                               module_genes)
            ),
            qtl_graph,
            mode='full'
        )
        if len(randomized_linksimvec) == 0:
            p_values.append(1)
        else:
            try: 
                p_values.append(sps.mannwhitneyu(real_linksimvec, randomized_linksimvec, alternative="two-sided")[1])
            except ValueError: # "all numbers are identical in mannwhitneyu"
                p_values.append(1)
        randomized_j_means.append(np.mean(randomized_linksimvec) 
                                  if len(randomized_linksimvec) != 0 else 0)

    real_avg_linksim = np.mean(real_linksimvec) if len(real_linksimvec) != 0 else 0
    significant = np.count_nonzero(np.array(p_values) < 0.05) > (0.1 * len(p_values))

    return real_avg_linksim, randomized_j_means, (p_values if return_pvalues else significant)

#----------------------------------------------------------------------------------------------------------------------#
#                                          LINKAGE SHARING ANALYZER (OBSOLETE)                                         #     
#----------------------------------------------------------------------------------------------------------------------#
class LinkageSharingStatistics:
    """Storage class for analysis results produced by the LinkageSharingAnalyzer class"""

    def __init__(self, analysis_type, interactions_type, module_name, module_type, qtl_type, qtl_df=None,
                 q_thresholds=None, test_results_df=None, module_genes=None):
        self.analysis_type = analysis_type
        self.interactions_type = interactions_type

        self.module_genes = module_genes
        self.module_name = module_name
        self.module_type = module_type

        self.qtl_df = qtl_df
        self.qtl_type = qtl_type

        self.test_results_df = test_results_df
        self.q_thresholds = q_thresholds

        self.destdir = "/".join(["./results", self.module_type, self.module_name, self.qtl_type])
        util.ensure_dir(self.destdir)

    def save(self):
        with open(self.destdir + "/" + self.analysis_type + "_test.pkl", 'wb') as pickle_file:
            pickle.dump(obj=self.__dict__, file=pickle_file)

    def load(self):
        with open(self.destdir + "/" + self.analysis_type + "_test.pkl", 'rb') as pickle_file:
            self.__dict__.update(pickle.load(file=pickle_file))

    def plot_test_results(self):
        plot_test_results(
            module_name=self.module_name,
            module_type=self.module_type,
            test_results_df=self.test_results_df,
            test_type=self.analysis_type,
            q_thresholds=self.q_thresholds,
            qtl_type=self.qtl_type
        )


def plot_analysis_results(expression_df, interactions_type, interactome_graph,
                          modules_dict, modules_type, qtl_type, qtl_df):
    module_results_dir = "./results/" + modules_type
    module_names = collections.defaultdict(list)
    for dirname in tqdm(os.listdir(module_results_dir), desc="subdirectories processed"):
        if os.path.isdir(module_results_dir + "/" + dirname):
            dirpath = '/'.join([module_results_dir, dirname, qtl_type])
            if os.path.exists(dirpath):
                module_names[qtl_type].append(dirname)
                for test_type in ["ppi", "pairwise"]:
                    test_stats = LinkageSharingStatistics(
                        analysis_type=test_type,
                        interactions_type=interactions_type,
                        module_name=dirname,
                        module_type=modules_type,
                        qtl_type=qtl_type
                    )
                    test_stats.load()
                    test_stats.plot_test_results()
                    del test_stats
                    gc.collect()

    LinkageSharingAnalyzer(
        expression_df=expression_df,
        interactions_type=interactions_type,
        interactome_graph=interactome_graph,
        modules_type=modules_type,
        modules_dict={module_name: modules_dict[module_name] for module_name in module_names[qtl_type]},
        qtl_type=qtl_type,
        qtl_df=qtl_df,
        q_thresholds=np.logspace(-5, -2, 10)
    ).plot_basic_module_info()
    plot_vignettes(module_type=modules_type, module_names=module_names[qtl_type])


class LinkageSharingAnalyzer:
    """
    This class was designed for bulk processing of functional modules.
    It's an early prototype of map-reduce classes used to generate figures.
    
    It is used to analyze to which extent interacting genes tend to share linkages.
    Given a datafame of linkages, a functional module graph and a list of q-value thresholds,
    this class perturbs the graph ~100 times and each time calculates averaged Jaccard coefficient
    for each threshold. Then, to ensure robustness of conclusions, simulation results are averaged and
    plotted altogether with 95% confidence intervals against actual data for the original graph.
    """

    # [-------------------------------------------------PUBLIC METHODS-------------------------------------------------]
    def __init__(self, expression_df, interactions_type, interactome_graph, modules_type, modules_dict,
                 qtl_type, qtl_df, q_thresholds, pairwise_test_iter=200, ppi_test_iter=1024):
        self.expression_df = expression_df
        self.interactions_type = interactions_type
        self.interactome = interactome_graph.simplify()
        self.interactome.vs.select(_degree=0).delete()
        self.qtl_type = qtl_type
        self.qtl_df = qtl_df
        self.modules_type = modules_type
        self.modules_dict = modules_dict
        self.qval_list = q_thresholds

        for module_name in self.modules_dict.keys():
            self.modules_dict[module_name] = np.intersect1d(self.modules_dict[module_name],
                                                            self.interactome.vs["name"])

        self.randomized_interactome = None
        self._default_row = np.array([0, 0, 1])

        self.pairwise_randiter_count = pairwise_test_iter
        self.ppi_randiter_count = ppi_test_iter

    def plot_basic_module_info(self):
        for module_name in tqdm(self.modules_dict.keys(), "graphs and hists plotted"):
            destdir = "/".join(["./results", self.modules_type, module_name, self.qtl_type])
            util.ensure_dir(destdir)

            module_graph = self.interactome.subgraph(
                set(self.interactome.vs["name"])
                & set(self.modules_dict[module_name])
            )

            plot_module_graph(
                destdir=destdir,
                module_graph=module_graph,
                qtl_df=self.qtl_df[self.qtl_df["gene"].isin(self.modules_dict[module_name])]
            )

            plot_q_hist(
                destdir=destdir,
                module_name=module_name,
                module_graph=module_graph,
                module_type=self.modules_type,
                qtl_df=self.qtl_df[self.qtl_df["gene"].isin(self.modules_dict[module_name])],
                qtl_type=self.qtl_type
            )

            del module_graph

    def pairwise_test(self):
        pool = mp.Pool(mp.cpu_count())
        real_means, randomized_means, empirical_p_values, lo_cis, hi_cis = \
            np.stack(
                pool.map(self._pairwise_process_threshold, self.qval_list),
                axis=1
            )
        pool.close()
        pool.join()

        results_df_dict = collections.defaultdict()
        for i, module_name in enumerate(self.modules_dict.keys()):
            results_df_dict[module_name] = pd.DataFrame(
                data=np.column_stack((
                    real_means[:, i],
                    randomized_means[:, i],
                    empirical_p_values[:, i],
                    empirical_p_values[:, i] < 0.05,
                    lo_cis[:, i],
                    hi_cis[:, i]
                )),
                columns=["real", "simulated", "p_value", "significant", "lo_ci", "hi_ci"]
            )

        for module_name in self.modules_dict.keys():
            pairwise_test_results_df = results_df_dict[module_name]

            LinkageSharingStatistics(
                analysis_type="pairwise",
                interactions_type=self.interactions_type,
                qtl_type=self.qtl_type,
                qtl_df=self.qtl_df[self.qtl_df["gene"].isin(self.modules_dict[module_name])],
                module_type=self.modules_type,
                module_name=module_name,
                module_genes=self.modules_dict[module_name],
                test_results_df=pairwise_test_results_df,
                q_thresholds=self.qval_list
            ).save()

    def ppi_test(self):
        """
        Interface function. Performs the simulation and saves the plot with results.
        return: Some measure of curve similarity between real and simulated data.
        """

        pool = mp.Pool(mp.cpu_count())
        raw_results_dict = pool.map(self._ppi_process_randiter, range(self.ppi_randiter_count))
        pool.close()
        pool.join()

        results_dict = {
            module_name: np.stack([rawd[module_name] for rawd in raw_results_dict])
            for module_name in self.modules_dict.keys()
        }

        ci_mx = {}
        for module_name in results_dict.keys():
            ci_mx[module_name] = np.apply_along_axis(
                func1d=lambda l: sps.rv_discrete(values=(l, np.full(len(l), 1 / len(l)))).interval(0.95),
                arr=results_dict[module_name][:, :, 1], axis=0).T

            j_columns = np.mean(results_dict[module_name], axis=0)
            significant = (
                    np.count_nonzero(results_dict[module_name][:, :, 2] < 0.05, axis=0)
                    >= self.ppi_randiter_count // 2
            )[np.newaxis].T

            results_dict[module_name] = np.column_stack((j_columns, significant))

        for module_name, module_vertices in self.modules_dict.items():
            ppi_test_results_df = pd.DataFrame(
                np.column_stack((results_dict[module_name], ci_mx[module_name])),
                columns=["real", "simulated", "p_value", "significant", "lo_ci", "hi_ci"]
            )

            LinkageSharingStatistics(
                analysis_type="ppi",
                interactions_type=self.interactions_type,
                qtl_type=self.qtl_type,
                qtl_df=self.qtl_df[self.qtl_df["gene"].isin(module_vertices)],
                module_type=self.modules_type,
                module_name=module_name,
                module_genes=self.modules_dict[module_name],
                test_results_df=ppi_test_results_df,
                q_thresholds=self.qval_list
            ).save()

    # [------------------------------------------------PRIVATE METHODS------------------------------------------------]
    def _pairwise_process_threshold(self, q_value_threshold):
        random_state = np.random.RandomState()
        gene_pool = self.expression_df["gene"].values
        qtl_graph = networks.graph_from_edges(
            self.qtl_df[self.qtl_df["q_value"] <= q_value_threshold][["SNP", "gene"]].values,
            directed=True
        )

        real_avg_link_sim = []
        random_avg_link_sim = []
        empirical_p_values = []
        lo_cis = []
        hi_cis = []

        for module_name, module_genes in self.modules_dict.items():
            full_g = ig.Graph.Full(len(module_genes))
            full_g.vs["name"] = module_genes
            real_link_sim = linkage_similarity(full_g, qtl_graph, mode='mean')
            real_avg_link_sim.append(real_link_sim)

            rand_link_sim = []
            for i in range(self.pairwise_randiter_count):
                full_g.vs["name"] = random_state.choice(gene_pool, len(module_genes), replace=False)
                rand_link_sim.append(linkage_similarity(full_g, qtl_graph, mode='mean'))

            rand_link_sim = np.array(rand_link_sim)
            rv = sps.rv_discrete(values=(rand_link_sim, np.full(len(rand_link_sim), 1. / len(rand_link_sim))))
            ci = rv.interval(0.95)

            lo_cis.append(ci[0])
            hi_cis.append(ci[1])

            empirical_p_values.append(np.mean(rand_link_sim >= real_link_sim))
            random_avg_link_sim.append(rv.median())

        return real_avg_link_sim, random_avg_link_sim, empirical_p_values, lo_cis, hi_cis

    def _ppi_process_threshold(self, q_value_threshold, module_name):
        qtl_graph = networks.graph_from_edges(
            edges=self.qtl_df[self.qtl_df["q_value"] <= q_value_threshold][["SNP", "gene"]].values,
            directed=True
        )
        genes_with_linkages = qtl_graph.vs.select(part=1)["name"]
        vertex_names = self.modules_dict[module_name]

        if set(genes_with_linkages).isdisjoint(set(vertex_names)):
            return self._default_row

        real_stats = linkage_similarity(
            module_graph=self.interactome.subgraph(
                set(self.interactome.vs["name"])
                & set(vertex_names)
            ),
            qtl_graph=qtl_graph,
            mode="full"
        )

        if len(real_stats) == 0 or all([sim_coeff == 0 for sim_coeff in real_stats]):
            return self._default_row

        randomized_stats = linkage_similarity(
            module_graph=self.randomized_interactome.subgraph(
                set(self.interactome.vs["name"]) &
                set(vertex_names)
            ),
            qtl_graph=qtl_graph,
            mode="full"
        )
        p_value = 0.
        if len(randomized_stats) != 0:
            try:
                _, p_value = sps.mannwhitneyu(real_stats, randomized_stats, alternative="two-sided")
            except ValueError:
                # prevents "all numbers are equal" error that is thrown when,
                # for example, [1, 1, 1] and [1, 1] are passed to MWU
                # it's really weird, but simply ignoring it seems the be best workaround so far
                pass

        return np.array([
            real_stats.mean(),
            randomized_stats.mean() if len(randomized_stats) != 0 else 0.,
            p_value
        ])

    def _ppi_process_randiter(self, iter_num):
        results_dict = collections.defaultdict(list)
        self.randomized_interactome = ig.Graph().Read_Pickle(
            "./data/randomized_interactome_copies/{0}/{1}.pkl"
                .format(self.interactions_type, iter_num)
        )
        for module_name, module_vertices in self.modules_dict.items():
            randiter_results = np.vstack(
                [self._ppi_process_threshold(q_thr, module_name)
                 for q_thr in self.qval_list]
            )
            if len(results_dict[module_name]) == 0:
                results_dict[module_name] = randiter_results
            else:
                results_dict[module_name] = np.dstack(
                    (results_dict[module_name], randiter_results)
                )
        return results_dict


#########################################################################################################################
#                                              JOINT EQTL/PQTL ANALYSIS                                                 #
#########################################################################################################################

import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

#----------------------------------------------------------------------------------------------------------------------#
#                                                    SIMPLE TESTS                                                      #     
#----------------------------------------------------------------------------------------------------------------------#

def qtl_overlap_test(eqtl_df, pqtl_df, gene_pool, modules_dict, randiter_count=100):
    linked_eQTLs, linked_pQTLs = set(), set()
    for module_name, module_genes in modules_dict.items():
        module_genes = np.intersect1d(module_genes, gene_pool)
        linked_eQTLs.update(set(linked_markers(eqtl_df, module_genes)))
        linked_pQTLs.update(set(linked_markers(pqtl_df, module_genes)))
    real_overlap = jaccard(linked_eQTLs, linked_pQTLs)

    random_state = np.random.RandomState()
    randomized_overlap = []
    for _ in range(randiter_count):
        linked_eQTLs.clear(), linked_pQTLs.clear()
        for module_name, module_genes in modules_dict.items():
            sample_genes = random_state.choice(gene_pool, len(module_genes), replace=False)
            linked_eQTLs.update(set(linked_markers(eqtl_df, sample_genes)))
            linked_pQTLs.update(set(linked_markers(pqtl_df, sample_genes)))
        randomized_overlap.append(jaccard(linked_eQTLs, linked_pQTLs))
    return real_overlap, np.mean(randomized_overlap)


def qtl_overlap_by_module_test(eqtl_df, pqtl_df, gene_pool, modules_dict, randiter_count=100):
    qtl_intersection_j = []
    for module_name, module_genes in modules_dict.items():
        module_genes = np.intersect1d(module_genes, gene_pool)
        linked_eQTLs = set(linked_markers(eqtl_df, module_genes))
        linked_pQTLs = set(linked_markers(pqtl_df, module_genes))
        qtl_intersection_j.append(jaccard(linked_eQTLs, linked_pQTLs))

    random_state = np.random.RandomState()
    randomized_qtl_intersection_j = []
    for module_name, module_genes in tqdm(modules_dict.items(), desc="random samples generated"):
        buffer = []
        for _ in range(randiter_count):
            sample_genes = random_state.choice(gene_pool, len(module_genes), replace=False)
            linked_eQTLs = set(linked_markers(eqtl_df, sample_genes))
            linked_pQTLs = set(linked_markers(pqtl_df, sample_genes))
            buffer.append(jaccard(linked_eQTLs, linked_pQTLs))
        randomized_qtl_intersection_j.append(np.median(buffer))
    return qtl_intersection_j, randomized_qtl_intersection_j

"""
WARNING: Inconsistent with the current directory structure!
"""
class PqtlPredictor:
    """
    Predicts pQTLs from eQTLs using PPI networks.
    Rationale:
        1. Good prediction quality indicates hidden relationship between eQTLs and pQTLs.
        2. Less candidates implies more pQTLs for the same FDR-cutoff in this context.
    Algorithm:
        1. Form a set of candidate linkages — (marker, gene) pairs.
            For each gene g, markers = eQTLs(g) U eQTLs(N(g)),
            where N(g) is a neighborhood of g in the PPI graph.
        2. Test for linkage naively, using MWU test and q-values (for FDR correction).
    """

    # [-------------------------------------------------PUBLIC METHODS-------------------------------------------------]
    def __init__(self, eqtls_df, pqtls_df,
                 pqtls_expression_df, pqtls_genotypes_df,
                 module_name, module_genes,
                 interactome_graph):
        self.pqtls_expr_df = pqtls_expression_df
        self.pqtls_gen_df = pqtls_genotypes_df

        # To parallelize p-value computation we must filter out genes
        # with no protein abundance data and markers not present in yeast
        self.possible_pQTL_markers = pqtls_genotypes_df["SNP"].values
        self.genes_with_protein_abundance = pqtls_expression_df["gene"].values

        self.eqtls_df = eqtls_df
        self.pqtls_df = pqtls_df

        # These are frequently queried during p-value calculations.
        # Conversion pandas.DataFrame to np.matrix allows to achieve twofold speed-up
        self.pqtls_expr_mx = \
            pqtls_expression_df.as_matrix(
                columns=pqtls_expression_df.columns[1:]
            )
        self.pqtls_gen_mx = \
            pqtls_genotypes_df.as_matrix(
                columns=pqtls_genotypes_df.columns[1:]
            )

        self.module_name = module_name.replace(' ', '_')
        self.module_genes = module_genes
        self.interactome_graph = interactome_graph

    def predict(self, fdr_cutoff=0.05):
        """
        Predicts pQTLs from eQTLs and returns analysis summary.

        :param fdr_cutoff: q-value threshold for predicted pQTLs
        :return: 5-tuple of a kind:
            (number of common_pQTLs,
            fraction of common pQTLs,
            number of initial pQTLs,
            number of predicted pQTLs,
            number of new pQTLs)
        """
        # 1. Construct a set of candidate linkages.
        results = [
            chunk for chunk in
            [self._get_possible_linkages(gene_name) for gene_name
             in np.intersect1d(self.module_genes, self.genes_with_protein_abundance)]
            if len(chunk) > 0
        ]

        if len(results) == 0:
            return (0, 0, self.pqtls_df.shape[0], 0, 0)

        # 2. Do QTL-mapping naively: MWU test with q-values.
        possible_linkages = np.vstack(results)
        p_values = np.array([self._compute_p_value(linkage) for linkage in possible_linkages])

        pd.DataFrame(p_values, columns=["pvalue"]).to_csv("./data/pQTLs/temp/pvalues_{}.csv".format(self.module_name),
                                                          sep='\t', index=False)
        path_to_rscript = "Rscript" if os.environ["HOME"] == "/home/ivanov_vv" else "/home/vvi/anaconda3/src/Rscript"
        subprocess.call([path_to_rscript, './src/p-values_to_q-values.R', self.module_name], shell=False)
        q_values = pd.read_table("./data/pQTLs/temp/qvalues_{}.csv".format(self.module_name)).values
        reject = q_values <= fdr_cutoff

        new_pqtls_df = pd.DataFrame(
            np.column_stack((possible_linkages, p_values, q_values, reject)),
            columns=["SNP", "gene", "pvalue", "qvalue", "reject"]
        ).query("reject == True")

        common = len(set(map(tuple, self.pqtls_df[["SNP", "gene"]].values))
                     & set(map(tuple, new_pqtls_df[["SNP", "gene"]].values)))

        old_pqtls_count = self.pqtls_df[self.pqtls_df["gene"].isin(self.module_genes)].shape[0]
        return (common, 100 * common / old_pqtls_count, old_pqtls_count,
                new_pqtls_df.shape[0], new_pqtls_df.shape[0] - common)

    # [------------------------------------------------PRIVATE METHODS------------------------------------------------]

    def _get_markers_linked_to_neighbors(self, gene_name):
        try:
            neighborhood = self.interactome_graph.vs[self.interactome_graph.neighborhood(gene_name)]["name"]
        except ValueError:
            neighborhood = []
        return linked_markers(qtl_df=self.eqtls_df, gene_list=np.append(np.array(gene_name), neighborhood))

    def _get_possible_linkages(self, gene_name):
        linked_markers = self._get_markers_linked_to_neighbors(gene_name)
        return np.column_stack((linked_markers, np.full(len(linked_markers), gene_name)))

    def _compute_p_value(self, candidate_linkage):
        marker_name, gene_name = candidate_linkage
        genotype_rowmask = self.pqtls_gen_df["SNP"] == marker_name
        genotype_row = self.pqtls_gen_mx[genotype_rowmask]
        expression_rowmask = self.pqtls_expr_df["gene"] == gene_name
        expression_row = self.pqtls_expr_mx[expression_rowmask]
        # Divide expression data by inherited marker variant
        from_BY = expression_row[genotype_row == 0]
        from_RM = expression_row[genotype_row == 2]
        _, p_value = stats.mannwhitneyu(from_BY, from_RM, alternative="two-sided")
        return p_value
