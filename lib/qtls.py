from __future__ import print_function

import collections
import multiprocessing as mp
import os
import subprocess

import igraph as ig
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.stats as stats
from scipy import ndimage
from scipy.spatial import distance

import networks
import util


########################################################################################################################
#                                                     UTILITIES                                                        #
########################################################################################################################


def linked_markers(qtl_df, gene_list):
    return np.unique(qtl_df[qtl_df["gene"].isin(gene_list)]["SNP"].values)


def linked_genes(qtl_df, marker_list):
    return np.unique(qtl_df[qtl_df["SNP"].isin(marker_list)]["gene"].values)


########################################################################################################################
#                                                  VISUALIZATIONS                                                      #
########################################################################################################################


'''
TODO: WHY NOT DOWNLOAD COORDINATES FOR ALL GENE FROM SOME DATABASE?
'''

def linkages2gencoords(qtl_df, full_genotypes_df):
    '''
    Given a dataframe of estimated QTL linkages and
    full strain genotype lookup table, plot the number
    of linkages against marker location in genome

    :param qtl_df: must contain columns "SNP" and "gene"
    :param full_genotypes_df: a list o genes sorted by genetic coordinate in ascending order
    :return: a list of linkages sorted in the same order
    '''
    qtl_graph = networks.graph_from_edges(
        edges=qtl_df[["SNP", "gene"]].values,
        directed=True
    )
    # Select only the marker vertices
    left_part = qtl_graph.vs.select(part=False)
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
    qtl_marker_to_linkages = \
        dict(
            zip(
                full_genotypes_df.iloc[:, 0],
                np.zeros(full_genotypes_df.shape[0])
            )
        )
    # Add the available linkage data
    qtl_marker_to_linkages.update(dict(
        zip(
            left_part["name"],
            left_part.outdegree()
        )
    ))
    # Sort the dictionary accordingly to the marker position on the chromosome
    # and unzip it to extract the ordered linkage list

    qtl_x, qtl_y = \
        map(list, zip(*sorted(qtl_marker_to_linkages.items(),
                              key=lambda p: marker_to_rownum[p[0]])))
    return qtl_x, qtl_y


def plot_module_graph(module_name, module_graph, qtl_type, qtl_df, dirname, format="png"):
    '''
    Plot functional module graph in .svg format with vertices colored in accordance to
    the number of linked QTLs. Vertices with no linkages are not shown.
    All other vertices have number of linked QTLs written in their labels.

    :param module_graph: interactome subgraph
    :param qtl_df: s/e
    :param dirname: s/e
    :param module_name: s/e
    :return: None
    '''
    num_linkages = [len(linked_markers(qtl_df, [gene_name])) for gene_name in module_graph.vs["name"]]
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
        target=dirname + '/' + qtl_type + "_" + module_name + "." + format,
        layout=module_graph.layout_fruchterman_reingold(maxiter=1000, area=n**3, repulserad=n**3),
        bbox=(1440, 720),
        vertex_label=[gene_name + ", " + str(num_linkages[i]) if num_linkages[i] != 0 else ''
                for i, gene_name in enumerate(module_graph.vs["name"])],
        margin=50,
    )


def plot_q_value_hist(destdir, module_graph, module_name, module_type, qtl_df, qtl_type):
    fig, ax = plt.subplots(figsize=(20, 10))
    plt.xscale("log")
    counts, bins, patches = plt.hist(
        qtl_df["q.value"],
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
                    qtl_df["q.value"].mean(),
                    qtl_df["q.value"].median(),
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

    plt.savefig(destdir + "/{}_q_value_distribution for {}; {}.png".format(
        qtl_type, module_name, module_type)
                )
    plt.close()


def plot_pairwise_test(linkage_sharing_df, module_name, module_type, q_thresholds, qtl_type):
    # for the plot to be informative, don't show the range
    # where there are no linkage with such q-values at all
    results_df = linkage_sharing_df.query("real_mean > 0")

    if results_df.shape[0] < 2:
        print("{} — no data!".format(module_name))
        return

    trimmed_qval_list = q_thresholds[-results_df.shape[0]:]

    fig, ax = plt.subplots(figsize=(20, 10))
    ax.set_title("{} linkage similarity for {}; {}".format(
        qtl_type, module_name, module_type
    ), fontsize=20)
    plt.xscale('log')

    ax.set_xlabel("linkage q-value threshold, lg-scale", fontsize=20)
    ax.set_ylabel("average linkage similarity", fontsize=20)
    ax.plot(trimmed_qval_list, results_df["real_mean"], label="functional module")
    ax.plot(trimmed_qval_list, results_df["simulated_mean"],
            label="random gene set\n")
    ax.legend(fontsize=15)

    ax.grid(linestyle='dotted')
    plt.setp(ax.get_xticklabels(), fontsize=15)
    plt.setp(ax.get_yticklabels(), fontsize=12)

    destdir = '/'.join(["./results", module_type, module_name, qtl_type])
    util.ensure_dir(destdir)
    plt.savefig(destdir + "/pairwise_{}_{}_{}.png".format(qtl_type, module_name, module_type))
    plt.close()


def plot_interactome_test(module_name, module_type, linkage_sharing_df,
                          robustness_score, q_thresholds, qtl_type):
    # for the plot to be informative, don't show the range
    # where there are no linkage with such q-values at all
    results_df = linkage_sharing_df[
        linkage_sharing_df["real_mean"] != 0
        ]

    if results_df.shape[0] < 2:
        print("{} — no data!".format(module_name))
        return

    trimmed_qval_list = q_thresholds[-results_df.shape[0]:]

    fig, ax1 = plt.subplots(figsize=(20, 10))
    ax1.set_title("{} linkage similarity for {}; {}".format(
        qtl_type, module_name, module_type
    ), fontsize=20)
    plt.xscale('log')
    ax2 = ax1.twinx()

    ax1.set_xlabel("linkage q-value threshold, lg-scale", fontsize=20)
    ax1.set_ylabel("average linkage similarity", fontsize=20)
    ax1.plot(trimmed_qval_list, results_df["real_mean"], label="real data")
    ax1.plot(trimmed_qval_list, results_df["simulated_mean"],
             label="synthesized data\n"
                   "curve similarity score (WMSQ): {0:.4f}".format(robustness_score)
             )

    ax2.plot(trimmed_qval_list, results_df["p_value"],
             color="r", linestyle="dashed", alpha=0.5,
             label="p-values\nmean p-value: {0:.4f}\nmedian p-value: {1:.4f}".format(
                 results_df["p_value"].mean(),
                 results_df["p_value"].median()
             )
    )
    ax2.set_ylabel("p-value of difference of edge scores", fontsize=20)

    lines, labels = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax2.legend(lines + lines2, labels + labels2, loc=2, fontsize=15)

    ax1.grid(linestyle='dotted')
    plt.setp(ax1.get_xticklabels(), fontsize=15)
    plt.setp(ax1.get_yticklabels(), fontsize=12)
    plt.setp(ax2.get_yticklabels(), fontsize=12)

    destdir = '/'.join(["./results", module_type, module_name, qtl_type])
    util.ensure_dir(destdir)
    plt.savefig(destdir + "/interactome_{}_{}_{}.png".format(
        qtl_type, module_name, module_type)
                )
    plt.close()

########################################################################################################################
#                                            LINKAGE SHARING ESTIMATION                                                #
########################################################################################################################

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
    :param interaction_graph: some subgraph of the interactome
    :param qtl_graph: bipartite graph of linkages
    :param mode: 'full' — return vector of Jaccard coefficients representing all edges
                 'mean' — return only average linkage similarity
    :return: some statistics about linkage similarity specified by "mode"
    """
    results = np.zeros(shape=module_graph.ecount())
    if mode == 'mean' and len(results) == 0:
        return 0.

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


class LinkageSharingStatistics:
    ''' Storage class for analysis results produced by the LinkageSharingAnalyzer class.'''

    def __init__(self, analysis_type, qtl_type, module_type, module_name, qtl_df=None, module_graph=None,
                 q_thresholds=None, linkage_sharing_df=None, robustness_score=None):
        self.analysis_type = analysis_type
        self.qtl_type = qtl_type
        self.module_type = module_type
        self.module_name = module_name

        self.qtl_df = qtl_df
        self.module_graph = module_graph
        self.q_thresholds = q_thresholds
        self.linkage_sharing_df = linkage_sharing_df
        self.robustness_score = robustness_score

        self.destdir = "/".join(["./results", self.analysis_type, self.module_type, self.module_name, self.qtl_type])
        util.ensure_dir(self.destdir)

    ''' IMPLEMENT PICKLING AND DEPICKLING '''
    def save(self):
        self.qtl_df.to_csv(
            "/".join([self.destdir, "qtl_df.csv"]),
            sep='\t', index=False
        )
        self.linkage_sharing_df.to_csv(
            "/".join([self.destdir, "linkage_sharing.csv"]),
            sep='\t', index=False
        )
        ig.write(
            graph=self.module_graph,
            filename="/".join([self.destdir, "module_graph.pkl"]),
            format="pickle"
        )
        np.savetxt("/".join([self.destdir, "q_thresholds.csv"]), self.q_thresholds, delimiter='\t')
        vardict = vars(self)
        pd.DataFrame.from_dict(
            data={
                key:vardict[key] for key in
                ["qtl_type", "module_type", "module_name", "robustness_score"]
            },
            orient="index"
        ).to_csv(
            "/".join([self.destdir, "statistics.csv"]),
            sep='\t', header=None
        )

    def load(self):
        self.qtl_df = pd.read_csv("/".join([self.destdir, "qtl_df.csv"]), sep='\t')
        self.linkage_sharing_df = pd.read_csv("/".join([self.destdir, "linkage_sharing.csv"]), sep='\t')
        self.module_graph = ig.Graph().Read_Pickle(fname="/".join([self.destdir, "module_graph.pkl"]))
        self.q_thresholds = np.recfromcsv(fname="/".join([self.destdir, "q_thresholds.csv"]), delimiter='\t')

        vardict = pd.read_csv("/".join([self.destdir, "statistics.csv"]), index_col=0, header=None, sep='\t').T\
            .to_dict(orient="records")[0]
        for key, val in vardict.items():
            if key == "robustness_score":
                self.robustness_score = float(val)
            else:
                self.__dict__[key] = val

    def features_dict(self):
        return collections.OrderedDict([
            ("robustness_score", self.robustness_score),
            ("mean_q_score", self.qtl_df["q.value"].mean()),
            ("median_q_score", self.qtl_df["q.value"].median()),
            ("mean_p_score", self.linkage_sharing_df["p_value"].mean()),
            ("median_p_score", self.linkage_sharing_df["p_value"].median())])

    def plot_module_graph(self, format):
        plot_module_graph(
            module_name=self.module_name,
            module_graph=self.module_graph,
            qtl_type=self.qtl_type,
            qtl_df=self.qtl_df,
            dirname=self.destdir,
            format=format
        )

    def plot_q_value_hist(self):
        plot_q_value_hist(
            destdir=self.destdir,
            module_graph=self.module_graph,
            module_name=self.module_name,
            module_type=self.module_type,
            qtl_df=self.qtl_df,
            qtl_type=self.qtl_type
        )


# Rewrite this, for god's sake
def plot_vignettes(module_type, module_names):

    for module_name in module_names:
        elplot, eqplot, epplot, egplot = None, None, None, None
        plplot, pqplot, ppplot, pgplot = None, None, None, None

        epath = "/".join(["./results", module_type, module_name, "eQTLs"])
        ppath = "/".join(["./results", module_type, module_name, "pQTLs"])

        if not os.path.exists(epath) and not os.path.exists(ppath):
            continue

        if os.path.exists(epath):
            elplot = ndimage.imread(epath + "/{}_{}_{}.png".format(
                "eQTLs", module_name, module_type), mode="RGBA")
            eqplot = ndimage.imread(epath + "/{}_q_value_distribution for {}; {}.png".format(
                "eQTLs", module_name, module_type), mode="RGBA")
            egplot = ndimage.imread(epath + "/{}_{}.png".format("eQTLs", module_name), mode="RGBA")
        if os.path.exists(ppath):
            plplot = ndimage.imread(ppath + "/{}_{}_{}.png".format(
                "pQTLs", module_name, module_type), mode="RGBA")
            pqplot = ndimage.imread(ppath + "/{}_q_value_distribution for {}; {}.png".format(
                "pQTLs", module_name, module_type), mode="RGBA")
            pgplot = ndimage.imread(ppath + "/{}_{}.png".format("pQTLs", module_name), mode="RGBA")

        ecplot, pcplot = None, None
        if elplot is not None:
            ecplot = np.vstack((eqplot, elplot, egplot))
        if plplot is not None:
            pcplot = np.vstack((pqplot, plplot, pgplot))

        if ecplot is not None and pcplot is not None:
            result = np.hstack((ecplot, pcplot))
        else:
            result = ecplot if ecplot is not None else pcplot
        plt.imsave("/".join(["./results", module_type, module_name + ".png"]), result)


class LinkageSharingAnalyzer:
    """
    Used to analyze to which extent interacting genes tend to share linkages.
    Given a datafame of linkages, a functional module graph and a list of q-value thresholds,
    this class perturbs the graph ~100 times and each time calculates averaged Jaccard coefficient
    for each threshold. Then, to ensure robustness of conclusions, simulation results are averaged and
    plotted altogether with 95% confidence intervals against actual data for the original graph.
    """

# [-------------------------------------------------PUBLIC METHODS-------------------------------------------------]
    def __init__(self, expression_df, interactome_graph, modules_type, module_dict, qtl_type, qtl_df, q_thresholds):
        self.expression_df = expression_df
        self.interactome = interactome_graph.simplify()
        self.interactome.vs.select(_degree=0).delete()
        self.qtl_type = qtl_type
        self.qtl_df = qtl_df
        self.modules_type = modules_type
        self.module_dict = module_dict
        self.qval_list = q_thresholds

        for module_name in self.module_dict.keys():
            self.module_dict[module_name] = np.intersect1d(self.module_dict[module_name],
                                                           self.interactome.vs["name"])

        self.statistics = collections.defaultdict()
        self.randomized_interactome = None

        self._default_row = np.zeros(3)

    def process_modules(self):
        self.all_pairs_test()
        self.interactome_test()

    def all_pairs_test(self, randiter_count=64):
        # Почему бы не сравнивать их просто со средним значением для всего графа?
        gene_pool = self.expression_df["gene"].values
        real_results, randomized_results = None, None

        for q_thr in self.qval_list:
            qtl_graph = networks.graph_from_edges(
                self.qtl_df[self.qtl_df["q.value"] <= q_thr][["SNP", "gene"]].values,
                directed=True
            )

            real_avg_link_sim = []
            random_avg_link_sim = []

            for module_name, module_genes in self.module_dict.items():
                full_g = ig.Graph.Full(len(module_genes))
                full_g.vs["name"] = module_genes
                real_avg_link_sim.append(linkage_similarity(full_g, qtl_graph, mode='mean'))
                randsum = 0.
                for i in range(randiter_count):
                    full_g.vs["name"] = np.random.choice(gene_pool, len(module_genes), replace=False)
                    randsum += linkage_similarity(full_g, qtl_graph, mode='mean')
                random_avg_link_sim.append(randsum / randiter_count)

            real_results = real_avg_link_sim if real_results is None else np.vstack((real_results, real_avg_link_sim))
            randomized_results = random_avg_link_sim if randomized_results is None \
                else np.vstack((randomized_results, random_avg_link_sim))

        results_df_dict = collections.defaultdict()
        for i, module_name in enumerate(self.module_dict.keys()):
            results_df_dict[module_name] = pd.DataFrame(
                data=np.column_stack((
                    real_results[:, i],
                    randomized_results[:, i]
                )),
                columns=["real_mean", "simulated_mean"]
            )

        for module_name in results_df_dict.keys():
            plot_pairwise_test(
                module_name=module_name,
                module_type="thecellmap",
                linkage_sharing_df=results_df_dict[module_name],
                q_thresholds=np.logspace(-5, -2, 100),
                qtl_type="eQTLs"
            )

        return results_df_dict

    # how to incapsulate implementation of the test into something nice?
    def interactome_test(self, randiter_count=8):
        '''
        Interface function. Performs the simulation and saves the plot with results.
        return: Some measure of curve similarity between real and simulated data.
        '''

        pool = mp.Pool(mp.cpu_count())
        raw_results_dict = pool.map(self._process_randiter, range(randiter_count))
        pool.close()
        pool.join()

        results_dict = {
            module_name:np.stack([rawd[module_name] for rawd in raw_results_dict])
            for module_name in self.module_dict.keys()
        }

        for module_name in results_dict.keys():
            results_dict[module_name] = results_dict[module_name].mean(axis=0)

        module_scores = collections.defaultdict(list)
        for module_name, module_vertices in self.module_dict.items():
            results_df = pd.DataFrame(results_dict[module_name], columns=["p_value", "real_mean", "simulated_mean"])
            results_df = results_df[results_df["real_mean"] != 0]
            if results_df.empty:
                print("{} — FATAL ERROR!".format(module_name))
                continue

            trimmed_qval_list = self.qval_list[-results_df.shape[0]:]

            curve_similarity = distance.sqeuclidean(
                results_df["real_mean"], results_df["simulated_mean"],
                -np.log10(trimmed_qval_list)
            )

            self.statistics[module_name] = LinkageSharingStatistics(
                qtl_type=self.qtl_type,
                qtl_df=self.qtl_df[self.qtl_df["gene"].isin(module_vertices)],
                module_type=self.modules_type,
                module_name=module_name,
                module_graph=self.interactome.subgraph(
                    set(self.interactome.vs["name"])
                    & set(module_vertices)
                ),
                linkage_sharing_df=results_df,
                robustness_score=curve_similarity,
                q_thresholds=self.qval_list
            )

            self.statistics[module_name].save()
            self.statistics[module_name].plot_module_graph()
            self.statistics[module_name].plot_q_value_hist()
            self.statistics[module_name].plot_linkage_sharing_comparison()

            module_scores[module_name] = self.statistics[module_name].features_dict()
        return module_scores

    # [-------------------------------------------------PRIVATE METHODS-------------------------------------------------]

    def _process_threshold(self, q_value_threshold, module_name):
        qtl_graph = networks.graph_from_edges(
            edges=self.qtl_df[self.qtl_df["q.value"] <= q_value_threshold][["SNP", "gene"]].values,
            directed=True
        )
        genes_with_linkages = qtl_graph.vs.select(part=1)["name"]
        vertex_names = self.module_dict[module_name]

        if set(genes_with_linkages).isdisjoint(set(vertex_names)):
            return self._default_row

        real_stats = linkage_similarity(
            module_graph=self.interactome.subgraph(
                set(genes_with_linkages) &
                set(vertex_names)
            ),
            qtl_graph=qtl_graph,
            mode="full"
        )

        if len(real_stats) == 0 or all([sim_coeff == 0 for sim_coeff in real_stats]):
            return self._default_row

        randomized_stats = linkage_similarity(
            module_graph=self.randomized_interactome.subgraph(
                set(genes_with_linkages) &
                set(vertex_names)
            ),
            qtl_graph=qtl_graph,
            mode="full"
        )
        p_value = 0.
        if len(randomized_stats) != 0:
            try:
                _, p_value = stats.mannwhitneyu(real_stats, randomized_stats, alternative="two-sided")
            except ValueError:
                # prevents "all numbers are equal" error that is thrown when,
                # for example, [1, 1, 1] and [1, 1] are passed to MWU
                # it's really weird, but simply ignoring it seems the be best workaround so far
                pass
        return np.array([p_value, real_stats.mean(), randomized_stats.mean() if len(randomized_stats) != 0 else 0.])

    def _process_randiter(self, iter_num):
        results_dict = collections.defaultdict(list)
        self.randomized_interactome = ig.Graph().Read_Pickle(
            "./data/randomized_interactome_copies/{}.pkl".format(iter_num)
        )
        for module_name, module_vertices in self.module_dict.items():
            randiter_results = np.vstack(
                [self._process_threshold(q_thr, module_name)
                 for q_thr in self.qval_list]
            )
            if len(results_dict[module_name]) == 0:
                results_dict[module_name] = randiter_results
            else:
                results_dict[module_name] = np.dstack(
                    (results_dict[module_name], randiter_results)
                )

        return results_dict


########################################################################################################################
#                                             PREDICTING PQTLS FROM EQTLS                                              #
########################################################################################################################

"""
TODO: 
- Refactor interface!
- Reimplement plotting 
- Add more statistics to returned string
"""


class PqtlPredictor:
# [-------------------------------------------------PUBLIC METHODS-------------------------------------------------]
    def __init__(self, eqtls_df, pqtls_df,
                 eqtls_expression_df, eqtls_genotypes_df,
                 pqtls_expression_df, pqtls_genotypes_df,
                 full_genotypes_df,
                 functional_module_name,
                 functional_module_graph):

        self.eqtls_df = eqtls_df
        self.pqtls_df = pqtls_df
        self.eqtls_expr_df = eqtls_expression_df
        self.eqtls_gen_df = eqtls_genotypes_df
        self.pqtls_expr_df = pqtls_expression_df
        self.pqtls_gen_df = pqtls_genotypes_df

        self.valid_markers = set(self.pqtls_gen_df["SNP"].values)
        self.valid_genes = set(self.pqtls_expr_df["gene"].values)

        self.pqtls_expr_mx = \
            self.pqtls_expr_df.as_matrix(
                columns=self.pqtls_expr_df.columns[1:]
            )
        self.pqtls_gen_mx = \
            self.pqtls_gen_df.as_matrix(
                columns=self.pqtls_gen_df.columns[1:]
            )

        self.full_gen_df = full_genotypes_df

        self.mname = functional_module_name
        self.mgraph = functional_module_graph

    # TO BE REFACTORED
    def predict(self):
        pool = mp.Pool(mp.cpu_count())
        results = pool.map(self._get_possible_linkages, self.valid_genes)
        possible_linkages = list(set(sum(results, [])))
        pool.close()
        pool.join()

        pool = mp.Pool(mp.cpu_count())
        p_values = pool.map(self._compute_p_value, possible_linkages)
        pool.close()
        pool.join()

        pd.DataFrame(p_values, columns=["pvalue"]).to_csv("./data/pQTLs/pvalues.csv", sep='\t', index=False)
        subprocess.check_call(['Rscript', './src/p-values_to_q-values.R'], shell=False)
        q_values = pd.read_table("./data/pQTLs/qvalues.csv").values
        reject = q_values <= 0.05

        new_pQTLs_list = []

        for i in range(len(possible_linkages)):
            marker_name, gene_name = possible_linkages[i]
            new_pQTLs_list.append((marker_name, gene_name, p_values[i], q_values[i], reject[i]))

        # Original and new results comparison via plots
        new_pQTLs_df = pd.DataFrame(new_pQTLs_list, columns=["SNP", "gene", "pvalue", "qvalue", "reject"])
        new_pQTLs_df = new_pQTLs_df[new_pQTLs_df["reject"] == True]
        new_pQTLs_df.to_csv("./data/pQTLs/new_results.csv", sep='\t', index=False, na_rep='NA')

        old_pQTL_x, old_pQTL_y = linkages2gencoords(self.pqtls_df, self.full_gen_df)
        new_pQTL_x, new_pQTL_y = linkages2gencoords(new_pQTLs_df, self.full_gen_df)

        plt.figure(figsize=(30, 10))
        plt.plot(old_pQTL_y, label="pQTLs found by brute force")
        plt.plot(new_pQTL_y, label="pQTLs predicted from interactions")
        plt.legend(fontsize=25)
        plt.xticks(fontsize=20)
        plt.yticks(fontsize=20)
        plt.xlabel("Genomic coordinates", fontsize=25)
        plt.ylabel("Number of linkages", fontsize=25)
        # plt.plot([-y for y in eQTL_y])
        plt.savefig("./img/linkages/pQTLs_old_and_new.png", format="png", dpi=300)

        new_pQTLs_df = pd.read_table("./data/pQTLs/new_results.csv", sep='\t')
        new_pQTLs_df = new_pQTLs_df[new_pQTLs_df["reject"] == "[ True]"]

        # approach quality evaluation:
        num_linked_to_gene = {gene_name: 0 for gene_name in self.pqtls_df["gene"].values}
        num_linked_to_marker = {marker_name: 0 for marker_name in self.pqtls_df["SNP"].values}
        for marker_name, gene_name in self.pqtls_df[["SNP", "gene"]].values:
            num_linked_to_marker[marker_name] += 1
            num_linked_to_gene[gene_name] += 1

        old_pQTL_linkage_pairs = set(
            (marker_name, gene_name)
            for marker_name, gene_name in self.pqtls_df[["SNP", "gene"]].values
        )

        common = 0
        for marker_name, gene_name in new_pQTLs_df[["SNP", "gene"]].values:
            linkage_pair = (marker_name, gene_name)
            common += linkage_pair in old_pQTL_linkage_pairs

        return "Common linkages: {}, {}%\nOld linkages, total: {}\nNew linkages, total: {}\nNew linkages found: {}". \
            format(common, 100 * common / self.pqtls_df.shape[0],
                   self.pqtls_df.shape[0],
                   new_pQTLs_df.shape[0],
                   new_pQTLs_df.shape[0] - common)

# [--------------------------------------------------PRIVATE METHODS--------------------------------------------------]

    def _get_interacting_genes(self, gene_name, bfs_depth=1):
        try:
            return set(self.mgraph.vs[
                           self.mgraph.neighborhood(gene_name, order=bfs_depth)
                       ]["name"])
        except ValueError:
            return set()

    def _get_linked_markers(self, gene_name, QTL_df):
        return set(QTL_df[QTL_df["gene"] == gene_name]["SNP"].values)

    def _get_linked_genes(self, marker_name, QTL_df):
        return set(QTL_df[QTL_df["SNP"] == marker_name]["gene"].values)

    def _get_linked_to_adjacent(self, gene_name):
        interacting_genes = set(gene_name) | self._get_interacting_genes(gene_name)
        linked = set()
        for neighbor in interacting_genes:
            linked |= self._get_linked_markers(neighbor, self.eqtls_df)
            if neighbor != gene_name:
                linked |= set(neighbor)
        return linked

    def _get_possible_linkages(self, gene_name):
        return [(marker_name, gene_name) for marker_name in
                self._get_linked_to_adjacent(gene_name) & self.valid_markers]

    def _compute_p_value(self, candidate_linkage):
        marker_name, gene_name = candidate_linkage
        genotype_rowmask = self.pqtls_gen_df["SNP"] == marker_name
        genotype_row = self.pqtls_gen_mx[genotype_rowmask]
        expression_rowmask = self.pqtls_expr_df["gene"] == gene_name
        expression_row = self.pqtls_expr_mx[expression_rowmask]
        from_BY = expression_row[genotype_row == 0]
        from_RM = expression_row[genotype_row == 2]
        _, p_value = stats.mannwhitneyu(from_BY, from_RM, alternative="two-sided")
        return p_value
