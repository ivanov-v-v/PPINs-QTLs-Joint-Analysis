from __future__ import print_function

import collections
import multiprocessing as mp
import pickle
import subprocess

import igraph as ig
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.stats as stats
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


def plot_q_hist(destdir, module_graph, module_name, module_type, qtl_df, qtl_type):
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

    plt.savefig(destdir + "/q_hist.png")
    plt.close()


def plot_pairwise_test(test_results_df, module_name, module_type, q_thresholds, qtl_type):
    # for the plot to be informative, don't show the range
    # where there are no linkage with such q-values at all
    fig, ax = plt.subplots(figsize=(20, 10))
    ax.set_title("{} pairwise test for {}; {}".format(
        qtl_type, module_name, module_type
    ), fontsize=20)
    plt.xscale('log')

    ax.set_xlabel("linkage q-value threshold, lg-scale", fontsize=20)
    ax.set_ylabel("average linkage similarity", fontsize=20)
    ax.plot(q_thresholds, test_results_df["real_mean"], label="functional module")
    ax.plot(q_thresholds, test_results_df["simulated_mean"],
            label="random gene set\n")
    ax.legend(fontsize=15)

    ax.grid(linestyle='dotted')
    plt.setp(ax.get_xticklabels(), fontsize=15)
    plt.setp(ax.get_yticklabels(), fontsize=12)

    destdir = '/'.join(["./results", module_type, module_name, qtl_type])
    util.ensure_dir(destdir)
    plt.savefig(destdir + "/pairwise_test.png")
    plt.close()


def plot_ppi_test(module_name, module_type, test_results_df,
                  robustness_score, q_thresholds, qtl_type):
    # for the plot to be informative, don't show the range
    # where there are no linkage with such q-values at all

    fig, ax1 = plt.subplots(figsize=(20, 10))
    ax1.set_title("{} PPIN test for {}; {}".format(
        qtl_type, module_name, module_type
    ), fontsize=20)
    plt.xscale('log')
    ax2 = ax1.twinx()

    ax1.set_xlabel("linkage q-value threshold, lg-scale", fontsize=20)
    ax1.set_ylabel("average linkage similarity", fontsize=20)
    ax1.plot(q_thresholds, test_results_df["real_mean"], label="real data")
    ax1.plot(q_thresholds, test_results_df["simulated_mean"],
             label="synthesized data\n"
                   "curve similarity score (WMSQ): {0:.4f}".format(robustness_score)
             )

    ax2.plot(
        q_thresholds, test_results_df["p_value"],
        color="r", linestyle="dashed", alpha=0.5,
        label="p-values\nmean p-value: {0:.4f}\nmedian p-value: {1:.4f}".format(
            test_results_df["p_value"].mean(),
            test_results_df["p_value"].median()
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
    plt.savefig(destdir + "/ppi_test.png")
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
    :param module_graph: simple graph of some functional module
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
    """Storage class for analysis results produced by the LinkageSharingAnalyzer class"""

    def __init__(self, analysis_type, module_name, module_type, qtl_type, qtl_df=None,
                 q_thresholds=None, test_results_df=None, module_scores=None):
        self.analysis_type = analysis_type

        self.module_name = module_name
        self.module_type = module_type

        self.qtl_df = qtl_df
        self.qtl_type = qtl_type

        self.test_results_df = test_results_df
        self.module_scores = module_scores
        self.q_thresholds = q_thresholds

        self.destdir = "/".join(["./results", self.module_type, self.module_name, self.qtl_type])
        util.ensure_dir(self.destdir)

    def save(self):
        with open(self.destdir + "/" + self.analysis_type + "_test.pkl", 'wb') as pickle_file:
            pickle.dump(obj=self.__dict__, file=pickle_file)

    def load(self):
        with open(self.destdir + "/" + self.analysis_type + "_test.pkl", 'rb') as pickle_file:
            self.__dict__.update(pickle.load(file=pickle_file))


def plot_vignettes(module_type, module_names):
    for module_name in module_names:
        print(module_name)
        plots = []
        for qtl_type in ["eQTLs", "pQTLs"]:
            path = "/".join(["./results", module_type, module_name, qtl_type])
            q_hist = util.try_imread(path + "/q_hist.png", mode="RGBA")
            pairwise_test = util.try_imread(path + "/pairwise_test.png", mode="RGBA")
            interactome_test = util.try_imread(path + "/ppi_test.png", mode="RGBA")
            network = util.try_imread(path + "/module_graph.png", mode="RGBA")

            vignette_side_parts = [img for img in [q_hist, pairwise_test, interactome_test, network] if img is not None]
            if len(vignette_side_parts) != 0:
                plots.append(np.vstack(vignette_side_parts))

        combined_plot = np.column_stack([img for img in plots if img is not None])
        plt.imsave("/".join(["./results", module_type, module_name + ".png"]), combined_plot)


class LinkageSharingAnalyzer:
    """
    Used to analyze to which extent interacting genes tend to share linkages.
    Given a datafame of linkages, a functional module graph and a list of q-value thresholds,
    this class perturbs the graph ~100 times and each time calculates averaged Jaccard coefficient
    for each threshold. Then, to ensure robustness of conclusions, simulation results are averaged and
    plotted altogether with 95% confidence intervals against actual data for the original graph.
    """

    # [-------------------------------------------------PUBLIC METHODS-------------------------------------------------]
    def __init__(self, expression_df, interactome_graph, modules_type, module_dict,
                 qtl_type, qtl_df, q_thresholds, pairwise_test_iter=8, ppi_test_iter=8):
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

        self.randomized_interactome = None
        self._default_row = np.zeros(3)

        self.pairwise_randiter_count = pairwise_test_iter
        self.ppi_randiter_count = ppi_test_iter

    def process_modules(self):
        for module_name in self.module_dict.keys():
            destdir = "/".join(["./results", self.modules_type, module_name, self.qtl_type])
            util.ensure_dir(destdir)

            module_graph = self.interactome.subgraph(
                set(self.interactome.vs["name"])
                & set(self.module_dict[module_name])
            )

            plot_module_graph(
                destdir=destdir,
                module_graph=module_graph,
                qtl_df=self.qtl_df[self.qtl_df["gene"].isin(self.module_dict[module_name])]
            )
            plot_q_hist(
                destdir=destdir,
                module_name=module_name,
                module_graph=module_graph,
                module_type=self.modules_type,
                qtl_df=self.qtl_df[self.qtl_df["gene"].isin(self.module_dict[module_name])],
                qtl_type=self.qtl_type
            )

        self.pairwise_test()
        self.ppi_test()

    def pairwise_test(self):
        # Почему бы не сравнивать их просто со средним значением для всего графа?

        pool = mp.Pool(mp.cpu_count())
        real_means, randomized_means = np.stack(
            pool.map(self._pairwise_process_threshold, self.qval_list),
            axis=1
        )
        pool.close()
        pool.join()

        results_df_dict = collections.defaultdict()
        for i, module_name in enumerate(self.module_dict.keys()):
            results_df_dict[module_name] = pd.DataFrame(
                data=np.column_stack((
                    real_means[:, i],
                    randomized_means[:, i]
                )),
                columns=["real_mean", "simulated_mean"]
            )

        module_scores = collections.defaultdict()
        for module_name in self.module_dict.keys():
            pairwise_test_results_df = results_df_dict[module_name]
            plot_pairwise_test(
                module_name=module_name,
                module_type=self.modules_type,
                test_results_df=pairwise_test_results_df,
                q_thresholds=self.qval_list,
                qtl_type=self.qtl_type
            )

            curve_similarity = distance.sqeuclidean(
                results_df_dict[module_name]["real_mean"],
                results_df_dict[module_name]["simulated_mean"],
                -np.log10(self.qval_list)
            )

            module_scores[module_name] = collections.OrderedDict([
                ("robustness_score", curve_similarity),
                ("mean_q_score", self.qtl_df[self.qtl_df["gene"].isin(
                    self.module_dict[module_name])]["q.value"].mean()),
                ("median_q_score", self.qtl_df[self.qtl_df["gene"].isin(
                    self.module_dict[module_name])]["q.value"].median()),
            ])

            LinkageSharingStatistics(
                analysis_type="pairwise",
                qtl_type=self.qtl_type,
                qtl_df=self.qtl_df[self.qtl_df["gene"].isin(self.module_dict[module_name])],
                module_type=self.modules_type,
                module_name=module_name,
                module_scores=module_scores[module_name],
                test_results_df=pairwise_test_results_df,
                q_thresholds=self.qval_list
            ).save()

        pd.DataFrame(module_scores).T.to_csv("/".join(
            ["./results/", self.modules_type, "pairwise_test_" + self.qtl_type + "_scores.csv"]),
            sep='\t'
        )
        return module_scores

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
            for module_name in self.module_dict.keys()
        }

        for module_name in results_dict.keys():
            results_dict[module_name] = results_dict[module_name].mean(axis=0)

        module_scores = collections.defaultdict(list)
        for module_name, module_vertices in self.module_dict.items():
            ppi_test_results_df = pd.DataFrame(results_dict[module_name],
                                               columns=["p_value", "real_mean", "simulated_mean"])

            curve_similarity = distance.sqeuclidean(
                ppi_test_results_df["real_mean"], ppi_test_results_df["simulated_mean"],
                -np.log10(self.qval_list)
            )

            plot_ppi_test(
                module_name=module_name,
                module_type=self.modules_type,
                test_results_df=ppi_test_results_df,
                robustness_score=curve_similarity,
                q_thresholds=self.qval_list,
                qtl_type=self.qtl_type
            )

            module_scores[module_name] = collections.OrderedDict([
                ("robustness_score", curve_similarity),
                ("mean_q_score", self.qtl_df[self.qtl_df["gene"].isin(
                    self.module_dict[module_name])]["q.value"].mean()),
                ("median_q_score", self.qtl_df[self.qtl_df["gene"].isin(
                    self.module_dict[module_name])]["q.value"].median()),
                ("mean_p_score", ppi_test_results_df["p_value"].mean()),
                ("median_p_score", ppi_test_results_df["p_value"].median())])

            LinkageSharingStatistics(
                analysis_type="ppi",
                qtl_type=self.qtl_type,
                qtl_df=self.qtl_df[self.qtl_df["gene"].isin(module_vertices)],
                module_type=self.modules_type,
                module_name=module_name,
                test_results_df=ppi_test_results_df,
                module_scores=module_scores[module_name],
                q_thresholds=self.qval_list
            ).save()

        pd.DataFrame(module_scores).T.to_csv("/".join(
            ["./results/", self.modules_type, "ppi_test_" + self.qtl_type + "_scores.csv"]),
            sep='\t'
        )
        return module_scores

    # [------------------------------------------------PRIVATE METHODS------------------------------------------------]

    def _pairwise_process_threshold(self, q_value_threshold):
        random_state = np.random.RandomState()
        gene_pool = self.expression_df["gene"].values
        qtl_graph = networks.graph_from_edges(
            self.qtl_df[self.qtl_df["q.value"] <= q_value_threshold][["SNP", "gene"]].values,
            directed=True
        )

        real_avg_link_sim = []
        random_avg_link_sim = []
        for module_name, module_genes in self.module_dict.items():
            full_g = ig.Graph.Full(len(module_genes))
            full_g.vs["name"] = module_genes
            real_avg_link_sim.append(linkage_similarity(full_g, qtl_graph, mode='mean'))
            randsum = 0.
            for i in range(self.pairwise_randiter_count):
                full_g.vs["name"] = random_state.choice(gene_pool, len(module_genes), replace=False)
                randsum += linkage_similarity(full_g, qtl_graph, mode='mean')
            random_avg_link_sim.append(randsum / self.pairwise_randiter_count)

        return real_avg_link_sim, random_avg_link_sim

    def _ppi_process_threshold(self, q_value_threshold, module_name):
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
                set(self.interactome.vs["name"]) &
                set(vertex_names)
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
                _, p_value = stats.mannwhitneyu(real_stats, randomized_stats, alternative="two-sided")
            except ValueError:
                # prevents "all numbers are equal" error that is thrown when,
                # for example, [1, 1, 1] and [1, 1] are passed to MWU
                # it's really weird, but simply ignoring it seems the be best workaround so far
                pass
        return np.array([p_value, real_stats.mean(), randomized_stats.mean() if len(randomized_stats) != 0 else 0.])

    def _ppi_process_randiter(self, iter_num):
        results_dict = collections.defaultdict(list)
        self.randomized_interactome = ig.Graph().Read_Pickle(
            "./data/randomized_interactome_copies/{}.pkl".format(iter_num)
        )
        for module_name, module_vertices in self.module_dict.items():
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


def process_ontologies(database_name, expression_dfs, qtl_dfs, qtl_types,
                       interactome_graph, ontology_to_genes, q_thresholds,
                       pairwise_test_iter=32, ppi_test_iter=1024):
    qtl_graphs = {
        qtl_type: networks.graph_from_edges(
            edges=qtl_df[qtl_df["q.value"] <= np.max(q_thresholds)][["SNP", "gene"]].values,
            directed=True
        )
        for qtl_type, qtl_df in list(zip(qtl_types, qtl_dfs))
    }

    for qtl_type, qtl_df, expression_df in list(zip(qtl_types, qtl_dfs, expression_dfs)):
        filtered_out_modules = []
        modules_dict = collections.defaultdict()
        for module_name, gene_list in ontology_to_genes.items():
            module_graph = interactome_graph.subgraph(
                set(ontology_to_genes[module_name])
                & set(interactome_graph.vs["name"])
            ).simplify()
            module_graph.vs.select(_degree=0).delete()

            if module_graph.vcount() < 4:
                filtered_out_modules.append((module_name, "less than 4 vertices"))
                continue

            genes_with_linkages = qtl_graphs[qtl_type].vs.select(part=1)["name"]

            if set(genes_with_linkages).isdisjoint(set(gene_list)):
                filtered_out_modules.append((module_name, "no linkages"))
                continue

            average_linkage_similarity = linkage_similarity(
                module_graph=interactome_graph.subgraph(
                    set(interactome_graph.vs["name"]) &
                    set(gene_list)
                ),
                qtl_graph=qtl_graphs[qtl_type],
                mode="mean"
            )
            if np.isclose(average_linkage_similarity, 0):
                filtered_out_modules.append((module_name, "linkage sets don't intersect"))
                continue

            modules_dict[module_name] = gene_list

        pd.DataFrame(filtered_out_modules, columns=["module_name", "reason"]).to_csv(
            "/".join(["./results", database_name, qtl_type + "_filtered_out_modules.csv"]), sep='\t', index=False
        )

        LinkageSharingAnalyzer(
            expression_df=expression_df,
            interactome_graph=interactome_graph,
            qtl_type=qtl_type,
            qtl_df=qtl_df,
            modules_type=database_name,
            module_dict=modules_dict,
            q_thresholds=q_thresholds,
            pairwise_test_iter=pairwise_test_iter,
            ppi_test_iter=ppi_test_iter
        ).process_modules()

        plot_vignettes(module_type=database_name, module_names=modules_dict.keys())


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
                 functional_module_name, functional_module_graph):

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
