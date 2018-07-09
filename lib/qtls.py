from __future__ import print_function

import collections
import multiprocessing as mp
import subprocess
from importlib import reload

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


def plot_module_svg(module_name, module_graph, qtl_type, qtl_df, dirname):
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

    edge_colors = [util.rgba2hex(rgba) for rgba in [[0, 0, 0, 0.05]] * module_graph.ecount()]
    module_graph.es["color"] = edge_colors

    ''' WRITE_SVG IS DEPRECATED, USE IG.PLOT INSTEAD '''
    module_graph.write_svg(
        fname=dirname + '/' + qtl_type + "_" + module_name + ".svg",
        layout=module_graph.layout_fruchterman_reingold(),
        width=1000,
        height=1000,
        labels=[gene_name + ", " + str(num_linkages[i]) if num_linkages[i] != 0 else ''
                for i, gene_name in enumerate(module_graph.vs["name"])],
        colors="color",
        shapes=[1] * module_graph.vcount(),
        vertex_size=8,
        font_size=25,
        edge_colors="color"
    )


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
    for i, edge in enumerate(module_graph.es):
        source = module_graph.vs[edge.source]
        target = module_graph.vs[edge.target]

        try:
            s_neigh = set(qtl_graph.neighbors(source["name"], mode="IN"))
            t_neigh = set(qtl_graph.neighbors(target["name"], mode="IN"))
            results[i] = jaccard(s_neigh, t_neigh)
        except ValueError:
            '''TODO: add description'''
            pass

    return results.mean() if mode == "mean" else results


class LinkageSharingStatistics:
    ''' Storage class for analysis results produced by the LinkageSharingAnalyzer class.'''
    def __init__(self, qtl_type, module_type, module_name, qtl_df=None, module_graph=None,
                 q_thresholds=None, linkage_sharing_df=None, robustness_score=None):
        self.qtl_type = qtl_type
        self.module_type = module_type
        self.module_name = module_name

        self.qtl_df = qtl_df
        self.module_graph = module_graph
        self.q_thresholds = q_thresholds
        self.linkage_sharing_df = linkage_sharing_df
        self.robustness_score = robustness_score

        self.destdir ="/".join(["./results", self.module_type, self.module_name, self.qtl_type])
        util.ensure_dir(self.destdir)

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
        self.module_graph.Read_Pickle(fname="/".join([self.destdir, "module_graph.pkl"]))
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
            ("vertex_q_score", self.qtl_df["q.value"].median()),
            ("edge_p_score", self.linkage_sharing_df["p_value"].median())])

    def plot_module_graph(self):
        plot_module_svg(
            module_name=self.module_name,
            module_graph=self.module_graph,
            qtl_type=self.qtl_type,
            qtl_df=self.qtl_df,
            dirname=self.destdir
        )

    def plot_q_value_hist(self):
        fig, ax = plt.subplots(figsize=(20, 10))
        counts, bins, patches = plt.hist(
            self.qtl_df["q.value"],
            bins='auto',
            alpha=0.6,
            color="xkcd:tangerine",
            edgecolor="black",
            linewidth=1.2,
            label="{} linkages in total\nmean q-value: {}\n|V| = {}\n|E| = {}"
                .format(self.qtl_df.shape[0], self.qtl_df["q.value"].mean(),
                        self.module_graph.vcount(), self.module_graph.ecount())
        )
        ax.grid(linestyle='dotted', alpha=0.8)
        ax.tick_params(labelsize=15)

        ax.set_xticks(bins)
        ax.tick_params(axis='x', labelrotation=90)
        ax.ticklabel_format(axis='x', style='sci', scilimits=(0, 0), useMathText=True)

        ax.set_title("{} linkage q-value distribution for {}; {}".format(
            self.qtl_type, self.module_name, self.module_type
        ), fontsize=20)
        ax.set_xlabel("q-value", fontsize=20)
        ax.set_ylabel("number of linkages", fontsize=20)
        ax.legend(fontsize=20)

        plt.savefig(self.destdir + "/{}_q_value_distribution for {}; {}.png".format(
            self.qtl_type, self.module_name, self.module_type)
        )
        plt.close()

    def plot_linkage_sharing_comparison(self):
        # for the plot to be informative, don't show the range
        # where there are no linkage with such q-values at all
        results_df = self.linkage_sharing_df[
            self.linkage_sharing_df["real_mean"] != 0
        ]
        trimmed_qval_list = self.q_thresholds[-results_df.shape[0]:]

        fig, ax1 = plt.subplots(figsize=(20, 10))
        ax1.set_title("{} linkage similarity for {}; {}".format(
            self.qtl_type, self.module_name, self.module_type
        ), fontsize=20)
        plt.xscale('log')
        ax2 = ax1.twinx()

        ax1.set_xlabel("linkage q-value threshold, log10-scale", fontsize=20)
        ax1.set_ylabel("average linkage similarity", fontsize=20)
        ax1.plot(trimmed_qval_list, results_df["real_mean"], label="real data")
        ax1.plot(trimmed_qval_list, results_df["simulated_mean"], label="synthesized data")

        ax2.plot(trimmed_qval_list, results_df["p_value"], color="r", linestyle="dashed", label="p-values", alpha=0.5)
        ax2.set_ylabel("p-value of difference of edge scores", fontsize=20)

        lines, labels = ax1.get_legend_handles_labels()
        lines2, labels2 = ax2.get_legend_handles_labels()
        ax2.legend(lines + lines2, labels + labels2, loc=0, fontsize=15)

        ax1.grid(linestyle='dotted')
        plt.setp(ax1.get_xticklabels(), fontsize=15)
        plt.setp(ax1.get_yticklabels(), fontsize=12)
        plt.setp(ax2.get_yticklabels(), fontsize=12)

        plt.savefig(self.destdir + "/{}_{}_{}.png".format(
            self.qtl_type, self.module_name, self.module_type)
        )
        plt.close()


class LinkageSharingAnalyzer:
    """
    Used to analyze to which extent interacting genes tend to share linkages.
    Given a datafame of linkages, a functional module graph and a list of q-value thresholds,
    this class perturbs the graph ~100 times and each time calculates averaged Jaccard coefficient
    for each threshold. Then, to ensure robustness of conclusions, simulation results are averaged and
    plotted altogether with 95% confidence intervals against actual data for the original graph.
    """

# [-------------------------------------------------PUBLIC METHODS-------------------------------------------------]
    def __init__(self, qtl_type, qtl_df,
                 module_type, module_name, module_graph,
                 q_thresholds):
        reload(util)
        self.qtl_type = qtl_type
        self.qtl_df = qtl_df

        self.mtype = module_type
        self.mname = module_name
        # Without this preprocessing the results will be __incorrect__,
        # inflated, as in linkage similarity calculations the same jaccard
        # coefficient will be calculated many times
        self.mgraph = module_graph.simplify()
        self.mgraph.vs.select(_degree=0).delete()

        self.qval_list = q_thresholds
        self.statistics = None

    def process_threshold(self, q_value_threshold, randiter=64):
        qtl_graph = networks.graph_from_edges(
            edges=self.qtl_df[self.qtl_df["q.value"] <= q_value_threshold][["SNP", "gene"]].values,
            directed=True
        )
        genes_with_linkages = qtl_graph.vs.select(part=1)["name"]
        vertex_names = self.mgraph.vs["name"]

        if set(genes_with_linkages).isdisjoint(set(vertex_names)):
            return np.zeros(4)

        real_stats = linkage_similarity(
            module_graph=self.mgraph.subgraph(
                set(genes_with_linkages) &
                set(vertex_names)
            ),
            qtl_graph=qtl_graph,
            mode="full"
        )

        if len(real_stats) == 0 or all([sim_coeff == 0 for sim_coeff in real_stats]):
            return np.zeros(4)

        simulation_means = np.zeros(randiter)
        simulation_stds = np.zeros(randiter)
        p_values = np.ones(randiter)

        for i in range(randiter):
            try:
                randomized_mgraph = self.mgraph.Degree_Sequence(
                    self.mgraph.degree(),
                    method="vl"
                )
            except RuntimeWarning: # create random edges in case there is only one graph with given degree sequence
                randedges = [
                    (self.mgraph.vs[i], self.mgraph.vs[j])
                    for i, j in util.sample_combinations(
                        (self.mgraph.vcount(), self.mgraph.vcount()),
                        self.mgraph.ecount()
                    )
                ]
                randomized_mgraph = networks.graph_from_edges(randedges)

            randomized_mgraph.vs["name"] = vertex_names
            randomized_stats = \
                linkage_similarity(
                    module_graph=randomized_mgraph.subgraph(
                        set(genes_with_linkages) &
                        set(vertex_names)
                    ),
                    qtl_graph=qtl_graph,
                    mode="full"
                )
            if len(randomized_stats) != 0:
                try:
                    _, p_values[i] = stats.mannwhitneyu(real_stats, randomized_stats, alternative="two-sided")
                except ValueError:
                    # prevents "all numbers are equal" error that is thrown when,
                    # for example, [1, 1, 1] and [1, 1] are passed to MWU
                    # it's really weird, but simply ignoring it seems the be best workaround so far
                    pass
                simulation_means[i] = randomized_stats.mean()
                simulation_stds[i] = randomized_stats.std()
        return np.array([p_values.mean(), real_stats.mean(), simulation_means.mean(), simulation_stds.mean()])

    def process_module(self, recompute=True):
        '''
        Interface function. Performs the simulation and saves the plot with results.
        return: Some measure of curve similarity between real and simulated data.
        '''

        if not recompute:
            if self.statistics is None:
                self.statistics = LinkageSharingStatistics(
                    qtl_type=self.qtl_type,
                    module_type=self.mtype,
                    module_name=self.mname
                )
            self.statistics.load()
            return self.statistics

        pool = mp.Pool(mp.cpu_count())
        results = np.vstack(pool.map(self.process_threshold, self.qval_list))
        pool.close()
        pool.join()

        results_df = pd.DataFrame(np.vstack(results),
                                  columns=["p_value", "real_mean", "simulated_mean", "simulated_std"])

        results_df = results_df[results_df["real_mean"] != 0]
        # filter those out __before__ calling this function
        # raise an error here
        if results_df.empty:
            print("{} — no data, no processing done!".format(self.mname))
            return 0

        trimmed_qval_list = self.qval_list[-results_df.shape[0]:]

        curve_similarity = distance.sqeuclidean(
            results_df["real_mean"], results_df["simulated_mean"],
            -np.log10(trimmed_qval_list)
        )

        self.results = LinkageSharingStatistics(
            qtl_type=self.qtl_type,
            qtl_df=self.qtl_df[self.qtl_df["gene"].isin(self.mgraph.vs["name"])],
            module_type=self.mtype,
            module_name=self.mname,
            module_graph=self.mgraph,
            linkage_sharing_df=results_df,
            robustness_score=curve_similarity,
            q_thresholds=self.qval_list
        )

        self.results.save()
        self.results.plot_module_graph()
        self.results.plot_q_value_hist()
        self.results.plot_linkage_sharing_comparison()

        return self.results


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
