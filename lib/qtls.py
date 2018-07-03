from __future__ import print_function

import subprocess
import sys
from multiprocessing import Pool, cpu_count

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.stats as stats

import networks


def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

# Given a dataframe of estimated QTL linkages and
# full strain genotype lookup table, plot the number
# of linkages against marker location in genome
def map_linkages_to_genome_location(QTL_df, full_genotypes_df):
    QTL_graph = networks.graph_from_edges(
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

''' TODO:   Стоит попробовать более совершенные метрики
            подобия, а также добиться лучшей скорости работы 
'''

def jaccard(set_A, set_B):
    '''
    :return: Jaccard coefficient for sets A and B
    '''
    return len(set_A & set_B) / len(set_A | set_B)


# Given a pair of graphs, representing gene-gene interactions
# and estimated QTL-linkages, calculate for each pair of interacting genes
# a Jaccard similarity coefficient, and then average it over all edges
def mean_linkage_similarity(interaction_graph, QTL_graph):
    genes_with_linkages = QTL_graph.vs.select(part=1)["name"]
    genes_with_interactions = interaction_graph.vs["name"]

    # Как поступать с петлями? Их учитывать или нет?
    subgraph_with_linkages = interaction_graph.subgraph(
        set(genes_with_linkages) & set(genes_with_interactions)
    )

    if not subgraph_with_linkages.ecount():
        return 0  # такой граф нет смысла проверять

    # Перебрать все рёбра и сопоставить каждому пару множеств —
    # eQTLs, которые линкуются с концами ребра —
    # а затем подсчитать для них долю общих элементов

    mean_jaccard = 0.
    for edge in subgraph_with_linkages.es:
        source = subgraph_with_linkages.vs[edge.source]
        target = subgraph_with_linkages.vs[edge.target]

        s_neigh = set(QTL_graph.neighbors(source["name"], mode="IN"))
        t_neigh = set(QTL_graph.neighbors(target["name"], mode="IN"))
        mean_jaccard += jaccard(s_neigh, t_neigh)

    return mean_jaccard / subgraph_with_linkages.ecount()


class LinkageSharingAnalyzer:
    '''
    Used to analyze to which extent interacting genes tend to share linkages.
    Given a datafame of linkages, a functional module graph and a list of q-value thresholds,
    this class perturbs the graph ~100 times and each time calculates averaged Jaccard coefficient
    for each threshold. Then, to ensure robustness of conclusions, simulation results are averaged and
    plotted altogether with 95% confidence intervals against actual data for the original graph.
    '''
# [--------------------------------------------------PUBLIC METHODS--------------------------------------------------]
    def __init__(self, qtl_type, qtl_df,
                 module_name, module_graph,
                 q_value_thresholds):
        self.qtl_type = qtl_type
        self.qtl_df = qtl_df
        self.mname = module_name
        self.mgraph = module_graph
        self.qval_list = q_value_thresholds

    def analyze_robustness_of_linkage_sharing(self, title, destination_folder):
        '''
        Interface function. Performs the simulation and saves the plot with results.
        '''
        actual_stats = self.actual_linkage_statistics()
        simulated_stats = self.simulated_linkage_statistics()
        self.plot_robustness_analysis(title, destination_folder,
                                      actual_stats, simulated_stats)

    def plot_robustness_analysis(self, title, destination_folder, actual_stats, simulation_stats):
        '''
        Plots simulation data (average among all the iterations and confidence intervals)
        against actual data obtained from the initial graph.

        :return: Saves the image self.qtl_type + _ + self.mname + ".png" in destination_folder
        '''
        plt.figure(figsize=(20, 10))
        plt.title(title, fontsize=25)
        plt.xscale('log')

        plt.plot(self.qval_list, actual_stats, label="edges")
        # Display also the confidence interval to estimate the statistical strength of conclusions
        plt.plot(self.qval_list, simulation_stats.mean(axis=0), label="random pairs")
        confidence_interval = \
            stats.t.interval(
                0.95,
                len(self.qval_list) - 1,
                loc=np.mean(simulation_stats, axis=0),
                scale=stats.sem(simulation_stats, axis=0)
            )
        plt.fill_between(self.qval_list,
                         confidence_interval[0], confidence_interval[1],
                         color="#FFD700", alpha=0.4)

        plt.xticks(fontsize=15)
        plt.yticks(fontsize=15)
        plt.grid(linestyle="dotted")

        plt.xlabel("linkage q-value threshold", fontsize=20)
        plt.ylabel("average linkage similarity", fontsize=20)
        plt.legend(fontsize=20)

        plt.savefig(destination_folder + self.qtl_type + '_' + self.mname + ".png", dpi=300)
        plt.close()

    def actual_linkage_statistics(self):
        return self._range_mean_jaccard(self.mgraph)

    def simulated_linkage_statistics(self, randomization_iterations=64):
        '''
        Runs self._randomize_and_recalc_stats in parallel using Pool.

        :param randomization_iterations: How many times to run the simulation.
        :return: Numpy matrix of simulation results (each row corresponds to one run)
        '''
        pool = Pool(processes=cpu_count())
        simulation_stats = np.vstack(
            pool.map(
                self._randomize_and_recalc_stats,
                range(randomization_iterations)
            ))
        pool.close()
        pool.join()
        return simulation_stats

# [--------------------------------------------------PRIVATE METHODS--------------------------------------------------]

    def _range_mean_jaccard(self, modified_mgraph):
        '''
        Given a dataframe of QTL linkages and a list of q-value thresholds, loops through them,
        filtering out the linkages with q-value above the threshold, and calculates
        the averaged linkage similarity (Jaccard coefficient) for QTL sets linked to interacting genes

        :param modified_mgraph: In this context "modified" means perturbed self.mgraph
                                for robustness testing preserving vertex degree
        :return: Numpy array of averaged Jaccard coefficients for every threshold given.

        '''
        num_thresholds = len(self.qval_list)
        jaccard = np.zeros(num_thresholds, dtype=np.float32)
        qtl_df = self.qtl_df
        for idx, q_threshold in enumerate(self.qval_list[::-1]):
            qtl_df = qtl_df[qtl_df['q.value'] <= q_threshold ]
            qtl_graph = networks.graph_from_edges(
                edges=qtl_df[["SNP", "gene"]].values,
                directed=True
            )
            jaccard[num_thresholds - idx - 1] = \
                mean_linkage_similarity(modified_mgraph, qtl_graph)
        return jaccard

    def _randomize_and_recalc_stats(self, __dummy_arg):
        '''
        Performs the simulation: perturbs the graph and calculates Jaccard coefficients.

        :param __dummy_arg: Dummy variable required for Pool usage.
        :return: Numpy array of Jaccard coefficients for each threshold from self.qval_list.
        '''
        try:
            randomized_mgraph = self.mgraph.Degree_Sequence(
                self.mgraph.degree(),
                method="vl"
            )
        except RuntimeWarning:
            randomized_mgraph = self.mgraph
            pass
        randomized_mgraph.vs["name"] = self.mgraph.vs["name"]
        return self._range_mean_jaccard(randomized_mgraph)

'''
TODO: 
- Refactor interface!
- Reimplement plotting 
- Add more statistics to returned string
'''
class PqtlPredictor:
# [--------------------------------------------------PUBLIC METHODS--------------------------------------------------]
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
        pool = Pool(cpu_count())
        results = pool.map(self._get_possible_linkages, self.valid_genes)
        possible_linkages = list(set(sum(results, [])))
        pool.close()
        pool.join()

        pool = Pool(cpu_count())
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

        old_pQTL_x, old_pQTL_y = map_linkages_to_genome_location(self.pqtls_df, self.full_gen_df)
        new_pQTL_x, new_pQTL_y = map_linkages_to_genome_location(new_pQTLs_df, self.full_gen_df)

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


