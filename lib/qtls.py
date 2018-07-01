from multiprocessing import Pool, cpu_count

import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as stats

import networks


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
    QTL_graph = networks.assemble_graph_of_interactions(
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


class FunctionalModuleAnalyzer:
    '''
    Used to analyze to which extent interacting genes tend to share linkages.
    Given a datafame of linkages, a functional module graph and a list of q-value thresholds,
    this class perturbs the graph ~100 times and each time calculates averaged Jaccard coefficient
    for each threshold. Then, to ensure robustness of conclusions, simulation results are averaged and
    plotted altogether with 95% confidence intervals against actual data for the original graph.
    '''
    def __init__(self, qtl_type, qtl_df,
                 module_name, module_graph,
                 q_value_thresholds):
        self.qtl_type = qtl_type
        self.qtl_df = qtl_df
        self.mname = module_name
        self.mgraph = module_graph
        self.qval_list = q_value_thresholds

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
            qtl_graph = networks.assemble_graph_of_interactions(
                edges=qtl_df[["SNP", "gene"]].values,
                directed=True
            )
            jaccard[num_thresholds - idx - 1] = \
                networks.mean_linkage_similarity(modified_mgraph, qtl_graph)
        return jaccard

    def _randomize_and_recalc_stats(self, __dummy_arg):
        '''
        Performs the simulation: perturbs the graph and calculates Jaccard coefficients.

        :param __dummy_arg: Dummy variable required for Pool usage.
        :return: Numpy array of Jaccard coefficients for each threshold from self.qval_list.
        '''
        randomized_mgraph = self.mgraph.Degree_Sequence(
            self.mgraph.degree(),
            method="vl"
        )
        randomized_mgraph.vs["name"] = self.mgraph.vs["name"]
        return self._range_mean_jaccard(randomized_mgraph)

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

    '''
    TODO: Add argument to specify the path where to save the image. 
    '''
    def plot_robustness_analysis(self, actual_stats, simulation_stats):
        '''
        Plots simulation data (average among all the iterations and confidence intervals)
        against actual data obtained from the initial graph.

        :return: Saves the image in ./img/interactions
        '''
        plt.figure(figsize=(20, 10))
        plt.title('{} average linkage similarity'.format(self.qtl_type), fontsize=35)
        plt.xscale('log')

        plt.plot(self.qval_list, actual_stats, label="original")
        # Display also the confidence interval to estimate the statistical strength of conclusions
        plt.plot(self.qval_list, simulation_stats.mean(axis=0), label="randomized")
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

        plt.xticks(fontsize=20)
        plt.yticks(fontsize=20)

        plt.xlabel("{}".format(self.mname), fontsize=25)
        plt.ylabel("Mean Jaccard coefficient", fontsize=25)
        plt.legend(fontsize=25)

        plt.savefig("./img/interactions/" + self.qtl_type + '_' + self.mname + ".png")
        plt.close()

    def analyze_robustness_of_linkage_sharing(self):
        '''
        Interface function. Performs the simulation and saves the plot with results.
        '''
        actual_stats = self.actual_linkage_statistics()
        simulated_stats = self.simulated_linkage_statistics()
        self.plot_robustness_analysis(actual_stats, simulated_stats)

