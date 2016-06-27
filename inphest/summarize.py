#! /usr/bin/env python
# -*- coding: utf-8 -*-

##############################################################################
##
##  Copyright 2015 Jeet Sukumaran and Mark T. Holder.
##  All rights reserved.
##
##  See "LICENSE.txt" for terms and conditions of usage.
##
##############################################################################

import collections
import itertools
import math
import dendropy
from dendropy.calculate import treecompare
from dendropy.calculate import profiledistance
from dendropy.utility import constants
from dendropy.calculate import statistics
from inphest import error

class SummaryStatsCalculator(object):

    @staticmethod
    def _euclidean_distance(v1, v2, is_weight_values_by_comparison_size=True):
        v1_size = len(v1)
        v2_size = len(v2)
        v1_idx = 0
        v2_idx = 0
        if v1_size > v2_size:
            v1_idx = v1_size - v2_size
            weight = float(v2_size)
        elif v2_size > v1_size:
            v2_idx = v2_size - v1_size
            weight = float(v1_size)
        else:
            weight = float(v1_size)
        if not is_weight_values_by_comparison_size:
            weight = 1.0
        ss = 0.0
        while v1_idx < v1_size and v2_idx < v2_size:
            ss += pow(v1[v1_idx]/weight - v2[v2_idx]/weight, 2)
            v1_idx += 1
            v2_idx += 1
        return math.sqrt(ss)

    @staticmethod
    def generate_induced_trees(
            tree,
            assemblage_leaf_sets,
            skip_null_assemblages=False,
            ):
        induced_trees = []
        for idx, assemblage_leaf_set in enumerate(assemblage_leaf_sets):
            if len(assemblage_leaf_set) == 0:
                if skip_null_assemblages:
                    continue
                raise error.PostTerminationFailedSimulationException("null assemblage set")
            node_filter_fn = lambda nd: nd in assemblage_leaf_set
            induced_tree = tree.extract_tree(
                               node_filter_fn=node_filter_fn,
                               is_apply_filter_to_leaf_nodes=True,
                               is_apply_filter_to_internal_nodes=False,
                               tree_factory=dendropy.Tree,
                               node_factory=dendropy.Node)
            induced_trees.append(induced_tree)
        return induced_trees

    def __init__(self, host_history, debug_mode):
        self.is_exchangeable_areas = True
        self.skip_null_symbiont_area_assemblages = True # If `False` requires all areas to have at least on symbiont lineage
        self.debug_mode = debug_mode
        self.ignore_incomplete_host_occupancies = False
        self.ignore_incomplete_area_occupancies = False
        self.bind_to_host_history(host_history)
        self.tree_shape_kernel = treecompare.TreeShapeKernel()
        self.num_profile_measurements = 6
        self.stat_name_delimiter = "."
        self.num_randomization_replicates = 100
        self.stat_name_prefix = "predictor"

    def get_profile_for_tree(self, tree):
        tree_profile = profiledistance.TreeProfile(
                tree=tree,
                is_measure_edge_lengths=True,
                is_measure_patristic_distances=True,
                is_measure_patristic_steps=True,
                is_measure_node_distances=True,
                is_measure_node_steps=True,
                is_measure_node_ages=True,
                is_measure_coalescence_intervals=False, # artificially imposes requirements that each area/host has at least 2 symbiont lineages
                is_normalize=True,
                ultrametricity_precision=constants.DEFAULT_ULTRAMETRICITY_PRECISION,
                tree_phylogenetic_distance_matrix=None,
                tree_node_distance_matrix=None,
                tree_id=None,
                is_skip_normalization_on_zero_division_error=True,
                )
        return tree_profile

    def bind_to_host_history(self, host_history):
        self.host_history = host_history
        self.host_tree = host_history.tree
        self.host_area_assemblage_trees = self.generate_induced_trees(
                tree=self.host_tree,
                assemblage_leaf_sets=self.host_history.area_assemblage_leaf_sets,
                skip_null_assemblages=False)
        self.host_tree_profile = self.get_profile_for_tree(self.host_tree)
        self.host_area_assemblage_tree_profiles = [self.get_profile_for_tree(t) for t in self.host_area_assemblage_trees]

    def calculate(self, symbiont_phylogeny, host_system, simulation_elapsed_time):
        old_taxon_namespace = self.preprocess_tree(symbiont_phylogeny)
        current_host_leaf_lineages = list(host_system.extant_host_lineages_at_current_time(simulation_elapsed_time))
        symbiont_phylogeny_leaf_sets_by_area = [set() for i in range(host_system.num_areas)]
        symbiont_phylogeny_leaf_sets_by_host = [set() for i in current_host_leaf_lineages]
        for leaf_idx, symbiont_lineage in enumerate(symbiont_phylogeny.leaf_node_iter()):
            for area in symbiont_lineage.area_iter():
                symbiont_phylogeny_leaf_sets_by_area[area.area_idx].add(symbiont_lineage.taxon)
            for host_idx, host_lineage in enumerate(current_host_leaf_lineages):
                # if symbiont_lineage.has_host(host_system.host_lineages_by_id[host.lineage_definition.lineage_id]):
                if symbiont_lineage.has_host(host_lineage):
                    symbiont_phylogeny_leaf_sets_by_host[host_idx].add(symbiont_lineage.taxon)
        if leaf_idx <= 2:
            raise error.InsufficientLineagesGenerated("Generated tree has too few lineages ({})".format(leaf_idx+1))
        if not self.ignore_incomplete_host_occupancies:
            if set() in symbiont_phylogeny_leaf_sets_by_host:
                raise error.IncompleteHostOccupancyException("incomplete host occupancy")
        if not self.ignore_incomplete_area_occupancies:
            if set() in symbiont_phylogeny_leaf_sets_by_area:
                raise error.IncompleteAreaOccupancyException("incomplete area occupancy")
        results = collections.OrderedDict()

        symbiont_pdm = symbiont_phylogeny.phylogenetic_distance_matrix()

        area_assemblage_descriptions = []
        for area_idx, areas in enumerate(symbiont_phylogeny_leaf_sets_by_area):
            regime = {
                "assemblage_basis_class_id": "area",
                "assemblage_basis_state_id": "state{}{}".format(self.stat_name_delimiter, area_idx),
            }
            area_assemblage_descriptions.append(regime)
        subresults = self._calc_community_ecology_stats(
            phylogenetic_distance_matrix=symbiont_pdm,
            assemblage_memberships=symbiont_phylogeny_leaf_sets_by_area,
            assemblage_descriptions=area_assemblage_descriptions,
            report_character_state_specific_results=False,
            report_character_class_wide_results=True,
            )
        if len(subresults) < 24:
            raise error.IncompleteAreaOccupancyException("Incomplete area occupancy")
        results.update(subresults)

        host_assemblage_descriptions = []
        for host_idx, hosts in enumerate(symbiont_phylogeny_leaf_sets_by_host):
            regime = {
                "assemblage_basis_class_id": "host",
                "assemblage_basis_state_id": "state{}{}".format(self.stat_name_delimiter, host_idx),
            }
            host_assemblage_descriptions.append(regime)
        subresults = self._calc_community_ecology_stats(
            phylogenetic_distance_matrix=symbiont_pdm,
            assemblage_memberships=symbiont_phylogeny_leaf_sets_by_host,
            assemblage_descriptions=host_assemblage_descriptions,
            report_character_state_specific_results=False,
            report_character_class_wide_results=True,
            )
        if len(subresults) < 24:
            raise error.IncompleteAreaOccupancyException("Incomplete host occupancy")
        results.update(subresults)

        self.restore_tree(symbiont_phylogeny, old_taxon_namespace)
        return results

    def _calc_community_ecology_stats(self,
            phylogenetic_distance_matrix,
            assemblage_memberships,
            assemblage_descriptions,
            report_character_state_specific_results=True,
            report_character_class_wide_results=True,
            ):

        assert len(assemblage_descriptions) == len(assemblage_memberships)

        summary_statistics_suite = collections.OrderedDict()
        results_by_character_class = collections.OrderedDict()
        stat_scores_to_be_harvested = ("obs", "z", "p",) # z = score, p = p-value (turns out this is quite informative)
        for sstbh in stat_scores_to_be_harvested:
            results_by_character_class[sstbh] = collections.defaultdict(list)

        for edge_weighted_desc in ("unweighted", "weighted"):
            if edge_weighted_desc:
                is_weighted_edge_distances = True
            else:
                is_weighted_edge_distances = False
            for underlying_statistic_type_desc in ("mpd", "mntd"):
                if underlying_statistic_type_desc == "mpd":
                    stat_fn_name = "standardized_effect_size_mean_pairwise_distance"
                else:
                    stat_fn_name = "standardized_effect_size_mean_nearest_taxon_distance"
                stat_fn = getattr(phylogenetic_distance_matrix, stat_fn_name)
                try:
                    results_group = stat_fn(
                        assemblage_memberships=assemblage_memberships,
                        is_weighted_edge_distances=is_weighted_edge_distances,
                        is_normalize_by_tree_size=True,
                        num_randomization_replicates=self.num_randomization_replicates,
                        )
                except dendropy.utility.error.SingleTaxonAssemblageException as e:
                    if not report_character_state_specific_results:
                        continue
                    else:
                        raise
                if not report_character_state_specific_results:
                    assert len(results_group) == len(assemblage_memberships)
                if len(results_group) == 0:
                    raise error.IncompleteStateSpaceOccupancyException
                for result, assemblage_desc in zip(results_group, assemblage_descriptions):
                    for ses_result_statistic in stat_scores_to_be_harvested:
                        character_class_statistic_prefix = self.stat_name_delimiter.join([
                            self.stat_name_prefix,
                            "community",
                            "by",
                            assemblage_desc["assemblage_basis_class_id"],
                            ])
                        statistic_subtype_desc = self.stat_name_delimiter.join([
                            edge_weighted_desc,
                            underlying_statistic_type_desc,
                            # assemblage_desc["assemblage_basis_state_id"],
                            ])
                        character_class_statistic_key = tuple([character_class_statistic_prefix, statistic_subtype_desc])
                        ses_result_statistic_value = getattr(result, ses_result_statistic)
                        if ses_result_statistic_value is None:
                            continue
                        if report_character_state_specific_results:
                            character_state_statistic_name = self.stat_name_delimiter.join([
                                character_class_statistic_prefix,
                                assemblage_desc["assemblage_basis_state_id"],
                                statistic_subtype_desc,
                                ses_result_statistic,
                                ])
                            assert character_state_statistic_name not in summary_statistics_suite
                            summary_statistics_suite[character_state_statistic_name] = ses_result_statistic_value
                        if report_character_class_wide_results:
                            results_by_character_class[ses_result_statistic][character_class_statistic_key].append(ses_result_statistic_value)
        if report_character_class_wide_results:
            for ses_result_statistic in results_by_character_class:
                if len(results_by_character_class[ses_result_statistic]) == 0:
                    continue
                for key in results_by_character_class[ses_result_statistic]:
                    character_class_statistic_prefix, statistic_subtype_desc = key
                    svalues = results_by_character_class[ses_result_statistic][key]
                    mean_var = statistics.mean_and_sample_variance(svalues)
                    for s, sdesc in zip( mean_var, ("mean", "var"), ):
                        sn_title = self.stat_name_delimiter.join([
                            character_class_statistic_prefix,
                            sdesc,
                            statistic_subtype_desc,
                            ses_result_statistic,
                            ])
                        assert sn_title not in summary_statistics_suite
                        summary_statistics_suite[sn_title] = s
        return summary_statistics_suite


    def calculate2(self, symbiont_phylogeny, host_system, simulation_elapsed_time):
        old_taxon_namespace = self.preprocess_tree(symbiont_phylogeny)

        current_host_leaf_lineages = list(host_system.extant_host_lineages_at_current_time(simulation_elapsed_time))
        symbiont_phylogeny_leaf_sets_by_area = [set() for i in range(host_system.num_areas)]
        symbiont_phylogeny_leaf_sets_by_host = [set() for i in current_host_leaf_lineages]
        for leaf_idx, symbiont_lineage in enumerate(symbiont_phylogeny.leaf_node_iter()):
            for area in symbiont_lineage.area_iter():
                symbiont_phylogeny_leaf_sets_by_area[area.area_idx].add(symbiont_lineage)
            for host_idx, host_lineage in enumerate(current_host_leaf_lineages):
                # if symbiont_lineage.has_host(host_system.host_lineages_by_id[host.lineage_definition.lineage_id]):
                if symbiont_lineage.has_host(host_lineage):
                    symbiont_phylogeny_leaf_sets_by_host[host_idx].add(symbiont_lineage)
        if leaf_idx <= 2:
            raise error.InsufficientLineagesGenerated("Generated tree has too few lineages ({})".format(leaf_idx+1))
        if not self.ignore_incomplete_host_occupancies:
            if set() in symbiont_phylogeny_leaf_sets_by_host:
                raise error.IncompleteHostOccupancyException("incomplete host occupancy")
        if not self.ignore_incomplete_area_occupancies:
            if set() in symbiont_phylogeny_leaf_sets_by_area:
                raise error.IncompleteAreaOccupancyException("incomplete area occupancy")

        symbiont_area_assemblage_trees = self.generate_induced_trees(
                tree=symbiont_phylogeny,
                assemblage_leaf_sets=symbiont_phylogeny_leaf_sets_by_area,
                skip_null_assemblages=False)
        symbiont_host_assemblage_trees = self.generate_induced_trees(
                tree=symbiont_phylogeny,
                assemblage_leaf_sets=symbiont_phylogeny_leaf_sets_by_host,
                skip_null_assemblages=False)

        results = collections.OrderedDict()
        expected_num_sum_stats = 0

        ## main trees kernel trick
        expected_num_sum_stats += 1
        results["predictor.primary.tree.tsktd"] = self.tree_shape_kernel(
                tree1=self.host_tree,
                tree2=symbiont_phylogeny)
        self.check_successful_subcalculation(expected_num_sum_stats, results, "predictor.primary.tree.tsktd")

        ## area trees kernel trick
        expected_num_sum_stats += ( min(len(self.host_area_assemblage_trees), len(symbiont_area_assemblage_trees)) )
        results.update(self.tree_shape_kernel_compare_trees(
            trees1=self.host_area_assemblage_trees,
            trees2=symbiont_area_assemblage_trees,
            fieldname_prefix="predictor.area.assemblage.tsktd",
            fieldname_suffix="",
            is_exchangeable_assemblage_classifications=True,
            default_value_for_missing_comparisons=False,
            ))
        for induced_tree in symbiont_area_assemblage_trees:
            self.tree_shape_kernel.remove_from_cache(induced_tree)
        self.check_successful_subcalculation(expected_num_sum_stats, results, "predictor.area.assemblage.tsktd")

        ## host trees kernel trick
        expected_num_sum_stats += ( min(len(self.host_area_assemblage_trees), len(symbiont_host_assemblage_trees)) )
        results.update(self.tree_shape_kernel_compare_trees(
            trees1=self.host_area_assemblage_trees,
            trees2=symbiont_host_assemblage_trees,
            fieldname_prefix="predictor.host.assemblage.tsktd.",
            fieldname_suffix="",
            is_exchangeable_assemblage_classifications=True,
            default_value_for_missing_comparisons=False,
            ))
        for induced_tree in symbiont_host_assemblage_trees:
            self.tree_shape_kernel.remove_from_cache(induced_tree)
        self.check_successful_subcalculation(expected_num_sum_stats, results, "predictor.host.assemblage.tsktd.")

        ## host trees vs. area trees kernel trick
        expected_num_sum_stats += ( min(len(self.host_area_assemblage_trees), len(symbiont_host_assemblage_trees)) )
        results.update(self.tree_shape_kernel_compare_trees(
            trees1=symbiont_area_assemblage_trees,
            trees2=symbiont_host_assemblage_trees,
            fieldname_prefix="predictor.host.vs.area.assemblage.tsktd",
            fieldname_suffix="",
            is_exchangeable_assemblage_classifications=True,
            default_value_for_missing_comparisons=False,
            ))
        for induced_tree in symbiont_host_assemblage_trees:
            self.tree_shape_kernel.remove_from_cache(induced_tree)
        self.check_successful_subcalculation(expected_num_sum_stats, results, "predictor.host.vs.area.assemblage.tsktd")

        ## main tree profile distance
        expected_num_sum_stats += self.num_profile_measurements
        symbiont_tree_profile = self.get_profile_for_tree(tree=symbiont_phylogeny)
        self.compare_profiles(
                profile1=self.host_tree_profile,
                profile2=symbiont_tree_profile,
                fieldname_prefix="predictor.profiledist.main.trees.",
                fieldname_suffix="",
                results=results)
        self.check_successful_subcalculation(expected_num_sum_stats, results, "predictor.profiledist.main.trees.")

        host_area_assemblage_profiles = [self.get_profile_for_tree(t) for t in self.host_area_assemblage_trees]
        symbiont_host_assemblage_profiles = [self.get_profile_for_tree(t) for t in symbiont_host_assemblage_trees]
        symbiont_area_assemblage_profiles = [self.get_profile_for_tree(t) for t in symbiont_area_assemblage_trees]

        expected_num_sum_stats += min( len(host_area_assemblage_profiles), len(symbiont_area_assemblage_profiles) ) * self.num_profile_measurements
        self.compare_multi_profiles(
                profiles1=host_area_assemblage_profiles,
                profiles2=symbiont_area_assemblage_profiles,
                fieldname_prefix="predictor.profiledist.area.assemblage.",
                fieldname_suffix="",
                results=results)
        self.check_successful_subcalculation(expected_num_sum_stats, results,"predictor.profiledist.area.assemblage." )

        expected_num_sum_stats += min( len(host_area_assemblage_profiles), len(symbiont_host_assemblage_profiles) ) * self.num_profile_measurements
        self.compare_multi_profiles(
                profiles1=host_area_assemblage_profiles,
                profiles2=symbiont_host_assemblage_profiles,
                fieldname_prefix="predictor.profiledist.host.assemblage.",
                fieldname_suffix="",
                results=results)
        self.check_successful_subcalculation(expected_num_sum_stats, results, "predictor.profiledist.host.assemblage.")

        expected_num_sum_stats += min( len(host_area_assemblage_profiles), len(symbiont_area_assemblage_profiles) ) * self.num_profile_measurements
        self.compare_multi_profiles(
                profiles1=host_area_assemblage_profiles,
                profiles2=symbiont_area_assemblage_profiles,
                fieldname_prefix="predictor.profiledist.host.vs.area.assemblage.",
                fieldname_suffix="",
                results=results)
        self.check_successful_subcalculation(expected_num_sum_stats, results, "predictor.profiledist.host.vs.area.assemblage.")

        self.restore_tree(symbiont_phylogeny, old_taxon_namespace)
        return results

    def check_successful_subcalculation(self, expected_num_sum_stats, results, message):
        if len(results) != expected_num_sum_stats:
            raise error.SummaryStatisticsCalculationFailure("Failed to calculate one or more summary statistics: expecting {} results but only found {} ({})".format(expected_num_sum_stats, len(results), message))

    def compare_profiles(self,
            profile1,
            profile2,
            fieldname_prefix,
            fieldname_suffix,
            results):
        d = profile1.measure_distances(profile2)
        for key in d:
            results["{}{}{}".format(fieldname_prefix, key, fieldname_suffix)] = d[key]
        return results

    def compare_multi_profiles(self,
            profiles1,
            profiles2,
            fieldname_prefix,
            fieldname_suffix,
            results,
            default_value_for_missing_comparisons=False,
            ):
        if len(profiles1) > len(profiles2):
            profiles2, profiles1 = profiles1, profiles2
        comparison_vectors = {}
        current_minimum_distances = {}
        current_joint_minimum_vectors = {}
        measurement_names = profiles1[0].measurement_names
        for name in measurement_names:
            comparison_vectors[name] = [0.0] * len(profiles1)
            current_minimum_distances[name] = None
            current_joint_minimum_vectors[name] = None
        for profile_permutation in itertools.permutations(profiles2, len(profiles1)):
            distances = {}
            for p2, p1 in zip(profile_permutation, profiles1):
                try:
                    pd = p1.measure_distances(p2)
                except ValueError:
                    ### TODO: better handling here: (1) specialized error raised by profiledistance + option of quitting or continuing
                    continue
                for name in measurement_names:
                    try:
                        distances[name].append(pd[name])
                    except KeyError:
                        distances[name] = [pd[name]]
            if not distances:
                ### TODO: ditto re: better handling
                continue
            for name in measurement_names:
                ed = self._euclidean_distance(distances[name], comparison_vectors[name])
                if current_minimum_distances[name] is None or ed < current_minimum_distances[name]:
                    current_minimum_distances[name] = ed
                    current_joint_minimum_vectors[name] = distances[name]
        if  not current_joint_minimum_vectors:
            raise error.InsufficientLineagesGenerated("Insufficient symbiont lineages in one or more hosts or areas")
        for name in measurement_names:
            for didx, d in enumerate(current_joint_minimum_vectors[name]):
                results["{}{}{}{}".format(fieldname_prefix, name, didx+1, fieldname_suffix, )] = d
            if default_value_for_missing_comparisons is not False:
                for didx in range(didx+1, len(profiles2)):
                    results["{}{}{}{}".format(fieldname_prefix, name, didx+1, fieldname_suffix, )] = default_value_for_missing_comparisons

    def tree_shape_kernel_compare_trees(self,
            trees1,
            trees2,
            fieldname_prefix,
            fieldname_suffix,
            is_exchangeable_assemblage_classifications,
            default_value_for_missing_comparisons,
            **kwargs
            ):
        score_table = collections.OrderedDict()
        if not is_exchangeable_assemblage_classifications:
            if len(trees1) != len(trees2):
                raise TypeError("Different numbers of induced trees not supported for non-exchangeable classifications: {} vs. {}".format(len(trees1), len(trees2)))
            for idx, (induced_tree1, induced_tree2) in enumerate(zip(trees1, trees2)):
                s = self.tree_shape_kernel(
                                tree1=induced_tree1,
                                tree2=induced_tree2,
                                is_tree1_cache_updated=kwargs.pop("is_tree1_cache_updated", True),
                                is_tree2_cache_updated=kwargs.pop("is_tree1_cache_updated", True),
                                )
                score_table["{}{}{}".format(fieldname_prefix, idx+1, fieldname_suffix, )] = s
        else:
            # ensure trees1 has the smaller number of elements
            if len(trees1) > len(trees2):
                trees2, trees1 = trees1, trees2
            comparison_vector = [0.0] * len(trees1)
            current_minimum_distance = None
            current_joint_minimum_vector = None
            for trees_permutation in itertools.permutations(trees2, len(trees1)):
                distances = []
                for t2, t1 in zip(trees_permutation, trees1):
                    distances.append(self.tree_shape_kernel(
                            tree1=t1,
                            tree2=t2,
                            is_tree1_cache_updated=kwargs.pop("is_tree1_cache_updated", True),
                            is_tree2_cache_updated=kwargs.pop("is_tree1_cache_updated", True),
                            ))
                euclidean_distance = self._euclidean_distance(distances, comparison_vector)
                if current_minimum_distance is None or euclidean_distance < current_minimum_distance:
                    current_minimum_distance = euclidean_distance
                    current_joint_minimum_vector = distances
            for didx, d in enumerate(current_joint_minimum_vector):
                score_table["{}{}{}".format(fieldname_prefix, didx+1, fieldname_suffix, )] = d
            if default_value_for_missing_comparisons is not False:
                for didx in range(didx+1, len(trees2)):
                    score_table["{}{}{}".format(fieldname_prefix, didx+1, fieldname_suffix, )] = default_value_for_missing_comparisons
        return score_table

    def preprocess_tree(self, tree):
        old_taxon_namespace = tree.taxon_namespace
        tree.taxon_namespace = dendropy.TaxonNamespace()
        for nd in tree:
            nd.taxon = tree.taxon_namespace.require_taxon(label="T{}".format(nd.index))
            nd.taxon.lineage = nd
        tree.is_rooted = True
        tree.encode_bipartitions()
        return old_taxon_namespace

    def restore_tree(self, tree, taxon_namespace):
        for nd in tree:
            nd.taxon = None
        tree.taxon_namespace = taxon_namespace

