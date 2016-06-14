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
from inphest import error

class AssemblageInducedTreeShapeKernel(treecompare.TreeShapeKernel, treecompare.AssemblageInducedTreeManager):

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

    def __init__(self, *args, **kwargs):
        self.exchangeable_assemblage_comparison_strategy = kwargs.pop("exchangeable_assemblage_comparison_strategy", "joint minimum")
        treecompare.TreeShapeKernel.__init__(self, *args, **kwargs)
        treecompare.AssemblageInducedTreeManager.__init__(self, *args, **kwargs)

    def remove_from_cache(self, tree):
        for induced_tree in self._tree_assemblage_induced_trees_map[tree]:
            treecompare.TreeShapeKernel.remove_from_cache(self, induced_tree)
        treecompare.TreeShapeKernel.remove_from_cache(self, tree)
        treecompare.AssemblageInducedTreeManager.remove_from_cache(self, tree)

    def update_assemblage_induced_tree_cache(self,
            tree,
            assemblage_leaf_sets):
        self.update_cache(tree=tree)
        induced_trees = self.generate_induced_trees(tree=tree,
                assemblage_leaf_sets=assemblage_leaf_sets)
        for induced_tree in induced_trees:
            self.update_cache(tree=induced_tree)

    def __call__(self,
            tree1,
            tree2,
            tree1_assemblage_leaf_sets,
            tree2_assemblage_leaf_sets,
            is_tree1_cache_updated=False,
            is_tree2_cache_updated=False,
            fieldname_prefix=None,
            ):
        if fieldname_prefix is None:
            fieldname_prefix = ""
        main_trees_score = treecompare.TreeShapeKernel.__call__(self,
                tree1=tree1,
                tree2=tree2,
                is_tree1_cache_updated=is_tree1_cache_updated,
                is_tree2_cache_updated=is_tree2_cache_updated,
                )
        if not is_tree1_cache_updated or tree1 not in self._tree_assemblage_induced_trees_map:
            if tree1_assemblage_leaf_sets is None:
                raise ValueError("Uncached tree requires specification of 'tree1_assemblage_leaf_sets'")
            self.update_assemblage_induced_tree_cache(
                    tree=tree1,
                    assemblage_leaf_sets=tree1_assemblage_leaf_sets)
        if not is_tree2_cache_updated or tree2 not in self._tree_assemblage_induced_trees_map:
            if tree2_assemblage_leaf_sets is None:
                raise ValueError("Uncached tree requires specification of 'tree2_assemblage_leaf_sets'")
            self.update_assemblage_induced_tree_cache(
                    tree=tree2,
                    assemblage_leaf_sets=tree2_assemblage_leaf_sets)
        ## ++ main tree score
        score_table = collections.OrderedDict()
        score_table["{}primary.tree.ktd".format(fieldname_prefix)] = main_trees_score
        induced_trees1 = self._tree_assemblage_induced_trees_map[tree1]
        induced_trees2 = self._tree_assemblage_induced_trees_map[tree2]
        # assert len(induced_trees1) == len(induced_trees2) == self._num_assemblage_classifications
        if not self.is_exchangeable_assemblage_classifications:
            if len(induced_trees1) != len(induced_trees2):
                raise TypeError("Different numbers of induced trees not supported for non-exchangeable classifications: {} vs. {}".format(len(induced_trees1), len(induced_trees2)))
            for idx, (induced_tree1, induced_tree2) in enumerate(zip(induced_trees1, induced_trees2)):
                s = treecompare.TreeShapeKernel.__call__(self,
                                tree1=induced_tree1,
                                tree2=induced_tree2,
                                is_tree1_cache_updated=True,
                                is_tree2_cache_updated=True,
                                )
                ## ++ raw scores direct comparisons of each of the induced trees
                score_table["{}induced.tree.{}.ktd".format(fieldname_prefix, idx+1)] = s
        else:
            if self.exchangeable_assemblage_comparison_strategy == "joint minimum":
                # if lengths are different, we want to fix the smaller set
                if len(induced_trees1) > len(induced_trees2):
                    induced_trees2, induced_trees1 = induced_trees1, induced_trees2
                comparison_vector = [0.0] * len(induced_trees1)
                current_minimum_distance = None
                current_joint_minimum_vector = None
                for induced_trees_permutation in itertools.permutations(induced_trees2, len(induced_trees1)):
                    distances = []
                    for t2, t1 in zip(induced_trees_permutation, induced_trees1):
                        distances.append(treecompare.TreeShapeKernel.__call__(self,
                                tree1=t1,
                                tree2=t2,
                                is_tree1_cache_updated=True,
                                is_tree2_cache_updated=True,))
                    euclidean_distance = self._euclidean_distance(distances, comparison_vector)
                    if current_minimum_distance is None or euclidean_distance < current_minimum_distance:
                        current_minimum_distance = euclidean_distance
                        current_joint_minimum_vector = distances
                for didx, d in enumerate(distances):
                    score_table["{}induced.tree.{}.ktd".format(fieldname_prefix, didx+1)] = d
                for didx in range(didx+1, self._num_assemblage_classifications):
                    score_table["{}induced.tree.{}.ktd".format(fieldname_prefix, didx+1)] = "NA"
            else:
                raise NotImplementedError()
        return score_table

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
                raise error.InsufficientFocalAreaLineagesSimulationException()
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

    def bind_to_host_history(self, host_history):
        self.host_history = host_history
        self.host_tree = host_history.tree
        self.host_area_assemblage_trees = self.generate_induced_trees(
                tree=self.host_tree,
                assemblage_leaf_sets=self.host_history.area_assemblage_leaf_sets,
                skip_null_assemblages=False)

    def calculate(self, symbiont_phylogeny, host_system, simulation_elapsed_time):
        old_taxon_namespace = self.preprocess_tree(symbiont_phylogeny)

        current_host_leaf_lineages = list(host_system.extant_host_lineages_at_current_time(simulation_elapsed_time))
        symbiont_phylogeny_leaf_sets_by_area = [set() for i in range(host_system.num_areas)]
        symbiont_phylogeny_leaf_sets_by_host = [set() for i in current_host_leaf_lineages]
        for symbiont_lineage in symbiont_phylogeny.leaf_node_iter():
            for area in symbiont_lineage.area_iter():
                symbiont_phylogeny_leaf_sets_by_area[area.area_idx].add(symbiont_lineage)
            for host_idx, host_lineage in enumerate(current_host_leaf_lineages):
                # if symbiont_lineage.has_host(host_system.host_lineages_by_id[host.lineage_definition.lineage_id]):
                if symbiont_lineage.has_host(host_lineage):
                    symbiont_phylogeny_leaf_sets_by_host[host_idx].add(symbiont_lineage)
        if not self.ignore_incomplete_host_occupancies:
            if set() in symbiont_phylogeny_leaf_sets_by_host:
                raise error.IncompleteHostOccupancyException("incomplete host occupancy")
        if not self.ignore_incomplete_area_occupancies:
            if set() in symbiont_phylogeny_leaf_sets_by_area:
                raise error.InsufficientFocalAreaLineagesSimulationException("incomplete area occupancy")

        ## main trees kernel trick
        results = collections.OrderedDict()
        results["predictor.primary.tree.tsktd"] = self.tree_shape_kernel(
                tree1=self.host_tree,
                tree2=symbiont_phylogeny)

        ## area trees kernel trick
        symbiont_area_assemblage_trees = self.generate_induced_trees(
                tree=symbiont_phylogeny,
                assemblage_leaf_sets=symbiont_phylogeny_leaf_sets_by_area,
                skip_null_assemblages=False)
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

        ## host trees kernel trick
        symbiont_host_assemblage_trees = self.generate_induced_trees(
                tree=symbiont_phylogeny,
                assemblage_leaf_sets=symbiont_phylogeny_leaf_sets_by_host,
                skip_null_assemblages=False)
        results.update(self.tree_shape_kernel_compare_trees(
            trees1=self.host_area_assemblage_trees,
            trees2=symbiont_host_assemblage_trees,
            fieldname_prefix="predictor.host.assemblage.tsktd",
            fieldname_suffix="",
            is_exchangeable_assemblage_classifications=True,
            default_value_for_missing_comparisons=False,
            ))
        for induced_tree in symbiont_area_assemblage_trees:
            self.tree_shape_kernel.remove_from_cache(induced_tree)

        self.restore_tree(symbiont_phylogeny, old_taxon_namespace)
        return results

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
            for didx, d in enumerate(distances):
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

    # def create_area_comparator(self):
    #     area_tree_shape_comparator = AssemblageInducedTreeShapeKernel(
    #             is_exchangeable_assemblage_classifications=self.is_exchangeable_areas,
    #             induced_tree_factory=dendropy.Tree,
    #             induced_tree_node_factory=dendropy.Node,
    #             )
    #     area_tree_shape_comparator.update_assemblage_induced_tree_cache(
    #             tree=self.host_history.tree,
    #             assemblage_leaf_sets=self.host_history.area_assemblage_leaf_sets,
    #             )
    #     return area_tree_shape_comparator

    # def create_host_comparator(self):
    #     host_tree_shape_comparator = AssemblageInducedTreeShapeKernel(
    #             is_exchangeable_assemblage_classifications=False,
    #             induced_tree_factory=dendropy.Tree,
    #             induced_tree_node_factory=dendropy.Node,
    #             skip_null_assemblages=True,
    #             )
    #     # Note this order will vary from run to run (because `host_system.extant_host_leaf_lineages` is a set).
    #     # Also note that we assume here that this calculation occurs at the end of the run; must use `extant_host_lineages_at_current_time` if not.
    #     self.extant_host_leaf_lineages = list(self.host_history.extant_leaf_lineages)
    #     ## basically, our induced subtrees area the entire tree for each host
    #     host_tree_shape_comparator.update_assemblage_induced_tree_cache(
    #             tree=self.host_history.tree,
    #             assemblage_leaf_sets=[set(self.extant_host_leaf_lineages) for i in self.extant_host_leaf_lineages],
    #             )
    #     return host_tree_shape_comparator

