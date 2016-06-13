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
        self.fieldname_prefix = kwargs.pop("fieldname_prefix", "")
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
            ):
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
        score_table["{}primary.tree.ktd".format(self.fieldname_prefix)] = main_trees_score
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
                score_table["{}induced.tree.{}.ktd".format(self.fieldname_prefix, idx+1)] = s
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
                    score_table["{}induced.tree.{}.ktd".format(self.fieldname_prefix, didx+1)] = d
                for didx in range(didx+1, self._num_assemblage_classifications):
                    score_table["{}induced.tree.{}.ktd".format(self.fieldname_prefix, didx+1)] = "NA"
            else:
                raise NotImplementedError()
        return score_table

class SummaryStatsCalculator(object):

    def __init__(self, host_history, debug_mode):
        self.is_exchangeable_areas = True
        self.host_history = host_history
        self.skip_null_symbiont_area_assemblages = True # If `False` requires all areas to have at least on symbiont lineage
        self.area_tree_shape_comparator = self.create_area_comparator()
        self.debug_mode = debug_mode

    def create_area_comparator(self):
        area_tree_shape_comparator = AssemblageInducedTreeShapeKernel(
                is_exchangeable_assemblage_classifications=self.is_exchangeable_areas,
                induced_tree_factory=dendropy.Tree,
                induced_tree_node_factory=dendropy.Node,
                fieldname_prefix="predictor.area.",
                )
        area_tree_shape_comparator.update_assemblage_induced_tree_cache(
                tree=self.host_history.tree,
                assemblage_leaf_sets=self.host_history.area_assemblage_leaf_sets,
                )
        return area_tree_shape_comparator

    def calculate(self, symbiont_phylogeny, host_system):
        old_taxon_namespace = self.preprocess_tree(symbiont_phylogeny)
        symbiont_phylogeny_leaf_sets = [set() for i in range(host_system.num_areas)]
        for area_idx, area in enumerate(host_system.areas):
            assert area_idx == area.area_idx
            for symbiont_lineage in area.symbiont_lineages:
                if self.debug_mode:
                    assert symbiont_lineage.has_area(area)
                symbiont_phylogeny_leaf_sets[area_idx].add(symbiont_lineage)
        results = self.area_tree_shape_comparator(
                tree1=self.host_history.tree,
                tree2=symbiont_phylogeny,
                tree1_assemblage_leaf_sets=None,
                tree2_assemblage_leaf_sets=symbiont_phylogeny_leaf_sets,
                is_tree1_cache_updated=True,
                is_tree2_cache_updated=False,
                )
        self.area_tree_shape_comparator.remove_from_cache(symbiont_phylogeny)
        self.restore_tree(symbiont_phylogeny, old_taxon_namespace)
        return results

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
