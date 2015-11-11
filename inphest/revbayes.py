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

"""
RevBayes biogeography analysis and parsing.
"""

import os
import sys
import csv
import collections
import re
import dendropy


class RevBayesBiogeographyParser(object):

    EVENT_PATTERN = re.compile(r"{(.*?)}")

    def __init__(self):
        ## information on each tree in sample
        self.tree_samples = []
        ## events organized by tree
        self.event_schedules_by_tree = {}
        ## information on each bipartition across all samples
        self.bipartition_samples = []
        ## events organized by bipartition
        self.event_schedules_by_edge = collections.defaultdict(dendropy.Bipartition)

    def parse(self, src, skip_first_row=True):
        if isinstance(src, str):
            src = open(src)
        if skip_first_row:
            next(src) # skip over header
        taxon_namespace = dendropy.TaxonNamespace()
        tree_idx = -1
        for row_idx, row in enumerate(src):
            row = row.strip("\n")
            if not row:
                continue

            tree_idx += 1
            parts = row.split("\t")
            iteration, posterior, likelihood, prior, tree_str = row.split("\t")

            tree_sample_entry = {}
            tree_sample_entry["tree_idx"] = tree_idx
            tree_sample_entry["iteration"] = float(iteration)
            tree_sample_entry["posterior"] = float(posterior)
            tree_sample_entry["ln_likelihood"] = float(likelihood)
            tree_sample_entry["prior"] = float(prior)

            self.tree_samples.append(tree_sample_entry)
            tree = dendropy.Tree.get(
                    data=tree_str,
                    schema="newick",
                    taxon_namespace=taxon_namespace,
                    terminating_semicolon_required=False,
                    rooting="force-rooted",
                    extract_comment_metadata=False,
                    )
            tree.encode_bipartitions()
            # print(tree.as_string("newick", suppress_annotations=True))
            tree.calc_node_ages(ultrametricity_precision=0.01)
            for nd in tree:
                edge_master_entry = {}
                edge_master_entry["tree_idx"] = tree_idx
                edge_master_entry["bipartition_idx"] = int(nd.edge.bipartition)
                edge_master_entry["split"] = nd.edge.bipartition.split_as_bitstring
                edge_master_entry["leafset"] = nd.edge.bipartition.leafset_as_bitstring
                edge_master_entry["duration"] = nd.edge.length
                if nd.parent_node:
                    edge_master_entry["starting_age"] = nd.parent_node.age
                else:
                    edge_master_entry["starting_age"] = None
                edge_master_entry["ending_age"] = nd.age
                if nd.is_leaf():
                    edge_master_entry["child0"] = None
                    edge_master_entry["child1"] = None
                else:
                    for ch_idx, ch in enumerate(nd.child_node_iter()):
                        edge_master_entry["child{}".format(ch_idx)] = int(ch.edge.bipartition)
                edge_metadata, edge_events = self._extract_comment_metadata(nd)
                edge_master_entry["revbayes_index"] = edge_metadata["index"]
                edge_master_entry["starting_state"] = edge_metadata["pa"]
                edge_master_entry["ending_state"] = edge_metadata["nd"]
                if "cs" in edge_metadata:
                    if edge_metadata["cs"] == "s":
                        edge_master_entry["cladogenetic_speciation_mode"] = "subset_sympatry"
                    elif edge_metadata["cs"] == "n":
                        edge_master_entry["cladogenetic_speciation_mode"] = "narrow_sympatry"
                    elif edge_metadata["cs"] == "w":
                        edge_master_entry["cladogenetic_speciation_mode"] = "widespread_sympatry"
                    elif edge_metadata["cs"] == "a":
                        edge_master_entry["cladogenetic_speciation_mode"] = "allopatry"
                    else:
                        raise ValueError("Unrecognized cladogenetic speciation mode event type: '{}'".format(edge_metadata["cs"]))
                else:
                    edge_master_entry["cladogenetic_speciation_mode"] = "NA"

    def _extract_comment_metadata(self, nd):

        #
        # Michael Landis, pers. comm., 2015-11-09:
        #
        #    Galago_senegalensis[&index=22;nd=0000000000000111000000000;pa=0100000000000101000000000;ev={{t:0.393107,a:10.3104,s:1,i:14},{t:0.70063,a:5.08596,s:0,i:1}}]:16.9889
        #
        #    which says the branch has index 22,
        #    it ends with state 0000000000000111000000000,
        #    it begins with state 0100000000000101000000000,
        #    with two intermediate events:
        #    area 14 transitions to state 1 at age 10.3104
        #    area 1 transitions to state 0 at age 5.08496
        #
        #    Basically, I take all these node metadata for all MCMC outputs and compile them into that dictionary in bg_parse. bg_parse.py has extremely poor documentation and might not see further development. I can dig into it more later, but here is an explanation from memory. In any case, if you are curious about these keys for a given node:
        #
        #    ['ch1', 'iteration', 'bn', 'nd', 'ch0', 'prior', 'posterior', 'pa', 'cs', 'ev', ' likelihood']
        #
        #    'bn' is the branch node index, I think...
        #    'nd' gives the node state preceding speciation
        #    'ch0','ch1' gives the daughter node states following speciation (left and right nodes)
        #    'cs' gives the cladogenic event type
        #    'pa' gives the ancestral state at the base of the branch
        #    'ev' gives an event vector of transition events
        #    'iteration', 'prior', 'posterior', 'likelihood' just report MCMC values
        #

        comment_str = ";".join(nd.comments)
        assert comment_str[0] == "&"
        comment_str = comment_str[1:]

        metadata_items_as_str = comment_str.split(";")
        edge_metadata = {}
        events_str = None
        for metadata_item in metadata_items_as_str:
            key, val = metadata_item.split("=")
            if key == "ev":
                events_str = val
            else:
                edge_metadata[key] = val
        edge_events = []
        if events_str:
            assert events_str[0] == "{" and events_str[-1] == "}", events_str
            events_str = events_str[1:-1] # strip outer braces
            if events_str:
                event_items_as_str = RevBayesBiogeographyParser.EVENT_PATTERN.findall(events_str)
                num_right_braces = events_str.count("}")
                # check that we pulled the correct number of items
                assert len(event_items_as_str) == num_right_braces, "{} != {}; '{}': {}".format(len(event_items_as_str), num_right_braces, event_items_as_str, events_str)
                for event_item_str in event_items_as_str:
                    event_item_parts = event_item_str.split(",")
                    assert len(event_item_parts) == 4, "{}:{}".format(event_item_str, event_item_parts)
                    event = {}
                    for event_item_part in event_item_parts:
                        key, val = event_item_part.split(":")
                        if key == "t":
                            event["time"] = float(val)
                        elif key == "a":
                            event["age"] = float(val)
                        elif key == "i":
                            event["area_idx"] = int(val)
                        elif key == "s":
                            event["to_state"] = val
                        else:
                            raise ValueError("Unrecognized entry key: '{}'".format(key))
                    assert len(event) == 4, event
                    assert "time" in event
                    assert "age" in event
                    assert "area_idx" in event
                    assert "to_state" in event
                    edge_events.append(event)
        return edge_metadata, edge_events

def main():
    rbp = RevBayesBiogeographyParser()
    rbp.parse(sys.argv[1])

if __name__ == "__main__":
    main()
