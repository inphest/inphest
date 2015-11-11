#! /usr/bin/env python

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
                # edge_master_entry["revbayes_index"] = nd.annotations.get_value("index")
                # edge_master_entry["starting_state"] = nd.annotations.get_value("pa")
                # edge_master_entry["ending_state"] = nd.annotations.get_value("nd")

    def _extract_comment_metadata(self, nd):
        comment_str = ";".join(nd.comments)
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
