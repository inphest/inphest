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
    TREE_SAMPLES_FIELDNAMES = ("tree_idx", "iteration", "posterior", "ln_likelihood", "prior")
    EDGE_SAMPLES_FIELDNAMES = (
            "tree_idx",
            "edge_id",
            "split_bitstring",
            "leafset_bitstring",
            "edge_starting_age",
            "edge_ending_age",
            "edge_duration",
            "child0_edge_id",
            "child1_edge_id",
            # "edge_revbayes_index",
            "edge_starting_state",
            "edge_ending_state",
            "edge_cladogenetic_speciation_mode",
            )
    EVENT_SAMPLES_FIELDNAMES = (
            "tree_idx",
            "edge_id",
            "age",
            "event_type",
            "event_subtype",
            "area_idx",
            "child0_edge_id",
            "child1_edge_id",
            )
    NULL_VALUE = "NA"

    def __init__(self):
        ## information on each tree in sample
        self.tree_entries = []
        ## events organized by tree
        self.event_schedules_by_tree = []
        ## information on each bipartition across all samples
        self.edge_entries = []
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

            tree_entry = {}
            tree_entry["tree_idx"] = tree_idx
            tree_entry["iteration"] = float(iteration)
            tree_entry["posterior"] = float(posterior)
            tree_entry["ln_likelihood"] = float(likelihood)
            tree_entry["prior"] = float(prior)

            self.tree_entries.append(tree_entry)
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
                edge_entry = {}
                edge_entry["tree_idx"] = tree_idx
                edge_entry["edge_id"] = int(nd.edge.bipartition)
                edge_entry["split_bitstring"] = nd.edge.bipartition.split_as_bitstring()
                edge_entry["leafset_bitstring"] = nd.edge.bipartition.leafset_as_bitstring()
                edge_entry["edge_duration"] = nd.edge.length
                edge_entry["edge_ending_age"] = nd.age
                if nd.parent_node:
                    edge_entry["edge_starting_age"] = nd.parent_node.age
                else:
                    # special case for root
                    edge_entry["edge_starting_age"] = 0.0
                    edge_entry["edge_duration"] = nd.age
                if nd.is_leaf():
                    edge_entry["child0_edge_id"] = RevBayesBiogeographyParser.NULL_VALUE
                    edge_entry["child1_edge_id"] = RevBayesBiogeographyParser.NULL_VALUE
                else:
                    for ch_idx, ch in enumerate(nd.child_node_iter()):
                        edge_entry["child{}_edge_id".format(ch_idx)] = int(ch.edge.bipartition)
                edge_metadata, edge_events = self._extract_comment_metadata(nd)
                # edge_entry["edge_revbayes_index"] = edge_metadata["index"]
                edge_entry["edge_starting_state"] = edge_metadata["pa"]
                edge_entry["edge_ending_state"] = edge_metadata["nd"]
                if "cs" in edge_metadata:
                    if edge_metadata["cs"] == "s":
                        edge_entry["edge_cladogenetic_speciation_mode"] = "subset_sympatry"
                    elif edge_metadata["cs"] == "n":
                        edge_entry["edge_cladogenetic_speciation_mode"] = "narrow_sympatry"
                    elif edge_metadata["cs"] == "w":
                        edge_entry["edge_cladogenetic_speciation_mode"] = "widespread_sympatry"
                    elif edge_metadata["cs"] == "a":
                        edge_entry["edge_cladogenetic_speciation_mode"] = "allopatry"
                    else:
                        raise ValueError("Unrecognized cladogenetic speciation mode event type: '{}'".format(edge_metadata["cs"]))
                else:
                    edge_entry["edge_cladogenetic_speciation_mode"] = RevBayesBiogeographyParser.NULL_VALUE
                self.edge_entries.append(edge_entry)
                for event in edge_events:
                    event_entry = {}
                    event_entry["tree_idx"] = edge_entry["tree_idx"]
                    event_entry["edge_id"] = edge_entry["edge_id"]
                    event_entry["age"] = event["age"]
                    # event_entry["time"] = event["time"]
                    event_entry["event_type"] = "anagenesis"
                    if event["to_state"] == "1":
                        event_entry["event_subtype"] = "area_gain"
                    elif event["to_state"] == "0":
                        event_entry["event_subtype"] = "area_loss"
                    else:
                        raise ValueError("Unexpected value for state: expecting '0' or '1' but found '{}'".format(event["to_state"]))
                    event_entry["area_idx"] = event["area_idx"]
                    # event_entry["to_state"] = event["to_state"]
                    self.event_schedules_by_tree.append(event_entry)
                ## handle splitting event
                if not nd.is_leaf():
                    self.event_schedules_by_tree.append({
                        "tree_idx": edge_entry["tree_idx"],
                        "edge_id": edge_entry["edge_id"],
                        "age": edge_entry["edge_ending_age"],
                        "event_type": "cladogenesis",
                        "event_subtype": edge_entry["edge_cladogenetic_speciation_mode"],
                        "child0_edge_id": edge_entry["child0_edge_id"],
                        "child1_edge_id": edge_entry["child1_edge_id"],
                            })

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

    def serialize_tables(self,
            output_prefix,
            delimiter="\t",
            output_suffix="tsv",
            ):
        # tree samples
        treef = open(output_prefix + ".tree-samples." + output_suffix, "w")
        writer = csv.DictWriter(treef,
                fieldnames=RevBayesBiogeographyParser.TREE_SAMPLES_FIELDNAMES,
                restval=RevBayesBiogeographyParser.NULL_VALUE,
                delimiter=delimiter,
                lineterminator=os.linesep,
                )
        writer.writeheader()
        writer.writerows(self.tree_entries)
        treef.close()

        # edge samples
        # events_by_treef = open(output_prefix + ".tree-events." + output_suffix, "w")
        edgef = open(output_prefix + ".edge-samples." + output_suffix, "w")
        writer = csv.DictWriter(edgef,
                fieldnames=RevBayesBiogeographyParser.EDGE_SAMPLES_FIELDNAMES,
                restval=RevBayesBiogeographyParser.NULL_VALUE,
                delimiter=delimiter,
                lineterminator=os.linesep,
                )
        writer.writeheader()
        writer.writerows(self.edge_entries)
        edgef.close()

        # event samples
        events_by_treef = open(output_prefix + ".event-samples." + output_suffix, "w")
        writer = csv.DictWriter(events_by_treef,
                fieldnames=RevBayesBiogeographyParser.EVENT_SAMPLES_FIELDNAMES,
                restval=RevBayesBiogeographyParser.NULL_VALUE,
                delimiter=delimiter,
                lineterminator=os.linesep,
                )
        writer.writeheader()
        writer.writerows(self.event_schedules_by_tree)
        events_by_treef.close()


def main():
    rbp = RevBayesBiogeographyParser()
    rbp.parse(sys.argv[1])
    rbp.serialize_tables(sys.argv[2])

if __name__ == "__main__":
    main()
