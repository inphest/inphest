#! /usr/bin/env python

try:
    from StringIO import StringIO # Python 2 legacy support: StringIO in this module is the one needed (not io)
except ImportError:
    from io import StringIO # Python 3
import sys
import os
import random
import collections
import argparse
import pprint
import copy
import json
from distutils.util import strtobool
import dendropy

from inphest import utility
from inphest import revbayes
from inphest import error

def weighted_choice(seq, weights, rng):
    """
    Selects an element out of seq, with probabilities of each element
    given by the list `weights` (which must be at least as long as the
    length of `seq` - 1).
    """
    if weights is None:
        weights = [1.0/len(seq) for count in range(len(seq))]
    else:
        weights = list(weights)
    if len(weights) < len(seq) - 1:
        raise Exception("Insufficient number of weights specified")
    sow = sum(weights)
    if len(weights) == len(seq) - 1:
        weights.append(1 - sow)
    return seq[weighted_index_choice(weights, sow, rng)]

def weighted_index_choice(weights, sum_of_weights, rng):
    """
    (From: http://eli.thegreenplace.net/2010/01/22/weighted-random-generation-in-python/)
    The following is a simple function to implement weighted random choice in
    Python. Given a list of weights, it returns an index randomly, according
    to these weights [1].
    For example, given [2, 3, 5] it returns 0 (the index of the first element)
    with probability 0.2, 1 with probability 0.3 and 2 with probability 0.5.
    The weights need not sum up to anything in particular, and can actually be
    arbitrary Python floating point numbers.
    If we manage to sort the weights in descending order before passing them
    to weighted_choice_sub, it will run even faster, since the random call
    returns a uniformly distributed value and larger chunks of the total
    weight will be skipped in the beginning.
    """
    rnd = rng.uniform(0, 1) * sum_of_weights
    for i, w in enumerate(weights):
        rnd -= w
        if rnd < 0:
            return i

class StatesVector(object):
    """
    A vector in which each element is an integer represents the state of a
    trait.

    E.g.,

        [1,0,1,2]

    is a 4-trait vector, where trait 0 is in state 1, trait 1 is in
    state 0, and so on.
    """

    def __init__(self,
            nchar,
            nstates=None,
            values=None,
            ):
        """
        Parameters
        ----------
        nchar : integer
            The number of traits to be tracked.
        nstates : list of integers
            The number of states for each trait. If not specified, defaults
            to binary (i.e., 2 states, 0 and 1). If specifed, must be a list of
            length `nchar`, with each element in the list being integer > 0.
        values : iterable of ints
            Vector of initial values. If not specified, defaults to all 0's.
        """
        self._nchar = nchar
        if nstates is not None:
            self._nstates = list(nstates)
        else:
            self._nstates = [2] * nchar
        if not values:
            self._states = [0] * nchar
        else:
            assert len(values) == nchar
            self._states = list(values)

    def clone(self):
        s = self.__class__(
                nchar=self._nchar,
                nstates=self._nstates,
            )
        s._states = list(self._states)
        return s

    @property
    def nchar(self):
        return len(self)

    def __len__(self):
        return self._nchar

    def __getitem__(self, idx):
        return self._states[idx]

    def __setitem__(self, idx, v):
        self._states[idx] = v

    def __repr__(self):
        return str(self._states)

class DistributionVector(StatesVector):

    @classmethod
    def from_string(cls, s):
        num_areas = len(s)
        values = [int(i) for i in s]
        return cls(num_areas=num_areas, values=values)

    def __init__(self, num_areas, values=None):
        StatesVector.__init__(self,
                nchar=num_areas,
                nstates=[2] * num_areas,
                values=values,
                )

    def presences(self):
        """
        Returns list of indexes in which lineage is present.
        """
        return [idx for idx, s in enumerate(self._states) if s == 1]

    def clone(self):
        s = self.__class__(num_areas=self._nchar)
        s._states = list(self._states)
        return s

class HostRegime(object):

    """
    The host history on which the symbiont history is conditioned.
    """

    HostLineage = collections.namedtuple("HostLineage", [
        # "tree_idx",                 #   identifer of tree from which this lineage has been sampled (same lineage, as given by split id will occur on different trees/histories)
        "lineage_id",               #   lineage (edge/split) id on which event occurs
        "lineage_start_distribution",  #   distribution/range (area set) at beginning of lineage
        "lineage_end_distribution",    #   distribution/range (area set) at end of lineage
    ])

    HostEvent = collections.namedtuple("HostEvent", [
        # "tree_idx",                 #   identifer of tree from which this event has been sampled
        "event_time",               #   time of event
        "weight",                   #   probability of event (1.0 if we take history as truth)
        "lineage_id",               #   lineage (edge/split) id on which event occurs
        "event_type",               #   type of event: anagenesis, cladogenesis
        "event_subtype",            #   if anagenesis: area_loss, area_gain; if cladogenesis: narrow sympatry etc.
        "area_idx",                 #   area involved in event (anagenetic)
        "child0_lineage_id",        #   split/edge id of first daughter (cladogenesis)
        "child1_lineage_id",        #   split/edge id of second daughter (cladogenesis)
        ])

    def __init__(self, taxon_namespace=None):
        if taxon_namespace is None:
            self.taxon_namespace = dendropy.TaxonNamespace()
        else:
            self.taxon_namespace = taxon_namespace
        self.next_event_index = 0
        self.events = [] # collection of HostEvent objects, sorted by time
        self.lineages = {} # keys: lineage_id (== int(Bipartition) == Bipartition.bitmask); values: HostLineage

    def compile(self):
        self.events.sort(key=lambda x: x.event_time)
        print(self.events)

class HostRegimeSamples(object):

    def __init__(self):
        self.host_regimes = []
        self.taxon_namespace = dendropy.TaxonNamespace()

    def parse_host_biogeography(self, src):
        """
        Reads the output of RevBayes biogeographical history.
        """
        rb = revbayes.RevBayesBiogeographyParser(taxon_namespace=self.taxon_namespace)
        rb.parse(src)

        # total_tree_ln_likelihoods = 0.0
        # for tree_entry in rb.tree_entries:
        #     total_tree_ln_likelihoods += tree_entry["posterior"]
        # for tree_entry in rb.tree_entries:
        #     self.tree_probabilities.append(tree_entry["posterior"]/total_tree_ln_likelihoods)

        tree_host_regimes = {}

        for edge_entry in rb.edge_entries:
            tree_idx = edge_entry["tree_idx"]
            if tree_idx not in tree_host_regimes:
                tree_host_regimes[tree_idx] = HostRegime(taxon_namespace=self.taxon_namespace)
            lineage_id = edge_entry["edge_id"]
            lineage = HostRegime.HostLineage(
                    # tree_idx=edge_entry["tree_idx"],
                    lineage_id=lineage_id,
                    lineage_start_distribution=DistributionVector.from_string(edge_entry["edge_starting_state"]),
                    lineage_end_distribution=DistributionVector.from_string(edge_entry["edge_ending_state"]),
                    )
            assert lineage.lineage_id not in tree_host_regimes[tree_idx].lineages
            tree_host_regimes[tree_idx].lineages[lineage_id] = lineage

        for event_entry in rb.event_schedules_by_tree:
            tree_idx = event_entry["tree_idx"]
            if tree_idx not in tree_host_regimes:
                tree_host_regimes[tree_idx] = HostRegime(taxon_namespace=self.taxon_namespace)
            event = HostRegime.HostEvent(
                # tree_idx=event_entry["tree_idx"],
                event_time=event_entry["age"],
                # weight=self.tree_probabilities[event_entry["tree_idx"]],
                weight=1.0,
                lineage_id=event_entry["edge_id"],
                event_type=event_entry["event_type"],
                event_subtype=event_entry["event_subtype"],
                area_idx=event_entry.get("area_idx", None),
                child0_lineage_id=event_entry.get("child0_edge_id", None),
                child1_lineage_id=event_entry.get("child1_edge_id", None),
                )
            assert event.lineage_id in tree_host_regimes[tree_idx].lineages
            tree_host_regimes[tree_idx].events.append(event)
        for host_regime in tree_host_regimes.values():
            host_regime.compile()
            self.host_regimes.append(host_regime)

if __name__ == "__main__":
    hh = HostRegimeSamples()
    rb_data = os.path.join(utility.TEST_DATA_PATH, "revbayes", "bg_large.events.txt")
    rb_data_src = open(rb_data, "r")
    hh.parse_host_biogeography(rb_data_src)
