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

class RateFunction(object):

    @classmethod
    def from_definition_dict(cls, rate_function_d):
        rf = cls()
        rf.parse_definition(rate_function_d)
        return rf

    def __init__(self,
            definition_type=None,
            definition_content=None,
            description=None,
            ):
        self.definition_type = definition_type # value, lambda, function, map
        self.definition_content = definition_content
        self.description = description
        self._compute_rate = None
        if self.definition_content is not None:
            self.compile_function()

    def __call__(self, lineage):
        return self._compute_rate(lineage)

    def parse_definition(self, rate_function_d):
        rate_function_d = dict(rate_function_d)
        self.definition_type = rate_function_d.pop("definition_type").replace("-", "_")
        self.definition_content = rate_function_d.pop("definition")
        self.description = rate_function_d.pop("description", "")
        if rate_function_d:
            raise TypeError("Unsupported function definition keywords: {}".format(rate_function_d))
        self.compile_function()

    def compile_function(self):
        self.definition_type = self.definition_type.replace("-", "_")
        if self.definition_type == "fixed_value":
            self.definition_content = float(self.definition_content)
            self._compute_rate = lambda lineage: self.definition_content
        elif self.definition_type == "lambda_definition":
            self._compute_rate = eval(self.definition_content)
        elif self.definition_type == "function_object":
            self._compute_rate = self.definition_content
        else:
            raise ValueError("Unrecognized function definition type: '{}'".format(self.definition_type))

    def as_definition(self):
        d = collections.OrderedDict()
        d["definition_type"] = self.definition_type
        if d["definition_type"] == "function_object":
            d["definition"] = str(self.definition_content)
        else:
            d["definition"] = self.definition_content
        d["description"] = self.description
        return d

class HostRegime(object):
    """
    A particular host history on which the symbiont history is conditioned.
    """

    HostRegimeLineageDefinition = collections.namedtuple("HostRegimeLineageDefinition", [
        # "tree_idx",                 #   identifer of tree from which this lineage has been sampled (same lineage, as given by split id will occur on different trees/histories)
        "lineage_id",               #   lineage (edge/split) id on which event occurs
        "leafset_bitstring",
        "split_bitstring",
        "rb_index",
        "lineage_start_time",          # time lineage appears
        "lineage_end_time",            # time lineage ends
        "lineage_start_distribution_bitstring",  #   distribution/range (area set) at beginning of lineage
        "lineage_end_distribution_bitstring",    #   distribution/range (area set) at end of lineage
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

    def __init__(self, taxon_namespace=None,):
        if taxon_namespace is None:
            self.taxon_namespace = dendropy.TaxonNamespace()
        else:
            self.taxon_namespace = taxon_namespace
        self.events = [] # collection of HostEvent objects, sorted by time
        self.lineages = {} # keys: lineage_id (== int(Bipartition) == Bipartition.bitmask); values: HostRegimeLineageDefinition
        self.start_time = None
        self.end_time = None

    def compile(self, start_time, end_time):
        self.start_time = start_time
        self.end_time = end_time
        self.events.sort(key=lambda x: x.event_time, reverse=False)
        for event in self.events:
            assert event.event_time >= self.start_time
            assert event.event_time <= self.end_time, "{} > {}".format(event.event_time, self.end_time)

    def validate(self):
        ## Basically, runs through the histories for each edge/lineage,
        ## ensuring the history correctly reproduces the end state given the start state
        lineage_events = collections.defaultdict(list)
        for event in self.events:
            lineage_events[event.lineage_id].append(event)
        for lineage_id in lineage_events:
            distribution_bitlist = list(self.lineages[lineage_id].lineage_start_distribution_bitstring)
            lineage_events[event.lineage_id].sort(key=lambda x: x.event_time, reverse=False)
            for event in lineage_events[lineage_id]:
                assert event.event_time >= self.lineages[lineage_id].lineage_start_time
                assert event.event_time <= self.lineages[lineage_id].lineage_end_time
                if event.event_type == "anagenesis" and event.event_subtype == "area_gain":
                    assert distribution_bitlist[event.area_idx] == "0", "Lineage {} at time {}: Trying to add area with index {} to distribution that already has area: {}".format(
                            lineage_id,
                            event.event_time,
                            event.area_idx,
                            "".join(distribution_bitlist))
                    distribution_bitlist[event.area_idx] == "1"
                elif event.event_type == "anagenesis" and event.event_subtype == "area_loss":
                    assert distribution_bitlist[event.area_idx] == "1", "Lineage {} at time {}: Trying to remove area with index {} from distribution that does not have area: {}".format(
                            lineage_id,
                            event.event_time,
                            event.area_idx,
                            "".join(distribution_bitlist))
                    distribution_bitlist[event.area_idx] == "0"
                elif event.event_type == "cladogenesis":
                    assert "".join(distribution_bitlist) == self.lineages[lineage_id].lineage_start_distribution_bitstring

class HostRegimeSamples(object):
    """
    A collection of host histories, one a single one of each a particular symbiont history will be conditioned.
    """

    def __init__(self):
        self.host_regimes = []
        self.taxon_namespace = dendropy.TaxonNamespace()

    def parse_host_biogeography(self,
            src,
            validate=True,
            ignore_validation_errors=False):
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
        # tree_root_heights = {}
        for edge_entry in rb.edge_entries:
            tree_idx = edge_entry["tree_idx"]
            if tree_idx not in tree_host_regimes:
                tree_host_regimes[tree_idx] = HostRegime(taxon_namespace=self.taxon_namespace)
            lineage_id = edge_entry["edge_id"]
            lineage = HostRegime.HostRegimeLineageDefinition(
                    # tree_idx=edge_entry["tree_idx"],
                    lineage_id=lineage_id,
                    leafset_bitstring=edge_entry["leafset_bitstring"],
                    split_bitstring=edge_entry["split_bitstring"],
                    rb_index=edge_entry["rb_index"],
                    lineage_start_time=edge_entry["edge_start_time"],
                    lineage_end_time=edge_entry["edge_end_time"],
                    lineage_start_distribution_bitstring=edge_entry["edge_starting_state"],
                    lineage_end_distribution_bitstring=edge_entry["edge_ending_state"],
                    )
            assert lineage.lineage_id not in tree_host_regimes[tree_idx].lineages
            assert lineage.lineage_start_time <= lineage.lineage_end_time, "{}, {}".format(lineage.lineage_start_time, lineage.lineage_end_time)
            tree_host_regimes[tree_idx].lineages[lineage_id] = lineage
            # try:
            #     tree_root_heights[tree_idx] = max(edge_entry["edge_ending_age"], tree_root_heights[tree_idx])
            # except KeyError:
            #     tree_root_heights[tree_idx] = edge_entry["edge_ending_age"]

        for event_entry in rb.event_schedules_across_all_trees:
            tree_idx = event_entry["tree_idx"]
            if tree_idx not in tree_host_regimes:
                tree_host_regimes[tree_idx] = HostRegime(taxon_namespace=self.taxon_namespace)
            event = HostRegime.HostEvent(
                # tree_idx=event_entry["tree_idx"],
                event_time=event_entry["time"],
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

        for tree_idx in tree_host_regimes:
            host_regime = tree_host_regimes[tree_idx]
            end_time = max(rb.tree_entries[tree_idx]["seed_node_age"], rb.max_event_times[tree_idx])
            host_regime.compile(
                    start_time=0.0,
                    end_time=end_time,
                    )
            if validate:
                host_regime.validate()
            self.host_regimes.append(host_regime)

class Area(object):

    """
    Manages the state of an area during a particular simulation replicate.
    """

    def __init__(self, area_idx):
        self.area_idx = area_idx
        self.host_lineages = set()
        self.symbiont_lineages = set()

    def __str__(self):
        return "Area{}".format(self.area_idx)

    def __repr__(self):
        return "<inphest.model.Area object at {} with index {}>".format(id(self), self.area_idx)

class HostLineage(object):
    """
    Manages the state of a host during a particular simulation replicate.
    """

    def __init__(self,
            host_regime_lineage_definition,
            host_system,
            debug_mode):
        self.host_regime_lineage_definition = host_regime_lineage_definition
        self.host_system = host_system
        self.lineage_id = host_regime_lineage_definition.lineage_id
        self.leafset_bitstring = host_regime_lineage_definition.leafset_bitstring
        self.split_bitstring = host_regime_lineage_definition.split_bitstring
        self.rb_index = host_regime_lineage_definition.rb_index
        self.start_time = host_regime_lineage_definition.lineage_start_time
        self.end_time = host_regime_lineage_definition.lineage_end_time
        self.start_distribution_bitstring = host_regime_lineage_definition.lineage_start_distribution_bitstring
        self.end_distribution_bitstring = host_regime_lineage_definition.lineage_end_distribution_bitstring
        self.current_distribution_bitstring = host_regime_lineage_definition.lineage_end_distribution_bitstring
        self._current_areas = set()
        assert len(self.host_system.areas) == len(self.start_distribution_bitstring)
        for (area_idx, area), (d_idx, presence) in zip(enumerate(self.host_system.areas), enumerate(self.start_distribution_bitstring)):
            assert area_idx == d_idx
            assert area.area_idx == area_idx
            if presence == "1":
                self.add_area(area)
            else:
                assert presence == "0"

    def add_area(self, area):
        assert area not in self._current_areas
        self._current_areas.add(area)
        area.host_lineages.add(self)

    def remove_area(self, area):
        assert area in self._current_areas
        self._current_areas.remove(area)
        area.host_lineages.remove(self)

    def has_area(self, area):
        return area in self._current_areas

    def current_area_iter(self):
        for area in self._current_areas:
            yield area

    def debug_check(self):
        for area in self.host_system.areas:
            if area in self._current_areas:
                assert self in area.host_lineages
            else:
                assert self not in area.host_lineages

    # def area_iter(self):
    #     """
    #     Iterate over areas in which this host occurs.
    #     """
    #     for area_idx in self.current_distribution_bitvector:
    #         if self.current_distribution_bitvector[area_idx] == 1:
    #             yield self.host_system.areas[area_idx]

    # def add_area(self, area_idx):
    #     self.current_distribution_bitvector[area_idx] = 1

    # def remove_area(self, area_idx):
    #     self.current_distribution_bitvector[area_idx] = 0

    # def has_area(self, area_idx):
    #     return self.current_distribution_bitvector[area_idx] == 1

class HostSystem(object):
    """
    Models the the collection of hosts for a particular simulation replicate,
    based on a HostRegime. The HostRegime provides the basic invariant
    definitions/rules that apply across all replicates, while the HostSystem
    tracks state etc. for a single replicate.
    """

    def __init__(self,
            host_regime,
            debug_mode=False,
            run_logger=None):
        self.debug_mode = debug_mode
        self.compile(host_regime, debug_mode=debug_mode)
        if self.debug_mode:
            ## setup for checking initial distributions ##
            self.check_lineage_distributions = {}
            for host_lineage in self.host_lineages:
                self.check_lineage_distributions[host_lineage] = list(host_lineage.start_distribution_bitstring)
                self.debug_check_initial_distribution(host_lineage)
        self._next_host_event = None

    def compile(self, host_regime, debug_mode=False):
        self.host_regime = host_regime
        self.start_time = self.host_regime.start_time
        self.end_time = self.host_regime.end_time
        num_areas = None

        # build areas
        for host_regime_lineage_id_definition in self.host_regime.lineages.values():
            if num_areas is None:
                num_areas = len(host_regime_lineage_id_definition.lineage_start_distribution_bitstring)
            assert num_areas == len(host_regime_lineage_id_definition.lineage_start_distribution_bitstring)
            assert num_areas == len(host_regime_lineage_id_definition.lineage_end_distribution_bitstring)
        self.num_areas = num_areas
        self.areas = []
        for area_idx in range(self.num_areas):
            area = Area(area_idx)
            self.areas.append(area)

        # self.area_host_symbiont_host_area_distribution = {}
        # for area in self.areas:
        #     self.area_host_symbiont_host_area_distribution[area] = {}
        #     for host_lineage in self.host_lineages_by_id.values():
        #         self.area_host_symbiont_host_area_distribution[area][host_lineage] = {}

        # compile lineages
        self.host_lineages = set()
        self.host_lineages_by_id = {}
        for host_regime_lineage_id_definition in self.host_regime.lineages.values():
            host = HostLineage(
                    host_regime_lineage_definition=host_regime_lineage_id_definition,
                    host_system=self,
                    debug_mode=debug_mode,
                    )
            self.host_lineages.add(host)
            self.host_lineages_by_id[host.lineage_id] = host

        # local copy of host events
        self.host_events = list(self.host_regime.events)

    def extant_host_lineages_at_current_time(self, current_time):
        ## TODO: if we hit this often, we need to construct a look-up table
        lineages = set()
        for host in self.host_lineages_by_id.values():
            if host.start_time <= current_time and host.end_time > current_time:
                lineages.add(host)
        return lineages

    def debug_check(self):
        for host_lineage in self.host_lineages:
            host_lineage.debug_check()
            for area in host_lineage._current_areas:
                assert area in self.areas
        for event in self.host_events:
            area_idx = event.area_idx
            if area_idx is None:
                continue
            assert self.areas[area_idx].area_idx == area_idx

    def debug_check_initial_distribution(self, host_lineage):
        for area_idx, presence in enumerate(self.check_lineage_distributions[host_lineage]):
            if presence == "1":
                assert self.areas[area_idx] in host_lineage._current_areas
            elif presence == "0":
                assert self.areas[area_idx] not in host_lineage._current_areas
            else:
                raise ValueError


class SymbiontHostAreaDistribution(object):
    """
    Manages the host-by-area distribution of a single symbiont lineage.
    """

    class NullDistributionException(Exception):
        pass

    def __init__(self, symbiont_lineage, host_system):
        self.symbiont_lineage = symbiont_lineage
        self.host_system = host_system

        ## distribution set: tuples of (host_lineage_id, area_idx)
        # self._host_area_distribution = set()

        self._host_area_distribution = {}
        for host_lineage in self.host_system.host_lineages:
            self._host_area_distribution[host_lineage] = {}
            for area in self.host_system.areas:
                self._host_area_distribution[host_lineage][area] = 0

        ## For quick look-up if host/area is infected
        self._infected_hosts = set()
        self._infected_areas = set()

    def add_host_in_area(self, host_lineage, area=None):
        """
        Adds a host to the distribution.
        If ``area`` is specified, then only the host in a specific area is infected.
        Otherwise, all hosts (of the given lineage) in all areas are infected.
        """
        if area is None:
            for area in host_lineage.current_area_iter():
                self._host_area_distribution[host_lineage][area] = 1
                area.symbiont_lineages.add(self.symbiont_lineage)
                self._infected_areas.add(area)
        else:
            assert host_lineage.has_area(area)
            self._host_area_distribution[host_lineage][area] = 1
            area.symbiont_lineages.add(self.symbiont_lineage)
        self._infected_hosts.add(host_lineage)

    def remove_host_in_area(self, host_lineage, area=None):
        """
        Removes a host from the distribution.
        If ``area_idx`` is specified, then only the host in that specific area is removed. Otherwise,
        Otherwise, all hosts (of the given lineage) of all areas are removed from the range.
        """
        if area is None:
            self.remove_host(host_lineage)
        else:
            assert host_lineage.has_area(area), "{} not in host area: {}".format(area, host_lineage._current_areas)
            self._host_area_distribution[host_lineage][area] = 0
            self.sync_area_cache(area)
            self.sync_host_cache(host_lineage)
            self.check_for_null_distribution()

    def remove_host(self, host_lineage):
        """
        Removes association with host from all areas.
        """
        for area in host_lineage.current_area_iter():
            if self._host_area_distribution[host_lineage][area] > 0:
                self._host_area_distribution[host_lineage][area] = 0
                self.sync_area_cache(area)
        self._infected_hosts.remove(host_lineage)
        self.check_for_null_distribution()

    def sync_area_cache(self, area, search_all_hosts=False):
        """
        Check if symbiont occurs in any host in the given area.
        If not, remove association with area.
        Otherwise, add it.
        If ``search_all_hosts`` is True, then all hosts in the system will be searched.
        Otherwise, only hosts known to be associated with this area will be searched.
        """
        if search_all_hosts:
            host_lineage_iter = iter(self.host_system.host_lineages)
        else:
            host_lineage_iter = iter(area.host_lineages)
        for host_lineage in area.host_lineages:
            if self._host_area_distribution[host_lineage][area] == 1:
                self._infected_areas.add(area)
                area.symbiont_lineages.add(self.symbiont_lineage)
                break
        else:
            self._infected_areas.remove(area)
            area.symbiont_lineages.remove(self.symbiont_lineage)

    def sync_host_cache(self, host_lineage, search_all_areas=False):
        """
        Check if symbiont occurs in the given host in any of its areas.
        If not, remove association with host.
        Otherwise, add it.
        If ``search_all_areas`` is True, then all areas in the system will be searched.
        Otherwise, only areas known to be associated with this host will be searched.
        """
        if search_all_areas:
            area_iter = iter(self.host_system.areas)
        else:
            area_iter = host_lineage.current_area_iter()
        for area in area_iter:
            if self._host_area_distribution[host_lineage][area] == 1:
                self._infected_hosts.add(host_lineage)
                break
        else:
            self._infected_hosts.remove(host_lineage)

    def check_for_null_distribution(self):
        """
        Ensures that lineage occurs at least in one host in one area.
        """
        if not self._infected_hosts or not self._infected_areas:
            raise SymbiontHostAreaDistribution.NullDistributionException()

    def has_host(self, host_lineage):
        """
        Returns True if host is infected in any of its areas.
        """
        return host_lineage in self._infected_hosts

    def has_host_in_area(self, host_lineage, area):
        """
        Returns True if host is infected in a particular area.
        """
        return self._host_area_distribution[host_lineage][area] == 1

    def host_iter(self):
        """
        Iterates over hosts in which lineage occurs.
        """
        for host in self._infected_hosts:
            yield host

    def area_iter(self):
        """
        Iterates over areas in which lineage occurs.
        """
        for area in self._infected_areas:
            yield area

    def areas_in_host_iter(self, host_lineage):
        """
        Iterates over areas in which lineage is associated with a particular host.
        """
        for area in self._host_area_distribution[host_lineage]:
            if self._host_area_distribution[host_lineage][area] == 1:
                yield area
            else:
                assert self._host_area_distribution[host_lineage][area] == 0

    def has_area(self, area):
        """
        Returns True if area is infected.
        """
        return area in self._infected_areas

    def clone(self, symbiont_lineage):
        """
        Returns deep-copy of this distribution.
        """
        d = self.__class__(
                symbiont_lineage=symbiont_lineage,
                host_system=self.host_system)
        for host_lineage in self.host_system.host_lineages:
            for area in self.host_system.areas:
                d._host_area_distribution[host_lineage][area] = d._host_area_distribution[host_lineage][area]


    def debug_check(self, simulation_elapsed_time=None):
        # check that, as an extant lineage, it occupies at least
        # one host/area
        infected_hosts = set()
        infected_areas = set()
        noninfected_hosts = set()
        noninfected_areas = set([area for area in self.host_system.areas])
        occurrences = 0
        for host_lineage in self.host_system.host_lineages:
            for area in self.host_system.areas:
                if self._host_area_distribution[host_lineage][area] == 1:
                    infected_hosts.add(host_lineage)
                    infected_areas.add(area)
                    assert host_lineage in area.host_lineages
                    assert self.symbiont_lineage in area.symbiont_lineages
                    noninfected_areas.discard(area)
                    noninfected_hosts.discard(host_lineage)
                    occurrences += 1
                elif self._host_area_distribution[host_lineage][area] != 0:
                    raise ValueError(self._host_area_distribution[host_lineage][area])
        assert occurrences > 0
        # check that the caches are in sync
        assert infected_hosts == self._infected_hosts
        assert infected_areas == self._infected_areas
        for area in noninfected_areas:
            assert self.symbiont_lineage not in area.symbiont_lineages
        # check that the infected hosts are supposed to exist at the current time
        if simulation_elapsed_time is not None:
            for host_lineage in self._infected_hosts:
                try:
                    assert simulation_elapsed_time >= host_lineage.start_time
                    assert simulation_elapsed_time < host_lineage.end_time
                except AssertionError:
                    # until we finish properly processing host events ...
                    print("!! IGNORING EXTINCT HOST CHECK ERROR: {} ({}): times = {} to {}, current simulation elapsed time = {}".format(
                        host_lineage.lineage_id,
                        host_lineage.leafset_bitstring,
                        host_lineage.start_time,
                        host_lineage.end_time,
                        simulation_elapsed_time))

class SymbiontLineage(dendropy.Node):

    def __init__(self, index, host_system):
        dendropy.Node.__init__(self)
        self.host_system = host_system
        self.host_area_distribution = SymbiontHostAreaDistribution(
                symbiont_lineage=self,
                host_system=self.host_system)
        self.host_area_distribution.lineage = self
        self.index = index
        self.is_extant = True
        self.edge.length = 0

class SymbiontPhylogeny(dendropy.Tree):

    def node_factory(cls, **kwargs):
        return SymbiontLineage(**kwargs)
    node_factory = classmethod(node_factory)

    def __init__(self, *args, **kwargs):
        if kwargs:
            self.model = kwargs.pop("model")
            self.model_id = self.model.model_id
            self.host_system = kwargs.pop("host_system")
            self.annotations.add_bound_attribute("model_id")
            self.rng = kwargs.pop("rng")
            self.debug_mode = kwargs.pop("debug_mode")
            self.run_logger = kwargs.pop("run_logger")
            self.lineage_indexer = utility.IndexGenerator(0)
            if "seed_node" not in kwargs:
                seed_node = self.node_factory(
                        index=next(self.lineage_indexer),
                        host_system=self.host_system,
                        )
                kwargs["seed_node"] = seed_node
            dendropy.Tree.__init__(self, *args, **kwargs)
            self.current_lineages = set([self.seed_node])
            extant_host_lineages = self.host_system.extant_host_lineages_at_current_time(0)
            for host_lineage in extant_host_lineages:
                self.seed_node.host_area_distribution.add_host_in_area(host_lineage=host_lineage)
        else:
            dendropy.Tree.__init__(self, *args, **kwargs)

    def __deepcopy__(self, memo=None):
        if memo is None:
            memo = {}
        memo[id(self.model)] = self.model
        memo[id(self.rng)] = None #self.rng
        memo[id(self.run_logger)] = self.run_logger
        memo[id(self.taxon_namespace)] = self.taxon_namespace
        return dendropy.Tree.__deepcopy__(self, memo)

    def current_lineage_iter(self):
        for lineage in self.current_lineages:
            yield lineage

    def split_lineage(self, lineage):
        pass
        # c1 = self.node_factory(
        #         index=next(self.lineage_indexer),
        #         host_area_distribution=dist1,
        #         )
        # c2 = self.node_factory(
        #         index=next(self.lineage_indexer),
        #         host_area_distribution=dist2,
        #         )

        # if self.debug_mode:
        #     self.run_logger.debug("Splitting {} with distribution {} under speciation mode {} to: {} (distribution: {}) and {} (distribution: {})".format(
        #         lineage,
        #         lineage.distribution_vector.presences(),
        #         speciation_mode,
        #         c1,
        #         dist1.presences(),
        #         c2,
        #         dist2.presences(),
        #         ))
        #     assert len(dist1.presences()) > 0
        #     assert len(dist2.presences()) > 0

        # lineage.is_extant = False
        # self.current_lineages.remove(lineage)
        # lineage.add_child(c1)
        # lineage.add_child(c2)
        # self.current_lineages.add(c1)
        # self.current_lineages.add(c2)

    def extinguish_lineage(self, lineage):
        self._make_lineage_extinct_on_phylogeny(lineage)

    def expand_lineage_host_set(self, lineage):
        pass

    def contract_lineage_host_set(self, lineage):
        pass

    def expand_lineage_area_set(self, lineage):
        pass

    def contract_lineage_area_set(self, lineage):
        pass

    def _make_lineage_extinct_on_phylogeny(self, lineage):
        if len(self.current_lineages) == 1:
            self.total_extinction_exception("no extant lineages remaining")
        lineage.is_extant = False
        self.current_lineages.remove(lineage)
        self.prune_subtree(lineage)

    def total_extinction_exception(self, msg):
        # self.run_logger.info("Total extinction: {}".format(msg))
        raise error.TotalExtinctionException(msg)

    # def evolve_trait(self, lineage, trait_idx, state_idx):
    #     lineage.traits_vector[trait_idx] = state_idx

    def disperse_lineage(self, lineage, dest_area_idx):
        lineage.distribution_vector[dest_area_idx] = 1

    # def focal_area_lineages(self):
    #     focal_area_lineages = set()
    #     for lineage in self.current_lineage_iter():
    #         for area_idx in self.model.geography.focal_area_indexes:
    #             if lineage.distribution_vector[area_idx] == 1:
    #                 focal_area_lineages.add(lineage)
    #                 break
    #     return focal_area_lineages

    # def num_focal_area_lineages(self):
    #     count = 0
    #     for lineage in self.current_lineage_iter():
    #         for area_idx in self.model.geography.focal_area_indexes:
    #             if lineage.distribution_vector[area_idx] == 1:
    #                 count += 1
    #                 break
    #     return count

    # def extract_focal_areas_tree(self):
    #     # tcopy = SymbiontPhylogeny(self)
    #     tcopy = copy.deepcopy(self)
    #     focal_area_lineages = tcopy.focal_area_lineages()
    #     if len(focal_area_lineages) < 2:
    #         raise error.InsufficientFocalAreaLineagesSimulationException("insufficient lineages in focal area at termination".format(len(focal_area_lineages)))
    #     try:
    #         tcopy.filter_leaf_nodes(filter_fn=lambda x: x in focal_area_lineages)
    #     except dendropy.SeedNodeDeletionException:
    #         raise error.InsufficientFocalAreaLineagesSimulationException("no extant lineages in focal area at termination".format(len(focal_area_lineages)))
    #     return tcopy

class InphestModel(object):

    _TRAITS_SEPARATOR = "."
    _LABEL_COMPONENTS_SEPARATOR = "^"
    _NULL_TRAITS = "NA"

    @classmethod
    def create(
            cls,
            host_regime,
            model_definition_source,
            model_definition_type,
            interpolate_missing_model_values=False,
            run_logger=None,
            ):
        """
        Create and return a model under which to run a simulation.

        Parameters
        ----------
        model_definition_source : object
            See 'model_definition_type' argument for values this can take.
        model_definition_type : str
            Whether 'model_definition_source' is:

                - 'python-dict' : a Python dictionary defining the model.
                - 'python-dict-str' : a string providing a Python dictionary
                  defining the model.
                - 'python-dict-filepath' : a path to a Python file to be evaluated;
                  the file should be a valid Python script containing nothing but a
                  dictionary defining the model.
                - 'json-filepath': a path to a JSON file containing a dictionary
                  defining the model.

        Returns
        -------
        m : ArchipelagoModel
            A fully-specified Archipelago model.

        Example
        -------

            hrs = HostRegimeSamples()
            rb_data = os.path.join(utility.TEST_DATA_PATH, "revbayes", "bg_large.events.txt")
            rb_data_src = open(rb_data, "r")
            hrs.parse_host_biogeography(rb_data_src)
            for host_regime in hrs.host_regimes:
                im = InphestModel.create(
                        host_regime=host_regime,
                        model_definition_source="model1.json",
                        model_definition_type="json",
                        )

        """
        if model_definition_type == "python-dict-filepath":
            src = open(model_definition_source, "r")
            model_definition = eval(src.read())
        elif model_definition_type == "python-dict-str":
            model_definition = eval(model_definition_source)
        elif model_definition_type == "python-dict":
            model_definition = model_definition_source
        elif model_definition_type == "json-filepath":
            src = open(model_definition_source, "r")
            model_definition = json.load(src)
        else:
            raise ValueError("Unrecognized model definition type: '{}'".format(model_definition_type))
        return cls.from_definition_dict(
                host_regime=host_regime,
                model_definition=model_definition,
                run_logger=run_logger,
                interpolate_missing_model_values=interpolate_missing_model_values)

    @classmethod
    def from_definition_dict(cls,
            host_regime,
            model_definition,
            interpolate_missing_model_values=False,
            run_logger=None):
        archipelago_model = cls()
        archipelago_model.parse_definition(
                model_definition=model_definition,
                host_regime=host_regime,
                interpolate_missing_model_values=interpolate_missing_model_values,
                run_logger=run_logger,
        )
        return archipelago_model

    @staticmethod
    def compose_encoded_label(symbiont_lineage):
        return "s{lineage_index}".format(lineage_index=symbiont_lineage.index)
        return encoding

    @staticmethod
    def decode_label(label):
        raise NotImplementedError
        # parts = label.split(ArchipelagoModel._LABEL_COMPONENTS_SEPARATOR)
        # traits_string = parts[1]
        # if not traits_string or traits_string == ArchipelagoModel._NULL_TRAITS:
        #     traits_vector = StatesVector(nchar=0)
        # else:
        #     traits_string_parts = traits_string.split(ArchipelagoModel._TRAITS_SEPARATOR)
        #     traits_vector = StatesVector(
        #             nchar=len(traits_string_parts),
        #             # The trait states need to be an integer if
        #             # archipelago-summarize.py coerces the user input to
        #             # integers
        #             # values=[int(i) for i in traits_string_parts],
        #             # The reason we do NOT want it parsed to an integer value
        #             # is to allow null traits 'NA', 'null', etc.
        #             values=[i for i in traits_string_parts],
        #             )
        # distribution_string = parts[2]
        # distribution_vector = DistributionVector(
        #         num_areas=len(distribution_string),
        #         values=[int(i) for i in distribution_string],)
        # return traits_vector, distribution_vector

    @staticmethod
    def set_lineage_data(
            tree,
            leaf_nodes_only=False,
            lineage_data_source="node",
            traits_filepath=None,
            areas_filepath=None,
            ):
        raise NotImplementedError
        # if lineage_data_source == "node":
        #     _decode = lambda x: ArchipelagoModel.decode_label(x.label)
        # elif lineage_data_source == "taxon":
        #     _decode = lambda x: ArchipelagoModel.decode_label(x.taxon.label)
        # else:
        #     raise ValueError("'lineage_data_source' must be 'node' or 'taxon'")
        # for nd in tree:
        #     if (not leaf_nodes_only or not nd._child_nodes) and (lineage_data_source == "node" or nd.taxon is not None):
        #         traits_vector, distribution_vector = _decode(nd)
        #         nd.traits_vector = traits_vector
        #         nd.distribution_vector = distribution_vector
        #     else:
        #         nd.traits_vector = None
        #         nd.distribution_vector = None

    def __init__(self):
        pass

    def parse_definition(self,
            model_definition,
            host_regime,
            run_logger=None,
            interpolate_missing_model_values=True):

        # initialize
        if model_definition is None:
            model_definition = {}
        else:
            model_definition = dict(model_definition)

        # model identification
        if "model_id" not in model_definition:
            model_definition["model_id"] = "Model1"
            if run_logger is not None:
                run_logger.warning("Model identifier not specified: defaulting to '{}'".format(model_definition["model_id"]))
        self.model_id = model_definition.pop("model_id", "Model1")
        if run_logger is not None:
            run_logger.info("Setting up model with identifier: '{}'".format(self.model_id))

        # host regime
        self.host_regime = host_regime

        # Diversification

        ## speciation
        diversification_d = dict(model_definition.pop("diversification", {}))
        if "symbiont_lineage_birth_rate" in diversification_d:
            self.symbiont_lineage_birth_rate_function = RateFunction.from_definition_dict(diversification_d.pop("symbiont_lineage_birth_rate"))
        else:
            self.symbiont_lineage_birth_rate_function = RateFunction(
                    definition_type="lambda_definition",
                    definition_content="lambda symbiont_lineage: 0.01",
                    description="fixed: 0.01",
                    )
        if run_logger is not None:
            run_logger.info("(DIVERSIFICATION) Setting symbiont lineage-specific birth rate function: {desc}".format(
                desc=self.symbiont_lineage_birth_rate_function.description,))

        ## extinction
        if "symbiont_lineage_death_rate" in diversification_d:
            self.symbiont_lineage_death_rate_function = RateFunction.from_definition_dict(diversification_d.pop("symbiont_lineage_death_rate"))
        else:
            self.symbiont_lineage_death_rate_function = RateFunction(
                    definition_type="lambda_definition",
                    definition_content="lambda symbiont_lineage: 0.0",
                    description="fixed: 0.0",
                    )
        if run_logger is not None:
            run_logger.info("(DIVERSIFICATION) Setting symbiont lineage-specific death rate function: {desc}".format(
                desc=self.symbiont_lineage_death_rate_function.description,))
        if diversification_d:
            raise TypeError("Unsupported diversification model keywords: {}".format(diversification_d))

        # Host Submodel

        ## Anagenetic Host Evolution Submodel

        anagenetic_host_range_evolution_d = dict(model_definition.pop("anagenetic_host_range_evolution", {}))

        ### Anagenetic Host Gain

        if "symbiont_lineage_host_gain_rate" in anagenetic_host_range_evolution_d:
            self.symbiont_lineage_host_gain_rate_function = RateFunction.from_definition_dict(anagenetic_host_range_evolution_d.pop("symbiont_lineage_host_gain_rate"))
        else:
            self.symbiont_lineage_host_gain_rate_function = RateFunction(
                    definition_type="lambda_definition",
                    definition_content="lambda symbiont_lineage: 0.01",
                    description="fixed: 0.01",
                    )
        if run_logger is not None:
            run_logger.info("(ANAGENETIC HOST RANGE EVOLUTION) Setting symbiont lineage-specific host gain weight function: {desc}".format(
                desc=self.symbiont_lineage_host_gain_rate_function.description,))

        ### Anagenetic Host Loss

        if "symbiont_lineage_host_loss_rate" in anagenetic_host_range_evolution_d:
            self.symbiont_lineage_host_loss_rate_function = RateFunction.from_definition_dict(anagenetic_host_range_evolution_d.pop("symbiont_lineage_host_loss_rate"))
        else:
            self.symbiont_lineage_host_loss_rate_function = RateFunction(
                    definition_type="lambda_definition",
                    definition_content="lambda symbiont_lineage: 0.0",
                    description="fixed: 0.0",
                    )
        if run_logger is not None:
            run_logger.info("(ANAGENETIC HOST RANGE EVOLUTION) Setting symbiont lineage-specific host loss weight function: {desc}".format(
                desc=self.symbiont_lineage_host_loss_rate_function.description,
                ))

        if anagenetic_host_range_evolution_d:
            raise TypeError("Unsupported keywords in anagenetic host range evolution submodel: {}".format(anagenetic_host_range_evolution_d))

        ## Cladogenetic Host Evolution Submodel

        cladogenetic_host_range_evolution = dict(model_definition.pop("cladogenetic_host_range_evolution", {}))
        self.cladogenesis_sympatric_subset_speciation_weight = float(cladogenetic_host_range_evolution.pop("sympatric_subset_speciation_weight", 1.0))
        self.cladogenesis_single_host_vicariance_speciation_weight = float(cladogenetic_host_range_evolution.pop("single_host_vicariance_speciation_weight", 1.0))
        self.cladogenesis_widespread_vicariance_speciation_weight = float(cladogenetic_host_range_evolution.pop("widespread_vicariance_speciation_weight", 1.0))
        self.cladogenesis_founder_event_speciation_weight = float(cladogenetic_host_range_evolution.pop("founder_event_speciation_weight", 0.0))
        if cladogenetic_host_range_evolution:
            raise TypeError("Unsupported keywords in cladogenetic range evolution submodel: {}".format(cladogenetic_host_range_evolution))
        if run_logger is not None:
            run_logger.info("(CLADOGENETIC HOST RANGE EVOLUTION) Base weight of sympatric subset speciation mode: {}".format(self.cladogenesis_sympatric_subset_speciation_weight))
            run_logger.info("(CLADOGENETIC HOST RANGE EVOLUTION) Base weight of single host vicariance speciation mode: {}".format(self.cladogenesis_single_host_vicariance_speciation_weight))
            run_logger.info("(CLADOGENETIC HOST RANGE EVOLUTION) Base weight of widespread vicariance speciation mode: {}".format(self.cladogenesis_widespread_vicariance_speciation_weight))
            run_logger.info("(CLADOGENETIC HOST RANGE EVOLUTION) Base weight of founder event speciation ('jump dispersal') mode: {} (note that the effective weight of this event for each lineage is actually the product of this and the lineage-specific host gain weight)".format(self.cladogenesis_founder_event_speciation_weight))

        # Geographical Range Evolution Submodel

        ## Anagenetic Geographical Range Evolution Submodel

        ### Anagenetic Geographical Area Gain

        anagenetic_geographical_range_evolution_d = dict(model_definition.pop("anagenetic_geographical_range_evolution", {}))
        if "lineage_area_gain_rate" in anagenetic_geographical_range_evolution_d:
            self.symbiont_lineage_area_gain_rate_function = RateFunction.from_definition_dict(anagenetic_geographical_range_evolution_d.pop("lineage_area_gain_rate"), self.trait_types)
        else:
            self.symbiont_lineage_area_gain_rate_function = RateFunction(
                    definition_type="lambda_definition",
                    definition_content="lambda lineage: 0.01",
                    description="fixed: 0.01",
                    )
        if run_logger is not None:
            run_logger.info("(ANAGENETIC GEOGRAPHICAL RANGE EVOLUTION) Setting symbiont lineage-specific area gain weight function: {desc}".format(
                desc=self.symbiont_lineage_area_gain_rate_function.description,))

        ### Anagenetic Geographical Area Gain

        if "lineage_area_loss_rate" in anagenetic_geographical_range_evolution_d:
            self.symbiont_lineage_area_loss_rate_function = RateFunction.from_definition_dict(anagenetic_geographical_range_evolution_d.pop("lineage_area_loss_rate"), self.trait_types)
        else:
            self.symbiont_lineage_area_loss_rate_function = RateFunction(
                    definition_type="lambda_definition",
                    definition_content="lambda lineage: 0.0",
                    description="fixed: 0.0",
                    )
        if run_logger is not None:
            run_logger.info("(ANAGENETIC GEOGRAPHICAL RANGE EVOLUTION) Setting symbiont lineage-specific area loss weight function: {desc}".format(
                desc=self.symbiont_lineage_area_loss_rate_function.description,
                ))

        if anagenetic_geographical_range_evolution_d:
            raise TypeError("Unsupported keywords in anagenetic geographical range evolution submodel: {}".format(anagenetic_geographical_range_evolution_d))

        ## Cladogenetic Geographical Evolution Submodel

        cladogenesis_geographical_range_evolution_d = dict(model_definition.pop("cladogenetic_geographical_range_evolution", {}))
        self.cladogenesis_sympatric_subset_speciation_weight = float(cladogenesis_geographical_range_evolution_d.pop("sympatric_subset_speciation_weight", 1.0))
        self.cladogenesis_single_area_vicariance_speciation_weight = float(cladogenesis_geographical_range_evolution_d.pop("single_area_vicariance_speciation_weight", 1.0))
        self.cladogenesis_widespread_vicariance_speciation_weight = float(cladogenesis_geographical_range_evolution_d.pop("widespread_vicariance_speciation_weight", 1.0))
        self.cladogenesis_founder_event_speciation_weight = float(cladogenesis_geographical_range_evolution_d.pop("founder_event_speciation_weight", 0.0))
        if cladogenesis_geographical_range_evolution_d:
            raise TypeError("Unsupported keywords in cladogenetic geographical range evolution submodel: {}".format(cladogenesis_geographical_range_evolution_d))
        if run_logger is not None:
            run_logger.info("(CLADOGENETIC GEOGRAPHICAL RANGE EVOLUTION) Base weight of sympatric subset speciation mode: {}".format(self.cladogenesis_sympatric_subset_speciation_weight))
            run_logger.info("(CLADOGENETIC GEOGRAPHICAL RANGE EVOLUTION) Base weight of single area vicariance speciation mode: {}".format(self.cladogenesis_single_area_vicariance_speciation_weight))
            run_logger.info("(CLADOGENETIC GEOGRAPHICAL RANGE EVOLUTION) Base weight of widespread vicariance speciation mode: {}".format(self.cladogenesis_widespread_vicariance_speciation_weight))
            run_logger.info("(CLADOGENETIC GEOGRAPHICAL RANGE EVOLUTION) Base weight of founder event speciation ('jump dispersal') mode: {} (note that the effective weight of this event for each lineage is actually the product of this and the lineage-specific area gain weight)".format(self.cladogenesis_founder_event_speciation_weight))

        termination_conditions_d = dict(model_definition.pop("termination_conditions", {}))
        self.max_time = termination_conditions_d.pop("max_time", self.host_regime.end_time)
        desc = "Simulation will terminate at time t = {}".format(self.max_time)
        if run_logger is not None:
            run_logger.info(desc)

        if model_definition:
            raise TypeError("Unsupported model keywords: {}".format(model_definition))

    def encode_lineage(self,
            symbiont_lineage,
            set_label=False,
            add_annotation=False,
            ):
        encoded_label = InphestModel.compose_encoded_label(symbiont_lineage=symbiont_lineage)
        if set_label:
            lineage.label =encoded_label
        if add_annotation:
            lineage.annotations.drop()
            lineage.annotations.add_new("traits_v", traits_v)
            lineage.annotations.add_new("distribution", areas_v)
            for trait_idx, trait in enumerate(self.trait_types):
                lineage.annotations.add_new(trait.label, lineage.traits_vector[trait_idx])
            area_list = []
            for area_idx, area in enumerate(self.geography.areas):
                if exclude_supplemental_areas and area.is_supplemental:
                    continue
                if lineage.distribution_vector[area_idx] == 1:
                    area_list.append(area.label)
            lineage.annotations.add_new("areas", area_list)
        return encoded_label

    def write_model(self, out):
        model_definition = collections.OrderedDict()
        model_definition["model_id"] = self.model_id
        model_definition["areas"] = self.geography.as_definition()
        model_definition["traits"] = self.trait_types.as_definition()
        model_definition["diversification"] = self.diversification_as_definition()
        model_definition["anagenetic_range_evolution"] = self.anagenetic_range_evolution_as_definition()
        model_definition["cladogenetic_range_evolution"] = self.cladogenetic_range_evolution_as_definition()
        model_definition["termination_conditions"] = self.termination_conditions_as_definition()
        json.dump(model_definition, out, indent=4, separators=(',', ': '))

    def diversification_as_definition(self):
        d = collections.OrderedDict()
        d["lineage_birth_rate"] = self.lineage_birth_rate_function.as_definition()
        d["lineage_death_rate"] = self.lineage_death_rate_function.as_definition()
        return d

    def anagenetic_range_evolution_as_definition(self):
        d = collections.OrderedDict()
        # if self.global_area_gain_rate is not None and self.mean_area_gain_rate is not None:
        #     raise TypeError("Both 'global_area_gain_rate' and 'mean_area_gain_rate' are populated")
        # elif self.global_area_gain_rate is None and self.mean_area_gain_rate is None:
        #     raise TypeError("Neither 'global_area_gain_rate' and 'mean_area_gain_rate' are populated")
        # elif self.global_area_gain_rate is not None:
        #     d["global_area_gain_rate"] = self.global_area_gain_rate
        # else:
        #     d["mean_area_gain_rate"] = self.mean_area_gain_rate
        d["lineage_area_gain_rate"] = self.lineage_area_gain_rate_function.as_definition()
        d["lineage_area_loss_rate"] = self.lineage_area_loss_rate_function.as_definition()
        return d

    def cladogenetic_range_evolution_as_definition(self):
        d = collections.OrderedDict()
        d["sympatric_subset_speciation_weight"] = self.cladogenesis_sympatric_subset_speciation_weight
        d["single_area_vicariance_speciation_weight"] = self.cladogenesis_single_area_vicariance_speciation_weight
        d["widespread_vicariance_speciation_weight"] = self.cladogenesis_widespread_vicariance_speciation_weight
        d["founder_event_speciation_weight"] = self.cladogenesis_founder_event_speciation_weight
        return d

    def termination_conditions_as_definition(self):
        d = collections.OrderedDict()
        d["target_focal_area_lineages"] = self.target_focal_area_lineages
        d["gsa_termination_focal_area_lineages"] = self.gsa_termination_focal_area_lineages
        d["max_time"] = self.max_time
        return d



if __name__ == "__main__":
    hrs = HostRegimeSamples()
    rb_data = os.path.join(utility.TEST_DATA_PATH, "revbayes", "bg_large.events.txt")
    rb_data_src = open(rb_data, "r")
    hrs.parse_host_biogeography(rb_data_src)
    for host_regime in hrs.host_regimes:
        im = InphestModel.create(
                host_regime=host_regime,
                model_definition_source={},
                model_definition_type="python-dict",
                )
