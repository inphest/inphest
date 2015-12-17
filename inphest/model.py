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

class HostRegime(object):
    """
    A particular host history on which the symbiont history is conditioned.
    """

    HostRegimeLineageDefinition = collections.namedtuple("HostRegimeLineageDefinition", [
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

    class HostDistributionVector(StatesVector):

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

    def __init__(self, taxon_namespace=None):
        if taxon_namespace is None:
            self.taxon_namespace = dendropy.TaxonNamespace()
        else:
            self.taxon_namespace = taxon_namespace
        self.next_event_index = 0
        self.events = [] # collection of HostEvent objects, sorted by time
        self.lineages = {} # keys: lineage_id (== int(Bipartition) == Bipartition.bitmask); values: HostRegimeLineageDefinition

    def compile(self):
        self.events.sort(key=lambda x: x.event_time)
        print(self.events)

class HostRegimeSamples(object):
    """
    A collection of host histories, one a single one of each a particular symbiont history will be conditioned.
    """

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
            lineage = HostRegime.HostRegimeLineageDefinition(
                    # tree_idx=edge_entry["tree_idx"],
                    lineage_id=lineage_id,
                    lineage_start_distribution=HostRegime.HostDistributionVector.from_string(edge_entry["edge_starting_state"]),
                    lineage_end_distribution=HostRegime.HostDistributionVector.from_string(edge_entry["edge_ending_state"]),
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

class HostSystem(object):
    """
    Models the host system, as defined by a HostRegime, for a particular simulation replicate.
    """

    class Host(object):
        """
        Manages the state of a host during a particular simulation replicate.
        """

        def __init__(self, host_regime_lineage_definition):
            self.host_regime_lineage_definition = host_regime_lineage_definition
            self.lineage_id = host_regime_lineage_definition.lineage_id
            self.start_distribution = host_regime_lineage_definition.lineage_start_distribution
            self.end_distribution = host_regime_lineage_definition.lineage_end_distribution

    def __init__(self, host_regime):
        self.compile(host_regime)

    def compile(self, host_regime):
        self.host_regime = host_regime
        self.host_lineages = {}
        num_areas = None
        for host_regime_lineage_id in self.host_regime.lineages.values():
            host = HostSystem.Host(host_regime_lineage_id_definition)
            self.host_lineages[host.lineage_id] = host
            if num_areas is None:
                num_areas = len(host.start_distribution)
            assert num_areas == len(host.start_distribution)
            assert num_areas == len(host.end_distribution)
        self.num_areas = num_areas

class Lineage(dendropy.Node):

    def __init__(self,
            index,
            distribution_vector=None,
            traits_vector=None,
            ):
        dendropy.Node.__init__(self)
        self.index = index
        self.distribution_vector = distribution_vector
        self.traits_vector = traits_vector
        self.is_extant = True
        self.edge.length = 0

class Phylogeny(dendropy.Tree):

    def node_factory(cls, **kwargs):
        return Lineage(**kwargs)
    node_factory = classmethod(node_factory)

    def __init__(self, *args, **kwargs):
        if kwargs:
            self.model = kwargs.pop("model")
            self.model_id = self.model.model_id
            self.annotations.add_bound_attribute("model_id")
            self.rng = kwargs.pop("rng")
            self.debug_mode = kwargs.pop("debug_mode")
            self.run_logger = kwargs.pop("run_logger")
            self.lineage_indexer = utility.IndexGenerator(0)
            if "seed_node" not in kwargs:
                seed_node = self.node_factory(
                        index=next(self.lineage_indexer),
                        distribution_vector=self.model.geography.new_distribution_vector(),
                        traits_vector=self.model.trait_types.new_traits_vector(),
                        )
                for trait_idx in range(len(self.model.trait_types)):
                    trait_states = [i for i in range(self.model.trait_types[trait_idx].nstates)]
                    seed_node.traits_vector[trait_idx] = self.rng.choice(trait_states)
                initial_area = self.rng.randint(0, len(seed_node.distribution_vector)-1)
                seed_node.distribution_vector[initial_area] = 1
                kwargs["seed_node"] = seed_node
            dendropy.Tree.__init__(self, *args, **kwargs)
            self.current_lineages = set([self.seed_node])
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

    def iterate_current_lineages(self):
        for lineage in self.current_lineages:
            yield lineage

    def split_lineage(self, lineage):
        dist1, dist2 = self._get_daughter_distributions(lineage)
        c1 = self.node_factory(
                index=next(self.lineage_indexer),
                distribution_vector=dist1,
                traits_vector=lineage.traits_vector.clone(),
                )
        c2 = self.node_factory(
                index=next(self.lineage_indexer),
                distribution_vector=dist2,
                traits_vector=lineage.traits_vector.clone(),
                )
        if self.debug_mode:
            self.run_logger.debug("Splitting {} with distribution {} under speciation mode {} to: {} (distribution: {}) and {} (distribution: {})".format(
                lineage,
                lineage.distribution_vector.presences(),
                speciation_mode,
                c1,
                dist1.presences(),
                c2,
                dist2.presences(),
                ))
            assert len(dist1.presences()) > 0
            assert len(dist2.presences()) > 0

        lineage.is_extant = False
        self.current_lineages.remove(lineage)
        lineage.add_child(c1)
        lineage.add_child(c2)
        self.current_lineages.add(c1)
        self.current_lineages.add(c2)

    def extinguish_lineage(self, lineage):
        self._make_lineage_extinct_on_phylogeny(lineage)

    def contract_lineage_range(self, lineage):
        presences = lineage.distribution_vector.presences()
        assert len(presences) > 0
        if len(presences) == 1:
            self._make_lineage_extinct_on_phylogeny(lineage)
        else:
            lineage.distribution_vector[ self.rng.choice(presences) ] = 0

    def _get_daughter_distributions(self, lineage):
        # speciation modes
        # 0:  single-area sympatric speciation
        #     -   ancestral range copied to both daughter species
        # 1:  sympatric subset: multi-area sympatric speciation
        #     -   d1: inherits complete range
        #     -   d2: inherits single area in ancestral range
        # 2:  (single-area) vicariance
        #     -   d1: single area
        #     -   d2: all other areas
        # 3:  (multi-area) vicariance
        #     -   ancestral range divided up unequally between two daughter
        #         species
        # 4:  founder-event jump dispersal
        #     -   single new area colonized
        presences = lineage.distribution_vector.presences()
        num_presences = len(presences)
        num_areas = len(self.model.geography.area_indexes)
        if num_presences <= 1:
            speciation_mode = 0
        else:
            if num_presences < num_areas:
                lineage_area_gain_rate = self.model.lineage_area_gain_rate_function(lineage)
                fes_weight = self.model.cladogenesis_founder_event_speciation_weight * lineage_area_gain_rate
            else:
                fes_weight = 0.0
            speciation_mode_weights = [
                self.model.cladogenesis_sympatric_subset_speciation_weight,
                self.model.cladogenesis_single_area_vicariance_speciation_weight,
                self.model.cladogenesis_widespread_vicariance_speciation_weight,
                fes_weight,
            ]
            sum_of_weights = sum(speciation_mode_weights)
            speciation_mode = 1 + weighted_index_choice(
                    weights=speciation_mode_weights,
                    sum_of_weights=sum_of_weights,
                    rng=self.rng)
        if speciation_mode == 0:
            # single-area sympatric speciation
            #     -   ancestral range copied to both daughter species
            dist1 = lineage.distribution_vector.clone()
            dist2 = lineage.distribution_vector.clone()
        elif speciation_mode == 1:
            # sympatric subset: multi-area sympatric speciation
            #     -   d1: inherits complete range
            #     -   d2: inherits single area in ancestral range
            dist1 = lineage.distribution_vector.clone()
            dist2 = self.model.geography.new_distribution_vector()
            # TODO: area diversity base speciation
            dist2[ self.rng.choice(presences) ] = 1
        elif speciation_mode == 2:
            # (single-area) allopatric vicariance
            #     -   d1: single area
            #     -   d2: all other areas
            dist1 = self.model.geography.new_distribution_vector()
            dist2 = self.model.geography.new_distribution_vector()
            self.rng.shuffle(presences)
            dist1[presences[0]] = 1
            for idx in presences[1:]:
                dist2[idx] = 1
        elif speciation_mode == 3:
            dist1 = self.model.geography.new_distribution_vector()
            dist2 = self.model.geography.new_distribution_vector()
            if num_presences == 2:
                dist1[presences[0]] = 1
                dist2[presences[1]] = 1
            else:
                n1 = self.rng.randint(1, num_presences-1)
                n2 = num_presences - n1
                if n2 == n1:
                    n1 += 1
                    n2 -= 1
                sample1 = set(self.rng.sample(presences, n1))
                for idx in self.model.geography.area_indexes:
                    if idx in sample1:
                        dist1[idx] = 1
                    else:
                        dist2[idx] = 1
        elif speciation_mode == 4:
            dist1 = lineage.distribution_vector.clone()
            dist2 = self.model.geography.new_distribution_vector()
            absences = [idx for idx in self.model.geography.area_indexes if idx not in presences]
            dist2[ self.rng.choice(absences) ] = 1
        else:
            raise ValueError(speciation_mode)
        return dist1, dist2

    def _make_lineage_extinct_on_phylogeny(self, lineage):
        if len(self.current_lineages) == 1:
            self.total_extinction_exception("no extant lineages remaining")
        lineage.is_extant = False
        self.current_lineages.remove(lineage)
        self.prune_subtree(lineage)

    def total_extinction_exception(self, msg):
        # self.run_logger.info("Total extinction: {}".format(msg))
        raise error.TotalExtinctionException(msg)

    def evolve_trait(self, lineage, trait_idx, state_idx):
        lineage.traits_vector[trait_idx] = state_idx

    def disperse_lineage(self, lineage, dest_area_idx):
        lineage.distribution_vector[dest_area_idx] = 1

    def focal_area_lineages(self):
        focal_area_lineages = set()
        for lineage in self.iterate_current_lineages():
            for area_idx in self.model.geography.focal_area_indexes:
                if lineage.distribution_vector[area_idx] == 1:
                    focal_area_lineages.add(lineage)
                    break
        return focal_area_lineages

    def num_focal_area_lineages(self):
        count = 0
        for lineage in self.iterate_current_lineages():
            for area_idx in self.model.geography.focal_area_indexes:
                if lineage.distribution_vector[area_idx] == 1:
                    count += 1
                    break
        return count

    def extract_focal_areas_tree(self):
        # tcopy = Phylogeny(self)
        tcopy = copy.deepcopy(self)
        focal_area_lineages = tcopy.focal_area_lineages()
        if len(focal_area_lineages) < 2:
            raise error.InsufficientFocalAreaLineagesSimulationException("insufficient lineages in focal area at termination".format(len(focal_area_lineages)))
        try:
            tcopy.filter_leaf_nodes(filter_fn=lambda x: x in focal_area_lineages)
        except dendropy.SeedNodeDeletionException:
            raise error.InsufficientFocalAreaLineagesSimulationException("no extant lineages in focal area at termination".format(len(focal_area_lineages)))
        return tcopy



if __name__ == "__main__":
    hrs = HostRegimeSamples()
    rb_data = os.path.join(utility.TEST_DATA_PATH, "revbayes", "bg_large.events.txt")
    rb_data_src = open(rb_data, "r")
    hrs.parse_host_biogeography(rb_data_src)
    hr = hrs.host_regimes[0]
    hs = HostSystem(hr)
