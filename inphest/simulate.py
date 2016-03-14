#! /usr/bin/env python

try:
    from StringIO import StringIO # Python 2 legacy support: StringIO in this module is the one needed (not io)
except ImportError:
    from io import StringIO # Python 3
import sys
import random
import collections
import argparse
import pprint
import copy
import json
import os
from distutils.util import strtobool

import dendropy
from dendropy.utility import textprocessing

import inphest
from inphest import model
from inphest import utility
from inphest import error

class InphestSimulator(object):

    @staticmethod
    def get_fixed_value_function(v, description):
        f = lambda x: v
        f.__doc__ = description
        return f

    @staticmethod
    def compose_trees_filepath(output_prefix):
        return output_prefix + ".trees"

    # @staticmethod
    # def compose_focal_areas_trees_filepath(output_prefix):
    #     return output_prefix + ".focal-areas.trees"

    # @staticmethod
    # def compose_all_areas_trees_filepath(output_prefix):
    #     return output_prefix + ".all-areas.trees"

    @staticmethod
    def simple_node_label_function(node):
        return "s{}".format(node.index)

    def __init__(self,
            inphest_model,
            config_d,
            is_verbose_setup):

        # configure
        self.elapsed_time = 0.0 # need to be here for logging
        config_d = dict(config_d) # make copy so we can pop items
        self.configure_simulator(config_d, verbose=is_verbose_setup)

        # set up model
        self.model = inphest_model

        # initialize host system
        self.host_system = model.HostSystem(host_history=self.model.host_history, run_logger=self.run_logger,)

        # track host events
        self.next_host_event = None
        self.processed_host_events = set()
        self.activated_host_lineages = set()
        self.deactivated_host_lineages = set()

        # initialize phylogeny
        self.phylogeny = model.SymbiontPhylogeny(
                model=self.model,
                host_system=self.host_system,
                rng=self.rng,
                debug_mode=self.debug_mode,
                run_logger=self.run_logger,
                )

        # begin logging generations
        self.run_logger.system = self

    def configure_simulator(self, config_d, verbose=True):

        self.name = config_d.pop("name", None)
        if self.name is None:
            self.name = str(id(self))
        self.output_prefix = config_d.pop("output_prefix", "inphest-{}".format(self.name))

        self.run_logger = config_d.pop("run_logger", None)
        if self.run_logger is None:
            self.run_logger = utility.RunLogger(
                    name="inphest",
                    stderr_logging_level=config_d.pop("standard_error_logging_level", "info"),
                    log_to_file=config_d.pop("log_to_file", True),
                    log_path=self.output_prefix + ".log",
                    file_logging_level=config_d.pop("file_logging_level", "info"),
                    )
        self.run_logger.system = self

        if verbose:
            self.run_logger.info("Configuring simulation '{}'".format(self.name))

        self.trees_file = config_d.pop("trees_file", None)
        if self.trees_file is None:
            self.trees_file = open(InphestSimulator.compose_trees_filepath(self.output_prefix), "w")
        if verbose:
            self.run_logger.info("Output trees filepath: {}".format(self.trees_file.name))

        if not self.trees_file:
            self.run_logger.warning("No trees will be stored!")

        self.is_suppress_internal_node_labels = config_d.pop("suppress_internal_node_labels", False)
        if verbose:
            self.run_logger.info("Internal node labels will{} be written on trees".format(" not" if self.is_suppress_internal_node_labels else ""))

        self.is_encode_nodes = config_d.pop("encode_nodes", True)
        if verbose:
            if self.is_encode_nodes:
                self.run_logger.info("Host associations and geographical ranges will be encoded on node labels")
            else:
                self.run_logger.info("Host associations and geographical ranges will NOT be encoded on node labels")

        self.is_annotate_nodes = config_d.pop("annotate_nodes", False)
        if verbose:
            if self.is_annotate_nodes:
                self.run_logger.info("Host associations and geographical ranges will be annotated on node labels")
            else:
                self.run_logger.info("Host associations and geographical ranges will NOT be annotated on node labels")

        self.rng = config_d.pop("rng", None)
        if self.rng is None:
            self.random_seed = config_d.pop("random_seed", None)
            if self.random_seed is None:
                self.random_seed = random.randint(0, sys.maxsize)
            if verbose:
                self.run_logger.info("Initializing with random seed {}".format(self.random_seed))
            self.rng = random.Random(self.random_seed)
        else:
            if "random_seed" in config_d:
                raise TypeError("Cannot specify both 'rng' and 'random_seed'")
            if verbose:
                self.run_logger.info("Using existing random number generator")

        self.log_frequency = config_d.pop("log_frequency", None)

        self.debug_mode = config_d.pop("debug_mode", False)
        if verbose and self.debug_mode:
            self.run_logger.info("Running in DEBUG mode")

        if config_d.pop("store_model_description", True):
            self.model_description_file = config_d.pop("model_description_file", None)
            if self.model_description_file is None:
                self.model_description_file = open(self.output_prefix + ".model.json", "w")
            if verbose:
                self.run_logger.info("Model description filepath: {}".format(self.model_description_file.name))
        else:
            self.model_description_file = None
            if verbose:
                self.run_logger.info("Model description will not be stored")

        if config_d:
            raise TypeError("Unsupported configuration keywords: {}".format(config_d))

    def activate_host_lineage(self, host_lineage):
        # if host_lineage.lineage_id == 8191:
        #     assert False
        #     print("\n\n\n\n\n\nactivating!\n\n\n\n\n")
        assert host_lineage not in self.activated_host_lineages, "Host lineage {} already activated".format(host_lineage.lineage_id)
        assert host_lineage not in self.deactivated_host_lineages
        host_lineage.activate(simulation_elapsed_time=self.elapsed_time, debug_mode=self.debug_mode)
        self.activated_host_lineages.add(host_lineage)

    def deactivate_host_lineage(self, host_lineage):
        assert host_lineage in self.activated_host_lineages
        assert host_lineage not in self.deactivated_host_lineages
        host_lineage.deactivate()
        self.activated_host_lineages.remove(host_lineage)
        self.deactivated_host_lineages.add(host_lineage)

    def run(self):

        ### Save model
        if self.model_description_file is not None:
            self.run_logger.warning("[WARNING] Model description output not implemented yet")
            # self.model.write_model(self.model_description_file)

        ### Initialize time
        self.elapsed_time = 0.0

        ### Initialize logging
        ### None: default logging, 0: no logging
        # if self.log_frequency is None:
        #     if self.model.target_focal_area_lineages:
        #         default_log_frequency = 1
        #     else:
        #         default_log_frequency = self.model.max_time/100
        # if self.log_frequency:
        #     if self.model.target_focal_area_lineages:
        #         last_logged_num_tips = 0
        #     else:
        #         last_logged_time = 0.0

        if self.log_frequency is None:
            default_log_frequency = self.model.max_time/100
        if self.log_frequency:
            last_logged_time = 0.0

        ### check system setup
        if self.debug_mode:
            self.host_system.debug_check(simulation_elapsed_time=None)

        ### Initialize run debugging
        if self.debug_mode:
            num_events = 0

        ### Initialize seed node distribution
        extant_host_lineages = self.host_system.extant_host_lineages_at_current_time(0)
        self.activate_host_lineage(self.host_system.seed_host_lineage)
        for lineage in self.phylogeny.current_lineages:
            lineage.add_host_in_area(host_lineage=self.host_system.seed_host_lineage)

        ### Initialize termination conditiong checking
        # ntips_in_focal_areas = self.phylogeny.num_focal_area_lineages()
        # ntips = len(self.phylogeny.current_lineages)

        while True:

            ### DEBUG
            if self.debug_mode:
                num_events += 1
                # self.run_logger.debug("Pre-event {}: debug check: {}".format(num_events, self.debug_compose_tree(self.phylogeny)))
                self.run_logger.debug("Post-event {}: debug check: current lineages = {}".format(num_events, len(self.phylogeny.current_lineages)))
                self.phylogeny._debug_check_tree()
                for lineage in self.phylogeny.current_lineages:
                    assert lineage.is_extant
                    lineage.debug_check(simulation_elapsed_time=self.elapsed_time)

            ### LOGGING
            if self.log_frequency:
                if self.model.target_focal_area_lineages:
                    if ntips_in_focal_areas - last_logged_num_tips >= self.log_frequency:
                        last_logged_num_tips = ntips_in_focal_areas
                        self.run_logger.info("{} lineages occurring in focal areas, {} lineages across all areas".format(ntips_in_focal_areas, ntips))
                else:
                    if self.elapsed_time - last_logged_time >= self.log_frequency:
                        last_logged_time = self.elapsed_time
                        self.run_logger.info("{} lineages occurring in focal areas, {} lineages across all areas".format(ntips_in_focal_areas, ntips))

            ### EVENT SCHEDULING
            event_calls, event_rates = self.schedule_events()
            sum_of_event_rates = sum(event_rates)
            if self.debug_mode:
                if sum_of_event_rates == 0:
                    self.run_logger.debug("Sum of event rates is 0: {}".format(event_rates))

            time_till_event = self.rng.expovariate(sum_of_event_rates)

            if self.next_host_event is None and self.host_system.host_events:
                try:
                    self.next_host_event = self.host_system.host_events.pop(0)
                except IndexError: # pop from empty list
                    self.next_host_event = None
            if self.next_host_event and self.next_host_event.event_time < (self.elapsed_time + time_till_event):
                time_till_event = self.next_host_event.event_time - self.elapsed_time
                if self.debug_mode:
                    self.run_logger.debug("Host Event {} of {}: {}".format(
                        len(self.host_system.host_history.events)-len(self.host_system.host_events),
                        len(self.host_system.host_history.events),
                        self.next_host_event))
                event_f = self.process_host_event
                event_kwargs = {"host_event": self.next_host_event}
                self.next_host_event = None
            else:
                event_idx = model.weighted_index_choice(
                        weights=event_rates,
                        sum_of_weights=sum_of_event_rates,
                        rng=self.rng)
                if self.debug_mode:
                    self.run_logger.debug("Symbiont Event {}: {}".format(num_events, event_calls[event_idx]))
                event_f = event_calls[event_idx][0]
                event_kwargs = event_calls[event_idx][1]

            self.elapsed_time += time_till_event
            if self.model.max_time and self.elapsed_time > self.model.max_time:
                self.elapsed_time = self.model.max_time
                assert len(self.processed_host_events) == len(self.host_system.host_history.events)
                assert len(self.host_system.host_events) == 0
                self.run_logger.info("Termination condition of t = {} reached: storing results and terminating".format(self.elapsed_time))
                self.store_sample(trees_file=self.trees_file)
                break
            for lineage in self.phylogeny.current_lineage_iter():
                lineage.edge.length += time_till_event

            ### EVENT EXECUTION
            try:
                event_f(**event_kwargs)
            except model.SymbiontLineage.NullDistributionException as lineage_null_distribution_exception:
                self.phylogeny.extinction_rate(symbiont_lineage=lineage_null_distribution_exception.lineage)

            ### DEBUG
            if self.debug_mode:
                # self.run_logger.debug("Post-event {}: debug check: {}".format(num_events, self.debug_compose_tree(self.phylogeny)))
                self.run_logger.debug("Post-event {}: debug check: current lineages = {}".format(num_events, len(self.phylogeny.current_lineages)))
                self.phylogeny._debug_check_tree()
                for lineage in self.phylogeny.current_lineages:
                    assert lineage.is_extant
                    lineage.debug_check(simulation_elapsed_time=self.elapsed_time)

            # if self.model.gsa_termination_focal_area_lineages and ntips_in_focal_areas >= self.model.gsa_termination_focal_area_lineages:
            #     # select/process one of the previously stored snapshots, write to final results file,
            #     # and then break
            #     raise NotImplementedError
            # elif self.model.gsa_termination_focal_area_lineages and ntips_in_focal_areas == self.target_focal_area_lineages:
            #     # store snapshot in log, but do not break
            #     raise NotImplementedError
            # elif self.model.target_focal_area_lineages and ntips_in_focal_areas >= self.model.target_focal_area_lineages:
            #     self.run_logger.info("Termination condition of {} lineages in focal areas reached at t = {}: storing results and terminating".format(self.model.target_focal_area_lineages, self.elapsed_time))
            #     self.store_sample(trees_file=self.trees_file)
            #     break

    def schedule_events(self):
        event_calls = []
        event_rates = []

        # if self.debug_mode:
        #     num_current_lineages = len(self.phylogeny.current_lineages)
        #     self.run_logger.debug("Scheduling events for {} current lineages".format(
        #         num_current_lineages))

        for lineage in self.phylogeny.current_lineage_iter():

            #---
            # Diversification Process: Birth (Speciation)
            speciation_rate = self.model.symbiont_lineage_birth_rate_function(symbiont_lineage=lineage, simulation_elapsed_time=self.elapsed_time)
            if speciation_rate:
                event_calls.append( (self.phylogeny.split_lineage, {"symbiont_lineage": lineage}) )
                event_rates.append(speciation_rate)

            #---
            # Diversification Process: Death (Extinction)
            extinction_rate = self.model.symbiont_lineage_death_rate_function(symbiont_lineage=lineage, simulation_elapsed_time=self.elapsed_time)
            if extinction_rate:
                event_calls.append( (self.phylogeny.extinguish_lineage, {"symbiont_lineage": lineage}) )
                event_rates.append(extinction_rate)

            #---
            # Anagenetic Host Set Evolution: Host Gain
            # ----------------------------------------
            # Submodel Design Objectives:
            # 1.    We want to be able to specify that the rate of gaining a
            #       particular host in a particular area as functions of:
            #       -   The area (e.g., the number of host lineages in the area; the
            #           number of symbiont lineages in the area; or the particular
            #           host or symbiont lineages present/absent from an area).
            #       -   The source and destination hosts (e.g., the phylogenetic
            #           distances between the source and destination hosts; the
            #           number of symbiont lineages in the destination host; etc.)
            # 2.    We want to be able to specify a mean per-lineage rate of
            #       transmission. Allows for estimating this from empirical data
            #       using some reasonable if simplified and low-fidelity model:
            #       e.g., as a per-lineage trait evolution rate, where the host set
            #       is a multistate character.
            # Objective (1) means that the rates must be calculated on a
            # per-source host per-destination host per area basis.
            # Objective (2) means that the transmission (host gain) rate
            # weights across all host gain events needs to sum to 1 (with the
            # actual rate obtained by multiplying with the system-wide mean
            # (per-lineage) host gain rate.)
            # NOTE: this can be optimized! Right now, more a naive
            # reference implementation of the logic.
            # NOTE: this needs to be tested.
            infected_hosts = {}
            uninfected_hosts = {}
            num_potential_new_host_infection_events = 0
            for area in lineage.area_iter():
                infected_hosts[area] = []
                uninfected_hosts[area] = []
                for host_lineage in area.host_lineages:
                    if self.debug_mode:
                        host_lineage.assert_correctly_extant(simulation_elapsed_time=self.elapsed_time)
                    if lineage.has_host_in_area(host_lineage, area):
                        infected_hosts[area].append( host_lineage )
                    else:
                        uninfected_hosts[area].append( host_lineage )
                    num_potential_new_host_infection_events += ( len(uninfected_hosts[area]) * len(infected_hosts[area]) )
            if num_potential_new_host_infection_events > 0:
                transmission_event_calls = []
                transmission_event_rates = []
                for area in lineage.area_iter():
                    for src_host in infected_hosts[area]:
                        for dest_host in uninfected_hosts[area]:
                            rate = self.model.symbiont_lineage_host_gain_rate_function(
                                    symbiont_lineage=lineage,
                                    from_host_lineage=src_host,
                                    to_host_lineage=dest_host,
                                    area=area,
                                    num_potential_new_host_infection_events=num_potential_new_host_infection_events,
                                    simulation_elapsed_time=self.elapsed_time)
                            transmission_event_calls.append( (lineage.add_host_in_area,  {"host_lineage": dest_host, "area": area,}) )
                            transmission_event_rates.append(rate)
                nf = float(sum(transmission_event_rates))
                transmission_event_rates = [tvr / nf for tvr in transmission_event_rates]
                event_calls.extend( transmission_event_calls )
                event_rates.extend( transmission_event_rates )

            #---
            # Anagenetic Host Set Evolution: Host Loss
            host_loss_rate = self.model.symbiont_lineage_host_loss_rate_function(symbiont_lineage=lineage, simulation_elapsed_time=self.elapsed_time)
            if host_loss_rate:
                for host_lineage in lineage.host_iter():
                    event_calls.append( (lineage.remove_host, {"host_lineage": host_lineage}) )
                    event_rates.append(host_loss_rate)

            #---
            # Anagenetic Area Set Evolution: Area Gain
            # (same notes apply as for "Anagenetic Host Set Evolution, Host Gain")
            occupied_areas = {}
            unoccupied_areas = {}
            num_potential_new_area_infection_events = 0
            for host_lineage in lineage.host_iter():
                if self.debug_mode:
                    host_lineage.assert_correctly_extant(simulation_elapsed_time=self.elapsed_time)
                occupied_areas[host_lineage] = []
                unoccupied_areas[host_lineage] = []
                for area in host_lineage.current_area_iter():
                    if lineage.has_host_in_area(host_lineage, area):
                        occupied_areas[host_lineage].append(area)
                    else:
                        unoccupied_areas[host_lineage].append(area)
                num_potential_new_area_infection_events += ( len(occupied_areas[host_lineage]) * len(unoccupied_areas[host_lineage]) )
            # if num_potential_new_area_infection_events > 0:
            #     per_host_area_gain_rate = self.model.symbiont_lineage_area_gain_rate_function(symbiont_lineage=lineage, simulation_elapsed_time=self.elapsed_time) / num_potential_new_area_infection_events
            #     for host_lineage in lineage.host_iter():
            #         for src_area in occupied_areas[host_lineage]:
            #             for dest_area in unoccupied_areas[host_lineage]:
            #                 event_calls.append( (lineage.add_host_in_area, {"host_lineage":host_lineage, "area": dest_area,}) )
            #                 event_rates.append(per_host_area_gain_rate)
            if num_potential_new_area_infection_events > 0:
                dispersal_event_calls = []
                dispersal_event_rates = []
                for host_lineage in lineage.host_iter():
                    if self.debug_mode:
                        host_lineage.assert_correctly_extant(simulation_elapsed_time=self.elapsed_time)
                    for src_area in occupied_areas[host_lineage]:
                        for dest_area in unoccupied_areas[host_lineage]:
                            rate = self.model.symbiont_lineage_area_gain_rate_function(
                                    symbiont_lineage=lineage,
                                    from_area_lineage=src_area,
                                    to_area_lineage=dest_area,
                                    host=host_lineage,
                                    num_potential_new_area_infection_events=num_potential_new_area_infection_events,
                                    simulation_elapsed_time=self.elapsed_time)
                            dispersal_event_calls.append( (lineage.add_host_in_area, {"host_lineage":host_lineage, "area": dest_area,}) )
                            dispersal_event_rates.append(rate)
                nf = float(sum(dispersal_event_rates))
                dispersal_event_rates = [tvr / nf for tvr in dispersal_event_rates]
                event_calls.extend( dispersal_event_calls )
                event_rates.extend( dispersal_event_rates )

            #---
            # Anagenetic Area Set Evolution: Area Loss
            area_loss_rate = self.model.symbiont_lineage_area_loss_rate_function(symbiont_lineage=lineage, simulation_elapsed_time=self.elapsed_time)
            if area_loss_rate:
                for host_lineage in lineage.host_iter():
                    for area in lineage.areas_in_host_iter():
                        event_calls.append( (lineage.remove_host_in_area, {"host_lineage": lineage, "area": area} ))
                        event_rates.append(area_loss_rate)

        # sum_of_event_rates = sum(event_rates)
        return event_calls, event_rates

    def process_host_event(self, host_event):
        # print(host_event)
        assert host_event not in self.processed_host_events
        host_lineage = self.host_system.host_lineages_by_id[host_event.lineage_id]
        assert host_lineage.start_time <= self.elapsed_time
        assert host_lineage.end_time >= self.elapsed_time
        if self.debug_mode:
            host_lineage.debug_check(simulation_elapsed_time=self.elapsed_time)
        if host_event.event_type == "anagenesis" and host_event.event_subtype == "area_gain":
            if self.debug_mode:
                assert host_lineage._current_distribution_check_bitlist[host_event.area_idx] == "0"
                host_lineage._current_distribution_check_bitlist[host_event.area_idx] = "1"
                self.run_logger.debug("Host lineage {}: anagenetic gain of area with index {}: {}".format(host_lineage.lineage_id, host_event.area_idx, host_lineage._current_distribution_check_bitlist))
            area = self.host_system.areas[host_event.area_idx]
            host_lineage.add_area(area)
        elif host_event.event_type == "anagenesis" and host_event.event_subtype == "area_loss":
            if self.debug_mode:
                assert host_lineage._current_distribution_check_bitlist[host_event.area_idx] == "1"
                host_lineage._current_distribution_check_bitlist[host_event.area_idx] = "0"
                self.run_logger.debug("Host lineage {}: anagenetic loss of area with index {}: {}".format(host_lineage.lineage_id, host_event.area_idx, host_lineage._current_distribution_check_bitlist))
            area = self.host_system.areas[host_event.area_idx]
            for symbiont_lineage in self.phylogeny.current_lineage_iter():
                if symbiont_lineage.has_host_in_area(host_lineage, area):
                    try:
                        symbiont_lineage.remove_host_in_area(host_lineage, area)
                    except model.SymbiontLineage.NullDistributionException:
                        self.phylogeny.extinguish_lineage(symbiont_lineage)
            host_lineage.remove_area(area)
        elif host_event.event_type == "cladogenesis":
            host_child0_lineage = self.host_system.host_lineages_by_id[host_event.child0_lineage_id]
            self.activate_host_lineage(host_child0_lineage)
            host_child1_lineage = self.host_system.host_lineages_by_id[host_event.child1_lineage_id]
            self.activate_host_lineage(host_child1_lineage)
            if self.debug_mode:
                assert host_child0_lineage.lineage_id == host_event.child0_lineage_id
                assert host_child1_lineage.lineage_id == host_event.child1_lineage_id
                self.run_logger.debug("Host lineage {} ({}): splitting into lineages {} ({}) and {} ({})".format(
                    host_lineage.lineage_id,
                    host_lineage._current_distribution_check_bitlist,
                    host_event.child0_lineage_id,
                    host_child0_lineage._current_distribution_check_bitlist,
                    host_event.child1_lineage_id,
                    host_child1_lineage._current_distribution_check_bitlist,
                    ))
                for ch_lineage in (host_child0_lineage, host_child1_lineage):
                    # assert ch_lineage.start_time <= self.elapsed_time
                    # assert ch_lineage.end_time >= self.elapsed_time
                    ch_lineage.debug_check(simulation_elapsed_time=self.elapsed_time)
            for symbiont_lineage in self.phylogeny.current_lineage_iter():
                if not symbiont_lineage.has_host(host_lineage):
                    continue
                lineage_areas_with_host = set(symbiont_lineage.areas_in_host_iter(host_lineage))
                hosts_in_areas_added = 0
                ## TODO: need to special case jump dispersal event subtype
                # The following scheme assigns a symbiont to a new daughter
                # host lineage in an area if the symbiont occupied the parent
                # in that area. This will mean in the case of a jump dispersal,
                # where the daughter lineage goes to a new area, no new
                # symbiont lineages will be assigned to the host at all.

                # if self.debug_mode:
                #     for area in lineage_areas_with_host:
                #         if not ( host_child0_lineage.has_area(area) or host_child1_lineage.has_area(area) ):
                #             print("Host {:5}, {:10}: {}".format(host_lineage.rb_index, host_lineage.leafset_bitstring, sorted( "{:03}".format(a.area_idx) for a in host_lineage.current_area_iter() )))
                #             print(" Ch0 {:5}, {:10}: {}".format(host_child0_lineage.rb_index, host_child0_lineage.leafset_bitstring, sorted( "{:03}".format(a.area_idx) for a in host_child0_lineage.current_area_iter() )))
                #             print(" Ch1 {:5}, {:10}: {}".format(host_child1_lineage.rb_index, host_child1_lineage.leafset_bitstring, sorted( "{:03}".format(a.area_idx) for a in host_child1_lineage.current_area_iter() )))
                #             print("Symb {}".format(sorted( "{:03}".format(a.area_idx) for a in symbiont_lineage.areas_in_host_iter(host_lineage) )))

                for ch_lineage in (host_child0_lineage, host_child1_lineage):
                    for area in ch_lineage.current_area_iter():
                        if area in lineage_areas_with_host:
                            symbiont_lineage.add_host_in_area(host_lineage=ch_lineage, area=area)
                            hosts_in_areas_added += 1
                assert hosts_in_areas_added > 0
                symbiont_lineage.remove_host(host_lineage)
            self.deactivate_host_lineage(host_lineage)

        if self.debug_mode:
            host_lineage.debug_check(simulation_elapsed_time=self.elapsed_time)

        self.processed_host_events.add(host_event)

    def store_sample(self, trees_file):
        self.write_tree(
                out=trees_file,
                tree=self.phylogeny,
                )

    def write_tree(self, out, tree):
        if self.is_encode_nodes:
            # host_lineages = self.host_system.extant_host_lineages_at_current_time(self.elapsed_time)
            # host_lineages = self.host_system.leaf_host_lineages
            host_lineages = self.host_system.host_lineages
            lineage_labels = {}
            for symbiont_lineage in tree:
                current_hosts = []
                for idx, host_lineage in enumerate(host_lineages):
                    if symbiont_lineage.has_host(host_lineage):
                        current_hosts.append("1")
                    else:
                        current_hosts.append("0")
                current_hosts_bitstring = "".join(current_hosts)
                label = "s{lineage_index}{sep}{host_occurrences}".format(
                        lineage_index=symbiont_lineage.index,
                        sep=model.InphestModel._LABEL_COMPONENTS_SEPARATOR,
                        host_occurrences=current_hosts_bitstring,
                        )
                lineage_labels[symbiont_lineage] = label
            labelf = lambda x: lineage_labels[x]
        else:
            labelf = InphestSimulator.simple_node_label_function
        tree.write_to_stream(
                out,
                schema="newick",
                suppress_annotations=False,
                node_label_compose_fn=labelf,
                suppress_internal_node_labels=self.is_suppress_internal_node_labels,
                )

def repeat_run(
        output_prefix,
        nreps,
        host_history_samples_path,
        model_definition_source,
        model_definition_type="python-dict",
        config_d=None,
        interpolate_missing_model_values=False,
        random_seed=None,
        stderr_logging_level="info",
        file_logging_level="debug",
        maximum_num_restarts_per_replicates=100,
        debug_mode=False):
    """
    Executes multiple runs of the Inphest simulator under identical
    parameters to produce the specified number of replicates, discarding failed
    runs.

    Parameters
    ----------
    nreps : integer
        Number of replicates to produce *per host biogeographical history sample*.
    config_d : dict
        Simulator configuration parameters as keyword-value pairs. To be
        re-used for each replicate.
    host_history_samples_path : object
        Path to host biogeographical history samples.
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

    interpolate_missing_model_values : bool
        Allow missing values in model to be populated by default values (inadvisable).
    random_seed : integer
        Random seed to be used (for single random number generator across all
        replicates).
    stderr_logging_level : string or None
        Message level threshold for screen logs; if 'none' or `None`, screen
        logs will be supprsed.
    file_logging_level : string or None
        Message level threshold for file logs; if 'none' or `None`, file
        logs will be supprsed.
    maximum_num_restarts_per_replicates : int
        A failed replicate (due to e.g., total extinction of all taxa) will be
        re-run. This limits the number of re-runs.
    """
    if output_prefix is None:
        output_prefix = config_d.pop("output_prefix", "inphest")
    if config_d is None:
        config_d = {}
    config_d["output_prefix"] = output_prefix
    if stderr_logging_level is None or stderr_logging_level.lower() == "none":
        log_to_stderr = False
    else:
        log_to_stderr = True
    if file_logging_level is None or file_logging_level.lower() == "none":
        log_to_file = False
    else:
        log_to_file = True
    if "run_logger" not in config_d:
        config_d["run_logger"] = utility.RunLogger(
                name="inphest",
                log_to_stderr=log_to_stderr,
                stderr_logging_level=stderr_logging_level,
                log_to_file=log_to_file,
                log_path=output_prefix + ".log",
                file_logging_level=file_logging_level,
                )
    config_d["debug_mode"] = debug_mode
    run_logger = config_d["run_logger"]
    run_logger.info("-inphest- Starting: {}".format(inphest.description()))
    if "rng" not in config_d:
        if random_seed is None:
            random_seed = config_d.pop("random_seed", None)
            if random_seed is None:
                random_seed = random.randint(0, sys.maxsize)
        else:
            random_seed = int(random_seed)
        run_logger.info("-inphest- Initializing with random seed: {}".format(random_seed))
        config_d["rng"] = random.Random(random_seed)
    else:
        run_logger.info("-inphest- Using existing RNG: {}".format(config_d["rng"]))
    if config_d.get("store_trees", True) and "trees_file" not in config_d:
        config_d["trees_file"] = open(InphestSimulator.compose_trees_filepath(output_prefix), "w")

    host_history_samples_path = os.path.normpath(host_history_samples_path)
    run_logger.info("-inphest- Using host biogeographical regime samples from: {}".format(host_history_samples_path))
    hrs = model.HostHistorySamples()
    rb_data_src = open(rb_data, "r")
    hrs.parse_host_biogeography(rb_data_src)
    # run_logger.info("-inphest- {} host biogeographical regime samples found in source '{}'".format(len(hrs.host_histories), host_history_samples_path))
    run_logger.info("-inphest- {} host biogeographical regime samples found in source".format(len(hrs.host_histories), host_history_samples_path))

    current_rep = 0
    while current_rep < nreps:
        for host_history_idx, host_history in enumerate(hrs.host_histories):
            simulation_name="Run_{}_{}".format(current_rep+1, host_history_idx+1)
            run_output_prefix = "{}.R{:04d}.H{:04d}".format(output_prefix, current_rep+1, host_history_idx+1)
            run_logger.info("-inphest- Replicate {} of {}, host regime {} of {}: Starting".format(current_rep+1, nreps, host_history_idx+1, len(hrs.host_histories)))
            num_restarts = 0
            while True:
                if num_restarts == 0 and current_rep == 0:
                    is_verbose_setup = True
                    model_setup_logger = run_logger
                else:
                    is_verbose_setup = False
                    model_setup_logger = None
                inphest_model = model.InphestModel.create(
                        host_history=host_history,
                        model_definition_source=model_definition_source,
                        model_definition_type=model_definition_type,
                        interpolate_missing_model_values=interpolate_missing_model_values,
                        run_logger=model_setup_logger,
                        )
                # inphest_model = model.InphestModel.create(
                #         model_definition_source=model_definition_source,
                #         model_definition_type=model_definition_type,
                #         interpolate_missing_model_values=interpolate_missing_model_values,
                #         run_logger=model_setup_logger,
                #         )
                inphest_simulator = InphestSimulator(
                    inphest_model=inphest_model,
                    config_d=config_d,
                    is_verbose_setup=is_verbose_setup)
                try:
                    inphest_simulator.run()
                    run_logger.system = None
                except error.InphestException as e:
                    run_logger.system = None
                    run_logger.info("-inphest- Replicate {} of {}, host regime {} of {}: Simulation failure before termination condition at t = {}: {}".format(current_rep+1, nreps, host_history_idx+1, len(hrs.host_histories), inphest_simulator.elapsed_time, e))
                    num_restarts += 1
                    if num_restarts > maximum_num_restarts_per_replicates:
                        run_logger.info("-inphest- Replicate {} of {}, host regime {} of {}: Maximum number of restarts exceeded: aborting".format(current_rep+1, nreps, host_history_idx+1, len(hrs.host_histories)))
                        break
                    else:
                        run_logger.info("-inphest- Replicate {} of {}, host regime {} of {}: Restarting replicate (number of restarts: {})".format(current_rep+1, nreps, host_history_idx+1, len(hrs.host_histories), num_restarts))
                else:
                    run_logger.system = None
                    run_logger.info("-inphest- Replicate {} of {}, host regime {} of {}: Completed to termination condition at t = {}".format(current_rep+1, nreps, host_history_idx+1, len(hrs.host_histories), inphest_simulator.elapsed_time))
                    num_restarts = 0
                    break
        current_rep += 1

if __name__ == "__main__":
    rb_data = os.path.join(utility.TEST_DATA_PATH, "revbayes", "bg_large.events.txt")
    if len(sys.argv) > 1:
        random_seed = int(sys.argv[1])
    else:
        random_seed = None
    repeat_run(
        output_prefix="test",
        nreps=1,
        host_history_samples_path=rb_data,
        model_definition_source={},
        model_definition_type="python-dict",
        config_d=None,
        interpolate_missing_model_values=False,
        random_seed=random_seed,
        stderr_logging_level="info",
        file_logging_level="debug",
        maximum_num_restarts_per_replicates=100,
        debug_mode=True)

    # hrs = model.HostHistorySamples()
    # rb_data = os.path.join(utility.TEST_DATA_PATH, "revbayes", "bg_large.events.txt")
    # rb_data_src = open(rb_data, "r")
    # hrs.parse_host_biogeography(rb_data_src)
    # for host_history in hrs.host_histories:
    #     im = model.InphestModel.create(
    #             host_history=host_history,
    #             model_definition_source={},
    #             model_definition_type="python-dict",
    #             )

