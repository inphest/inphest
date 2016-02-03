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

        # initialize phylogeny
        self.phylogeny = model.Phylogeny(
                model=self.model,
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
                self.run_logger.info("Trait states and ranges will be encoded on node labels")
            else:
                self.run_logger.info("Trait states and ranges will NOT be encoded on node labels")

        self.is_annotate_nodes = config_d.pop("annotate_nodes", False)
        if verbose:
            if self.is_annotate_nodes:
                self.run_logger.info("Trait states and ranges will be annotated on node labels")
            else:
                self.run_logger.info("Trait states and ranges will NOT be annotated on node labels")

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

        ### Initialize debugging
        if self.debug_mode:
            num_events = 0

        ### Initialize termination conditiong checking
        ntips_in_focal_areas = self.phylogeny.num_focal_area_lineages()
        ntips = len(self.phylogeny.current_lineages)

        while True:

            ### DEBUG
            if self.debug_mode:
                num_events += 1
                # self.run_logger.debug("Pre-event {}: debug check: {}".format(num_events, self.debug_compose_tree(self.phylogeny)))
                self.run_logger.debug("Post-event {}: debug check: current lineages = {}".format(num_events, len(self.phylogeny.current_lineages)))
                self.phylogeny._debug_check_tree()
                for lineage in self.phylogeny.current_lineages:
                    assert lineage.is_extant
                    assert len(lineage.distribution_vector.presences()) > 0

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
            self.elapsed_time += time_till_event
            if self.model.max_time and self.elapsed_time > self.model.max_time:
                self.elapsed_time = self.model.max_time
                self.run_logger.info("Termination condition of t = {} reached: storing results and terminating".format(self.elapsed_time))
                self.store_sample(trees_file=self.trees_file)
                break
            for lineage in self.phylogeny.iterate_current_lineages():
                lineage.edge.length += time_till_event

            ### EVENT SELECTION AND EXECUTION
            event_idx = model.weighted_index_choice(
                    weights=event_rates,
                    sum_of_weights=sum_of_event_rates,
                    rng=self.rng)
            if self.debug_mode:
                self.run_logger.debug("Event {}: {}".format(num_events, event_calls[event_idx]))
            event_calls[event_idx][0](*event_calls[event_idx][1:])

            ### DEBUG
            if self.debug_mode:
                # self.run_logger.debug("Post-event {}: debug check: {}".format(num_events, self.debug_compose_tree(self.phylogeny)))
                self.run_logger.debug("Post-event {}: debug check: current lineages = {}".format(num_events, len(self.phylogeny.current_lineages)))
                self.phylogeny._debug_check_tree()
                for lineage in self.phylogeny.current_lineages:
                    assert lineage.is_extant
                    assert len(lineage.distribution_vector.presences()) > 0


            ntips_in_focal_areas = self.phylogeny.num_focal_area_lineages()
            ntips = len(self.phylogeny.current_lineages)

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

        for lineage in self.phylogeny.iterate_current_lineages():
            # speciation
            speciation_rate = self.model.lineage_birth_rate_function(lineage)
            if speciation_rate:
                event_calls.append( (self.phylogeny.split_lineage, lineage) )
                event_rates.append(speciation_rate)
            # extinction
            extinction_rate = self.model.lineage_death_rate_function(lineage)
            if extinction_rate:
                event_calls.append( (self.phylogeny.extinguish_lineage, lineage) )
                event_rates.append(extinction_rate)
            # extinction
            area_loss_rate = self.model.lineage_area_loss_rate_function(lineage)
            if area_loss_rate:
                event_calls.append( (self.phylogeny.contract_lineage_range, lineage) )
                event_rates.append(area_loss_rate)
            # trait evolution
            for trait_idx, current_state_idx in enumerate(lineage.traits_vector):
                ## normalized
                for proposed_state_idx in range(self.model.trait_types[trait_idx].nstates):
                    if proposed_state_idx == current_state_idx:
                        continue
                    trait_transition_rate = self.model.trait_types[trait_idx].transition_rate_matrix[current_state_idx][proposed_state_idx]
                    if trait_transition_rate:
                        event_calls.append( (self.phylogeny.evolve_trait, lineage, trait_idx, proposed_state_idx) )
                        event_rates.append(trait_transition_rate)
            # dispersal
            lineage_area_gain_rate = self.model.lineage_area_gain_rate_function(lineage)
            if not lineage_area_gain_rate:
                continue
            for dest_area_idx in self.model.geography.area_indexes:
                if lineage.distribution_vector[dest_area_idx]:
                    # already occurs here: do we model it or not?
                    continue
                sum_of_area_connection_weights_to_dest = 0.0
                for src_area_idx, occurs in enumerate(lineage.distribution_vector):
                    if not occurs:
                        continue
                    if dest_area_idx == src_area_idx:
                        continue
                    sum_of_area_connection_weights_to_dest += lineage_area_gain_rate * self.model.geography.area_connection_weights[src_area_idx][dest_area_idx]
                if sum_of_area_connection_weights_to_dest:
                    event_calls.append( (self.phylogeny.disperse_lineage, lineage, dest_area_idx) )
                    event_rates.append(sum_of_area_connection_weights_to_dest)
        # sum_of_event_rates = sum(event_rates)
        return event_calls, event_rates

    def store_sample(self, focal_areas_tree_out, all_areas_tree_out):
        if focal_areas_tree_out is not None:
            focal_areas_tree = self.phylogeny.extract_focal_areas_tree()
            n = len(focal_areas_tree.seed_node._child_nodes)
            if n < 2:
                raise error.InsufficientFocalAreaLineagesSimulationException("Insufficient lineages in focal area: {}".format(n))
            self.write_focal_areas_tree(
                    out=focal_areas_tree_out,
                    tree=focal_areas_tree,
                    )
        if all_areas_tree_out is not None:
            self.write_all_areas_tree(
                    out=all_areas_tree_out,
                    tree=self.phylogeny,
                    )

    def write_focal_areas_tree(self, out, tree):
        if self.is_encode_nodes:
            labelf = lambda x: self.model.encode_lineage(x,
                    set_label=False,
                    add_annotation=self.is_annotate_nodes,
                    exclude_supplemental_areas=True)
        else:
            labelf = InphestSimulator.simple_node_label_function
        tree.write_to_stream(
                out,
                schema="newick",
                suppress_annotations=False,
                node_label_compose_fn=labelf,
                suppress_internal_node_labels=self.is_suppress_internal_node_labels,
                )

    def write_all_areas_tree(self, out, tree):
        if self.is_encode_nodes:
            labelf = lambda x: self.model.encode_lineage(x,
                    set_label=False,
                    add_annotation=self.is_annotate_nodes,
                    exclude_supplemental_areas=False)
        else:
            labelf = InphestSimulator.simple_node_label_function
        tree.write_to_stream(
                out,
                schema="newick",
                suppress_annotations=False,
                node_label_compose_fn=labelf,
                suppress_internal_node_labels=self.is_suppress_internal_node_labels,
                )

    def debug_compose_tree(self, tree):
        labelf = lambda x: self.model.encode_lineage(x,
                set_label=False,
                add_annotation=False,
                exclude_supplemental_areas=False)
        s = tree.as_string(
                "newick",
                node_label_compose_fn=self.model.encode_all_areas_lineage,
                suppress_edge_lengths=True)
        return s.replace("\n", "")


def repeat_run(
        output_prefix,
        nreps,
        host_regime_samples_path,
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
    host_regime_samples_path : object
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

    host_regime_samples_path = os.path.normpath(host_regime_samples_path)
    run_logger.info("-inphest- Using host biogeographical regime samples from: {}".format(host_regime_samples_path))
    hrs = model.HostRegimeSamples()
    rb_data_src = open(rb_data, "r")
    hrs.parse_host_biogeography(rb_data_src)
    # run_logger.info("-inphest- {} host biogeographical regime samples found in source '{}'".format(len(hrs.host_regimes), host_regime_samples_path))
    run_logger.info("-inphest- {} host biogeographical regime samples found in source".format(len(hrs.host_regimes), host_regime_samples_path))

    current_rep = 0
    while current_rep < nreps:
        for host_regime_idx, host_regime in enumerate(hrs.host_regimes):
            simulation_name="Run_{}_{}".format(current_rep+1, host_regime_idx+1)
            run_output_prefix = "{}.R{:04d}.H{:04d}".format(output_prefix, current_rep+1, host_regime_idx+1)
            run_logger.info("-inphest- Replicate {} of {}, host regime {} of {}: Starting".format(current_rep+1, nreps, host_regime_idx+1, len(hrs.host_regimes)))
            num_restarts = 0
            while True:
                if num_restarts == 0 and current_rep == 0:
                    is_verbose_setup = True
                    model_setup_logger = run_logger
                else:
                    is_verbose_setup = False
                    model_setup_logger = None
                inphest_model = model.InphestModel.create(
                        host_regime=host_regime,
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
                    run_logger.info("-inphest- Replicate {} of {}, host regime {} of {}: Simulation failure before termination condition at t = {}: {}".format(current_rep+1, nreps, host_regime_idx+1, len(hrs.host_regimes), inphest_simulator.elapsed_time, e))
                    num_restarts += 1
                    if num_restarts > maximum_num_restarts_per_replicates:
                        run_logger.info("-inphest- Replicate {} of {}, host regime {} of {}: Maximum number of restarts exceeded: aborting".format(current_rep+1, nreps, host_regime_idx+1, len(hrs.host_regimes)))
                        break
                    else:
                        run_logger.info("-inphest- Replicate {} of {}, host regime {} of {}: Restarting replicate (number of restarts: {})".format(current_rep+1, nreps, host_regime_idx+1, len(hrs.host_regimes), num_restarts))
                else:
                    run_logger.system = None
                    run_logger.info("-inphest- Replicate {} of {}, host regime {} of {}: Completed to termination condition at t = {}".format(current_rep+1, nreps, host_regime_idx+1, len(hrs.host_regimes), inphest_simulator.elapsed_time))
                    num_restarts = 0
                    break
            current_rep += 1

if __name__ == "__main__":
    rb_data = os.path.join(utility.TEST_DATA_PATH, "revbayes", "bg_large.events.txt")
    repeat_run(
        output_prefix="test",
        nreps=1,
        host_regime_samples_path=rb_data,
        model_definition_source={},
        model_definition_type="python-dict",
        config_d=None,
        interpolate_missing_model_values=False,
        random_seed=None,
        stderr_logging_level="info",
        file_logging_level="debug",
        maximum_num_restarts_per_replicates=100,
        debug_mode=False)

    # hrs = model.HostRegimeSamples()
    # rb_data = os.path.join(utility.TEST_DATA_PATH, "revbayes", "bg_large.events.txt")
    # rb_data_src = open(rb_data, "r")
    # hrs.parse_host_biogeography(rb_data_src)
    # for host_regime in hrs.host_regimes:
    #     im = model.InphestModel.create(
    #             host_regime=host_regime,
    #             model_definition_source={},
    #             model_definition_type="python-dict",
    #             )
