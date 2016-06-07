from __future__ import print_function, division
import itertools
import json
import sim.metadata as sm
import sim.methods as smethods
from pyutils import fs
import paths


class Estimators(object):
    def __init__(self, estimators_json):
        self.estimators = []
        for e_def in estimators_json:
            self.estimators.append(Estimators.create_estimator_from_json(e_def))

    def __iter__(self):
        return self.estimators.__iter__()

    @classmethod
    def create_estimator_from_json(cls, json_entry, debug=False):
        method_name = json_entry['method']; del json_entry['method']
        if debug: print('creating estimator', method_name)
        return smethods.find_method(method_name)(**json_entry)


class Simulations(object):
    def __init__(self, simulation_names):
        self.simulations = []
        for sim_name in simulation_names:
            self.simulations.append(sm.Simulation(sim_name))

    def __iter__(self):
        return self.simulations.__iter__()


class Experiment(object):
    def __init__(self, name):
        self.name = name
        self.__dict__.update(json.load(
            open(paths.experiments + name + '.json')))
        self.estimators = Estimators(self.estimators)
        self.truth = Estimators.create_estimator_from_json(self.truth)
        self.simulations = Simulations(self.simulations)

        if not hasattr(self, 'exclude'):
            self.exclude = []
        self.estimators = list(filter(
            lambda e: e.params.pretty_name not in self.exclude,
            self.estimators))

    def results_folder(self, create=True):
        path = paths.results + self.name + '/'
        if create:
            fs.makedir(path)
        return path

    @property
    def purpose_filename(self):
        return self.results_folder() + 'purpose.txt'

    @property
    def estimators_and_truth(self):
        return itertools.chain(self.estimators, [self.truth])

if __name__ == '__main__':
    exp = Experiment('2016.04.13.1_testmethod')

    for s in exp.simulations:
        for e in exp.estimators_and_truth:
            print(str(e), str(s.name))
