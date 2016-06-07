from __future__ import print_function, division
import argparse
import pandas as pd
from estimator import Estimator
import gprim.annotation as pa
import paths


class TestMethodTypeA(Estimator):
    parser = argparse.ArgumentParser(add_help=False, parents=[Estimator.parser])
    parser.add_argument('--myparam', type=float, required=True,
            help='the value of the magical thing!')
    parser.add_argument('--sannot_chr', type=str, required=True,
            help='path to an annotation file')

    def fsid(self):
        return 'tmA-myp{}'.format(
                self.params.myparam)

    def required_files(self, s):
        a = pa.Annotation(paths.annotations + self.params.sannot_chr)
        return [a.filestem(c)+'.testprocess' for c in s.chromosomes]

    def preprocess(self, s):
        print('TestMethod is preprocessing', s.name,
                'with refpanel=', self.params.refpanel)

    def run(self, s, beta_num):
        print('TestMethodA is running', s.name, 'on beta', beta_num,
                'with refpanel=', self.params.refpanel)
        print(self.params)
        return pd.DataFrame(columns=['ESTIMATE', 'STDERR'])

class TestMethodTypeB(TestMethodTypeA):
    def fsid(self):
        return 'tmB-myp{}'.format(
                self.params.myparam)

    def run(self, s, beta_num):
        print('TestMethodB is running', s.name, 'on beta', beta_num,
                'with refpanel=', self.params.refpanel)
        print(self.params)
        return pd.DataFrame(columns=['ESTIMATE', 'STDERR'])
