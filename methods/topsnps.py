from __future__ import print_function, division
import argparse
import numpy as np
import scipy.stats
import pandas as pd
from estimator import Estimator
import gprim.annotation as pa
from pyutils import memo
import paths

class TopSnps(Estimator):
    parser = argparse.ArgumentParser(add_help=False, parents=[Estimator.parser])
    parser.add_argument('--coeff', type=str, required=True,
            help='the coefficient to report')
    parser.add_argument('--annot_chr', type=str, required=True,
            help='path to annotgz files, not incl chromosome number and extension')
    parser.add_argument('--vthresh', type=float, required=True,
            help='the threshold to use for absolute value of v')
    parser.add_argument('--pthresh', type=float, default=5e-8,
            help='the threshold to use for p-values')

    @property
    @memo.memoized
    def annotation(self):
        return pa.Annotation(paths.annotations + self.params.annot_chr)

    def fsid(self):
        return 'topsnps.coeff={},A={}'.format(
                self.params.coeff,
                self.params.annot_chr.replace('/','_')) + \
                ',pt={:0.1e}'.format(self.params.pthresh) + \
                ',vt={:0.1e}'.format(self.params.vthresh)

    def required_files(self, s):
        return []

    def preprocess(self, s):
        pass

    def run(self, s, beta_num):
        print('topsnps is running', s.name, 'on beta', beta_num,
                'with refpanel=', self.params.refpanel)
        print(self.params)

        print('reading v')
        v = pd.concat([
            self.annotation.sannot_df(c) for c in s.chromosomes], axis=0)
        v['VSIG'] = np.abs(v[self.params.coeff]) > self.params.vthresh

        print('reading sumstats')
        alphahat = pd.read_csv(s.sumstats_file(beta_num), delim_whitespace=True)

        print('computing pvalues')
        alphahat['P'] = scipy.stats.chi2.sf(alphahat['Z']**2, 1)
        alphahat['PSIG'] = alphahat['P'] < self.params.pthresh

        print('tallying counts')
        alphahat['BOTHSIG'] = v['VSIG'].values & alphahat['PSIG'].values
        alphahat['VSIGONLY'] = v['VSIG'].values & (~alphahat['PSIG'].values)
        alphahat['NOTPSIG'] = ~alphahat['PSIG'].values

        print('computing result')
        n1 = np.sum(alphahat['PSIG'])
        X1 = np.sum(alphahat['BOTHSIG'])
        n2 = np.sum(alphahat['NOTPSIG'])
        X2 = np.sum(alphahat['VSIGONLY'])
        if n1 == 0: # no top snps ==> make sure we get a null result
            X1 = X2
            n1 = n2
        print('X1={}, n1={}\nX2={}, n2={}'.format(X1, n1, X2, n2))
        p1hat = X1 / n1
        p2hat = X2 / n2
        phat = (X1 + X2) / (n1 + n2)
        stat = p1hat - p2hat
        var = phat*(1-phat)*(1/n1+1/n2)
        print([stat, np.sqrt(var), scipy.stats.norm.sf(stat/np.sqrt(var))])

        return pd.DataFrame(columns=['ESTIMATE', 'STDERR', 'P'],
                data=[[stat, np.sqrt(var), scipy.stats.norm.sf(stat/np.sqrt(var))]])


if __name__ == '__main__':
    import sim.metadata as sm
    s = sm.Simulation('mock_10xenrichednm')
    # est = TruthRE(refpanel='1000G3.wim5u', coeff='ANNOT', output='total')
    # est = TruthFE(refpanel='GERAimp.wim5u', coeff='ANNOT', annot_chr='GERAimp.wim5u/mock/')
    # print(est.run(s, 1))
    est = TopSnps(refpanel='GERAimp.wim5unm', coeff='ANNOT',
        annot_chr='GERAimp.wim5unm/mock/', vthresh=0.1, pthresh=5e-8)
    print(est.run(s, 1))

