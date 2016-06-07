from __future__ import print_function, division
import argparse
import numpy as np
import pandas as pd
import os
from estimator import Estimator
import gprim.annotation as pa
from pyutils import memo, bsub
import sim.fcor.paths as paths

class TruthRE(Estimator):
    parser = argparse.ArgumentParser(add_help=False, parents=[Estimator.parser])
    parser.add_argument('--coeff', type=str, required=True,
            help='the coefficient to report')
    parser.add_argument('--annot_chr', type=str, required=True,
            help='path to annotgz files, not incl chromosome number and extension')
    parser.add_argument('--print', type=str, default='all',
            help='all means print normal estimate. num means print numerator. denom ' + \
                    'means print denominator')

    @property
    @memo.memoized
    def annotation(self):
        return pa.Annotation(paths.annotations + self.params.annot_chr)

    def fsid(self):
        return 'truthRE-{}'.format(
                self.params.coeff,
                self.params.annot_chr.replace('/','_')) + \
                    (',num' if self.params.print == 'num' else '') + \
                    (',den' if self.params.print == 'denom' else '')

    def required_files(self, s):
        if self.params.print == 'num' or self.params.print == 'denom':
            return [self.annotation.conv_filename(c, full=True)
                    for c in s.chromosomes]
        else:
            return []

    def preprocess(self, s):
        print('TruthRE is preprocessing', s.name,
                'with refpanel=', self.params.refpanel)
        print(self.params)

        print('preprocessing', self.annotation.filestem())
        for c in s.chromosomes:
            if not os.path.exists(self.annotation.conv_filename(c, full=True)):
                conv_command = [
                    'python', '-u', paths.code + 'acor/acor.py',
                    '--annot-chr', self.annotation.stem_chr,
                    '--bfile-chr', self.refpanel.bfile_chr,
                    '-fullconv',
                    'conv',
                    '--chroms', str(c)]
                print(' '.join(conv_command))
                outfilepath = self.annotation.filestem(c) + '.convbsub_out'
                bsub.submit(
                        conv_command,
                        outfilepath,
                        jobname=self.params.annot_chr.replace('/','_') + \
                                ',conv,chr='+str(c))

    def run(self, s, beta_num):
        print('TruthRE is running', s.name, 'on beta', beta_num,
                'with refpanel=', self.params.refpanel)
        print(self.params)

        mean_effects, var_effects = s.architecture.params()
        mu = mean_effects[self.params.coeff]
        if self.params.print == 'all':
            print('result:', mu)
            return pd.DataFrame(columns=['ESTIMATE', 'STDERR', 'P'],
                    data=[[mu, 0, 0]])
        else:
            print('reading Rv')
            Rv = pd.concat(
                [pd.read_csv(self.annotation.conv_filename(
                    chrom, full=True), sep = "\t")
                    for chrom in s.chromosomes],
                axis = 0)[self.params.coeff + '.conv1']
            # print('reading ldscores')
            # l2 = pd.concat(
            #     [pd.read_csv(self.refpanel.bfile(
            #         chrom) + '.l2.ldscore.gz', sep='\t')
            #         for chrom in s.chromosomes],
                # axis=0)['L2']

            print('computing result')
            result = np.dot(Rv, Rv)
            if self.params.print == 'num':
                print('numerator:', mu*result)
                return pd.DataFrame(columns=['ESTIMATE', 'STDERR', 'P'],
                        data=[[mu*result, 0, 0]])
            else:
                print('denominator:', result)
                return pd.DataFrame(columns=['ESTIMATE', 'STDERR', 'P'],
                        data=[[result, 0, 0]])


class TruthFE(Estimator):
    parser = argparse.ArgumentParser(add_help=False, parents=[Estimator.parser])
    parser.add_argument('--coeff', type=str, required=True,
            help='the coefficient to report')
    parser.add_argument('--annot_chr', type=str, required=True,
            help='path to annotgz files, not incl chromosome number and extension')
    parser.add_argument('--fullconv', type=int, default=0,
            help='0 means use ld blocks in computing conv files, 1 means dont')

    @property
    @memo.memoized
    def annotation(self):
        return pa.Annotation(paths.annotations + self.params.annot_chr)

    def fsid(self):
        return 'truthFE-{}-A{}'.format(
                self.params.coeff,
                self.params.annot_chr.replace('/','_')) + \
                    (',fc' if self.params.fullconv else '')

    def required_files(self, s):
        return [self.annotation.conv_filename(c, full=self.params.fullconv)
                for c in s.chromosomes]

    def preprocess(self, s):
        print('TruthFE is preprocessing', s.name,
                'with refpanel=', self.params.refpanel)
        print(self.params)

        print('preprocessing', self.annotation.filestem())
        for c in s.chromosomes:
            if not os.path.exists(self.annotation.conv_filename(c, full=self.params.fullconv)):
                conv_command = [
                    'python', '-u', paths.code + 'acor/acor.py',
                    '--annot-chr', self.annotation.stem_chr,
                    '--bfile-chr', self.refpanel.bfile_chr] + \
                    (['-fullconv'] if self.params.fullconv else []) + \
                    ['conv',
                    '--chroms', str(c)]
                print(' '.join(conv_command))
                outfilepath = self.annotation.filestem(c) + '.convbsub_out'
                bsub.submit(
                        conv_command,
                        outfilepath,
                        jobname=self.params.annot_chr.replace('/','_') + \
                                ',conv,chr='+str(c))

    def run(self, s, beta_num):
        print('TruthFE is running', s.name, 'on beta', beta_num,
                'with refpanel=', self.params.refpanel)
        print(self.params)

        print('reading beta')
        beta = pd.concat(
            [pd.read_csv(s.beta_file(beta_num, chrom), sep = "\t")
                for chrom in s.chromosomes],
            axis = 0)['BETA']
        print('reading Rv')
        Rv = pd.concat(
            [pd.read_csv(self.annotation.conv_filename(
                chrom, full=self.params.fullconv), sep = "\t")
                for chrom in s.chromosomes],
            axis = 0)[self.params.coeff + '.conv1']
        print('reading v')
        v = pd.concat(
            [self.annotation.sannot_df(chrom)
                for chrom in s.chromosomes],
            axis = 0)[self.params.coeff]
        print('computing result')
        result = np.dot(beta, Rv) / np.dot(v, Rv)
        return pd.DataFrame(columns=['ESTIMATE', 'STDERR', 'P'],
                data=[[result, 0, 0]])

class TruthFENoLD(Estimator):
    parser = argparse.ArgumentParser(add_help=False, parents=[Estimator.parser])
    parser.add_argument('--coeff', type=str, required=True,
            help='the coefficient to report')
    parser.add_argument('--annot_chr', type=str, required=True,
            help='path to sannotgz files, not incl chromosome number and extension')

    @property
    @memo.memoized
    def annotation(self):
        return pa.Annotation(paths.annotations + self.params.annot_chr)

    def fsid(self):
        return 'truthFENoLD-{}-A{}'.format(
                self.params.coeff,
                self.params.annot_chr.replace('/','_'))

    def required_files(self, s):
        return []

    def preprocess(self, s):
        pass

    def run(self, s, beta_num):
        print('TruthFE is running', s.name, 'on beta', beta_num,
                'with refpanel=', self.params.refpanel)
        print(self.params)

        print('reading beta')
        beta = pd.concat(
            [pd.read_csv(s.beta_file(beta_num, chrom), sep = "\t")
                for chrom in s.chromosomes],
            axis = 0)['BETA']
        print('reading v')
        v = pd.concat(
            [self.annotation.sannot_df(chrom)
                for chrom in s.chromosomes],
            axis = 0)[self.params.coeff]
        print('computing result')
        result = np.dot(beta, v) / np.dot(v, v)
        return pd.DataFrame(columns=['ESTIMATE', 'STDERR', 'P'],
                data=[[result, 0, 0]])


if __name__ == '__main__':
    import sim.metadata as sm
    s = sm.Simulation('mock_nullnm')
    # est = TruthRE(refpanel='1000G3.wim5u', coeff='ANNOT', output='total')
    # est = TruthFE(refpanel='GERAimp.wim5u', coeff='ANNOT', annot_chr='GERAimp.wim5u/mock/')
    # print(est.run(s, 1))
    est = TruthFENoLD(refpanel='GERAimp.wim5unm', coeff='ANNOT', annot_chr='GERAimp.wim5unm/mock/')
    print(est.run(s, 1))

