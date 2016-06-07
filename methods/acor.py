from __future__ import print_function, division
import subprocess, os
import argparse
import numpy as np
import pandas as pd
from pyutils import bsub
from estimator import Estimator
import gprim.annotation as pa
import sim.fcor.paths as paths


class Acor(Estimator):
    parser = argparse.ArgumentParser(add_help=False, parents=[Estimator.parser])
    parser.add_argument('--annot_chr', type=str, required=True,
            help='comma-delimited set of paths to annotgz files, ' + \
                    'not including chromosome number and extension')
    parser.add_argument('--coeff', type=str, required=True,
            help='the name of the coefficient to report')
    parser.add_argument('--kind', type=str, required=True,
            help='re for random effects, fe for fixed')
    parser.add_argument('--fullconv', type=int, default=0,
            help='0 means use ld blocks in computing conv files, 1 means dont')
    parser.add_argument('--weights', type=int, default=1,
            help='1 means weight the RE regression, 0 means dont')
    parser.add_argument('--maf_thresh', type=float, default=0,
            help='threshold for which snps to include in RE regression')
    parser.add_argument('--biascorrect', type=int, default=0,
            help='1 means bias correct the denominator of the RE regression')
    parser.add_argument('--regvar', type=int, default=0,
            help='1 means report a std. error based on a random-beta model even for acorfe')
    parser.add_argument('--print', type=str, default='all',
            help='all means print normal estimate. num means print numerator. denom ' + \
                    'means print denominator')

    def fsid(self):
        return 'Acor.{},coeff={},A={}'.format(
            self.params.kind,
            self.params.coeff,
            self.params.annot_chr.replace('/','_')) + \
                (',fc' if self.params.fullconv else '') + \
                (',nowt' if not self.params.weights else '') + \
                (',num' if self.params.print == 'num' else '') + \
                (',den' if self.params.print == 'denom' else '') + \
                (',bc' if self.params.biascorrect else '') + \
                ('maf{:0.3f}'.format(self.params.maf_thresh)
                        if self.params.maf_thresh>0 else '') + \
                (',rv' if self.params.regvar else '')

    # NOTE: the function below assumes that the ldscores for the reference panel have
    # already been computed
    def required_files(self, s):
        annots = [pa.Annotation(paths.annotations + aname)
                for aname in self.params.annot_chr.split(',')]
        return [a.conv_filename(c, full=self.params.fullconv)
                for c in s.chromosomes for a in annots]

    def preprocess(self, s):
        print('Acor is preprocessing', s.name,
                'with refpanel=', self.params.refpanel)
        print(self.params)

        annots = [pa.Annotation(paths.annotations + aname)
                for aname in self.params.annot_chr.split(',')]

        for a in annots:
            print('preprocessing', a.filestem())
            for c in s.chromosomes:
                if not os.path.exists(a.conv_filename(c, full=self.params.fullconv)):
                    conv_command = [
                        'python', '-u', paths.code + 'acor/acor.py',
                        '--annot-chr', a.stem_chr,
                        '--bfile-chr', self.refpanel.bfile_chr] + \
                        (['-fullconv'] if self.params.fullconv else []) + \
                        ['conv',
                        '--chroms', str(c)]
                    print(' '.join(conv_command))
                    outfilepath = a.filestem(c) + '.' + \
                            ('full' if self.params.fullconv else '') + \
                            'convbsub_out'
                    bsub.submit(
                            conv_command,
                            outfilepath,
                            jobname=self.params.annot_chr.replace('/','_') + \
                                    ',conv,chr='+str(c))

    def run(self, s, beta_num):
        print('Acor is running', s.name, 'on beta', beta_num,
                'with refpanel=', self.params.refpanel)
        print(self.params)
        annots = [pa.Annotation(paths.annotations + aname)
                for aname in self.params.annot_chr.split(',')]

        cmd = [
                'python', '-u', paths.code + 'acor/acor.py',
                '--annot-chr', ' '.join([a.stem_chr for a in annots]),
                '--bfile-chr', self.refpanel.bfile_chr] + \
                (['-fullconv'] if self.params.fullconv else []) + \
                ['cor',
                '--ldscores-chr', self.refpanel.bfile_chr,
                '--sumstats', s.sumstats_filename(beta_num),
                '--out', self.result_filename(s, beta_num),
                self.params.kind] + \
                (['-reg-var'] if self.params.regvar else []) + \
                (['-noweights'] if not self.params.weights else []) + \
                (['-biascorrect'] if self.params.biascorrect else []) + \
                (['--maf-thresh', str(self.params.maf_thresh)] if self.params.maf_thresh > 0
                        else []) + \
                ['--chroms'] + [str(c) for c in s.chromosomes]
        print(' '.join(cmd))
        subprocess.call(cmd)

        acorresults = pd.read_csv(self.result_filename(s, beta_num)+'.results',
                delim_whitespace=True,
                header=0)
        rowindex = np.where(np.concatenate([
            a.names(s.chromosomes[-1]) for a in annots]) == self.params.coeff)[0][0]
        if self.params.print == 'all':
            estimate = acorresults['MU_EST'][rowindex]
        elif self.params.print == 'num':
            estimate = acorresults['TOP'][rowindex]
        elif self.params.print == 'denom':
            estimate = acorresults['BOTTOM'][rowindex]
        stderr = acorresults['MU_STDERR'][rowindex]
        pval = acorresults['MU_P'][rowindex]

        return pd.DataFrame(columns=['ESTIMATE', 'STDERR', 'P'],
                data=[[estimate, stderr, pval]])


if __name__ == '__main__':
    import sim.metadata as sm
    est = Acor(refpanel='1000G3.wim5u',
            annot_chr='1000G3.wim5u/mock/',
            kind='re',
            coeff='mock')
    s = sm.Simulation('test')
    # est.preprocess(s)
    est.run(s, 1)

