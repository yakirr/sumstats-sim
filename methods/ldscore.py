from __future__ import print_function, division
import subprocess, os
import argparse
import numpy as np
import pandas as pd
from pyutils import bsub
from estimator import Estimator
import primitives.annotation as pa
import paths


class LDSC(Estimator):
    parser = argparse.ArgumentParser(add_help=False, parents=[Estimator.parser])
    parser.add_argument('--ld_window', type=int, required=True,
            help='the maximal ld window to consider, in units corresponding to ' + \
                    'the refpanels bim files')
    parser.add_argument('--baseline', type=int, required=True,
            help='1 to use the baseline model, 0 to not')
    parser.add_argument('--constrain_intercept', type=int, required=False, default=1,
            help='1 to constrain the intercept, 0 to not')
    parser.add_argument('--weights_chr', type=str, required=True,
            help='path to a set of l2ldscoregz files with weights for the regression')
    parser.add_argument('--annot_chr', type=str, required=True,
            help='path to a set of annotgz files, not including chromosome number')
    parser.add_argument('--coeff', type=str, required=True,
            help='the name of the coefficient to report')

    def fsid(self):
        return 'LDSC,b={},c={},coeff={},A={},w={},ref={},ldw={}'.format(
                self.params.baseline,
                self.params.constrain_intercept,
                self.params.coeff,
                self.params.annot_chr.replace('/','_'),
                self.params.weights_chr.replace('/','_'),
                self.params.refpanel,
                self.params.ld_window)

    # NOTE: the function below assumes that the baseline model ldscores have already
    # been computed if they are necessary
    def required_files(self, s):
        annots = [pa.Annotation(paths.annotations + aname)
                for aname in self.params.annot_chr.split(',')]
        return [a.ldscores_filename(c) for c in range(1,23) for a in annots]

    def preprocess(self, s):
        print('LDSC is preprocessing', s.name,
                'with refpanel=', self.params.refpanel)
        print(self.params)

        for annot_chr in self.params.annot_chr.split(','):
            a = pa.Annotation(paths.annotations + annot_chr)
            for c in range(1,23):
                if not os.path.exists(a.ldscores_filename(c)):
                    ldscores_command = [
                            'python', '-u', paths.foreign + 'ldsc/ldsc.py',
                            '--l2',
                            '--ld-wind-cm', str(self.params.ld_window),
                            '--bfile', self.refpanel.bfile(c),
                            '--annot', a.annot_filename(c),
                            '--out', a.filestem(c)]
                    print(' '.join(ldscores_command))
                    outfilepath = a.filestem(c) + '.ldscoresbsub_out'
                    bsub.submit(
                            ldscores_command,
                            outfilepath,
                            jobname=self.params.annot_chr.replace('/','_') + \
                                    ',ldcores,chr='+str(c))

    #TODO: implement baseline model
    def run(self, s, beta_num):
        print('LDSC is running', s.name, 'on beta', beta_num,
                'with refpanel=', self.params.refpanel)
        print(self.params)
        annots = [pa.Annotation(paths.annotations + aname)
                for aname in self.params.annot_chr.split(',')]

        if self.params.constrain_intercept:
            extra_args = ['--no-intercept']
        else:
            extra_args = []
        cmd = [
                'python', '-u', paths.foreign + 'ldsc/ldsc.py',
                '--h2', s.sumstats_filename(beta_num),
                '--ref-ld-chr', ','.join([a.stem_chr for a in annots]),
                '--w-ld-chr', self.refpanel.bfile_chr,
                '--overlap-annot',
                '--print-coefficients',
                '--chisq-max', '999999',
                '--frqfile-chr', self.refpanel.bfile_chr,
                '--out', self.result_filename(s, beta_num)] + extra_args
        print(' '.join(cmd))
        subprocess.call(cmd)

        ldscresults = pd.read_csv(self.result_filename(s, beta_num)+'.results',
                delim_whitespace=True,
                header=0)
        rowindex = np.where(np.concatenate([
            a.names(22) for a in annots]) == self.params.coeff)[0][0]
        estimate = ldscresults['Prop._h2'][rowindex]
        stderr = ldscresults['Prop._h2_std_error'][rowindex]
        pval = ldscresults['Enrichment_p'][rowindex]

        return pd.DataFrame(columns=['ESTIMATE', 'STDERR', 'P'],
                data=[[estimate, stderr, pval]])

    # # returns tuple of (ld scores (M x n_annot), weights (M x 1), number of snps (1 x n_annot)
    # def ld_score_info(self):
    #     cols = [self.params.region+'L2'] + (['OTHERL2'] if not self.params.baseline else [])
    #     ldscores = []
    #     M_annot = np.zeros((len(cols) +
    #         (0 if not self.params.baseline else len(LDSC.baseline_model_regions)),))

    #     # ld scores and number of snps in each category
    #     for chrnum in self.refpanel.chromosomes():
    #         df = pd.read_csv(self.annotation_l2_filestem(chrnum) + '.l2.ldscore.gz', sep='\t',
    #                 usecols=cols)
    #         M_annot_chr = np.loadtxt(self.annotation_l2_filestem(chrnum)+'.l2.M')[:len(cols)]

    #         if self.params.baseline:
    #             df_baseline = pd.read_csv(
    #                     self.baseline_l2_filestem(chrnum) + '.l2.ldscore.gz', sep='\t',
    #                     usecols=[r+'L2' for r in LDSC.baseline_model_regions])
    #             M_annot_baseline_chr = np.loadtxt(
    #                     self.baseline_l2_filestem(chrnum)+'.l2.M')
    #             M_annot_chr = np.concatenate([M_annot_chr, M_annot_baseline_chr])
    #             df = pd.concat([df, df_baseline], axis=1)

    #         ldscores.append(np.array(df))
    #         M_annot += M_annot_chr
    #     ldscores = np.concatenate(ldscores)

    #     # weights
    #     if self.params.baseline:
    #         w_ld = ldscores[:,1].reshape((-1, 1))
    #     else:
    #         w_ld = np.sum(ldscores, axis=1).reshape((-1,1))

    #     return ldscores, w_ld, M_annot.reshape((1,-1))

    # def overlap_vector(self):
    #     counts = np.zeros((
    #         2 + (len(LDSC.baseline_model_regions)-1 if self.params.baseline else 0),))
    #     for chrnum in self.refpanel.chromosomes():
    #         annot = pd.read_csv(self.annotation_filename(chrnum), delim_whitespace=True,
    #                 compression='gzip')
    #         annot = annot.ix[:,4:]
    #         if self.params.baseline:
    #             annot = annot.ix[:,0]
    #             annot_baseline = pd.read_csv(self.baseline_filename(chrnum),
    #                     delim_whitespace=True, compression='gzip')
    #             annot_baseline = annot_baseline.ix[:,4:]
    #             annot = pd.concat([annot, annot_baseline], axis=1)
    #         counts += np.dot(annot.ix[:,0], annot)
    #     return counts

    # def preprocess(self):
    #     if self.params.baseline and not self.baseline_preprocessing_in_progress():
    #         print('baseline model not found. creating...')
    #         self.declare_baseline_preprocessing_in_progress()
    #         self.create_baseline_model()

    #     print('submitting ld score jobs for annotation of interest')
    #     gs = GenomicSubset(self.params.region)

    #     # create the annotation file
    #     for chrnum in self.refpanel.chromosomes():
    #         d = Dataset(self.params.refpanel, chrnum=chrnum)
    #         ss = SnpSubset(d, gs.restricted_to_chrom_bedtool(chrnum))
    #         SnpSubset.print_subsets(self.annotation_filename(chrnum),
    #                 [ss], [self.params.region], add_other=True)

    #     # create the ldscores file
    #     for chrnum in self.refpanel.chromosomes():
    #         d = Dataset(self.params.refpanel, chrnum=chrnum)
    #         ldscores_command = [
    #                 'python', '-u', paths.foreign + 'ldsc/ldsc.py',
    #                 '--l2',
    #                 '--ld-wind-cm', str(self.params.ld_window / 1000.),
    #                 '--bfile', d.genotypes_bedfile.filename,
    #                 '--annot', self.annotation_filename(chrnum),
    #                 '--out', self.annotation_l2_filestem(chrnum)]
    #         print(' '.join(ldscores_command))
    #         outfilepath = self.annotation_l2_filestem(chrnum) + '.bsub_out'
    #         bsub.submit(
    #                 ldscores_command,
    #                 outfilepath,
    #                 jobname=self.preprocessing_foldername()+',chr='+str(chrnum))

    # def run(self, beta_num, sim):
    #     print('loading data set and region info')
    #     d = Dataset(sim.dataset)
    #     gs = GenomicSubset(self.params.region)
    #     ss = SnpSubset(d, bedtool=gs.bedtool)

    #     print('loading ld score info')
    #     ref_ldscores, w_ld, M_annot = self.ld_score_info()
    #     N = np.ones((d.M, 1)) * d.N

    #     print(('ref_ldscores shape:{}\nw_ld shape:{}\nN shape:{}\n' + \
    #             'M_annot shape:{}').format(
    #                 ref_ldscores.shape,
    #                 w_ld.shape,
    #                 N.shape,
    #                 M_annot.shape))

    #     overlaps = self.overlap_vector()
    #     print('num snps overlapping with each category:', overlaps)
    #     results = []
    #     variances = []
    #     for alphahat in sim.sumstats_files(beta_num):
    #         alphahat = d.N * alphahat ** 2
    #         if self.params.constrain_intercept:
    #             hsqhat = ldsc.ldscore.regressions.Hsq(
    #                     alphahat.reshape((d.M,1)),
    #                     ref_ldscores,
    #                     w_ld,
    #                     N,
    #                     M_annot,
    #                     intercept=1)
    #         else:
    #             hsqhat = ldsc.ldscore.regressions.Hsq(
    #                     alphahat.reshape((d.M,1)),
    #                     ref_ldscores,
    #                     w_ld,
    #                     N,
    #                     M_annot)
    #         results.append(hsqhat.coef.dot(overlaps))
    #         variances.append(overlaps.dot(hsqhat.coef_cov).dot(overlaps))
    #         print('intercept:', hsqhat.intercept)
    #         print(len(results), results[-1], variances[-1])

    #     return np.concatenate([np.array([results]).T, np.array([variances]).T], axis=1)


if __name__ == '__main__':
    import sim.metadata as sm
    est = LDSC(refpanel='1000G3.wim5u', ld_window=1, baseline=0,
            annot_chr='1000G3.wim5u/mock/,1000G3.wim5u/all/',
            weights_chr='1000G3.wim5u/all/', coeff='mock')
    s = sm.Simulation('test')
    # est.preprocess(s)
    est.run(s, 1)

