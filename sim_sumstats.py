from __future__ import print_function, division
import argparse, subprocess, os
import pandas as pd
import scipy.stats as stats
import numpy as np
from pyutils import bsub, pretty
import metadata as sm

# create beta and noiseless phenotype for each chromosome
def create_beta_and_profiles(s, beta_num):
    for chrnum in s.chromosomes:
        # create beta, normalize by maf since genotypes will be centered but not standardized
        print('creating beta for chr', chrnum)
        beta = s.architecture.draw_beta(chrnum)
        beta.to_csv(s.beta_file(beta_num, chrnum, mode='w'), index=False, sep='\t')
        maf = s.dataset.frq_df(chrnum)['MAF'].values
        beta['BETA'] /= np.sqrt(2 * maf * (1-maf)) # convert beta to plink's per-allele scale

        # write sparse beta as well (plink works faster if beta is explicitly sparse)
        sparsebeta = beta.loc[beta['BETA'] != 0]
        sparsebeta.to_csv(s.sparsebeta_file(beta_num,chrnum,mode='w'),
                index=False, sep='\t')
        print('norm of beta{}*sqrt(2maf(1-maf)) equals {}'.format(
            chrnum, np.linalg.norm(beta['BETA'] * np.sqrt(2*maf * (1-maf)))**2))
        print('beta{} is length {} and has support of size {}'.format(
            chrnum, len(beta), len(sparsebeta)))

        # call plink to create profile file
        cmd = ['plink',
                '--bfile', s.dataset.bfile(chrnum),
                '--allow-no-sex',
                '--score', s.sparsebeta_filename(beta_num, chrnum), '1', '2', '4',
                'sum',
                'center',
                '--out', s.chr_filestem(beta_num, chrnum)]

        print('executing', ' '.join(cmd))
        subprocess.call(cmd)

# sum the noiseless phenotypes across chromosomes
def make_noiseless_pheno(s, beta_num):
    def get_profile(fname):
        df = pd.read_csv(fname, header=0, delim_whitespace=True)
        df.drop(['PHENO', 'CNT', 'CNT2'], axis=1, inplace=True)
        return df

    print('merging phenotypes. chr', s.chromosomes[0])
    phenotype = get_profile(s.noiselessYchr_filename(beta_num, s.chromosomes[0]))
    phenotype.rename(columns={'SCORESUM' : 'PHENO'}, inplace=True)
    for chrnum in s.chromosomes[1:]:
        print('merging phenotypes. chr', chrnum)
        profile = get_profile(s.noiselessYchr_filename(beta_num, chrnum))
        phenotype = pd.merge(phenotype, profile, how='inner', on=['FID', 'IID'])
        phenotype['PHENO'] += phenotype['SCORESUM']
        phenotype.drop(['SCORESUM'], axis=1, inplace=True)

    phenotype.to_csv(s.noiselessY_filename(beta_num),
            index=False,
            sep='\t')
    return phenotype

# add noise to resulting phenotype
def add_noise_and_save(s, beta_num, phenotype):
    sigma2e = 1-s.h2g
    print('adding noise. sigma2e =', sigma2e)
    phenotype['PHENO'] += np.sqrt(sigma2e) * np.random.randn(len(phenotype))
    phenotype.to_csv(s.noisyY_filename(beta_num),
            index=False,
            sep='\t')

# call plink to compute sumstats
def make_qassoc(s, beta_num):
    # compute one set of sumstats per chromosome
    for chrnum in s.chromosomes:
        print('computing sumstats for chr', chrnum)
        try:
            cmd = ['plink',
                    '--bfile', s.dataset.bfile(chrnum),
                    '--pheno', s.noisyY_filename(beta_num),
                    '--allow-no-sex',
                    '--linear', 'hide-covar',
                    '--covar', s.covariate_file] + s.covariate_plinkinfo + \
                    ['--out', s.chr_filestem(beta_num, chrnum)]
        except AttributeError: # means there are no covariates specified
            cmd = ['plink',
                    '--bfile', s.dataset.bfile(chrnum),
                    '--pheno', s.noisyY_filename(beta_num),
                    '--allow-no-sex',
                    '--assoc',
                    '--out', s.chr_filestem(beta_num, chrnum)]
        print('executing', ' '.join(cmd))
        subprocess.call(cmd)

def make_sumstats(s, beta_num):
    # create one large df from the qassoc files.
    print('merging chromosomes of sumstats')
    sumstats = []
    for chrnum in s.chromosomes:
        print(chrnum)
        if os.path.exists(s.chr_filestem(beta_num, chrnum)+'.assoc.linear'):
            assoc = pd.read_csv(s.chr_filestem(beta_num, chrnum)+'.assoc.linear',
                header=0,
                delim_whitespace=True)
            assoc['A2'] = s.dataset.bim_df(chrnum).A2 # A1 is already in there
            sumstats.append(assoc[['SNP', 'A1', 'A2', 'STAT', 'NMISS']].rename(
                columns={'STAT':'T'}))
        else: # means there were no covariates
            qassoc = pd.read_csv(s.chr_filestem(beta_num, chrnum)+'.qassoc',
                header=0,
                delim_whitespace=True)
            qassoc['A1'] = s.dataset.bim_df(chrnum)['A1']
            qassoc['A2'] = s.dataset.bim_df(chrnum)['A2']
            sumstats.append(qassoc[['SNP', 'A1', 'A2', 'T', 'NMISS']])

    sumstats = pd.concat(sumstats, axis=0)

    # filter down to the snps that we want
    if s.print_snps != 'none':
        print('filtering to', s.print_snps)
        rsids = pd.read_csv(s.print_snps, header=None, names=['SNP'])
        sumstats = pd.merge(sumstats, rsids, on='SNP', how='inner')
        print(len(sumstats), 'snps left')

    # since sample size is large, we assume that the T statistic and the Z score are
    # approximately equal.
    sumstats.rename(columns={'T':'Z', 'NMISS':'N'}, inplace=True)

    if np.sum(np.isnan(sumstats['Z'])) > 0 or np.sum(np.isinf(sumstats['Z'])) > 0:
        print('ERROR: some summary statistics were either inf or nan. aborting')
        print(np.where(np.isnan(sumstats['Z']))[0])
        print(np.where(np.isinf(sumstats['Z']))[0])
        return

    # save the big sumstats file
    print('saving result to', s.sumstats_filename(beta_num))
    sumstats.to_csv(s.sumstats_file(beta_num, mode='w'), sep='\t', index=False)

    # delete intermediate files
    if s.low_space:
        import gzip, shutil
        print('compressing/deleting intermediate files')
        for chrnum in s.chromosomes:
            fname = s.chr_filestem(beta_num, chrnum)+'.betanz'
            with open(fname, 'rb') as f_in, gzip.open(fname+'.gz', 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)

            for suff in ['beta.gz','betanz','log','nosex','profile','qassoc','assoc.linear']:
                fname = s.chr_filestem(beta_num, chrnum)+'.'+suff
                print('removing', fname)
                if os.path.exists(fname): os.remove(fname)
        fname = s.noiselessY_filename(beta_num)
        print('removing', fname)
        if os.path.exists(fname): os.remove(fname)


def main(args):
    np.random.seed(args.beta_num)
    s = sm.Simulation(args.sim_name)
    create_beta_and_profiles(s, args.beta_num)
    phenotype = make_noiseless_pheno(s, args.beta_num)
    add_noise_and_save(s, args.beta_num, phenotype)
    make_qassoc(s, args.beta_num)
    make_sumstats(s, args.beta_num)

def submit(args):
    s = sm.Simulation(args.sim_name)
    my_args = ['--sim-name', args.sim_name] + \
            ['main',
            '--beta-num', '$LSB_JOBINDEX']
    outfilename = s.beta_folder('%I', create=False) + 'sim_sumstats.out'
    bsub.submit(['python', '-u', __file__] + my_args,
            outfilename,
            # queue='medium', time_in_hours=40,
            queue='short', time_in_hours=12,
            jobname=args.sim_name+'.simsumstats[1-{}]'.format(s.num_betas),
            debug=args.debug)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--sim-name', type=str, required=True,
            help='name of the set of the json file of simulation parameters, without ext')

    main_parser, submit_parser = bsub.add_main_and_submit(parser, main, submit)

    main_parser.add_argument('--beta-num', type=int, required=True,
            help='index of the beta to simulate. 1-indexed!')
    submit_parser.add_argument('-debug', action='store_true', default=False)

    bsub.choose_parser_and_run(parser)


