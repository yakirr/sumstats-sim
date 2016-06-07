from __future__ import print_function
import numpy as np
import pandas as pd
from experiment import Experiment
import argparse

def get_power(exp, sig_levels, exclude=[]):
    def tsvfilename(s):
        return '{}{}.power.tsv'.format(exp.results_folder(), s.name)
    def latexfilename(s):
        return '{}{}.power.tex'.format(exp.results_folder(), s.name)

    for s in exp.simulations:
        print(s.name)
        power_array = []
        for est in exp.estimators:
            if est.params.pretty_name in exclude:
                continue
            results_df = est.results(s)
            row = [est.params.pretty_name]
            for sig_level in sig_levels:
                row.append(np.mean(results_df['P'] < sig_level))
            row.append(np.median(results_df['P']))
            power_array.append(row)
        df = pd.DataFrame(
                columns=['Estimators'] + \
                        ['P(reject), $\\alpha='+str(sl)+'$' for sl in sig_levels] + \
                        ['median p-val'],
                data=power_array)
        print(df)
        print()
        df.to_csv(tsvfilename(s), sep='\t', index=False)
        with open(latexfilename(s), 'w') as f:
            f.write(df.to_latex(index=False, escape=False))

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--exp-name', help='name of the json experiment file')
    parser.add_argument('--sig-levels', nargs='+', type=float, help='significance level')
    parser.add_argument('--exclude', nargs='+', default=[],
            help='prettynames of estimators to ignore')
    args = parser.parse_args()

    exp = Experiment(args.exp_name)
    get_power(exp, args.sig_levels, exclude=args.exclude)
