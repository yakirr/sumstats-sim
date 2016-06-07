from __future__ import print_function, division
import numpy as np
import pandas as pd
import argparse
import matplotlib.pyplot as plt
from experiment import Experiment


def create_plot(exp, s, to_plot, all_results, plotestimate=False, pretty=False,
        drop_mses=True):
    print('inside create_plot')
    mses = {e:'{:.1e}'.format(np.mean(all_results[e]['BIAS']**2)) for e in to_plot}

    def label(e):
        if pretty:
            stem = e.params.pretty_name.replace(',','\n')
        else:
            stem = e.fsid().replace(',','\n')
        suffix = '' if drop_mses else ('\nMSE:' + mses[e])
        return stem + suffix
    def filename(stdaxes=False):
        return '{}{}.biasplot'.format(exp.results_folder(), s.name) + \
                ('.axes' if stdaxes else '') + '.png'

    scorename = ('ESTIMATE' if plotestimate else 'BIAS')
    horizontal_line = np.mean(all_results.values()[0]['TRUTH']) if plotestimate else 0

    # create box plots
    plt.figure()
    plt.boxplot([all_results[e][scorename] for e in to_plot],
        labels=[label(e) for e in to_plot],
        sym='',
        widths=0.75)

    # add individual points with jitter
    for i, e in enumerate(to_plot):
        # x = np.random.normal(1+i, 0.06, size=len(error_lists[e]))
        delta = 0.6/len(all_results[e])
        x = (1 + i) + np.arange(-0.3, 0.3, delta) + delta/2
        plt.plot(x, all_results[e][scorename], 'r.', alpha=1)
        plt.plot([1+i], np.mean(all_results[e][scorename]), 'k.')

    # formatting
    plt.xticks(rotation=35)
    plt.axhline(y=horizontal_line)
    if plotestimate:
        plt.ylabel('estimate')
    else:
        plt.ylabel('error')
    if not pretty:
        plt.title(s.name)
    fig = plt.gcf()
    if plotestimate and horizontal_line != 0:
        fig.gca().set_ylim(bottom=0, top=2*horizontal_line)
    fig.set_size_inches(2 + len(to_plot) * 1.5, 10)
    fig.subplots_adjust(bottom=0.25)

    print('saving figure')
    fig.savefig(filename(), dpi=300)
    print('showing figure')
    plt.show()
    # fig.gca().set_ylim([-0.1, 0.1])
    # print('saving with standardized axes')
    # fig.savefig(filename(stdaxes=True), dpi=300)

def write_results(exp, s, all_results):
    def filename(est):
        return '{}{}.{}.biases.tsv'.format(exp.results_folder(), s.name, est.fsid())

    for est in all_results.keys():
        all_results[est].to_csv(filename(est), sep='\t', index=False)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--exp-name',
            help='the name of the json experiment file whose results we wish to plot')
    parser.add_argument('-pretty', action='store_true', default=False,
            help='print short pretty names for methods')
    parser.add_argument('-drop-mses', action='store_true', default=False,
            help='dont print mses')
    parser.add_argument('--exclude', nargs='+', default=[],
            help='list of pretty names of methods to exclude')
    parser.add_argument('-plotestimate', action='store_true', default=False,
            help='rather than plotting errors relative to 0, plot estimates ' + \
                    'relative to average value of truth.')
    args = parser.parse_args()

    print('loading experiment', args.exp_name)
    exp = Experiment(args.exp_name)
    if not hasattr(exp, 'exclude'):
        exp.exclude = []

    if not args.pretty:
        plt.rcParams.update({'axes.titlesize': 'medium'})
        plt.rcParams.update({'font.size': 7})
    else:
        plt.rcParams.update({'axes.titlesize': 'medium'})
        plt.rcParams.update({'font.size': 9})

    print('plotting results')
    def plot_results_for(s):
        true_results = exp.truth.results(s)
        true_results.rename(columns={'ESTIMATE':'TRUTH'}, inplace=True)
        true_results.drop(['P','STDERR'], axis=1, inplace=True)
        all_results = {}
        to_plot = []
        for est in exp.estimators:
            if est.params.pretty_name in args.exclude:
                continue
            to_plot.append(est)
            results = est.results(s)
            merged = pd.merge(results, true_results, how='inner', on='BETA_NUM')
            merged['BIAS'] = merged['ESTIMATE'] - merged['TRUTH']
            merged = merged[['BETA_NUM', 'ESTIMATE', 'TRUTH', 'BIAS']]
            all_results[est] = merged
        write_results(exp, s, all_results)
        create_plot(exp, s, to_plot, all_results,
                plotestimate=args.plotestimate, pretty=args.pretty,
                drop_mses=args.drop_mses)
    map(plot_results_for, exp.simulations)

