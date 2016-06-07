from __future__ import print_function, division
import argparse
from experiment import Experiment


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--exp-name',
            help='the name of the json experiment file whose results we wish to verify')
    parser.add_argument('-v', action='store_true', default=False,
            help='print the actual filenames of the missing result files')
    parser.add_argument('-s', action='store_true', default=False,
            help='suppress information corresponding to sim/estimator pairs with no problem')
    parser.add_argument('-l', action='store_true', default=False,
            help='print exactly 1 line per missing file, so one can pipe to wc')
    args = parser.parse_args()

    if not args.l: print('loading experiment', args.exp_name)
    exp = Experiment(args.exp_name)

    if not args.l: print('verifying results')
    def verify(s):
        if not args.l: print(s.name)
        for est in exp.estimators_and_truth:
            if est.missing_results(s):
                if not args.l: print('XXX', str(est))
                if args.v or args.l:
                    print('\n\t'.join(est.missing_results(s)))
                if not args.l: print('')
            else:
                if not args.s and not args.l:
                    print('vvv', str(est))
                    print('')
        if not args.l: print('')
    map(verify, exp.simulations)

