from __future__ import print_function, division
import argparse
from experiment import Experiment


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--exp-name', type=str, required=True,
            help='the name of the json experiment file to run')
    parser.add_argument('-overwrite', action='store_true', default=False,
            help='forces submission of the relevant jobs even if results \
                    already appear to exist.')
    parser.add_argument('-debug', action='store_true', default=False,
            help='prevents actual submission of the jobs')
    args = parser.parse_args()

    print('loading experiment', args.exp_name)
    exp = Experiment(args.exp_name)

    print('submitting runs')
    for s in exp.simulations:
        for e in exp.estimators_and_truth:
            e.submit_runs(s, overwrite=args.overwrite, debug=args.debug)
