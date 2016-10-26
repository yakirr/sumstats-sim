from __future__ import print_function, division
import numpy as np
import pandas as pd
import json
from collections import defaultdict
from gprim.annotation import Annotation
import paths


class EffectSizeDist(object):
    def __init__(self, name):
        self.name = name
        self.__dict__.update(json.load(
            open(paths.architectures + 'effect_distributions/' + name + '.json')))

    # returns an array of shape (num_snps,)
    def draw_effect_sizes(self, num_snps, per_snp_variance):
        if self.causal_dist == 'normal':
            result = np.zeros(num_snps)
            causal_mask = np.random.binomial(
                    1,
                    self.p_causal,
                    size=num_snps).astype(bool)
            num_causal = np.sum(causal_mask)
            result[causal_mask] = \
                    np.sqrt(per_snp_variance / self.p_causal) * np.random.randn(num_causal)
            return result
        else:
            print('error, unsupported causal effect-size distribution')
            return None

class Architecture(object):
    def __init__(self, name):
        self.name = name
        self.__dict__.update(
                json.load(open(paths.architectures + name + '.json')))
        self.annotations = {}
        for annot_file in self.annot_files:
            self.annotations[annot_file] = Annotation(paths.annotations+annot_file)

    def set_pheno_var(self, h2g, chroms):
        self.sqnorms = pd.DataFrame()
        self.sizes = pd.DataFrame()
        for a in self.annotations.values():
            print('loading size/norm info on', a.stem_chr)
            if a.names(chroms[0]) in self.sqnorms.columns.values:
                print('WARNING: duplicate annotation column names!')
            self.sqnorms[a.names(chroms[0])] = a.total_sqnorms(chroms)
            self.sizes[a.names(chroms[0])] = a.total_sizes(chroms)
        print('done')
        self.h2g = h2g

    # returns a tuple consisting of a dict of mean effects and a dict of variance effects
    def params(self):
        return (defaultdict(float,
            { n: np.sqrt(e['h2g_explained']/self.sqnorms[n][0] * self.h2g) * e['sign']
                for n, e in self.mean_effects.items()}),
            defaultdict(float,
            { n: e['h2g_explained'] / self.sizes[n][0] * self.h2g
                for n,e in self.variance_effects.items()}))


    def __find_annot(self, n, chrnum, signed=True):
        matches = [a for a in self.annotations.values() if n in a.names(chrnum)]
        if signed:
            return matches[0].sannot_df(chrnum)[n].values # there should only be one match
        else:
            return matches[0].annot_df(chrnum)[n].values # there should only be one match

    # assumes that set_pheno_var has already been called
    def draw_beta(self, chrnum):
        result = self.annotations.values()[0].sannot_df(chrnum).copy()
        result = result[['SNP', 'A1', 'A2']]
        result['BETA'] = 0

        mean_effects, var_effects = self.params()

        # add variance effects
        for n, norm_var in var_effects.items():
            print('size of {} is {}'.format(n, self.sizes[n][0]))
            print('per-SNP variance is', norm_var)
            print('contributing variance of', norm_var * self.sizes[n][0], 'over all chr')
            esd = EffectSizeDist(self.variance_effects[n]['effectdist'])
            v = self.__find_annot(n, chrnum)
            result.ix[np.flatnonzero(v), 'BETA'] += esd.draw_effect_sizes(
                    np.count_nonzero(v), norm_var)

        # add mean effects
        for n, norm_mu in mean_effects.items():
            print('norm of {} is {}'.format(n, self.sqnorms[n][0]))
            print('per-SNP mu is', norm_mu)
            print('contributing variance of', norm_mu**2 * self.sqnorms[n][0],
                    'over all chr')
            result['BETA'] += self.__find_annot(n, chrnum) * norm_mu

        return result

if __name__ == '__main__':
    np.random.seed(0)
    esd = EffectSizeDist('sparse')
    print('below should be a sparse length-100 vector with large non-zero values')
    print(esd.draw_effect_sizes(100, 100))
    print()

    a = Architecture('mock_null_varenriched')
    a.set_pheno_var(0.5, [1, 22])
    beta = a.draw_beta(22)
    print(beta.columns)
