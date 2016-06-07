from __future__ import print_function, division
import gzip
import json
import pandas as pd
from pyutils import fs, memo
from phenotype import Architecture
import primitives.dataset as prd
import paths

class Dataset(prd.Dataset):
    def __init__(self, name, datasets_dict_path=paths.metadata+'datasets.json'):
        self.name=name
        self.__dict__.update(
                json.load(open(datasets_dict_path))[name])

class Simulation(object):
    def __init__(self, name, path=paths.simulations):
        self.name = name
        self.__dict__.update(
                json.load(open(path+name+'.json')))
        self.dataset = Dataset(self.dataset)

    @property
    @memo.memoized
    def architecture(self):
        self.__architecture = Architecture(self.architecture_name)
        self.__architecture.set_pheno_var(self.h2g, self.chromosomes)
        return self.__architecture

    def root_folder(self, create=True):
        folder = '{}{}/{}/'.format(paths.simsumstats,
                self.dataset.name, self.name)
        if create:
            fs.makedir(folder)
        return folder
    def beta_folder(self, beta_num, create=True):
        folder = '{}{}/'.format(self.root_folder(), beta_num)
        if create:
            fs.makedir(folder)
        return folder
    def beta_filename(self, beta_num, chrnum):
        return '{}{}.beta.gz'.format(
            self.beta_folder(beta_num), chrnum)
    def beta_file(self, beta_num, chrnum, mode='r'):
        return gzip.open(self.beta_filename(beta_num, chrnum), mode=mode)
    def sparsebeta_filename(self, beta_num, chrnum, mode='r'):
        return '{}{}.betanz'.format(
            self.beta_folder(beta_num), chrnum)
    def sparsebeta_file(self, beta_num, chrnum, mode='r'):
        return open(self.sparsebeta_filename(beta_num, chrnum), mode=mode)

    def chr_filestem(self, beta_num, chrnum):
        return '{}{}'.format(self.beta_folder(beta_num), chrnum)
    def noiselessYchr_filename(self, beta_num, chrnum):
        return self.chr_filestem(beta_num, chrnum) + '.profile'
    def noiselessY_filename(self, beta_num):
        return self.beta_folder(beta_num) + 'noiseless.pheno'
    def noisyY_filename(self, beta_num):
        return self.beta_folder(beta_num) + 'noisy.pheno'
    def sumstatschr_filename(self, beta_num, chrnum):
        return self.chr_filestem(beta_num, chrnum) + '.qassoc'
    def sumstats_filename(self, beta_num):
        return self.beta_folder(beta_num) + 'all.sumstats.gz'
    def sumstats_file(self, beta_num, mode='r'):
        return gzip.open(self.sumstats_filename(beta_num), mode=mode)


if __name__ == '__main__':
    d = Dataset('1000G3.wim5u')
    print(d.path)
    print(d.bfile_chr)
    print(d.bfile(22))

    s = Simulation('smalltest')
    print(s.architecture)
