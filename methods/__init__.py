import inspect
import testmethod
import ldscore
import acor
import truth
import topsnps

allmethods = inspect.getmembers(testmethod, inspect.isclass) + \
        inspect.getmembers(ldscore, inspect.isclass) + \
        inspect.getmembers(truth, inspect.isclass) + \
        inspect.getmembers(topsnps, inspect.isclass) + \
        inspect.getmembers(acor, inspect.isclass)

def find_method(name):
    for n, m in allmethods:
        if n == name:
            return m
    return None
