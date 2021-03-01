import time
from ..constants import logger
from .cache import cached
from typing import List


pluck = lambda dict, *args: (dict[arg] for arg in args)


class catch_time():
    def __init__(self, name: str):
        self.name = name

    def __enter__(self):
        self.t = time.time()
        return self

    def __exit__(self, type, value, traceback):
        self.t = time.time() - self.t
        logger.info(f'{self.name} took {self.t:.3f}s')

def chunkify(l: List, n: int):
    for i in range(0, len(l), n):
        yield l[i:i+n]


# Utility function parallizing array operations
def parallelized(func, l):
    from joblib import Parallel, delayed
    import multiprocessing
    num_cores = multiprocessing.cpu_count()
    return Parallel(n_jobs=num_cores)(delayed(func)(m) for m in l)