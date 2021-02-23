import time
from ..constants import logger
from .cache import cached


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
