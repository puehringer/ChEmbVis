import time
from ..constants import logger
import inspect

def _is_class_method(f):
    spec = inspect.getargspec(f)
    return spec.args and spec.args[0] == 'self'

_waiting = set()
_cache = {}

def cached(f):
    is_class_method = _is_class_method(f)
    def wrapper(*args, **kw):
        # TODO: Change args[1:] to args[0:] to include all parameters
        key = f'{id(f)}_{str(args[1:] if is_class_method else args[0:])}_{str(kw)}'
        if key not in _cache:
            # TODO: Find better way using Promise-like syntax
            # If it is already computing, wait for it...
            if key in _waiting:
                logger.info('Cache miss, waiting for value')
                # Wait for n-seconds and check if the value already exists in-between
                i = 0
                while key in _waiting and i < 600:
                    time.sleep(1)
                    i += 1
            # Compute the value if it still isn't in the cache
            if key not in _cache:
                logger.info('Cache miss, computing value')
                # Add the key to the waiting set
                _waiting.add(key)
                try:
                    # Compute the value
                    _cache[key] = f(*args, **kw)
                finally:
                    # Remove the key from the waiting set (i.e. "notify" other threads)
                    _waiting.remove(key)
        else:
            logger.info('Cache hit')
        return _cache[key]
    return wrapper
