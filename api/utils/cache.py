from ..constants import logger

_cache = {}


def cached(f):
    def wrapper(*args, **kw):
        key = f'{id(f)}_{str(args[1:])}_{str(kw)}'
        if key not in _cache:
            logger.info('Cache miss')
            _cache[key] = f(*args, **kw)
        else:
            logger.info('Cache hit')
        return _cache[key]
    return wrapper
