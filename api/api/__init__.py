from ..utils import logger
try:
    from .embedding import *  # NOQA: F401
except Exception:
    logger.exception('Error importing embedding API')

try:
    from .molecule import MoleculeImageAPI  # NOQA: F401
except Exception:
    logger.exception('Error importing molecule API')

try:
    from .projection import *  # NOQA: F401
except Exception:
    logger.exception('Error importing projection API')

try:
    from .pso import PSOAPI  # NOQA: F401
except Exception:
    logger.exception('Error importing pso API')

try:
    from .interpolation import *  # NOQA: F401
except Exception:
    logger.exception('Error importing interpolation API')

try:
    from .mmp import *  # NOQA: F401
except Exception:
    logger.exception('Error importing mmp API')

try:
    from .stoned_selfies import *  # NOQA: F401
except Exception:
    logger.exception('Error importing stoned_selfies API')

try:
    from .registry import *  # NOQA: F401
except Exception:
    logger.exception('Error importing registry API')
