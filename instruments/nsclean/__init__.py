from .NSClean import NSClean
from .NSCleanSubarray import NSCleanSubarray
from .util import make_subarray_filter
from .util import chsuf
from .util import make_lowpass_filter
from .util import png2fits

# Package constants
from .config import __version__
from .config import NY
from .config import NX
from .config import RB

# Say if we are using GPU acceleration
import os
if os.getenv('NSCLEAN_USE_CUPY') == 'YES':
    print('NSClean is using GPU acceleration')
    
