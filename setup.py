#!/usr/bin/env python

import distutils.sysconfig, os, sys

from   distutils.core      import setup
from   distutils.extension import Extension
from   Cython.Distutils    import build_ext

INCLUDE_DIRS = [os.path.join(os.environ['SAGE_ROOT'], 'local/lib/python2.5/site-packages/', 'numpy/core/include/numpy/')]

setup( name        = 'pyca',
       version     = '0.2.1',
       description = 'SAGE/Python Cellular Automata Toolkit',
       author      = 'Iztok Jeras',
       author_email= 'iztok.jeras@rattus.info',
       url         = 'http://code.google.com/p/cellular-automata-sage-toolkit/',
       py_modules  = [ 'pyca', 'ca_vizual' ],
       ext_modules = [ Extension ('ca1d', ['ca1d.pyx'], include_dirs=INCLUDE_DIRS),
                       Extension ('ca2d', ['ca2d.pyx'], include_dirs=INCLUDE_DIRS) ],
       cmdclass    = {'build_ext': build_ext}
     )
