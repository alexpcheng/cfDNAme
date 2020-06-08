#!/usr/bin/env python

#License: BSD

#Copyright (c) 2008 Li Charlie Xia
#All rights reserved.
#
#Redistribution and use in source and binary forms, with or without
#modification, are permitted provided that the following conditions
#are met:
#1. Redistributions of source code must retain the above copyright
#   notice, this list of conditions and the following disclaimer.
#2. Redistributions in binary form must reproduce the above copyright
#   notice, this list of conditions and the following disclaimer in the
#   documentation and/or other materials provided with the distribution.
#3. The name of the author may not be used to endorse or promote products
#   derived from this software without specific prior written permission.
#
#THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR
#IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
#OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
#IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT, INDIRECT,
#INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT
#NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
#DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
#THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
#(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
#THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

"""GRAMMY: GENOME Relative Abundance based on Mixture Model theorY.

This C++ python module provide tools for metagenomic abundance estimation
based on read probability partitions (alignment or composition based).

TODO:
"""

from setuptools import setup, find_packages
from distutils.core import Extension
from distutils.command import build
import os, sys

doclines=__doc__.splitlines()

if os.path.exists('MANIFEST'): os.remove('MANIFEST')

os.environ['CC'] = 'g++'  #temporary measure to trick distutils use g++, need update to distutils2

class my_build(build.build):
# different order: build_ext *before* build_py, so that 
# build_py can use ctypes! 
  sub_commands = [('build_ext', build.build.has_ext_modules),
                  ('build_py', build.build.has_pure_modules),
                  ('build_clib', build.build.has_c_libraries),
                  ('build_scripts', build.build.has_scripts), ]

setup(name="grammy",
  version="1.0.0",
  description=doclines[0],
  long_description="\n".join(doclines[2:]),
  author="Li Charlie Xia",
  author_email="lxia@usc.edu",
  url="http://meta.usc.edu/softs/grammy",
  license="BSD",
  platforms=["Linux"],
  packages=find_packages(exclude=['ez_setup', 'test', 'doc', 'dist', 'gem.egg-info']),
  include_package_data=True,
  zip_safe=False,
  install_requires=['python >= 2.7','scipy >= 0.6','numpy >= 1.0'],
  provides=['grammy'],
  ext_modules = [ Extension( 'grammy._gemcore', 
                      sources = ['grammy/gemcore_wrap.cpp', 'grammy/gemcore.cpp', 'grammy/mmio.c'],
                      depends=['grammy/gemcore.hpp', 'grammy/mmio.h'],
                      language ='c++',
                      libraries = ['stdc++'],  #trick to seduce distutils use g++ instead of gcc
                      #include_dirs = [ os.environ['HOME']+'/usr/include' ],
                      #library_dirs = [ os.environ['HOME']+'/usr/lib' ],
                  ),
                ],
  py_modules = ['grammy.gemcore', 'grammy.gemlib', 'grammy.gemaux', 'grammy.gemmath', 'grammy.gemutil'],
  cmdclass = {'build': my_build},
  #keywords = ('python', 'grammy', 'EM', 'Mixture Models'),
  #classifiers = [ 'Development Status :: 5 - Production/Stable',
  #                    'Environment :: Console',
  #                    'License :: OSI Approved :: BSD License',
  #                    'Programming Language :: Python',
  #                    'Programming Language :: C++',
  #              ],
  entry_points = { 'console_scripts': [
                      'grammy_rdt = grammy.grammy_rdt:main',
                      'grammy_gdt = grammy.grammy_gdt:main',
                      'grammy_pre = grammy.grammy_pre:main',
                      'grammy_em = grammy.grammy_em:main',
                      'grammy_post = grammy.grammy_post:main'
    ] 
  },
)
