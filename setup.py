#!/usr/bin/env python3

#    setup.py: script to install TisOpt
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.

from setuptools import setup

setup(
    name="TisOpt",
    version="0.1",
    packages=['TisOpt'],
    package_dir={'TisOpt': 'TisOpt'},
    include_package_data=True,
    package_data={'TisOpt': ['data/*.csv'],'': ['LICENSE', '*.md']},
    python_requires=">=3.6",
    install_requires=[
        'numpy',
        'pandas',
        'RNA'
    ],
    author="Xavier Hernandez-Alias",
    author_email="xavier.hernandez@crg.eu",
    description="Codon optimizer for tissue-specific expression",
    url="https://github.com/hexavier/TisOpt",
)
