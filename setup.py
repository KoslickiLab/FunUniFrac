import os
from setuptools import setup, find_packages

# Utility function to read the README file.
# Used for the long_description.  It's nice, because now 1) we have a top level
# README file and 2) it's easier to type in the README file than to put a raw
# string in below ...

# SCRIPTS = []
# SCRIPTS.extend([os.path.join("scripts", script)
# 				for script in os.listdir(os.path.join(os.path.dirname(__file__), "scripts"))
# 				if script.endswith(".py")])

HERE = os.path.abspath(os.path.dirname(__file__))
with open(os.path.join(HERE, 'README.md'), 'r') as fid:
	LONG_DESCRIPTION = fid.read()

setup(
	name="fununifrac",
	version="0.0.1",
	author="David Koslicki",
	author_email="dmkoslicki@gmail.com",
	description=("A repository to implement UniFrac, but on functional profiles of metagenomic data."),
	long_description=LONG_DESCRIPTION,
	#license="BSD-3-Clause",  # see classifiers
	keywords="unifrac kegg emd genomics metagenomics",
	url="https://github.com/KoslickiLab/FunUniFrac",
	packages=find_packages(),
	install_requires=[
	    'blist',
        'scipy==1.8.0',
        'networkx==2.8.4',
        'numpy==1.23.2',
        'pandas==1.4.3',
        'pyemd==0.5.1',
        'sparse',
        'requests',
        'seaborn'
    ],
	zip_safe=False,
	# package_data={'CMash': ['data/*.fna', 'tests/Organisms/*.fna.gz']},
	# scripts=SCRIPTS,
	entry_points={
        'console_scripts': [ 
            'compute_fununifrac.py = fununifrac.compute_fununifrac:main',
	    	'compute_edges.py = fununifrac.compute_edges:main',
		    'create_edge_matrix.py = fununifrac.create_edge_matrix:main'
        ]
    },
	classifiers=[
		"Development Status :: 3 - Alpha",
		"Topic :: Scientific/Engineering :: Bio-Informatics",
		"Topic :: Scientific/Engineering :: Mathematics",
		"License :: OSI Approved :: BSD License",
		"Intended Audience :: Science/Research",
		"Programming Language :: Python :: 3.6",
		"Programming Language :: Python :: 3.7",
		"Programming Language :: Python :: 3.8",
		"Programming Language :: Python :: 3.9",
		"Natural Language :: English",
		"Operating System :: MacOS :: MacOS X",
		"Operating System :: POSIX :: Linux",
	]
)