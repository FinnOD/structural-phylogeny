from setuptools import setup
from pathlib import Path

this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text(encoding='utf-8')

setup(
	name='structphy',
	version='0.0.1',
	description="Create phylogenetic trees from protein structure",
	long_description=long_description,
	long_description_content_type='text/markdown',
	author="Finn O'Donoghue",
	author_email='finnodonoghue@gmail.com',
	url='https://github.com/finnod/structural-phylogeny',
	
	keywords=[
		'protein',
		'informatics',
		'bioinformatics',
	],
	classifiers=[
		'Development Status :: 4 - Beta',
		'Environment :: Console',
		'Topic :: Scientific/Engineering :: Bio-Informatics',
		'Topic :: Scientific/Engineering :: Information Analysis',
		'Topic :: Scientific/Engineering :: Mathematics',
		'Intended Audience :: End Users/Desktop',
		'Operating System :: MacOS :: MacOS X',
		'Operating System :: Microsoft :: Windows',
		'Operating System :: POSIX',
		'Programming Language :: Python',
		'Programming Language :: Python :: 3 :: Only',
		'Programming Language :: Python :: 3'
	],

	packages=['structphy'],
	entry_points={
		'console_scripts': [
			'structphy=structphy.__main__:main',
		],
	},
	package_data={'structphy': [
		'example_data/*'
		]
	},

	install_requires=[
        'click',
        'requests',
        'numpy',
        'pandas',
        'tqdm',
	],
)