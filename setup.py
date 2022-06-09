from setuptools import setup

setup(
    name='recentrifuge',
    version='1.8.2',
    packages=['recentrifuge'],
    url='http://www.recentrifuge.org',
    license='AGPL except krona.js, with its own license by BNBI',
    author='Jose Manuel Martí',
    author_email='jse.mnl@gmail.com',
    description='Robust comparative analysis and contamination removal for metagenomics',
    long_description="""
**Robust comparative analysis and contamination removal for metagenomics**

[![Retest](https://github.com/khyox/Recentrifuge/actions/workflows/retest.yaml/badge.svg?branch=v1.8.2)](https://github.com/khyox/recentrifuge/actions/workflows/retest.yaml)

With Recentrifuge, researchers can interactively explore what organisms are in their samples and at which level of confidence, thus enabling a robust comparative analysis of multiple samples in any metagenomic study.

 * Removes diverse contaminants, including crossovers, using a novel **robust contamination removal** algorithm.
 * Provides a **confidence level for every result**, since the calculated score propagates to all the downstream analysis and comparisons.
 * Unveils the generalities and specificities in the metagenomic samples, thanks to a new **comparative analysis engine**.
 
Recentrifuge's novel approach combines robust statistics, arithmetic of scored taxonomic trees, and parallel computational algorithms.

Recentrifuge is especially useful when a more **reliable detection of minority organisms** is needed (e.g. in the case of low microbial biomass metagenomic studies) in clinical, environmental, or forensic analysis. Beyond the standard confidence levels, Recentrifuge implements others devoted to variable length reads, very convenient for complex datasets generated by **nanopore sequencers**.

____

To play with an example of a webpage generated by Recentrifuge, click on the next screenshot: 

<p align="center">
  <a href="https://rawgit.com/khyox/rcf-aux/master/TEST.rcf.html?dataset=5&node=0&collapse=false&color=true&depth=30&font=12&key=true" target="_blank">
    <img src="https://raw.githubusercontent.com/khyox/rcf-aux/master/RCF_screenshot_750.png" alt="Recentrifuge test screenshot" width="750px"/></a></p>
<p align="center">
  <a href="https://rawgit.com/khyox/rcf-aux/master/TEST.rcf.html?dataset=5&node=0&collapse=false&color=true&depth=30&font=12&key=true" target="_blank">Recentrifuge webpage example</a><p align="center">
""",
    long_description_content_type='text/markdown',
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Environment :: Console',
        'Environment :: Web Environment',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU Affero General Public License v3 or later (AGPLv3+)',
        'Natural Language :: English',
        'Operating System :: OS Independent',
        'Programming Language :: JavaScript',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Scientific/Engineering :: Medical Science Apps.'
    ],
    keywords='metagenomics comparative-analysis contamination-removal',
    project_urls={
        'Documentation': 'https://github.com/khyox/recentrifuge/wiki',
        'Source': 'https://github.com/khyox/recentrifuge',
        'Tracker': 'https://github.com/khyox/recentrifuge/issues',
    },
    scripts=['rcf', 'rextract', 'retaxdump', 'remock', 'retest'],
    install_requires=['biopython'],
    python_requires='>=3.6',
    include_package_data=True,
    package_data={
        'recentrifuge': ['img/*.uri', 'test/*.mck', 'test/*.xlsx'],
    },
)
