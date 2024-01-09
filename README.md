<p align="center"><a href="http://www.recentrifuge.org" target="_blank">
<img src="https://raw.githubusercontent.com/khyox/rcf-aux/master/RCFheader.png" alt="Recentrifuge" width="900px"/></a></p><hr>
<p align="center"><b>Robust comparative analysis and contamination removal for metagenomics</b>
</p> 


____
[![Retest](https://github.com/khyox/Recentrifuge/actions/workflows/retest.yaml/badge.svg?branch=v1.13.1)](https://github.com/khyox/recentrifuge/actions/workflows/retest.yaml)
[![](https://img.shields.io/maintenance/yes/2024.svg)](http://www.recentrifuge.org)
[![](https://img.shields.io/github/languages/top/khyox/recentrifuge.svg)](https://pypi.org/project/recentrifuge/)
[![](https://img.shields.io/pypi/pyversions/recentrifuge.svg)](https://pypi.org/project/recentrifuge/)
[![](https://img.shields.io/pypi/v/recentrifuge.svg)](https://pypi.org/project/recentrifuge/)
[![](https://img.shields.io/pypi/wheel/recentrifuge.svg)](https://pypi.org/project/recentrifuge/)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/recentrifuge/badges/version.svg)](https://anaconda.org/bioconda/recentrifuge)

[![](https://img.shields.io/badge/platforms-linux%20%7C%20macos%20%7C%20win-lightgrey.svg)](http://www.recentrifuge.org)
[![](https://img.shields.io/github/languages/count/khyox/recentrifuge.svg)](http://www.recentrifuge.org)
[![](https://img.shields.io/website-up-down-green-red/http/www.recentrifuge.org.svg?label=recentrifuge.org)](http://www.recentrifuge.org)
[![](https://img.shields.io/badge/Publication-PLOS_Computational_Biology-violet.svg)](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1006967)


### Robust comparative analysis

With Recentrifuge, researchers can interactively explore what organisms are in their samples and at which level of confidence, thus enabling a robust comparative analysis of multiple samples in any metagenomic study.

Recentrifuge enables researchers to analyze results from [Centrifuge](http://www.ccb.jhu.edu/software/centrifuge/), [LMAT](https://computation.llnl.gov/projects/livermore-metagenomics-analysis-toolkit), [CLARK](http://clark.cs.ucr.edu/), [Kraken](http://ccb.jhu.edu/software/kraken/), and many other taxonomic classifiers using interactive pie charts, by placing great emphasis on the confidence level (score) of the taxonomic classifications. 

The arithmetic of scored taxonomic trees of Recentrifuge supports the 48 taxonomic ranks of the [NCBI Taxonomy](https://www.ncbi.nlm.nih.gov/taxonomy), including several infraspecific levels such as strain or isolate. In addition to the scored charts for the original samples, different sets of scored plots are generated for each taxonomic level of interest, like shared taxa and exclusive taxa per sample.

### Robust contamination removal

If there are one or more negative control samples in the study, Recentrifuge will also generate additional control-subtracted interactive plots with the contamination removed from every sample and from the shared taxa specimen. 

The novel and robust contamination removal algorithm of Recentrifuge detects and selectively removes various types of pollutants, including crossovers.

### Reliable detection of minority organisms in complex datasets

Recentrifuge quickly analyzes complex metagenomic datasets using parallel computational algorithms. It is especially useful in the case of low microbial biomass studies and when a more reliable detection of minority organisms is needed, like in clinical, environmental, and forensic applications. 

Beyond the standard confidence levels, Recentrifuge implements others devoted to variable length reads, very convenient for datasets generated by nanopore sequencers.

____
For usage and docs, please, see [the Recentrifuge wiki](https://github.com/khyox/recentrifuge/wiki) or, depending on your classifier, [installation](https://github.com/khyox/recentrifuge/wiki/Installation) and:
 * [Recentrifuge for Centrifuge users](https://github.com/khyox/recentrifuge/wiki/Running-recentrifuge-for-Centrifuge)
 * [Recentrifuge for LMAT users](https://github.com/khyox/recentrifuge/wiki/Running-recentrifuge-for-LMAT)
 * [Recentrifuge for CLARK users](https://github.com/khyox/recentrifuge/wiki/Running-recentrifuge-for-CLARK)
 * [Recentrifuge for Kraken users](https://github.com/khyox/recentrifuge/wiki/Running-recentrifuge-for-Kraken)
 * [Recentrifuge for users of other classifiers](https://github.com/khyox/recentrifuge/wiki/Running-recentrifuge-for-a-generic-classifier)
____
To play with an example of webpage generated by Recentrifuge, click on the screenshot: 

<p align="center">
  <a href="https://rawgit.com/khyox/rcf-aux/master/TEST.rcf.html?dataset=5&node=0&collapse=false&color=true&depth=30&font=12&key=true" target="_blank">
    <img src="https://raw.githubusercontent.com/khyox/rcf-aux/master/RCF_screenshot_750.png" alt="Recentrifuge test screenshot" width="750px"/></a></p>
<p align="center">
  <a href="https://rawgit.com/khyox/rcf-aux/master/TEST.rcf.html?dataset=5&node=0&collapse=false&color=true&depth=30&font=12&key=true" target="_blank">Recentrifuge webpage example</a><p align="center">

____
Recentrifuge is in a stable development status. For further details, please check the [PLOS CB article](https://doi.org/10.1371/journal.pcbi.1006967).

If you use Recentrifuge in your research, please cite the paper. Thanks!

Martí JM (2019) **Recentrifuge: Robust comparative analysis and contamination removal for metagenomics**. _PLOS Computational Biology_ 15(4): e1006967.  https://doi.org/10.1371/journal.pcbi.1006967
____
