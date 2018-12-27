# Recentrifuge
**Robust comparative analysis and contamination removal for metagenomics**
____
[![Build Status](https://travis-ci.org/khyox/recentrifuge.svg?branch=master)](https://travis-ci.org/khyox/recentrifuge)
[![](https://img.shields.io/github/languages/top/khyox/recentrifuge.svg)](https://pypi.org/project/recentrifuge/)
[![](https://img.shields.io/pypi/pyversions/recentrifuge.svg)](https://pypi.org/project/recentrifuge/)
[![](https://img.shields.io/pypi/v/recentrifuge.svg)](https://pypi.org/project/recentrifuge/)
[![](https://img.shields.io/pypi/wheel/recentrifuge.svg)](https://pypi.org/project/recentrifuge/)

[![](https://img.shields.io/maintenance/yes/2019.svg)](http://www.recentrifuge.org)
[![](https://img.shields.io/badge/platforms-linux%20%7C%20macos%20%7C%20win-lightgrey.svg)](http://www.recentrifuge.org)
[![](https://img.shields.io/github/languages/count/khyox/recentrifuge.svg)](http://www.recentrifuge.org)
[![](https://img.shields.io/website-up-down-green-red/http/www.recentrifuge.org.svg?label=recentrifuge.org)](http://www.recentrifuge.org)



### Robust comparative analysis

With Recentrifuge, researchers can interactively explore what organisms are in their samples and at which level of confidence, thus enabling a robust comparative analysis of multiple samples in any metagenomic study.

Recentrifuge enables researchers to analyze results from [Centrifuge](http://www.ccb.jhu.edu/software/centrifuge/), [LMAT](https://computation.llnl.gov/projects/livermore-metagenomics-analysis-toolkit), [CLARK](http://clark.cs.ucr.edu/), [Kraken](http://ccb.jhu.edu/software/kraken/), and many other taxonomic classifiers using interactive pie charts, by placing great emphasis on the confidence level (score) of the taxonomic classifications.

In addition to the scored charts for the original samples, different sets of scored plots are generated for each taxonomic level of interest, like shared taxa and exclusive taxa per sample.

### Robust contamination removal

If there are one or more negative control samples in the study, Recentrifuge will also generate additional control-subtracted interactive plots with the contamination removed from every sample and from the shared taxa one. 

The novel and robust contamination removal algorithm of Recentrifuge detects and selectively removes various types of pollutants, including crossovers.

### Low biomass or reliable minority detection

Recentrifuge can analyze any metagenomic dataset. However, it is especially useful in the case of low microbial biomass studies and when a more reliable detection of minority organisms is needed, like in clinical, environmental and forensic applications.

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
For further details, please check the [bioRxiv pre-print](https://doi.org/10.1101/190934).
____
