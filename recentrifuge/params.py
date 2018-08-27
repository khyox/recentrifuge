"""
This module contains the parameters of the robust contamination removal method.

"""
#
# NOTE: See the online documentation and/or the article for information
#   about the robust contamination algorithm and these parameters.
#
# Epsilon for floating point comparisons with zero
EPS: float = 1e-14
# Min num of non-ctrl samples to enable advCtrl
ROBUST_MIN_SAMPLES: int = 1
# Cutoff for crossover outlier test, typically in [3, 5]
ROBUST_XOVER_OUTLIER = 5
# Relfreq order of magnitude dif in crossover test, typically in [2, 3]
ROBUST_XOVER_ORD_MAG = 3
# Min rel frequency of severe contaminant
SEVR_CONTM_MIN_RELFREQ: float = 0.01
# Min rel frequency of mild contaminant
MILD_CONTM_MIN_RELFREQ: float = 0.001
