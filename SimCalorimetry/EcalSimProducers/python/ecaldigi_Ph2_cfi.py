#raise RuntimeError, "Do not import obsolete file ecalDigi_cfi.py. If you need parameters, use 'from SimGeneral.MixingModule.ecalDigitizer_cfi import *'"
import FWCore.ParameterSet.Config as cms
print "[ecaldigi_Ph2_cfi.py] including SimGeneral.MixingModule.ecalDigitizer_Ph2_cfi"
from SimGeneral.MixingModule.ecalDigitizer_Ph2_cfi import *
