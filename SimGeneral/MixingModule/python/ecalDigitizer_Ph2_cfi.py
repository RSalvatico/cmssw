import FWCore.ParameterSet.Config as cms
print "[ecalDigitizer_Ph2_cfi] including SimCalorimetry.EcalSimProducers.ecalDigiParameters_Ph2_cff and SimCalorimetry.EcalSimProducers.ecalElectronicsSim_Ph2_cff"
from SimCalorimetry.EcalSimProducers.ecalDigiParameters_Ph2_cff import *
from SimCalorimetry.EcalSimProducers.apdSimParameters_cff import *
from SimCalorimetry.EcalSimProducers.ecalSimParameterMap_cff import *
from SimCalorimetry.EcalSimProducers.ecalElectronicsSim_Ph2_cff import *
from SimCalorimetry.EcalSimProducers.ecalNotContainmentSim_cff import *
from SimCalorimetry.EcalSimProducers.ecalCosmicsSim_cff import *

print "[ecalDigitizer_Ph2_cfi]:"
print "ecalDigitizer_Ph2 = cms.PSet("
print "ecal_digi_parameters,"
print "apd_sim_parameters"
print "ecal_electronics_sim,"
print "ecal_cosmics_sim,"
print "ecal_sim_parameter_map,"
print "ecal_notCont_sim,"
print "hitsProducer = cms.string(\'g4SimHits\'),"
print "accumulatorType = cms.string(\"EcalDigiProducer_Ph2\"),"
print "makeDigiSimLinks = cms.untracked.bool(False)" 
print ")"
ecalDigitizer_Ph2 = cms.PSet(
    ecal_digi_parameters,
    apd_sim_parameters,
    ecal_electronics_sim,
    ecal_cosmics_sim,
    ecal_sim_parameter_map,
    ecal_notCont_sim,
    hitsProducer = cms.string('g4SimHits'),
    accumulatorType = cms.string("EcalDigiProducer_Ph2"),
    makeDigiSimLinks = cms.untracked.bool(False)
)

from Configuration.Eras.Modifier_fastSim_cff import fastSim
fastSim.toModify(ecalDigitizer_Ph2, hitsProducer = "fastSimProducer")

print "[ecalDigitizer_Ph2_cfi]: ecalDigitizer_Ph2.doEB = cms.bool(True)"    
ecalDigitizer_Ph2.doEB = cms.bool(True)

