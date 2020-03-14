import FWCore.ParameterSet.Config as cms
print "loading mix simcastor digi"
print "In aliases_cfi many modifications without duplicate it"
simCastorDigis = cms.EDAlias(
    mix = cms.VPSet(
      cms.PSet(type = cms.string('CastorDataFramesSorted'))
    )
)
print "loading mix simEcalUnsuppressedDigi"
simEcalUnsuppressedDigis = cms.EDAlias(
    mix = cms.VPSet(
#        cms.PSet(type = cms.string('EBDigiCollectionPh2')),
        cms.PSet(type = cms.string('EBDigiCollection')),
        cms.PSet(type = cms.string('EEDigiCollection')),
        cms.PSet(type = cms.string('ESDigiCollection'))
    )
)
print "loading mix simHcalUnsuppressedDigi"
simHcalUnsuppressedDigis = cms.EDAlias(
    mix = cms.VPSet(
      cms.PSet(type = cms.string('HBHEDataFramesSorted')),
      cms.PSet(type = cms.string('HFDataFramesSorted')),
      cms.PSet(type = cms.string('HODataFramesSorted')),
      cms.PSet(type = cms.string('ZDCDataFramesSorted')),
      cms.PSet(type = cms.string('QIE10DataFrameHcalDataFrameContainer')),
      cms.PSet(type = cms.string('QIE11DataFrameHcalDataFrameContainer'))
    )
)
print "loading pixel common"
_pixelCommon = cms.VPSet(
    cms.PSet(type = cms.string('PixelDigiedmDetSetVector')),
    cms.PSet(type = cms.string('PixelDigiSimLinkedmDetSetVector'))
)
print "loading mix pixel digis"
simSiPixelDigis = cms.EDAlias(
    mix = _pixelCommon
) 
print "loading mix strp digis"
simSiStripDigis = cms.EDAlias(
    mix = cms.VPSet(
      cms.PSet(type = cms.string('SiStripDigiedmDetSetVector')),
      cms.PSet(type = cms.string('SiStripRawDigiedmDetSetVector')),
      cms.PSet(type = cms.string('StripDigiSimLinkedmDetSetVector'))
    )
)
print "loading mix hgcal digis"
simHGCalUnsuppressedDigis = cms.EDAlias(
    mix = cms.VPSet(
        cms.PSet(
            type = cms.string("DetIdHGCSampleHGCDataFramesSorted"),
            fromProductInstance = cms.string("HGCDigisEE"),
            toProductInstance = cms.string("EE"),
        ),
        cms.PSet(
            type = cms.string("DetIdHGCSampleHGCDataFramesSorted"),
            fromProductInstance = cms.string("HGCDigisHEfront"),
            toProductInstance = cms.string("HEfront"),
        ),
        cms.PSet(
            type = cms.string("DetIdHGCSampleHGCDataFramesSorted"),
            fromProductInstance = cms.string("HGCDigisHEback"),
            toProductInstance = cms.string("HEback"),
        ),
    )
)
print "loading mix HFNose"
simHFNoseUnsuppressedDigis = cms.EDAlias(
    mix = cms.VPSet(
        cms.PSet(
            type = cms.string("DetIdHGCSampleHGCDataFramesSorted"),
            fromProductInstance = cms.string("HFNoseDigis"),
            toProductInstance = cms.string("HFNose"),
        ),
    )
)

#print "loading mix simAPV saturation"
#simAPVsaturation = cms.EDAlias(
#    mix = cms.VPSet(
#        cms.PSet(type = cms.string('bool'))
#    )
#)

simAPVsaturation = cms.EDAlias(
    mixData = cms.VPSet(
        cms.PSet(
            type = cms.string('bool'),
            fromProductInstance = cms.string('siStripDigisDMSimulatedAPVDynamicGain'),
            toProductInstance = cms.string('SimulatedAPVDynamicGain'),
            )
    )
)




print "imprt run3 common"
from Configuration.Eras.Modifier_run3_common_cff import run3_common
run3_common.toModify(simCastorDigis, mix = None)
print "imprt phase2 hgcal"
from Configuration.Eras.Modifier_phase2_hgcal_cff import phase2_hgcal
(~phase2_hgcal).toModify(simHGCalUnsuppressedDigis, mix = None)

from Configuration.ProcessModifiers.premix_stage1_cff import premix_stage1
(premix_stage1 & phase2_hgcal).toModify(simHGCalUnsuppressedDigis,
    mix = {
        0 : dict(type = "PHGCSimAccumulator"),
        1 : dict(type = "PHGCSimAccumulator"),
        2 : dict(type = "PHGCSimAccumulator"),
    }
)
print "import phase2 nose"
from Configuration.Eras.Modifier_phase2_hfnose_cff import phase2_hfnose
(~phase2_hfnose).toModify(simHFNoseUnsuppressedDigis, mix = None)
print "import phase1 pixel"
from Configuration.Eras.Modifier_phase1Pixel_cff import phase1Pixel
phase1Pixel.toModify(simSiPixelDigis, mix = _pixelCommon + [cms.PSet(type = cms.string('PixelFEDChanneledmNewDetSetVector'))])

# no castor,pixel,strip digis in fastsim
print "fast sim"
from Configuration.Eras.Modifier_fastSim_cff import fastSim
fastSim.toModify(simCastorDigis, mix = None)
fastSim.toModify(simSiPixelDigis, mix = None)
fastSim.toModify(simSiStripDigis, mix = None)
fastSim.toModify(simAPVsaturation, mix = None)
