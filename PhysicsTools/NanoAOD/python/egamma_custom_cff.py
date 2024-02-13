import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.common_cff import *
from PhysicsTools.NanoAOD.electrons_cff import *
from PhysicsTools.NanoAOD.photons_cff import *
from PhysicsTools.NanoAOD.NanoAODEDMEventContent_cff import *
from PhysicsTools.NanoAOD.nanoDQM_cfi import nanoDQM
from PhysicsTools.NanoAOD.nanoDQM_cff import _Photon_extra_plots, _Electron_extra_plots

def addExtraEGammaVarsCustomize(process):
    #photon
    process.nanoTableTaskCommon.remove(process.photonTablesTask)
    process.nanoTableTaskCommon.add(process.photonTablesExtraTask)
    process.nanoDQM.vplots.Photon.plots = _Photon_extra_plots
    #electron
    process.nanoTableTaskCommon.remove(process.electronTablesTask)
    process.nanoTableTaskCommon.add(process.electronTablesExtraTask)
    process.nanoDQM.vplots.Electron.plots = _Electron_extra_plots
    return process
