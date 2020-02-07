#ifndef EcalSimAlgos_EcalElectronicsSim_Ph2_h
#define EcalSimAlgos_EcalElectronicsSim_Ph2_h 1

#include "CalibFormats/CaloObjects/interface/CaloTSamples_Ph2.h"
#include "SimCalorimetry/CaloSimAlgos/interface/CaloVNoiseSignalGenerator.h"
#include "CondFormats/EcalObjects/interface/EcalConstants.h"


class EcalCoder_Ph2;
class EcalDataFrame;
class EcalSimParameterMap_Ph2;

namespace CLHEP {
  class HepRandomEngine;
}

/* \class EcalElectronicsSim_Ph2
 * \brief Converts CaloDataFrame in CaloTimeSample and vice versa.
 * 
 */

class EcalElectronicsSim_Ph2 {
public:
  typedef CaloTSamples_Ph2<float, ecalPh2::sampleSize> EcalSamples;
  //typedef CaloTSamples_Ph2<float, 10> EcalSamples;
  EcalElectronicsSim_Ph2(const EcalSimParameterMap_Ph2* parameterMap,
                     EcalCoder_Ph2* coder,
                     bool applyConstantTerm,
                     double rmsConstantTerm);

  ~EcalElectronicsSim_Ph2();

  /// from EcalSamples to EcalDataFrame
  void analogToDigital(CLHEP::HepRandomEngine*, EcalSamples& clf, EcalDataFrame& df) const;

  void newEvent() {}

  void setNoiseSignalGenerator(const CaloVNoiseSignalGenerator* noiseSignalGenerator) {
    theNoiseSignalGenerator = noiseSignalGenerator;
  }

private:
  /// input signal is in pe.  Converted in GeV
  void amplify(EcalSamples& clf, CLHEP::HepRandomEngine*) const;

  /// map of parameters

  const EcalSimParameterMap_Ph2* m_simMap;

  const CaloVNoiseSignalGenerator* theNoiseSignalGenerator;

  EcalCoder_Ph2* m_theCoder;

  const double m_thisCT;
  const bool m_applyConstantTerm;
};

#endif
