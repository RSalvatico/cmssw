#ifndef EcalSimAlgos_EcalDigitizerTraits_Ph2_h
#define EcalSimAlgos_EcalDigitizerTraits_Ph2_h

#include "DataFormats/EcalDigi/interface/EcalDigiCollections_Ph2.h"
#include "SimCalorimetry/EcalSimAlgos/interface/EcalElectronicsSim_Ph2.h"
#include "CalibFormats/CaloObjects/interface/CaloTSamples_Ph2.h"
#include "DataFormats/EcalDigi/interface/EcalDataFrame.h"
#include "DataFormats/EcalDigi/interface/EBDataFrame.h"
#include "CondFormats/EcalObjects/interface/EcalConstants.h"
class EcalHitResponse ;

class EBDigitizerTraits_Ph2 {
   public:
      /// the digis collection
  typedef EBDigiCollectionPh2 DigiCollection;
      /// the dataframes
      typedef EBDataFrame Digi;
      /// the electronics simulation
      typedef EcalElectronicsSim_Ph2 ElectronicsSim_Ph2;

      typedef CaloTSamples_Ph2<float,ecalPh2::sampleSize> EcalSamples ;

      static void fix( Digi& digi, edm::DataFrame df ) {};
};

#endif

