#ifndef DIGIECAL_ECALDIGICOLLECTION_PH2_H
#define DIGIECAL_ECALDIGICOLLECTION_PH2_H

#include "DataFormats/EcalDigi/interface/EBDataFrame.h"
#include "DataFormats/EcalDigi/interface/EEDataFrame.h"
#include "DataFormats/EcalDigi/interface/ESDataFrame.h"
#include "DataFormats/EcalDigi/interface/EcalTimeDigi.h"
#include "DataFormats/EcalDigi/interface/EcalTriggerPrimitiveDigi.h"
#include "DataFormats/EcalDigi/interface/EcalEBTriggerPrimitiveDigi.h"
#include "DataFormats/EcalDigi/interface/EcalTrigPrimCompactColl.h"
#include "DataFormats/EcalDigi/interface/EcalPseudoStripInputDigi.h"
#include "DataFormats/EcalDigi/interface/EBSrFlag.h"
#include "DataFormats/EcalDigi/interface/EESrFlag.h"
#include "DataFormats/EcalDigi/interface/EcalPnDiodeDigi.h"
#include "DataFormats/EcalDigi/interface/EcalMatacqDigi.h"
#include "DataFormats/Common/interface/SortedCollection.h"

#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/EcalDetId/interface/EcalSubdetector.h"
#include "DataFormats/Common/interface/DataFrameContainer.h"
#include "CondFormats/EcalObjects/interface/EcalConstants.h"

class EcalDigiCollectionPh2 : public edm::DataFrameContainer {
public:
  typedef edm::DataFrameContainer::size_type size_type;
  static const size_type MAXSAMPLES = ecalPh2::sampleSize;
  explicit EcalDigiCollectionPh2(size_type istride=MAXSAMPLES, int isubdet=0)  : 
    edm::DataFrameContainer(istride, isubdet){}
  void swap(DataFrameContainer& other) {this->DataFrameContainer::swap(other);}
  
};

// make edm (and ecal client) happy
class EBDigiCollectionPh2 : public  EcalDigiCollectionPh2 {
public:
  typedef edm::DataFrameContainer::size_type size_type;
  typedef EBDataFrame Digi;
  typedef Digi::key_type DetId;

  EBDigiCollectionPh2(size_type istride=MAXSAMPLES) : 
    EcalDigiCollectionPh2(istride, EcalBarrel){}
  void swap(EBDigiCollectionPh2& other) {this->EcalDigiCollectionPh2::swap(other);}
  void push_back(const Digi& digi){ DataFrameContainer::push_back(digi.id(), digi.frame().begin()); }
  void push_back(id_type iid){DataFrameContainer::push_back(iid);}
  void push_back(id_type iid,data_type const* idata){DataFrameContainer::push_back(iid,idata);}
  
};



// Free swap functions
inline
void swap(EcalDigiCollectionPh2& lhs, EcalDigiCollectionPh2& rhs) {
  lhs.swap(rhs);
}

inline
void swap(EBDigiCollectionPh2& lhs, EBDigiCollectionPh2& rhs) {
  lhs.swap(rhs);
}


typedef edm::SortedCollection<EcalTimeDigi> EcalTimeDigiCollection;
//DA ELIMINARE???
typedef edm::SortedCollection<EcalTriggerPrimitiveDigi> EcalTrigPrimDigiCollection;
typedef edm::SortedCollection<EcalEBTriggerPrimitiveDigi> EcalEBTrigPrimDigiCollection;

typedef edm::SortedCollection<EcalPseudoStripInputDigi> EcalPSInputDigiCollection;
typedef edm::SortedCollection<EBSrFlag> EBSrFlagCollection;
typedef edm::SortedCollection<EcalPnDiodeDigi> EcalPnDiodeDigiCollection;
typedef edm::SortedCollection<EcalMatacqDigi> EcalMatacqDigiCollection;

#endif
