#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/EcalDigi/interface/EcalDigiCollections.h"
#include "SimDataFormats/CaloHit/interface/PCaloHit.h"
#include "TH1.h"
#include "TH2.h"
#include "SimDataFormats/Track/interface/SimTrack.h"
#include <fstream>
#include "CondFormats/EcalObjects/interface/EcalConstants.h"

using namespace std;
//
// class declaration
//

class PhaseIAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
 public:
  explicit PhaseIAnalyzer(const edm::ParameterSet&);
  ~PhaseIAnalyzer();     

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


 private:
  virtual void beginJob() ;
  void analyze(const edm::Event&, const edm::EventSetup&) override ;
  virtual void endJob() ;

  edm::InputTag digiTagEB_;
  
  edm::EDGetTokenT<EBDigiCollection> digiTokenEB_; 
  

  edm::InputTag hitTagEB_;

  
  edm::EDGetTokenT<vector<PCaloHit>> hitTokenEB_; 


  edm::InputTag trackTag_;
  edm::EDGetTokenT<vector<SimTrack>> trackToken_; 
  
  //Histograms
  TH1I *EBEnergyHisto[ecalPh2::sampleSize];
TH1I *EBGainHisto[ecalPh2::sampleSize];


   TH2D* meEBDigiOccupancy_;

   TH1D* meEBDigiMultiplicity_;

  TH1D*  meEBDigiADCGlobal_;

   TH1I* SingleChannelE;

  TH1D*  meEBDigiADCAnalog_[ecalPh2::sampleSize];
  TH1D*  meEBDigiADCgS_[ecalPh2::sampleSize];
  TH1D*  meEBDigiGain_[ecalPh2::sampleSize]; 


 TH1D*  meEBPedestal_;
 
TH1D*   meEBMaximumgt100ADC_; 
 
 TH1D* meEBMaximumgt10ADC_; 
 
 TH1D*  meEBnADCafterSwitch_;
 








  
};
