// -*- C++ -*-
//
// Package:    PhaseI/PhaseIAnalyzer
// Class:      PhaseIAnalyzer
// 
/**\class PhaseIAnalyzer PhaseIAnalyzer.cc PhaseI/PhaseIAnalyzer/plugins/PhaseIAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Dario Soldi
//         Created:  Tue, 16 Jan 2018 11:56:05 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

//My includes
#include "DataFormats/EcalDigi/interface/EcalDigiCollections.h"
#include "CondFormats/EcalObjects/interface/EcalPedestals.h"
#include "CondFormats/DataRecord/interface/EcalPedestalsRcd.h"
#include "CondFormats/DataRecord/interface/EcalIntercalibConstantsMCRcd.h"
#include "CondFormats/EcalObjects/interface/EcalIntercalibConstantsMC.h"

#include "PhaseIAnalyzer.h"

#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "SimDataFormats/CaloHit/interface/PCaloHit.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "SimDataFormats/Track/interface/SimTrack.h"
#include "CondFormats/EcalObjects/interface/EcalConstants.h"

using namespace std;
using namespace edm;

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
PhaseIAnalyzer::PhaseIAnalyzer(const edm::ParameterSet& iConfig)
{
//now do what ever initialization is needed
usesResource("TFileService");

 digiTagEB_= iConfig.getParameter<edm::InputTag>("BarrelDigis");

 digiTokenEB_ = consumes<EBDigiCollection>(digiTagEB_);
 
 //Files:
 edm::Service<TFileService> fs;
 //Histograms

 for(int isample=0;isample<16;isample++) {
   EBEnergyHisto[isample] = fs->make<TH1I>(Form("EnergyEB_%d", isample), Form("Energy sample %d  Barrel;ADC",isample), 950 , 100 , 2000 );
   EBGainHisto[isample]   = fs->make<TH1I>(Form("GainEB_%d", isample), Form("Gain Barrel sample %d;Gain",isample) , 5 , 0 , 4 );
 }
  SingleChannelE = fs->make<TH1I>("SingleChannelE", "Energy single channel Barrel;ADC", 17 , -0.5 , 16.5 );
  SingleChannelELow = fs->make<TH1I>("SingleChannelELow", "Energy single channel Barrel;ADC", 17 , -0.5 , 16.5 );

  Char_t histo[200];
 
  sprintf (histo, "EcalDigiTaskBarrelOccupancy" ) ;
  meEBDigiOccupancy_ = fs->make<TH2D>(histo, histo, 360, 0., 360., 170, -85., 85.);

  sprintf (histo, "EcalDigiTaskBarrelOccupancyHigh" ) ;
  meEBDigiOccupancyHigh_ = fs->make<TH2D>(histo, histo, 360, 0., 360., 170, -85., 85.);

  sprintf (histo, "EcalDigiTaskBarrelOccupancyLow" ) ;
  meEBDigiOccupancyLow_ = fs->make<TH2D>(histo, histo, 360, 0., 360., 170, -85., 85.);
  
  sprintf (histo, "EcalDigiTaskBarrelOccupancyMid" ) ;
  meEBDigiOccupancyMid_ = fs->make<TH2D>(histo, histo, 360, 0., 360., 170, -85., 85.);

  sprintf (histo, "EcalDigiTaskBarrelDigisMultiplicity" ) ;
  meEBDigiMultiplicity_ = fs->make<TH1D>(histo, histo, 612, 0., 61200);
  
    
  for (int i = 0; i < ecalPh2::sampleSize ; i++ ) {

    sprintf (histo, "EcalDigiTaskBarrelAnalogPulse%02d", i+1) ;
    meEBDigiADCAnalog_[i] = fs->make<TH1D>(histo, histo, 4000, 0., 400.);

    sprintf (histo, "EcalDigiTaskBarrelADCPulse%02dGain0_Saturated", i+1) ;
    meEBDigiADCgS_[i] = fs->make<TH1D>(histo, histo, 4096, -0.5, 4095.5);
    sprintf (histo, "EcalDigiTaskBarrelGainPulse%02d", i+1) ;
    meEBDigiGain_[i] = fs->make<TH1D>(histo, histo, 4, 0, 4);

  }
    
  sprintf (histo, "EcalDigiTaskBarrelPedestalForPreSample" ) ;
  meEBPedestal_ = fs->make<TH1D>(histo, histo, 4096, -0.5, 4095.5) ;

  sprintf (histo, "EcalDigiTaskBarrelMaximumPositionGt100ADC" ) ;
  meEBMaximumgt100ADC_ = fs->make<TH1D>(histo, histo, 10, 0., 10.) ;

  sprintf (histo, "EcalDigiTaskBarrelMaximumPositionGt10ADC" ) ;
  meEBMaximumgt10ADC_ = fs->make<TH1D>(histo, histo, 10, 0., 10.) ;

  //   sprintf (histo, "EcalDigiTask Barrel ADC counts after gain switch" ) ;
  //meEBnADCafterSwitch_ = fs->make<TH1D>histo, histo, 10, 0., 10.) ;


 
}
  

PhaseIAnalyzer::~PhaseIAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
PhaseIAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  
  //LogInfo("PhaseI") << "new event ";

   Handle<EBDigiCollection> pDigiEB;
   iEvent.getByToken(digiTokenEB_,pDigiEB);

   // edm::ESHandle<EcalIntercalibConstantsMC> ical;
   // iSetup.get<EcalIntercalibConstantsMCRcd>().get(ical);
   // const EcalIntercalibConstantsMC* ical_map = ical.product();
   // EcalIntercalibConstantsMC::const_iterator itical = ical_map->getMap().find(2);
   // cout << "intercalib: " << (*itical) << endl;

   const int MAXSAMPLES=ecalPh1::sampleSize; 

   std::vector<double> ebAnalogSignal ;
   std::vector<double> ebADCCounts ;
   std::vector<double> ebADCGains_temp ;
   std::vector<double> ebADCGains ;


   ebAnalogSignal.reserve(EBDataFrame::MAXSAMPLES);
   ebADCCounts.reserve(EBDataFrame::MAXSAMPLES);
   ebADCGains_temp.reserve(EBDataFrame::MAXSAMPLES);
   ebADCGains.reserve(EBDataFrame::MAXSAMPLES);

   //Take Pedestals:
  //  edm::ESHandle<EcalLiteDTUPedestals> peds;
  //  iSetup.get<EcalLiteDTUPedestalsRcd>().get(peds);
  //  const EcalLiteDTUPedestals* myped = peds.product();
  //   int cnt=0;
  // for( EcalLiteDTUPedestals::const_iterator it = myped->barrelItems().begin(); it != myped->barrelItems().end(); ++it)
  //   {
  //     std::cout << "EcalPedestal: " << " BARREL " << cnt << " "
  //               << "  mean:  " <<(*it).mean(0) << " rms: " << (*it).rms(0);                
  //     std::cout << std::endl;
  //           ++cnt;
  //   }
  
   int nDigis=0;
   int LowGain(0),MidGain(0),HighGain(0);
   int maxADCValue=0;

   for (EBDigiCollection::const_iterator pDigi=pDigiEB->begin(); pDigi!=pDigiEB->end(); ++pDigi) 
    {
      LowGain=0;
      MidGain=0;
      HighGain=0;
      maxADCValue=0;
      
      EBDataFrame digi( *pDigi );
       int nrSamples = digi.size();
       cout<<"NSamples found: "<<nrSamples<<endl;
       EBDetId ebid = digi.id () ;
       cout<<" Crystall ID "<<ebid<<endl;
       nDigis++;
       cout<<" nDigis aaaaaaa "<<nDigis<<endl;
       double Emax = 0. ;
       int Pmax = 0 ;
       double pedestalPreSample = 0.;
       double pedestalPreSampleAnalog = 0.;
       int countsAfterGainSwitch = -1;
     
       for (int sample = 0 ; sample < nrSamples; ++sample) 
        {
          ebAnalogSignal[sample] = 0.;
          ebADCCounts[sample] = 0.;
          ebADCGains_temp[sample] = 0.;
          ebADCGains[sample] = 0.;
        }

       double  gainConv_[4]={0,12,6,1};
      // saturated channels
      double barrelADCtoGeV_ = 0.048; //GeV

      // EcalIntercalibConstantsMC* ical = 

      // EcalIntercalibConstantMCMap &icalMap = 


      
      for (int sample = 0 ; sample < nrSamples; ++sample) 
       {

	 
	  int thisSample = digi[sample];
	  
          ebADCCounts[sample] = (thisSample&0xFFF);
          ebADCGains[sample]  = (thisSample&(0x3<<12))>>12;
          ebAnalogSignal[sample] = (ebADCCounts[sample]*gainConv_[(int)ebADCGains[sample]]*barrelADCtoGeV_);

	  if(ebADCCounts[sample] > maxADCValue) maxADCValue = ebADCCounts[sample];

	  if( ebid.iphi()==333 and ebid.ieta()==83) SingleChannelE -> SetBinContent(sample,ebADCCounts[sample]);

	  if( ebid.iphi()==334 and ebid.ieta()==83) {
	    SingleChannelELow -> SetBinContent(sample,ebADCCounts[sample]);
	  }
	  
	  if( ebADCGains[sample]==3) LowGain=1;
	  if( ebADCGains[sample]==2) MidGain=1;
	  if( ebADCGains[sample]==1) HighGain=1;
	  
	   
	  if(ebADCCounts[sample]>250) {
	    cout<<"Channel: "<<ebid<<endl; 
          cout<<"Sample "<<sample<<endl;
          cout<<"		Full data "<<thisSample<<endl;
          cout<<"		ebADCCounts "<<ebADCCounts[sample]<<endl;
          cout<<"		ebADCGains "<<ebADCGains[sample]<<endl;
          cout<<"		gainConv_ "<<gainConv_[(int)ebADCGains[sample]]<<endl;
          cout<<"		barrelADCtoGeV_ "<<barrelADCtoGeV_<<endl;
          cout<<"		ebAnalogSignal "<<ebAnalogSignal[sample]<<endl;
	  }
      
          if(Emax < ebAnalogSignal[sample] ) 
           {
              Emax = ebAnalogSignal[sample] ;
              Pmax = sample ;          
           }

      }
      
    if(1==1)cout<<"P max "<<Pmax<<endl;
    cout<<"High Gain: "<<HighGain<<" MidGain "<<MidGain<<" LowGain "<<LowGain<<endl;
      cout<<"maxADCValue "<<maxADCValue<<endl;
      if(maxADCValue > 250) {
	meEBDigiOccupancy_->SetBinContent( ebid.iphi(), ebid.ieta(),maxADCValue );
	if(HighGain==1) meEBDigiOccupancyHigh_->SetBinContent( ebid.iphi(), ebid.ieta(),maxADCValue );
	if(MidGain==1) 	meEBDigiOccupancyMid_->SetBinContent( ebid.iphi(), ebid.ieta(),maxADCValue );
	if(LowGain==1) 	meEBDigiOccupancyLow_->SetBinContent( ebid.iphi(), ebid.ieta(),maxADCValue );

      }
      
  }//end digi
   
   
}


// ------------ method called once each job just before starting event loop  ------------
void 
PhaseIAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
PhaseIAnalyzer::endJob() 
{
}


// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
PhaseIAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(PhaseIAnalyzer);
