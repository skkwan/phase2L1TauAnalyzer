// -*- C++ -*-
//
// Package:    L1Trigger/phase2L1TauAnalyzer
// Class:      phase2L1TauAnalyzer
// 
/**\class phase2L1TauAnalyzer phase2L1TauAnalyzer.cc L1Trigger/phase2L1TauAnalyzer/plugins/phase2L1TauAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Isobel Ojalvo
//         Created:  Fri, 26 May 2017 15:44:30 GMT
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
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"

#include "L1Trigger/phase2Demonstrator/interface/triggerGeometryTools.hh"
#include "DataFormats/L1CaloTrigger/interface/L1CaloCollections.h"

#include "L1Trigger/phase2L1TauAnalyzer/plugins/helpers.h"
#include "DataFormats/Math/interface/deltaR.h"

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class phase2L1TauAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit phase2L1TauAnalyzer(const edm::ParameterSet&);
      ~phase2L1TauAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      // ----------member data ---------------------------
  typedef std::vector<reco::GenParticle> GenParticleCollectionType;

  struct genVisTau{
    reco::Candidate::LorentzVector p4;
    int decayMode;
  };

  edm::EDGetTokenT< L1CaloClusterCollection > L1ClustersToken_;
  edm::EDGetTokenT<std::vector<reco::GenParticle> > genToken_;
  edm::InputTag genSrc_;

  double genPt, genEta, genPhi;
  int decayMode, run, lumi, event;


};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
phase2L1TauAnalyzer::phase2L1TauAnalyzer(const edm::ParameterSet& cfg):
  L1ClustersToken_( consumes< L1CaloClusterCollection >(cfg.getParameter<edm::InputTag>("L1Clusters"))),
  genSrc_ ((cfg.getParameter<edm::InputTag>( "genParticles")))
{
   //now do what ever initialization is needed
   usesResource("TFileService");
   genToken_ =     consumes<std::vector<reco::GenParticle> >(genSrc_);

}


phase2L1TauAnalyzer::~phase2L1TauAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)


}


//
// member functions
//

// ------------ method called for each event  ------------
void
phase2L1TauAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

  edm::Handle< std::vector<L1CaloCluster> > l1Clusters;
  iEvent.getByToken( L1ClustersToken_, l1Clusters);

   // Get genParticles
   edm::Handle<GenParticleCollectionType> genParticleHandle;
   if(!iEvent.getByToken(genToken_,genParticleHandle))
     std::cout<<"No gen Particles Found "<<std::endl;
   else
     std::cout<<"Gen Particles size "<<genParticleHandle->size()<<std::endl;

   std::vector<reco::GenParticle> genTaus;
   std::vector<reco::GenParticle> genPiZeros;
   std::vector<reco::GenParticle> genParticles;

   //reco::GenParticle* genTau;
   for(unsigned int i = 0; i< genParticleHandle->size(); i++){
     edm::Ptr<reco::GenParticle> ptr(genParticleHandle, i);
     genParticles.push_back(*ptr);

     if(abs(ptr->pdgId())==111 && abs(ptr->eta()<1.74)){
       genPiZeros.push_back(*ptr);
       std::cout<<"Found PiZero PDGID 111 pt: "<<ptr->pt()<<" eta: "<<ptr->eta()<<" phi: "<<ptr->phi()<<std::endl;
     }
     if(abs(ptr->pdgId())==15){
       genTaus.push_back(*ptr);
     }
   }


   for(auto genPiZero: genPiZeros){

     reco::Candidate::LorentzVector visGenPiZero = getVisMomentum(&genPiZero, &genParticles);
     std::cout<<"pt: "<<visGenPiZero.pt()<< " eta: "<<visGenPiZero.eta()<<" phi: "<<visGenPiZero.phi()<<std::endl;
     std::cout<<"------printing matched clusters-------"<<std::endl;
     for(int i = 0; i < l1Clusters->size(); i++){
       if( reco::deltaR(l1Clusters->at(i).p4().Eta(), l1Clusters->at(i).p4().Phi(),genPiZero.p4().eta(),genPiZero.p4().phi()) < 0.15){
	 std::cout<<l1Clusters->at(i)<<std::endl;
       }
     }
   }


   std::vector<genVisTau> GenOneProngTaus;
   std::vector<genVisTau> GenOneProngPi0Taus;
   std::vector<genVisTau> GenThreeProngTaus;
   //Find and Sort the 1 Prong, 1 Prong + pi0 and 3 Prong Taus

   for(auto genTau: genTaus){
     reco::Candidate::LorentzVector visGenTau= getVisMomentum(&genTau, &genParticles);
     genVisTau Temp;
     genPt = visGenTau.pt();
     genEta = visGenTau.eta();
     genPhi = visGenTau.phi();
     decayMode = GetDecayMode(&genTau);
     Temp.p4 = visGenTau;
     Temp.decayMode = decayMode;
     std::cout<<"Tau Decay Mode "<<decayMode<<std::endl;
     std::cout<<"tau vis pt: "<<genPt<<" genEta: "<<genEta<<" genPhi: "<<genPhi<<std::endl;

     if(decayMode >21 ){
       std::cout<<"found 3 prong tau: "<<decayMode<<std::endl;
       GenThreeProngTaus.push_back(Temp);
     }

     if(decayMode == 10 ){
       GenOneProngTaus.push_back(Temp);
     }

     if(decayMode > 10 && decayMode < 20 ){
       GenOneProngPi0Taus.push_back(Temp);
     }
   }

  
}


// ------------ method called once each job just before starting event loop  ------------
void 
phase2L1TauAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
phase2L1TauAnalyzer::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
phase2L1TauAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(phase2L1TauAnalyzer);
