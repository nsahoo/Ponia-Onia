// -*- C++ -*-
//
// Package:    Onia2MuMuRootupler
// Class:      Onia2MuMuRootupler
// 
// Description: Dump  Onia(mu+ mu-)  decays
//
// Author:  Alberto Sanchez Hernandez
//

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TLorentzVector.h"
#include "TTree.h"

//
// class declaration
//

class Onia2MuMuRootupler:public edm::EDAnalyzer {
      public:
	explicit Onia2MuMuRootupler(const edm::ParameterSet &);
	~Onia2MuMuRootupler();

	static void fillDescriptions(edm::ConfigurationDescriptions & descriptions);

      private:
        std::vector<UInt_t> getTriggerBits(const edm::Event &);
        bool   isAncestor(const reco::Candidate *, const reco::Candidate *);
        const  reco::Candidate* GetAncestor(const reco::Candidate *);

	virtual void beginJob();
	virtual void analyze(const edm::Event &, const edm::EventSetup &);
	virtual void endJob();

	virtual void beginRun(edm::Run const &, edm::EventSetup const &);
	virtual void endRun(edm::Run const &, edm::EventSetup const &);
	virtual void beginLuminosityBlock(edm::LuminosityBlock const &, edm::EventSetup const &);
	virtual void endLuminosityBlock(edm::LuminosityBlock const &, edm::EventSetup const &);

	// ----------member data ---------------------------
	std::string file_name;
	edm::EDGetTokenT<pat::CompositeCandidateCollection> dimuon_Label;
        edm::EDGetTokenT<reco::VertexCollection> primaryVertices_Label;
        edm::EDGetTokenT<edm::TriggerResults> triggerResults_Label;
        int  pdgid_;
        std::vector<double> OniaMassCuts_;
	bool isMC_;
        bool OnlyBest_;
        bool OnlyGen_;

	UInt_t run;
	UInt_t event;
        Int_t  irank;
        UInt_t trigger;
        std::vector<UInt_t> trigvec;
  UInt_t passbit0, passbit1, passbit2, passbit3, passbit4, passbit5, passbit6, passbit7, passbit8, passbit9, passbit10, passbit11, passbit12, passbit13, passbit14, passbit15, passbit16, passbit17, passbit18;
        Int_t  charge; 

	TLorentzVector dimuon_p4;
	TLorentzVector muonP_p4;
	TLorentzVector muonN_p4;

        Float_t MassErr;
        Float_t vProb;
        Float_t DCA;
        Float_t ppdlPV;
        Float_t ppdlErrPV;
        Float_t ppdlBS;
        Float_t ppdlErrBS;
        Float_t cosAlpha;
        Float_t lxyPV;
        Float_t lxyBS;

	Int_t numPrimaryVertices;

	TTree *onia_tree;

        Int_t mother_pdgId;
        Int_t dimuon_pdgId;
	TLorentzVector gen_dimuon_p4;
	TLorentzVector gen_muonP_p4;
        TLorentzVector gen_muonM_p4;

          
        edm::EDGetTokenT<reco::GenParticleCollection> genCands_;
        edm::EDGetTokenT<pat::PackedGenParticleCollection> packCands_;

};

//
// constructors and destructor
//

Onia2MuMuRootupler::Onia2MuMuRootupler(const edm::ParameterSet & iConfig):
dimuon_Label(consumes<pat::CompositeCandidateCollection>(iConfig.getParameter< edm::InputTag>("dimuons"))),
primaryVertices_Label(consumes<reco::VertexCollection>(iConfig.getParameter< edm::InputTag>("primaryVertices"))),
triggerResults_Label(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("TriggerResults"))),
pdgid_(iConfig.getParameter<uint32_t>("onia_pdgid")),
OniaMassCuts_(iConfig.getParameter<std::vector<double>>("onia_mass_cuts")),
isMC_(iConfig.getParameter<bool>("isMC")),
OnlyBest_(iConfig.getParameter<bool>("OnlyBest")),
OnlyGen_(iConfig.getParameter<bool>("OnlyGen"))
{
  edm::Service < TFileService > fs;
  onia_tree = fs->make < TTree > ("oniaTree", "Tree of Onia2MuMu");

  if (!OnlyGen_) {
    onia_tree->Branch("run",       &run,       "run/i");
    onia_tree->Branch("event",     &event,     "event/i");
    onia_tree->Branch("irank",     &irank,     "irank/I");
    onia_tree->Branch("trigger",   &trigger,   "trigger/i");
    onia_tree->Branch("passbit0",  &passbit0,  "passbit0/i");
    onia_tree->Branch("passbit1",  &passbit1,  "passbit1/i");
    onia_tree->Branch("passbit2",  &passbit2,  "passbit2/i");
    onia_tree->Branch("passbit3",  &passbit3,  "passbit3/i");
    onia_tree->Branch("passbit4",  &passbit4,  "passbit4/i");
    onia_tree->Branch("passbit5",  &passbit5,  "passbit5/i");
    onia_tree->Branch("passbit6",  &passbit6,  "passbit6/i");
    onia_tree->Branch("passbit7",  &passbit7,  "passbit7/i");
   
    onia_tree->Branch("passbit8",  &passbit8,  "passbit8/i");
    onia_tree->Branch("passbit9",  &passbit9,  "passbit9/i");

    onia_tree->Branch("passbit10", &passbit10, "passbit10/i");
    onia_tree->Branch("passbit11", &passbit11, "passbit11/i");
    onia_tree->Branch("passbit12", &passbit12, "passbit12/i");
    onia_tree->Branch("passbit13", &passbit13, "passbit13/i");
    onia_tree->Branch("passbit14", &passbit14, "passbit14/i");
    onia_tree->Branch("passbit15", &passbit15, "passbit15/i");
    onia_tree->Branch("passbit16", &passbit16, "passbit16/i");
    onia_tree->Branch("passbit17", &passbit17, "passbit17/i");
    onia_tree->Branch("passbit18", &passbit18, "passbit18/i");

    onia_tree->Branch("charge",    &charge,    "charge/I");

    onia_tree->Branch("dimuon_p4", "TLorentzVector", &dimuon_p4);
    onia_tree->Branch("muonP_p4",  "TLorentzVector", &muonP_p4);
    onia_tree->Branch("muonN_p4",  "TLorentzVector", &muonN_p4);

    onia_tree->Branch("MassErr",   &MassErr,    "MassErr/F");
    onia_tree->Branch("vProb",     &vProb,      "vProb/F");
    onia_tree->Branch("DCA",       &DCA,        "DCA/F");
    onia_tree->Branch("ppdlPV",    &ppdlPV,     "ppdlPV/F");
    onia_tree->Branch("ppdlErrPV", &ppdlErrPV,  "ppdlErrPV/F");
    onia_tree->Branch("ppdlBS",    &ppdlBS,     "ppdlBS/F");
    onia_tree->Branch("ppdlErrBS", &ppdlErrBS,  "ppdlErrBS/F");
    onia_tree->Branch("cosAlpha",  &cosAlpha,   "cosAlpha/F");
    onia_tree->Branch("lxyPV",     &lxyPV,      "lxyPV/F");
    onia_tree->Branch("lxyBS",     &lxyBS,      "lxyBS/F");

    onia_tree->Branch("numPrimaryVertices", &numPrimaryVertices, "numPrimaryVertices/I");
  }

  if (isMC_ || OnlyGen_) {
     onia_tree->Branch("mother_pdgId",  &mother_pdgId,     "mother_pdgId/I");
     onia_tree->Branch("dimuon_pdgId",  &dimuon_pdgId,     "dimuon_pdgId/I");
     onia_tree->Branch("gen_dimuon_p4", "TLorentzVector",  &gen_dimuon_p4);
     onia_tree->Branch("gen_muonP_p4",  "TLorentzVector",  &gen_muonP_p4);
     onia_tree->Branch("gen_muonN_p4",  "TLorentzVector",  &gen_muonM_p4);
  }
  genCands_ = consumes<reco::GenParticleCollection>((edm::InputTag)"prunedGenParticles");
  packCands_ = consumes<pat::PackedGenParticleCollection>((edm::InputTag)"packedGenParticles");
}

Onia2MuMuRootupler::~Onia2MuMuRootupler() {}

//
// member functions
//

const reco::Candidate* Onia2MuMuRootupler::GetAncestor(const reco::Candidate* p) {
   if (p->numberOfMothers()) {
      if  ((p->mother(0))->pdgId() == p->pdgId()) return GetAncestor(p->mother(0));
      else return p->mother(0);
   }
   //std::cout << "GetAncestor: Inconsistet ancestor, particle does not have a mother " << std::endl;
   return p;
}

//Check recursively if any ancestor of particle is the given one
bool Onia2MuMuRootupler::isAncestor(const reco::Candidate* ancestor, const reco::Candidate * particle) {
   if (ancestor == particle ) return true;
   for (size_t i=0; i< particle->numberOfMothers(); i++) {
      if (isAncestor(ancestor, particle->mother(i))) return true;
   }
   return false;
}

/* Grab Trigger information. Save it in variable trigger, trigger is an int between 0 and 127, in binary it is:
   (pass 2)(pass 1)(pass 0)
   ex. 7 = pass 0, 1 and 2
   ex. 6 = pass 2, 3
   ex. 1 = pass 0
*/

std::vector <UInt_t> Onia2MuMuRootupler::getTriggerBits(const edm::Event& iEvent ) {
   UInt_t itrigger = 0;
   std::vector <UInt_t> trgbitsN;
   edm::Handle<edm::TriggerResults> triggerResults_handle;
   iEvent.getByToken(triggerResults_Label, triggerResults_handle);
   if ( triggerResults_handle.isValid() ) { 
      const edm::TriggerNames & TheTriggerNames = iEvent.triggerNames(*triggerResults_handle);
      std::vector <unsigned int> bits_0, bits_1, bits_2, bits_3, bits_4, bits_5, bits_6, bits_7, bits_8, bits_9, bits_10, bits_11, bits_12, bits_13, bits_14, bits_15, bits_16, bits_17, bits_18;
      //std::vector <unsigned int>  bits_0, bits_1, bits_2, bits_3, bits_4, bits_5, bits_6, bits_7, bits_8, bits_9 ;
     
      for ( int version = 1; version<3; version ++ ) {
	std::stringstream ss0,ss1,ss2,ss3,ss4,ss5,ss6,ss7,ss8,ss9,ss10,ss11,ss12,ss13,ss14,ss15,ss16,ss17,ss18;
	//std::stringstream ss0,ss1,ss2,ss3,ss4,ss5,ss6,ss7,ss8,ss9 ;
	//quarkonia paths


	ss0<<"HLT_Dimuon16_Jpsi_v"<<version;
	bits_0.push_back(TheTriggerNames.triggerIndex( edm::InputTag(ss0.str()).label().c_str()));       
	ss1<<"HLT_Dimuon20_Jpsi_v"<<version;
        bits_1.push_back(TheTriggerNames.triggerIndex( edm::InputTag(ss1.str()).label().c_str()));

	ss2<<"HLT_Dimuon13_PsiPrime_v"<<version;
	bits_2.push_back(TheTriggerNames.triggerIndex( edm::InputTag(ss2.str()).label().c_str()));
	ss3<<"HLT_Dimuon13_Upsilon_v"<<version;
	bits_3.push_back(TheTriggerNames.triggerIndex( edm::InputTag(ss3.str()).label().c_str()));

	ss4<<"HLT_Dimuon10_Jpsi_Barrel_v"<<version;
	bits_4.push_back(TheTriggerNames.triggerIndex( edm::InputTag(ss4.str()).label().c_str()));
	ss5<<"HLT_Dimuon8_PsiPrime_Barrel_v"<<version;
	bits_5.push_back(TheTriggerNames.triggerIndex( edm::InputTag(ss5.str()).label().c_str()));
	ss6<<"HLT_Dimuon8_Upsilon_Barrel_v"<<version;
	bits_6.push_back(TheTriggerNames.triggerIndex( edm::InputTag(ss6.str()).label().c_str()));
	ss7<<"HLT_Dimuon0_Phi_Barrel_v"<<version;
	bits_7.push_back(TheTriggerNames.triggerIndex( edm::InputTag(ss7.str()).label().c_str()));

	ss8<<"HLT_Mu25_TkMu0_dEta18_Onia_v"<<version;
	bits_8.push_back(TheTriggerNames.triggerIndex( edm::InputTag(ss8.str()).label().c_str()));	
	ss9<<"HLT_Mu16_TkMu0_dEta18_Onia_v"<<version;
	bits_9.push_back(TheTriggerNames.triggerIndex( edm::InputTag(ss9.str()).label().c_str()));

	 //TnP paths
	 ss10<<"HLT_Mu7p5_L2Mu2_Jpsi_v"<<version;
	 bits_10.push_back(TheTriggerNames.triggerIndex( edm::InputTag(ss10.str()).label().c_str()));
	 ss11<<"HLT_Mu7p5_L2Mu2_Upsilon_v"<<version;
	 bits_11.push_back(TheTriggerNames.triggerIndex( edm::InputTag(ss11.str()).label().c_str()));

	 ss12<<"HLT_Mu7p5_Track2_Jpsi_v"<<version;
	 bits_12.push_back(TheTriggerNames.triggerIndex( edm::InputTag(ss12.str()).label().c_str()));
	 ss13<<"HLT_Mu7p5_Track3p5_Jpsi_v"<<version;
	 bits_13.push_back(TheTriggerNames.triggerIndex( edm::InputTag(ss13.str()).label().c_str()));
	 ss14<<"HLT_Mu7p5_Track7_Jpsi_v"<<version;
	 bits_14.push_back(TheTriggerNames.triggerIndex( edm::InputTag(ss14.str()).label().c_str()));
	 ss15<<"HLT_Mu7p5_Track2_Upsilon_v"<<version;
	 bits_15.push_back(TheTriggerNames.triggerIndex( edm::InputTag(ss15.str()).label().c_str()));
	 ss16<<"HLT_Mu7p5_Track3p5_Upsilon_v"<<version;
	 bits_16.push_back(TheTriggerNames.triggerIndex( edm::InputTag(ss16.str()).label().c_str()));
	 ss17<<"HLT_Mu7p5_Track7_Upsilon_v"<<version;
	 bits_17.push_back(TheTriggerNames.triggerIndex( edm::InputTag(ss17.str()).label().c_str()));

	 ss18<<"HLT_Mu16_TkMu0_dEta18_Phi_v2"<<version;
	 bits_18.push_back(TheTriggerNames.triggerIndex( edm::InputTag(ss17.str()).label().c_str()));

      }

      //std::cout << "trg handle result: " << triggerResults_handle->size() << std::endl;
      
      for (unsigned int i=0; i<bits_0.size(); i++) {
         unsigned int bit = bits_0[i];
         if ( bit < triggerResults_handle->size() ){
	   if ( triggerResults_handle->accept( bit ) && !triggerResults_handle->error( bit ) ) {
             itrigger += 1;
	     trgbitsN.push_back(1);
             break;
	   } else {
	     trgbitsN.push_back(0);
	   //printf("INFO: Bit 0 size is more than handle size.\n");
	   }
	 }
      }

      for (unsigned int i=0; i<bits_1.size(); i++) {
         unsigned int bit = bits_1[i];
         if ( bit < triggerResults_handle->size() ){
           if ( triggerResults_handle->accept( bit ) && !triggerResults_handle->error( bit ) ) {
             itrigger += 2;
	     trgbitsN.push_back(1);
             break;
           }  else{
	     trgbitsN.push_back(0);
	   //printf("INFO: Bit 1 size is more than handle size.\n");
	   }
	 }
      }

      for (unsigned int i=0; i<bits_2.size(); i++) {
         unsigned int bit = bits_2[i];
         if ( bit < triggerResults_handle->size() ){
           if ( triggerResults_handle->accept( bit ) && !triggerResults_handle->error( bit ) ) {
             itrigger += 4;
	     trgbitsN.push_back(1);
             break;
           }  else{
	     trgbitsN.push_back(0);
	   //printf("INFO: Bit 2 size is more than handle size.\n");
	   }
	 }
      }

      for (unsigned int i=0; i<bits_3.size(); i++) {
         unsigned int bit = bits_3[i];
         if ( bit < triggerResults_handle->size() ){
           if ( triggerResults_handle->accept( bit ) && !triggerResults_handle->error( bit ) ) {
             itrigger += 8;
	     trgbitsN.push_back(1);
             break;
           }  else{
	   trgbitsN.push_back(0);
	   //printf("INFO: Bit 3 size is more than handle size.\n");
	   }
	 }
      }

      for (unsigned int i=0; i<bits_4.size(); i++) {
         unsigned int bit = bits_4[i];
         if ( bit < triggerResults_handle->size() ){
           if ( triggerResults_handle->accept( bit ) && !triggerResults_handle->error( bit ) ) {
             itrigger += 16;
	     trgbitsN.push_back(1);
             break;
           } else{
	   trgbitsN.push_back(0);
	   //printf("INFO: Bit 4 size is more than handle size.\n");
	   }
	 }
      }

      for (unsigned int i=0; i<bits_5.size(); i++) {
         unsigned int bit = bits_5[i];
         if ( bit < triggerResults_handle->size() ){
           if ( triggerResults_handle->accept( bit ) && !triggerResults_handle->error( bit ) ) {
             itrigger += 32;
	     trgbitsN.push_back(1);
             break;
           } else{
	   trgbitsN.push_back(0);
	   //printf("INFO: Bit 5 size is more than handle size.\n");
	   }
	 }
      }

      for (unsigned int i=0; i<bits_6.size(); i++) {
         unsigned int bit = bits_6[i];
         if ( bit < triggerResults_handle->size() ){
           if ( triggerResults_handle->accept( bit ) && !triggerResults_handle->error( bit ) ) {
             itrigger += 64;
	     trgbitsN.push_back(1);
             break;
           } else{
	   trgbitsN.push_back(0);
	   //printf("INFO: Bit 6 size is more than handle size.\n");
	   }
	 }
      }

      for (unsigned int i=0; i<bits_7.size(); i++) {
         unsigned int bit = bits_7[i];
         if ( bit < triggerResults_handle->size() ){
           if ( triggerResults_handle->accept( bit ) && !triggerResults_handle->error( bit ) ) {
             itrigger += 128;
	     trgbitsN.push_back(1);
             break;
           } else{
	   trgbitsN.push_back(0);
	   //printf("INFO: Bit 7 size is more than handle size.\n");
	   }
	 }
      }

      //
      
      for (unsigned int i=0; i<bits_8.size(); i++) {
	unsigned int bit = bits_8[i];
	if ( bit < triggerResults_handle->size() ){
	  if ( triggerResults_handle->accept( bit ) && !triggerResults_handle->error( bit ) ) {
	    itrigger += 256;
	    trgbitsN.push_back(1);
	    break;
	  } else{
	  trgbitsN.push_back(0);
	  //printf("INFO: Bit 8 size is more than handle size.\n");
	  }
	}
      }
      

      for (unsigned int i=0; i<bits_9.size(); i++) {
	unsigned int bit = bits_9[i];
	if ( bit < triggerResults_handle->size() ){
	  if ( triggerResults_handle->accept( bit ) && !triggerResults_handle->error( bit ) ) {
	    itrigger += 512;
	    trgbitsN.push_back(1);
	    break;
	  } else{
	  trgbitsN.push_back(0);
	  //printf("INFO: Bit 9 size is more than handle size.\n");
	  }
	}
      }

      //
      for (unsigned int i=0; i<bits_10.size(); i++) {
	unsigned int bit = bits_10[i];
	if ( bit < triggerResults_handle->size() ){
	  if ( triggerResults_handle->accept( bit ) && !triggerResults_handle->error( bit ) ) {
	    itrigger += pow(2,10);
	    trgbitsN.push_back(1);
	    break;
	  } else{
	  trgbitsN.push_back(0);
	  //printf("INFO: Bit 10 size is more than handle size.\n");
	  }
	}
      }


      for (unsigned int i=0; i<bits_11.size(); i++) {
	unsigned int bit = bits_11[i];
	if ( bit < triggerResults_handle->size() ){
	  if ( triggerResults_handle->accept( bit ) && !triggerResults_handle->error( bit ) ) {
	    itrigger += pow(2,11);
	    trgbitsN.push_back(1);
	    break;
	  } else{
	  trgbitsN.push_back(0);
	  //printf("INFO: Bit 11 size is more than handle size.\n");
	  }
	}
      }

      for (unsigned int i=0; i<bits_12.size(); i++) {
	unsigned int bit = bits_12[i];
	if ( bit < triggerResults_handle->size() ){
	  if ( triggerResults_handle->accept( bit ) && !triggerResults_handle->error( bit ) ) {
	    itrigger += pow(2,12);
	    trgbitsN.push_back(1);
	    break;
	  } else{
	  trgbitsN.push_back(0);
	  //printf("INFO: Bit 12 size is more than handle size.\n");
	  }
	}
      }

      for (unsigned int i=0; i<bits_13.size(); i++) {
	unsigned int bit = bits_13[i];
	if ( bit < triggerResults_handle->size() ){
	  if ( triggerResults_handle->accept( bit ) && !triggerResults_handle->error( bit ) ) {
	    itrigger += pow(2,13);
	    trgbitsN.push_back(1);
	    break;
	  } else{
	  trgbitsN.push_back(0);
	  //printf("INFO: Bit 13 size is more than handle size.\n");
	  }
	}
      }

      for (unsigned int i=0; i<bits_14.size(); i++) {
	unsigned int bit = bits_14[i];
	if ( bit < triggerResults_handle->size() ){
	  if ( triggerResults_handle->accept( bit ) && !triggerResults_handle->error( bit ) ) {
	    itrigger += pow(2,14);
	    trgbitsN.push_back(1);
	    break;
	  } else{
	  trgbitsN.push_back(0);
	  //printf("INFO: Bit 14 size is more than handle size.\n");
	  }
	}
      }

      for (unsigned int i=0; i<bits_15.size(); i++) {
	unsigned int bit = bits_15[i];
	if ( bit < triggerResults_handle->size() ){
	  if ( triggerResults_handle->accept( bit ) && !triggerResults_handle->error( bit ) ) {
	    itrigger += pow(2,15);
	    trgbitsN.push_back(1);
	    break;
	  } else{
	  trgbitsN.push_back(0);
	  //printf("INFO: Bit 15 size is more than handle size.\n");
	  }
	}
      }

      for (unsigned int i=0; i<bits_16.size(); i++) {
	unsigned int bit = bits_16[i];
	if ( bit < triggerResults_handle->size() ){
	  if ( triggerResults_handle->accept( bit ) && !triggerResults_handle->error( bit ) ) {
	    itrigger += pow(2,16);
	    trgbitsN.push_back(1);
	    break;
	  } else{
	  trgbitsN.push_back(0);
	  //printf("INFO: Bit 16 size is more than handle size.\n");
	  }
	}
      }

      for (unsigned int i=0; i<bits_17.size(); i++) {
        unsigned int bit = bits_17[i];
        if ( bit < triggerResults_handle->size() ){
          if ( triggerResults_handle->accept( bit ) && !triggerResults_handle->error( bit ) ) {
            itrigger += pow(2,17);
            trgbitsN.push_back(1);
            break;
          } else{
          trgbitsN.push_back(0);
          //printf("INFO: Bit 17 size is more than handle size.\n");
	  }
	}
      }

      for (unsigned int i=0; i<bits_18.size(); i++) {
        unsigned int bit = bits_18[i];
        if ( bit < triggerResults_handle->size() ){
          if ( triggerResults_handle->accept( bit ) && !triggerResults_handle->error( bit ) ) {
            itrigger += pow(2,18);
            trgbitsN.push_back(1);
            break;
          } else{
	    trgbitsN.push_back(0);
	    //printf("INFO: Bit 18 size is more than handle size.\n");                                                                                                                                     
          }
        }
      } 


   }

   trgbitsN.push_back(itrigger);
   //std::cout << "size: " << trgbitsN.size() << std::endl;
   return trgbitsN ;
}

// ------------ method called for each event  ------------
void Onia2MuMuRootupler::analyze(const edm::Event & iEvent, const edm::EventSetup & iSetup) {

  edm::Handle<pat::CompositeCandidateCollection> dimuons;
  iEvent.getByToken(dimuon_Label,dimuons);

  edm::Handle<reco::VertexCollection> primaryVertices_handle;
  iEvent.getByToken(primaryVertices_Label, primaryVertices_handle);
 
  if (!OnlyGen_) {
    numPrimaryVertices = -1;
    if (primaryVertices_handle.isValid()) numPrimaryVertices = (int) primaryVertices_handle->size();


    //trigger   = getTriggerBits(iEvent).at(18);
    /*
    trigger   = getTriggerBits(iEvent).at(18);
    passbit0  = getTriggerBits(iEvent).at(0);
    passbit1  = getTriggerBits(iEvent).at(1);
    passbit2  = getTriggerBits(iEvent).at(2);
    passbit3  = getTriggerBits(iEvent).at(3);
    passbit4  = getTriggerBits(iEvent).at(4);
    passbit5  = getTriggerBits(iEvent).at(5);
    passbit6  = getTriggerBits(iEvent).at(6);
    passbit7  = getTriggerBits(iEvent).at(7);
    passbit8   = getTriggerBits(iEvent).at(8);
    passbit9   = getTriggerBits(iEvent).at(9);
    passbit10  = getTriggerBits(iEvent).at(10);
    passbit11  = getTriggerBits(iEvent).at(11);
    passbit12  = getTriggerBits(iEvent).at(12);
    passbit13  = getTriggerBits(iEvent).at(13);
    passbit14  = getTriggerBits(iEvent).at(14);
    passbit15  = getTriggerBits(iEvent).at(15);
    passbit16  = getTriggerBits(iEvent).at(16);
    passbit17  = getTriggerBits(iEvent).at(17);
    */


    trigvec   = getTriggerBits(iEvent);

    passbit0  = trigvec.at(0);
    passbit1  = trigvec.at(1);
    passbit2  = trigvec.at(2);
    passbit3  = trigvec.at(3);
    passbit4  = trigvec.at(4);
    passbit5  = trigvec.at(5);
    passbit6  = trigvec.at(6);
    passbit7  = trigvec.at(7);
    passbit8  = trigvec.at(8);
    passbit9  = trigvec.at(9);
    passbit10 = trigvec.at(10);
    passbit11 = trigvec.at(11);
    passbit12 = trigvec.at(12);
    passbit13 = trigvec.at(13);
    passbit14 = trigvec.at(14);
    passbit15 = trigvec.at(15);
    passbit16 = trigvec.at(16);
    passbit17 = trigvec.at(17);
    passbit18 = trigvec.at(18);

    trigger   = trigvec.at(19);



    std::cout << "trigvec size: " << trigvec.size() << std::endl;
    if (trigvec.size() != 20) printf("ERROR: size is more.\n");

    run     = iEvent.id().run();
    event   = iEvent.id().event();

    std::cout << "event: " << event << std::endl;
    //printf("total triggered events: %d \n", trigger);
  }

  dimuon_pdgId = 0;
  mother_pdgId = 0;
  irank = 0;

  // Pruned particles are the one containing "important" stuff
  edm::Handle<reco::GenParticleCollection> pruned;
  iEvent.getByToken(genCands_, pruned);

  // Packed particles are all the status 1, so usable to remake jets
  // The navigation from status 1 to pruned is possible (the other direction should be made by hand)
  edm::Handle<pat::PackedGenParticleCollection> packed;
  iEvent.getByToken(packCands_,  packed);

  if ((isMC_ || OnlyGen_) && packed.isValid() && pruned.isValid()) {
     dimuon_pdgId  = 0;
     gen_dimuon_p4.SetPtEtaPhiM(0.,0.,0.,0.);
     int foundit   = 0;

     for (size_t i=0; i<pruned->size(); i++) {
         int p_id = abs((*pruned)[i].pdgId());
         const reco::Candidate *aonia = &(*pruned)[i];
         if (( p_id == pdgid_ ) && (aonia->status() == 2)) {
            dimuon_pdgId = p_id;
            foundit++;
            for (size_t j=0; j<packed->size(); j++) { //get the pointer to the first survied ancestor of a given packed GenParticle in the prunedCollection
               const reco::Candidate * motherInPrunedCollection = (*packed)[j].mother(0);
               const reco::Candidate * d = &(*packed)[j];
               if ( motherInPrunedCollection != nullptr && (d->pdgId() == 13 ) && isAncestor(aonia , motherInPrunedCollection) ){
                  gen_muonM_p4.SetPtEtaPhiM(d->pt(),d->eta(),d->phi(),d->mass());
                  foundit++;
               } 
               if ( motherInPrunedCollection != nullptr && (d->pdgId() == -13 ) && isAncestor(aonia , motherInPrunedCollection) ) {
                  gen_muonP_p4.SetPtEtaPhiM(d->pt(),d->eta(),d->phi(),d->mass());
                  foundit++;
               }
               if ( foundit == 3 ) break;               
            }
            if ( foundit == 3 ) {
               gen_dimuon_p4 = gen_muonM_p4 + gen_muonP_p4;   // this should take into account FSR
               mother_pdgId  = GetAncestor(aonia)->pdgId();
               break;
            } else {
               foundit = 0;
               dimuon_pdgId = 0;
               mother_pdgId = 0;
               gen_dimuon_p4.SetPtEtaPhiM(0.,0.,0.,0.);
            }            
         }  // if ( p_id
     } // for (size

     // sanity check
     //if ( ! dimuon_pdgId ) std::cout << "Onia2MuMuRootupler: does not found the given decay " << run << "," << event << std::endl;
  }  // end if isMC

  float OniaMassMax_ = OniaMassCuts_[1];
  float OniaMassMin_ = OniaMassCuts_[0];

  if ( ! OnlyGen_ && dimuons.isValid() && dimuons->size() > 0) {
     for (pat::CompositeCandidateCollection::const_iterator  dimuonCand = dimuons->begin(); dimuonCand!= dimuons->end(); ++dimuonCand) {
        if (dimuonCand->mass() > OniaMassMin_ && dimuonCand->mass() < OniaMassMax_ && dimuonCand->charge() == 0) {
           dimuon_p4.SetPtEtaPhiM(dimuonCand->pt(), dimuonCand->eta(), dimuonCand->phi(), dimuonCand->mass());
           reco::Candidate::LorentzVector vP = dimuonCand->daughter("muon1")->p4();
           reco::Candidate::LorentzVector vM = dimuonCand->daughter("muon2")->p4();
           if ( dimuonCand->daughter("muon1")->charge() < 0) {
              vP = dimuonCand->daughter("muon2")->p4();
              vM = dimuonCand->daughter("muon1")->p4();
           }
           muonP_p4.SetPtEtaPhiM(vP.pt(), vP.eta(), vP.phi(), vP.mass());
           muonN_p4.SetPtEtaPhiM(vM.pt(), vM.eta(), vM.phi(), vM.mass());
           MassErr = dimuonCand->userFloat("MassErr");
           vProb = dimuonCand->userFloat("vProb");
           DCA = dimuonCand->userFloat("DCA");
           ppdlPV = dimuonCand->userFloat("ppdlPV");
           ppdlErrPV = dimuonCand->userFloat("ppdlErrPV");
           ppdlBS = dimuonCand->userFloat("ppdlBS");
           ppdlErrBS = dimuonCand->userFloat("ppdlErrBS");
           cosAlpha = dimuonCand->userFloat("cosAlpha");
           charge = dimuonCand->charge(); 
           TVector3 pperp(dimuonCand->px(), dimuonCand->py(), 0);
           lxyPV = ppdlPV * pperp.Perp() / dimuonCand->mass();
           lxyBS = ppdlBS * pperp.Perp() / dimuonCand->mass();
           irank++;
           onia_tree->Fill();
           if (OnlyBest_) break;
        } 
     }
  } else {
     if (dimuon_pdgId && (OnlyGen_ || isMC_)) onia_tree->Fill();
     //else std::cout << "Onia2MuMuRootupler: does not find a valid dimuon combination " << run << "," << event << std::endl;
  }

  dimuon_p4.SetPtEtaPhiM(0.,0.,0.,0.);
  gen_dimuon_p4.SetPtEtaPhiM(0.,0.,0.,0.);
  dimuon_pdgId = 0;
  mother_pdgId = 0;
  vProb=-1;
}

// ------------ method called once each job just before starting event loop  ------------
void Onia2MuMuRootupler::beginJob() {}

// ------------ method called once each job just after ending the event loop  ------------
void Onia2MuMuRootupler::endJob() {}

// ------------ method called when starting to processes a run  ------------
void Onia2MuMuRootupler::beginRun(edm::Run const &, edm::EventSetup const &) {}

// ------------ method called when ending the processing of a run  ------------
void Onia2MuMuRootupler::endRun(edm::Run const &, edm::EventSetup const &) {}

// ------------ method called when starting to processes a luminosity block  ------------
void Onia2MuMuRootupler::beginLuminosityBlock(edm::LuminosityBlock const &, edm::EventSetup const &) {}

// ------------ method called when ending the processing of a luminosity block  ------------
void Onia2MuMuRootupler::endLuminosityBlock(edm::LuminosityBlock const &, edm::EventSetup const &) {}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void Onia2MuMuRootupler::fillDescriptions(edm::ConfigurationDescriptions & descriptions) {
	//The following says we do not know what parameters are allowed so do no validation
	// Please change this to state exactly what you do use, even if it is no parameters
	edm::ParameterSetDescription desc;
	desc.setUnknown();
	descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(Onia2MuMuRootupler);
