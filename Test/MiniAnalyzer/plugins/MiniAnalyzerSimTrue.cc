// -*- C++ -*-
//
// Package:    Test/MiniAnalyzer
// Class:      MiniAnalyzer
// 
/**\class MiniAnalyzer MiniAnalyzer.cc Test/MiniAnalyzer/plugins/MiniAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  
//         Created:  Thu, 01 Feb 2024 02:16:09 GMT


// system include files
#include <memory>
#include <vector>
#include "TFile.h"
#include "TH1D.h"
#include <fstream>
#include <iostream>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

//classes to extract Muon information
#include "DataFormats/PatCandidates/interface/Muon.h"

#include "DataFormats/Math/interface/LorentzVector.h"

#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"


// class declaration

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

// edm::one::WatchRuns /////after SharedResources, before >

class MiniAnalyzerSimTrue : public edm::one::EDAnalyzer<edm::one::SharedResources> {
   public:
      explicit MiniAnalyzerSimTrue(const edm::ParameterSet&);
      ~MiniAnalyzerSimTrue();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      edm::EDGetTokenT<std::vector<reco::GenParticle>> GenParticleToken_;
      //edm::EDGetTokenT<std::vector<reco::GenParticle>> GenParticleToken_;
      edm::EDGetTokenT<GenEventInfoProduct> weightToken_;


      TH1D *simh_particle_pt;
      TH1D *simh_particle_eta;
      TH1D *simh_particle_phi;
      TH1D *simh_particle_energy;
      TH1D *simh_particle_mass;
      TH1D *simh_particle_leading;
      TH1D *simh_particle_subleading;
      TH1D *simh_Z_pt;
      TH1D *simh_Z_eta;
      TH1D *simh_Z_phi;
      TH1D *simh_Z_energy;
      TH1D *simh_Z_mass;
      TH1D *simh_Z_mass_fine;
      TFile *fs;

      double weight_sum; 
      double xsec = 6422; // 6422 for first, 20480 for second
      double lumi = 16494;
};

// constants, enums and typedefs
// static data member definitions
// constructors and destructor


MiniAnalyzerSimTrue::MiniAnalyzerSimTrue(const edm::ParameterSet& iConfig):
      GenParticleToken_(consumes<std::vector<reco::GenParticle>>(iConfig.getUntrackedParameter<edm::InputTag>("GenParticle"))),
      weightToken_(consumes<GenEventInfoProduct>(iConfig.getUntrackedParameter<edm::InputTag>("GenEventInfo")))

{
   usesResource("TFileService");
   simh_particle_pt = new TH1D("simh_particle_pt", "Particle PT", 100, 0, 150);
   simh_particle_eta = new TH1D("simh_particle_eta", "Particle ETA", 100, -2.5, 2.5);
   simh_particle_phi = new TH1D("simh_particle_phi", "Particle PHI", 100, -3.14, 3.14);
   simh_particle_energy = new TH1D("simh_particle_energy", "Particle ENERGY", 100, 0, 150);
   simh_particle_mass = new TH1D("simh_particle_mass", "Particle MASS", 100, 0.1055, 0.1059);
   simh_particle_leading = new TH1D("simh_particle_leading", "Particle LEADING", 100, 0, 150);
   simh_particle_subleading = new TH1D("simh_particle_subleading", "Particle SUBLEADING", 100, 0, 150);

   simh_Z_pt = new TH1D("simh_Z_pt", "Z Boson PT", 100, 0, 300);
   simh_Z_eta = new TH1D("simh_Z_eta", "Z Boson ETA", 100, -2.5, 2.5);
   simh_Z_phi = new TH1D("simh_Z_phi", "Z Boson PHI", 100, -3.14, 3.14);
   simh_Z_energy = new TH1D("simh_Z_energy", "Z Boson ENERGY", 100, 0, 500);
   simh_Z_mass = new TH1D("simh_Z_mass", "Z Boson MASS", 1000, 30, 1000);
   simh_Z_mass_fine = new TH1D("simh_Z_mass_fine", "Z BOSON MASS", 150, 70, 110);
}


MiniAnalyzerSimTrue::~MiniAnalyzerSimTrue()
{
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
}


void
MiniAnalyzerSimTrue::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

   using namespace edm;

   edm::Handle<GenEventInfoProduct> weightHandle;
   iEvent.getByToken(weightToken_, weightHandle);
   double event_weight = weightHandle.isValid() ? weightHandle->weight() : 1.0;
   double norm_weight = event_weight / std::abs(event_weight);
   double weight = norm_weight * xsec * lumi / weight_sum;
   // std::cout << "weight: " << weight << std::endl;

   edm::Handle<std::vector<reco::GenParticle>> particles;
   iEvent.getByToken(GenParticleToken_, particles);

   if(particles.isValid() && particles->size() >= 2){
      std::vector<reco::GenParticle> selectedparticles;
      for (const auto& genParticle : *particles){
         if (abs(genParticle.pdgId()) == 13 && genParticle.fromHardProcessFinalState() == true) { //&& genParticle.status() == 1  parCand.fromHardProcessFinalState()
            selectedparticles.push_back(genParticle);
         }
      }
      // std::cout << "there are " <<  selectedparticles.size()<<" particles" << std::endl;

      if (selectedparticles.size() == 2){
         
         double pt1 = selectedparticles[0].pt();
         double eta1 = selectedparticles[0].eta();
         double phi1 = selectedparticles[0].phi();
         double energy1 = selectedparticles[0].energy();
         double mass1 = selectedparticles[0].mass();

         double pt2 = selectedparticles[1].pt();
         double eta2 = selectedparticles[1].eta();
         double phi2 = selectedparticles[1].phi();
         double energy2 = selectedparticles[1].energy();
         double mass2 = selectedparticles[1].mass();

         simh_particle_pt->Fill(pt1, weight);
         simh_particle_pt->Fill(pt2, weight);
         simh_particle_eta->Fill(eta1, weight);
         simh_particle_eta->Fill(eta2, weight);
         simh_particle_phi->Fill(phi1, weight);
         simh_particle_phi->Fill(phi2, weight);
         simh_particle_energy->Fill(energy1, weight);
         simh_particle_energy->Fill(energy2, weight);
         simh_particle_mass->Fill(mass1, weight);
         simh_particle_mass->Fill(mass2, weight);
         simh_particle_leading->Fill(pt1, weight);
         simh_particle_subleading->Fill(pt2, weight);
         
         // std::cout << "Particle 1: pt=" << pt1 << ", eta=" << eta1 << ", phi=" << phi1 << ", energy=" << energy1 << ", mass=" << mass1 << std::endl;
         // std::cout << "Particle 2: pt=" << pt2 << ", eta=" << eta2 << ", phi=" << phi2 << ", energy=" << energy2 << ", mass=" << mass2 << std::endl;

         math::PtEtaPhiELorentzVector particle1P4(pt1, eta1, phi1, energy1);
         math::PtEtaPhiELorentzVector particle2P4(pt2, eta2, phi2, energy2);
         auto ZbosonP4 = particle1P4 + particle2P4;

         double Zboson_pt = ZbosonP4.pt();
         double Zboson_eta = ZbosonP4.eta();
         double Zboson_phi = ZbosonP4.phi();
         double Zboson_energy = ZbosonP4.energy();
         double Zboson_mass = ZbosonP4.mass ();

         simh_Z_pt->Fill(Zboson_pt, weight);
         simh_Z_eta->Fill(Zboson_eta, weight);
         simh_Z_phi->Fill(Zboson_phi, weight);
         simh_Z_energy->Fill(Zboson_energy, weight);
         simh_Z_mass->Fill(Zboson_mass, weight);
         simh_Z_mass_fine->Fill(Zboson_mass, weight);

         // std::cout << "Z boson: pt=" << Zboson_pt << ", eta=" << Zboson_eta << ", phi=" << Zboson_phi << ", energy=" << Zboson_energy << ", mass=" << Zboson_mass<< std::endl;
      }
   }
}
   

/** 
#ifdef THIS_IS_AN_EVENT_EXAMPLE
   Handle<ExampleData> pIn;
   iEvent.getByLabel("example",pIn);
#endif
   
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
#endif
*/

// ------------ method called once each job just before starting event loop  ------------
void 
MiniAnalyzerSimTrue::beginJob()
{
   fs = new TFile("rootoutputs/simoutputtrue1.root","RECREATE"); //simoutputtrue1.root for first, simoutputtrue2.root for second

   std::ifstream inFile("weightsums/weight_sum1.txt"); //weight_sum1.txt for first, weight_sum2.txt for second
   if (inFile.is_open())
   {
      inFile >> weight_sum;
      inFile.close();
   }
}

// ------------ method called once each job just after ending the event loop  ------------
void 
MiniAnalyzerSimTrue::endJob() 
{
   fs->cd();
   simh_particle_pt->Write();
   simh_particle_eta->Write();
   simh_particle_phi->Write();
   simh_particle_energy->Write();
   simh_particle_mass->Write();
   simh_particle_leading->Write();
   simh_particle_subleading->Write();
   simh_Z_pt->Write();
   simh_Z_eta->Write();
   simh_Z_phi->Write();
   simh_Z_energy->Write();
   simh_Z_mass->Write();
   simh_Z_mass_fine->Write();
   fs->Close();

}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
MiniAnalyzerSimTrue::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
   //The following says we do not know what parameters are allowed so do no validation
   // Please change this to state exactly what you do use, even if it is no parameters
   edm::ParameterSetDescription desc;
   desc.setUnknown();
   descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(MiniAnalyzerSimTrue);
