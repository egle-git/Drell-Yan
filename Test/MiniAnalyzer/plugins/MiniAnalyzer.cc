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

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

//classes to extract Muon information
#include "DataFormats/PatCandidates/interface/Muon.h"

#include "DataFormats/Math/interface/LorentzVector.h"

// class declaration

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.


class MiniAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources> {
   public:
      explicit MiniAnalyzer(const edm::ParameterSet&);
      ~MiniAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      // ----------member data ---------------------------

      edm::EDGetTokenT<std::vector<pat::Muon>> muonToken_;
      TH1D *h_muon_pt;
      TH1D *h_muon_eta;
      TH1D *h_muon_phi;
      TH1D *h_muon_energy;
      TH1D *h_muon_mass;
      TH1D *h_muon_leading;
      TH1D *h_muon_subleading;
      TH1D *h_Z_pt;
      TH1D *h_Z_eta;
      TH1D *h_Z_phi;
      TH1D *h_Z_energy;
      TH1D *h_Z_mass;
      TH1D *h_Z_mass_eq;
      TH1D *h_Z_mass_fine;
      TFile *fs;
};

// constants, enums and typedefs
// static data member definitions
// constructors and destructor


MiniAnalyzer::MiniAnalyzer(const edm::ParameterSet& iConfig): muonToken_(consumes<std::vector<pat::Muon>>(iConfig.getParameter<edm::InputTag>("muons")))

{
   usesResource("TFileService");
   h_muon_pt = new TH1D("h_muon_pt", "Muon PT", 100, 0, 150);
   h_muon_eta = new TH1D("h_muon_eta", "Muon ETA", 100, -2.5, 2.5);
   h_muon_phi = new TH1D("h_muon_phi", "Muon PHI", 100, -3.14, 3.14);
   h_muon_energy = new TH1D("h_muon_energy", "Muon ENERGY", 100, 0, 150);
   h_muon_mass = new TH1D("h_muon_mass", "Muon MASS", 100, 0.1055, 0.1059);
   h_muon_leading = new TH1D("h_muon_leading", "Muon LEADING", 100, 0, 150);
   h_muon_subleading = new TH1D("h_muon_subleading", "Muon SUBLEADING", 100, 0, 150);

   double bins[37]= {40,45,50,55,60,64,68,72,76,81,86,91,96,101,106,110,115,120,126,133,141,150,160,171,185,200,220,243,273,320,380,440,510,600,700,830,1000};
   int nbins = 36;

   h_Z_pt = new TH1D("h_Z_pt", "Z Boson PT", 100, 0, 300);
   h_Z_eta = new TH1D("h_Z_eta", "Z Boson ETA", 100, -2.5, 2.5);
   h_Z_phi = new TH1D("h_Z_phi", "Z Boson PHI", 100, -3.14, 3.14);
   h_Z_energy = new TH1D("h_Z_energy", "Z Boson ENERGY", 100, 0, 500);
   h_Z_mass = new TH1D("h_Z_mass", "Z Boson MASS", nbins, bins);
   h_Z_mass_eq = new TH1D("h_Z_mass_eq", "Z Boson MASS", 1000, 30, 1000);
   h_Z_mass_fine = new TH1D("h_Z_mass_fine", "Z BOSON MASS", 150, 70, 110);
}


MiniAnalyzer::~MiniAnalyzer()
{
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
}

// member functions
// ------------ method called for each event  ------------

void
MiniAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   edm::Handle<std::vector<pat::Muon>> muons;
   iEvent.getByToken(muonToken_, muons);

   if(muons.isValid() && muons->size() > 2){
      const pat::Muon* muon1 = nullptr;
      const pat::Muon* muon2 = nullptr;
      for (const auto& muon : *muons){
         double isoSum = muon.pfIsolationR03().sumChargedHadronPt + muon.pfIsolationR03().sumNeutralHadronEt + muon.pfIsolationR03().sumPhotonEt;
         double muonPt = muon.pt();
         double relIso = isoSum / muonPt;

         if (isoSum < 100 && relIso < 0.15){
            if (!muon1 || muon.pt() > muon1->pt()) {
               muon2 = muon1;
               muon1 = &muon;
            }
            else if (!muon2 || muon.pt() > muon2->pt()) {
               muon2 = &muon;
            }
         }
      }
      if (muon1 && muon2){
         
         double pt1 = muon1->pt();
         double eta1 = muon1->eta();
         double phi1 = muon1->phi();
         double energy1 = muon1->energy();
         double mass1 = muon1->mass();

         double pt2 = muon2->pt();
         double eta2 = muon2->eta();
         double phi2 = muon2->phi();
         double energy2 = muon2->energy();
         double mass2 = muon2->mass();

         if (pt1>=20 && pt2>=12) {
            h_muon_pt->Fill(pt1);
            h_muon_pt->Fill(pt2);
            h_muon_eta->Fill(eta1);
            h_muon_eta->Fill(eta2);
            h_muon_phi->Fill(phi1);
            h_muon_phi->Fill(phi2);
            h_muon_energy->Fill(energy1);
            h_muon_energy->Fill(energy2);
            h_muon_mass->Fill(mass1);
            h_muon_mass->Fill(mass2);
            h_muon_leading->Fill(pt1);
            h_muon_subleading->Fill(pt2);
         
            // std::cout << "Muon 1: pt=" << pt1 << ", eta=" << eta1 << ", phi=" << phi1 << ", energy=" << energy1 << ", mass=" << mass1 << std::endl;
            // std::cout << "Muon 2: pt=" << pt2 << ", eta=" << eta2 << ", phi=" << phi2 << ", energy=" << energy2 << ", mass=" << mass2 << std::endl;

            math::PtEtaPhiELorentzVector muon1P4(pt1, eta1, phi1, energy1);
            math::PtEtaPhiELorentzVector muon2P4(pt2, eta2, phi2, energy2);
            auto ZbosonP4 = muon1P4 + muon2P4;

            double Zboson_pt = ZbosonP4.pt();
            double Zboson_eta = ZbosonP4.eta();
            double Zboson_phi = ZbosonP4.phi();
            double Zboson_energy = ZbosonP4.energy();
            double Zboson_mass = ZbosonP4.mass ();

            h_Z_pt->Fill(Zboson_pt);
            h_Z_eta->Fill(Zboson_eta);
            h_Z_phi->Fill(Zboson_phi);
            h_Z_energy->Fill(Zboson_energy);
            h_Z_mass->Fill(Zboson_mass);
            h_Z_mass_eq->Fill(Zboson_mass);
            h_Z_mass_fine->Fill(Zboson_mass);

            // std::cout << "Z boson: pt=" << Zboson_pt << ", eta=" << Zboson_eta << ", phi=" << Zboson_phi << ", energy=" << Zboson_energy << ", mass=" << Zboson_mass<< std::endl;
         }
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
}
*/

// ------------ method called once each job just before starting event loop  ------------
void 
MiniAnalyzer::beginJob()
{
   fs = new TFile("rootoutput/output.root","RECREATE");
}

// ------------ method called once each job just after ending the event loop  ------------
void 
MiniAnalyzer::endJob() 
{
   fs->cd();
   h_muon_pt->Write();
   h_muon_eta->Write();
   h_muon_phi->Write();
   h_muon_energy->Write();
   h_muon_mass->Write();
   h_muon_leading->Write();
   h_muon_subleading->Write();
   h_Z_pt->Write();
   h_Z_eta->Write();
   h_Z_phi->Write();
   h_Z_energy->Write();
   h_Z_mass->Write();
   h_Z_mass_eq->Write();
   h_Z_mass_fine->Write();
   fs->Close();
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
MiniAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(MiniAnalyzer);
