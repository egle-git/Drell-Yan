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
#include "TTree.h"

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


class MiniAnalyzerTest : public edm::one::EDAnalyzer<edm::one::SharedResources> {
    public:
        explicit MiniAnalyzerTest(const edm::ParameterSet&);
        ~MiniAnalyzerTest();

        static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


    private:
        virtual void beginJob() override;
        virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
        virtual void endJob() override;

        // ----------member data ---------------------------

        edm::EDGetTokenT<std::vector<pat::Muon>> muonToken_;
    
        TFile *fs;
        TTree *tree;

        float muon_pt1, muon_eta1, muon_phi1, muon_energy1, muon_mass1, muon_p1, muon_px1, muon_py1, muon_pz1, muon_charge1;
        float muon_pt2, muon_eta2, muon_phi2, muon_energy2, muon_mass2, muon_p2, muon_px2, muon_py2, muon_pz2, muon_charge2;
        float Z_pt, Z_eta, Z_phi, Z_energy, Z_mass, Z_px, Z_py, Z_pz;
};

// constants, enums and typedefs
// static data member definitions
// constructors and destructor


MiniAnalyzerTest::MiniAnalyzerTest(const edm::ParameterSet& iConfig): muonToken_(consumes<std::vector<pat::Muon>>(iConfig.getParameter<edm::InputTag>("muons")))

{
    usesResource("TFileService");
}


MiniAnalyzerTest::~MiniAnalyzerTest()
{
    // do anything here that needs to be done at desctruction time
    // (e.g. close files, deallocate resources etc.)
}

// member functions
// ------------ method called for each event  ------------

void MiniAnalyzerTest::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
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
         
            muon_pt1 = muon1->pt();
            muon_eta1 = muon1->eta();
            muon_phi1 = muon1->phi();
            muon_energy1 = muon1->energy();
            muon_mass1 = muon1->mass();
            muon_p1 = muon1->p();
            muon_px1 = muon1->px();
            muon_py1 = muon1->py();
            muon_pz1 = muon1->pz();
            muon_charge1 = muon1->charge();

            muon_pt2 = muon2->pt();
            muon_eta2 = muon2->eta();
            muon_phi2 = muon2->phi();
            muon_energy2 = muon2->energy();
            muon_mass2 = muon2->mass();
            muon_p2 = muon1->p();
            muon_px2 = muon1->px();
            muon_py2 = muon1->py();
            muon_pz2 = muon1->pz();
            muon_charge2 = muon1->charge();

            if (muon_pt1>=20 && muon_pt2>=12) {
                math::PtEtaPhiELorentzVector muon1P4(muon_pt1, muon_eta1, muon_phi1, muon_energy1);
                math::PtEtaPhiELorentzVector muon2P4(muon_pt2, muon_eta2, muon_phi2, muon_energy2);
                auto ZbosonP4 = muon1P4 + muon2P4;

                Z_pt = ZbosonP4.pt();
                Z_eta = ZbosonP4.eta();
                Z_phi = ZbosonP4.phi();
                Z_energy = ZbosonP4.energy();
                Z_mass = ZbosonP4.mass();
                Z_px = ZbosonP4.px();
                Z_py = ZbosonP4.py();
                Z_pz = ZbosonP4.pz();

                tree->Fill();
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
void MiniAnalyzerTest::beginJob()
{
    fs = new TFile("rootoutputs/output_tree_real.root","RECREATE");
    tree = new TTree("Events", "");

    tree->Branch("muon_pt1", &muon_pt1);
    tree->Branch("muon_eta1", &muon_eta1);
    tree->Branch("muon_phi1", &muon_phi1);
    tree->Branch("muon_energy1", &muon_energy1);
    tree->Branch("muon_mass1", &muon_mass1);
    tree->Branch("muon_p1", &muon_p1);
    tree->Branch("muon_px1", &muon_px1);
    tree->Branch("muon_py1", &muon_py1);
    tree->Branch("muon_pz1", &muon_pz1);
    tree->Branch("muon_charge1", &muon_charge1);

    tree->Branch("muon_pt2", &muon_pt2);
    tree->Branch("muon_eta2", &muon_eta2);
    tree->Branch("muon_phi2", &muon_phi2);
    tree->Branch("muon_energy2", &muon_energy2);
    tree->Branch("muon_mass2", &muon_mass2);
    tree->Branch("muon_p2", &muon_p2);
    tree->Branch("muon_px2", &muon_px2);
    tree->Branch("muon_py2", &muon_py2);
    tree->Branch("muon_pz2", &muon_pz2);
    tree->Branch("muon_charge2", &muon_charge2);

    tree->Branch("Z_pt", &Z_pt);
    tree->Branch("Z_eta", &Z_eta);
    tree->Branch("Z_phi", &Z_phi);
    tree->Branch("Z_energy", &Z_energy);
    tree->Branch("Z_mass", &Z_mass);
    tree->Branch("Z_px", &Z_px);
    tree->Branch("Z_py", &Z_py);
    tree->Branch("Z_pz", &Z_pz);
}

// ------------ method called once each job just after ending the event loop  ------------
void MiniAnalyzerTest::endJob() 
{
    fs->cd();
    tree->Write();
    fs->Close();
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void MiniAnalyzerTest::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    //The following says we do not know what parameters are allowed so do no validation
    // Please change this to state exactly what you do use, even if it is no parameters
    edm::ParameterSetDescription desc;
    desc.setUnknown();
    descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(MiniAnalyzerTest);
