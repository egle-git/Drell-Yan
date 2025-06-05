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
#include <fstream>
#include <cmath>
#include <bits/stdc++.h>

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

#include "DataFormats/VertexReco/interface/Vertex.h"

// class declaration

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

// edm::one::WatchRuns /////after SharedResources, before >

class MiniAnalyzerTreeSim12 : public edm::one::EDAnalyzer<edm::one::SharedResources> {
    public:
        explicit MiniAnalyzerTreeSim12(const edm::ParameterSet&);
        ~MiniAnalyzerTreeSim12();

        static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


    private:
        virtual void beginJob() override;
        virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
        virtual void endJob() override;


        edm::EDGetTokenT<std::vector<pat::Muon>> muonToken_;
        edm::EDGetTokenT<std::vector<reco::GenParticle>> GenParticleToken_;
        edm::EDGetTokenT<GenEventInfoProduct> weightToken_;
        edm::EDGetTokenT<std::vector<reco::Vertex>> vertexToken_;

        TFile *fs;
        TTree *tree_muon;
        TTree *tree_tau;

        float muon_pt1, muon_eta1, muon_phi1, muon_energy1, muon_mass1, muon_p1, muon_px1, muon_py1, muon_pz1, muon_charge1;
        float muon_pt2, muon_eta2, muon_phi2, muon_energy2, muon_mass2, muon_p2, muon_px2, muon_py2, muon_pz2, muon_charge2;
        bool muon_loose1, muon_medium1, muon_tight1, muon_loose2, muon_medium2, muon_tight2;
        float Z_pt, Z_eta, Z_phi, Z_energy, Z_mass, Z_px, Z_py, Z_pz;
        float muon_isoSum1, muon_isoSum2, muon_isoSumCorr1, muon_isoSumCorr2, muon_relIso1, muon_relIso2, muon_event_weight, muon_norm_weight;

        float tau_pt1, tau_eta1, tau_phi1, tau_energy1, tau_mass1, tau_p1, tau_px1, tau_py1, tau_pz1, tau_charge1;
        float tau_pt2, tau_eta2, tau_phi2, tau_energy2, tau_mass2, tau_p2, tau_px2, tau_py2, tau_pz2, tau_charge2;
        float tauZ_pt, tauZ_eta, tauZ_phi, tauZ_energy, tauZ_mass, tauZ_px, tauZ_py, tauZ_pz;
        float tau_isoSum1, tau_isoSum2, tau_isoSumCorr1, tau_isoSumCorr2, tau_relIso1, tau_relIso2, tau_event_weight, tau_norm_weight;

        std::string mcProcess_;
};

// constants, enums and typedefs
// static data member definitions
// constructors and destructor


MiniAnalyzerTreeSim12::MiniAnalyzerTreeSim12(const edm::ParameterSet& iConfig):
        muonToken_(consumes<std::vector<pat::Muon>>(iConfig.getParameter<edm::InputTag>("muons"))),
        GenParticleToken_(consumes<std::vector<reco::GenParticle>>(iConfig.getUntrackedParameter<edm::InputTag>("GenParticle"))),
        weightToken_(consumes<GenEventInfoProduct>(iConfig.getUntrackedParameter<edm::InputTag>("GenEventInfo"))),
        vertexToken_(consumes<std::vector<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("vertices"))),
        mcProcess_(iConfig.getParameter<std::string>("mcProcess"))

{
    usesResource("TFileService");
}


MiniAnalyzerTreeSim12::~MiniAnalyzerTreeSim12()
{
    // do anything here that needs to be done at desctruction time
    // (e.g. close files, deallocate resources etc.)
}


void
MiniAnalyzerTreeSim12::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

    using namespace edm;
    using namespace std;

    edm::Handle<GenEventInfoProduct> weightHandle;
    edm::Handle<std::vector<pat::Muon>> muons;
    edm::Handle<std::vector<reco::Vertex>> vertices;
    iEvent.getByToken(weightToken_, weightHandle);
    iEvent.getByToken(muonToken_, muons);
    iEvent.getByToken(vertexToken_, vertices);

    double event_weight = weightHandle.isValid() ? weightHandle->weight() : 1.0;
    double norm_weight = event_weight / std::abs(event_weight);

    if(muons.isValid() && muons->size() > 2){
        const pat::Muon* muon1 = nullptr;
        const pat::Muon* muon2 = nullptr;
        const reco::Vertex& primaryVertex = (*vertices)[0];
        for (const auto& muon : *muons){
            isoSum = muon.pfIsolationR03().sumChargedHadronPt + muon.pfIsolationR03().sumNeutralHadronEt + muon.pfIsolationR03().sumPhotonEt;
            isoSumCorr = muon.pfIsolationR03().sumChargedHadronPt + max(0., muon.pfIsolationR03().sumNeutralHadronEt + muon.pfIsolationR03().sumPhotonEt - 0.5 * muon.pfIsolationR03().sumPUPt);
            double muonPt = muon.pt();
            relIso = isoSum / muonPt;

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

            double pt2 = muon2->pt();
            double eta2 = muon2->eta();
            double phi2 = muon2->phi();
            double energy2 = muon2->energy();

            math::PtEtaPhiELorentzVector muon1P4(pt1, eta1, phi1, energy1);
            math::PtEtaPhiELorentzVector muon2P4(pt2, eta2, phi2, energy2);
            auto ZbosonP4 = muon1P4 + muon2P4;

            if (pt1>=20 && pt2>=12) {
                edm::Handle<std::vector<reco::GenParticle>> genparticles;
                iEvent.getByToken(GenParticleToken_, genparticles);

                bool MuonsFinalState=true, TauFinalState=false;

                if(genparticles.isValid() && genparticles->size() >= 2){
                    std::vector<reco::GenParticle> selectedparticles;
                    for (const auto& genParticle : *genparticles){
                        if (genParticle.isHardProcess() && (std::abs(genParticle.pdgId()) == 13 || std::abs(genParticle.pdgId()) == 15)) { //&& genParticle.status() == 1  parCand.fromHardProcessFinalState()
                            selectedparticles.push_back(genParticle);
                        }
                    }
                    if (selectedparticles.size() == 2 && std::abs(selectedparticles[0].pdgId()) == 13 && std::abs(selectedparticles[1].pdgId()) == 13) {
                        MuonsFinalState = true;
                        TauFinalState = false;
                    }
                    else if (selectedparticles.size() == 2 && std::abs(selectedparticles[0].pdgId()) == 15 && std::abs(selectedparticles[1].pdgId()) == 15) {
                        MuonsFinalState = false;
                        TauFinalState = true;
                    }
                }

                if (MuonsFinalState) {
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
                    muon_loose1 = muon1->isLooseMuon();
                    muon_medium1 = muon1->isMediumMuon();
                    muon_tight1 = muon1->isTightMuon(primaryVertex);
                    muon_isoSum1 = muon1.pfIsolationR03().sumChargedHadronPt + muon1.pfIsolationR03().sumNeutralHadronEt + muon1.pfIsolationR03().sumPhotonEt;
                    muon_isoSumCorr1 = muon1.pfIsolationR03().sumChargedHadronPt + max(0., muon1.pfIsolationR03().sumNeutralHadronEt + muon1.pfIsolationR03().sumPhotonEt - 0.5 * muon1.pfIsolationR03().sumPUPt);
                    muon_relIso1 = muon_isoSum1 / muon_pt1;

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
                    muon_loose2 = muon2->isLooseMuon();
                    muon_medium2 = muon2->isMediumMuon();
                    muon_tight2 = muon2->isTightMuon(primaryVertex);
                    muon_isoSum2 = muon2.pfIsolationR03().sumChargedHadronPt + muon2.pfIsolationR03().sumNeutralHadronEt + muon2.pfIsolationR03().sumPhotonEt;
                    muon_isoSumCorr2 = muon2.pfIsolationR03().sumChargedHadronPt + max(0., muon2.pfIsolationR03().sumNeutralHadronEt + muon2.pfIsolationR03().sumPhotonEt - 0.5 * muon2.pfIsolationR03().sumPUPt);
                    muon_relIso2 = muon_isoSum2 / muon_pt2;

                    muon_event_weight = event_weight;
                    muon_norm_weight = norm_weight;

                    Z_pt = ZbosonP4.pt();
                    Z_eta = ZbosonP4.eta();
                    Z_phi = ZbosonP4.phi();
                    Z_energy = ZbosonP4.energy();
                    Z_mass = ZbosonP4.mass();
                    Z_px = ZbosonP4.px();
                    Z_py = ZbosonP4.py();
                    Z_pz = ZbosonP4.pz();

                    tree_muon->Fill();
                }
                else if (TauFinalState) {
                    tau_pt1 = muon1->pt();
                    tau_eta1 = muon1->eta();
                    tau_phi1 = muon1->phi();
                    tau_energy1 = muon1->energy();
                    tau_mass1 = muon1->mass();
                    tau_p1 = muon1->p();
                    tau_px1 = muon1->px();
                    tau_py1 = muon1->py();
                    tau_pz1 = muon1->pz();
                    tau_charge1 = muon1->charge();
                    tau_isoSum1 = muon1.pfIsolationR03().sumChargedHadronPt + muon1.pfIsolationR03().sumNeutralHadronEt + muon1.pfIsolationR03().sumPhotonEt;
                    tau_isoSumCorr1 = muon1.pfIsolationR03().sumChargedHadronPt + max(0., muon1.pfIsolationR03().sumNeutralHadronEt + muon1.pfIsolationR03().sumPhotonEt - 0.5 * muon1.pfIsolationR03().sumPUPt);
                    tau_relIso1 = tau_isoSum1 / tau_pt1;

                    tau_pt2 = muon2->pt();
                    tau_eta2 = muon2->eta();
                    tau_phi2 = muon2->phi();
                    tau_energy2 = muon2->energy();
                    tau_mass2 = muon2->mass();
                    tau_p2 = muon1->p();
                    tau_px2 = muon1->px();
                    tau_py2 = muon1->py();
                    tau_pz2 = muon1->pz();
                    tau_charge2 = muon1->charge();
                    tau_isoSum1 = muon2.pfIsolationR03().sumChargedHadronPt + muon2.pfIsolationR03().sumNeutralHadronEt + muon2.pfIsolationR03().sumPhotonEt;
                    tau_isoSumCorr1 = muon2.pfIsolationR03().sumChargedHadronPt + max(0., muon2.pfIsolationR03().sumNeutralHadronEt + muon2.pfIsolationR03().sumPhotonEt - 0.5 * muon2.pfIsolationR03().sumPUPt);
                    tau_relIso1 = tau_isoSum2 / tau_pt2;

                    tau_event_weight = event_weight;
                    tau_norm_weight = norm_weight;

                    tauZ_pt = ZbosonP4.pt();
                    tauZ_eta = ZbosonP4.eta();
                    tauZ_phi = ZbosonP4.phi();
                    tauZ_energy = ZbosonP4.energy();
                    tauZ_mass = ZbosonP4.mass();
                    tauZ_px = ZbosonP4.px();
                    tauZ_py = ZbosonP4.py();
                    tauZ_pz = ZbosonP4.pz();

                    tree_tau->Fill();
                }
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
*/

// ------------ method called once each job just before starting event loop  ------------
void 
MiniAnalyzerTreeSim12::beginJob()
{

    std::string outputfile;
    if (mcProcess_ == "sim1")
        outputfile = "rootoutputs/treeout1Test.root";
    else if (mcProcess_ == "sim2")
        outputfile = "rootoutputs/treeout2Test.root";
    fs = new TFile(outputfile.c_str(), "RECREATE");
    tree_muon = new TTree("Muon", "");
    tree_tau = new TTree("Tau", "");

    tree_muon->Branch("muon_pt1", &muon_pt1);
    tree_muon->Branch("muon_eta1", &muon_eta1);
    tree_muon->Branch("muon_phi1", &muon_phi1);
    tree_muon->Branch("muon_energy1", &muon_energy1);
    tree_muon->Branch("muon_mass1", &muon_mass1);
    tree_muon->Branch("muon_p1", &muon_p1);
    tree_muon->Branch("muon_px1", &muon_px1);
    tree_muon->Branch("muon_py1", &muon_py1);
    tree_muon->Branch("muon_pz1", &muon_pz1);
    tree_muon->Branch("muon_charge1", &muon_charge1);
    tree_muon->Branch("muon_loose1", &muon_loose1);
    tree_muon->Branch("muon_medium1", &muon_medium1);
    tree_muon->Branch("muon_tight1", &muon_tight1);
    tree_muon->Branch("muon_isoSum1", &muon_isoSum1);
    tree_muon->Branch("muon_isoSumCorr1", &muon_isoSumCorr1);
    tree_muon->Branch("muon_relIso1", &muon_relIso1);

    tree_muon->Branch("muon_pt2", &muon_pt2);
    tree_muon->Branch("muon_eta2", &muon_eta2);
    tree_muon->Branch("muon_phi2", &muon_phi2);
    tree_muon->Branch("muon_energy2", &muon_energy2);
    tree_muon->Branch("muon_mass2", &muon_mass2);
    tree_muon->Branch("muon_p2", &muon_p2);
    tree_muon->Branch("muon_px2", &muon_px2);
    tree_muon->Branch("muon_py2", &muon_py2);
    tree_muon->Branch("muon_pz2", &muon_pz2);
    tree_muon->Branch("muon_charge2", &muon_charge2);
    tree_muon->Branch("muon_loose2", &muon_loose2);
    tree_muon->Branch("muon_medium2", &muon_medium2);
    tree_muon->Branch("muon_tight2", &muon_tight2);
    tree_muon->Branch("muon_isoSum2", &muon_isoSum2);
    tree_muon->Branch("muon_isoSumCorr2", &muon_isoSumCorr2);
    tree_muon->Branch("muon_relIso2", &muon_relIso2);

    tree_muon->Branch("Z_pt", &Z_pt);
    tree_muon->Branch("Z_eta", &Z_eta);
    tree_muon->Branch("Z_phi", &Z_phi);
    tree_muon->Branch("Z_energy", &Z_energy);
    tree_muon->Branch("Z_mass", &Z_mass);
    tree_muon->Branch("Z_px", &Z_px);
    tree_muon->Branch("Z_py", &Z_py);
    tree_muon->Branch("Z_pz", &Z_pz);

    tree_muon->Branch("muon_event_weight", &muon_event_weight);
    tree_muon->Branch("muon_norm_weight", &muon_norm_weight);



    tree_tau->Branch("tau_pt1", &tau_pt1);
    tree_tau->Branch("tau_eta1", &tau_eta1);
    tree_tau->Branch("tau_phi1", &tau_phi1);
    tree_tau->Branch("tau_energy1", &tau_energy1);
    tree_tau->Branch("tau_mass1", &tau_mass1);
    tree_tau->Branch("tau_p1", &tau_p1);
    tree_tau->Branch("tau_px1", &tau_px1);
    tree_tau->Branch("tau_py1", &tau_py1);
    tree_tau->Branch("tau_pz1", &tau_pz1);
    tree_tau->Branch("tau_charge1", &tau_charge1);
    tree_tau->Branch("tau_isoSum1", &tau_isoSum1);
    tree_tau->Branch("tau_isoSumCorr1", &tau_isoSumCorr1);
    tree_tau->Branch("tau_relIso1", &tau_relIso1);

    tree_tau->Branch("tau_pt2", &tau_pt2);
    tree_tau->Branch("tau_eta2", &tau_eta2);
    tree_tau->Branch("tau_phi2", &tau_phi2);
    tree_tau->Branch("tau_energy2", &tau_energy2);
    tree_tau->Branch("tau_mass2", &tau_mass2);
    tree_tau->Branch("tau_p2", &tau_p2);
    tree_tau->Branch("tau_px2", &tau_px2);
    tree_tau->Branch("tau_py2", &tau_py2);
    tree_tau->Branch("tau_pz2", &tau_pz2);
    tree_tau->Branch("tau_charge2", &tau_charge2);
    tree_tau->Branch("tau_isoSum2", &tau_isoSum2);
    tree_tau->Branch("tau_isoSumCorr2", &tau_isoSumCorr2);
    tree_tau->Branch("tau_relIso2", &tau_relIso2);

    tree_tau->Branch("tauZ_pt", &tauZ_pt);
    tree_tau->Branch("tauZ_eta", &tauZ_eta);
    tree_tau->Branch("tauZ_phi", &tauZ_phi);
    tree_tau->Branch("tauZ_energy", &tauZ_energy);
    tree_tau->Branch("tauZ_mass", &tauZ_mass);
    tree_tau->Branch("tauZ_px", &tauZ_px);
    tree_tau->Branch("tauZ_py", &tauZ_py);
    tree_tau->Branch("tauZ_pz", &tauZ_pz);

    tree_tau->Branch("tau_event_weight", &tau_event_weight);
    tree_tau->Branch("tau_norm_weight", &tau_norm_weight);



    // jei atskirai weight sum
    // std::ifstream inFile("weight_sum1.txt");
    // // weight_sum1.txt  -  sim1
    // // weight_sum2.txt  -  sim2  
    // if (inFile.is_open())
    // {
    //    inFile >> weight_sum;
    //    inFile.close();
    // }

}

// ------------ method called once each job just after ending the event loop  ------------
void 
MiniAnalyzerTreeSim12::endJob() 
{
    fs->cd();
    tree_muon->Write();
    tree_tau->Write();
    fs->Close();

//    std::string weightFilename;
//    if (mcProcess_ == "sim1")
//       weightFilename = "weightsums/TESTweight_sum1.txt";
//    else if (mcProcess_ == "sim2")
//       weightFilename = "weightsums/TESTweight_sum2.txt";
//    std::ofstream outFile(weightFilename.c_str());
//    if (outFile.is_open())
//    {
//       outFile << weight_sum << std::endl;
//       outFile.close();
//    }
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
MiniAnalyzerTreeSim12::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    //The following says we do not know what parameters are allowed so do no validation
    // Please change this to state exactly what you do use, even if it is no parameters
    edm::ParameterSetDescription desc;
    desc.setUnknown();
    descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(MiniAnalyzerTreeSim12);
