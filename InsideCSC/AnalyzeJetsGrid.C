#include "TFile.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TObjArray.h"
#include "TH1.h"
#include "TProfile.h"
#include "TH2.h"
#include "TH3.h"
#include "TParticle.h"
#include "TParticlePDG.h"
#include "TSystem.h"
#include "TROOT.h"
#include "Riostream.h"

#include "AliRunLoader.h"
#include "AliStack.h"
#include "AliHeader.h"
#include "AliGenPythiaEventHeader.h"
#include "AliFOCALCluster.h"

#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequenceArea.hh"
#ifndef __FJCORE__
#include "fastjet/GhostedAreaSpec.hh" // for area support
#endif                                // __FJCORE__

using std::cout;
using std::endl;

Float_t eta(TParticle *); // declaration, implementation at end of file
Float_t phi(TParticle *);

void GetMomentum(TLorentzVector &p, const Float_t *vertex, AliFOCALCluster *foClust, Float_t mass);
void GetGeometricalMatchingLevel(fastjet::PseudoJet &jet1, fastjet::PseudoJet &jet2, Double_t &d);
void DoJetMatching(std::vector<fastjet::PseudoJet> &jetArray1, std::vector<fastjet::PseudoJet> &jetArray2);

class JetMatchingParams : public fastjet::PseudoJet::UserInfoBase
{
public:
    // default ctor
    //  - index1 jet1 index
    //  - index2 jet2 index
    //  - d1 distance to jet1
    //  - d2 distance to jet2
    JetMatchingParams() = default;
    JetMatchingParams(const int &index1, const int &index2, const double &distance1, const double &distance2) : mIndex1(index1), mIndex2(index2), mDistance1(distance1), mDistance2(distance2) {}
    JetMatchingParams(JetMatchingParams &params) : mIndex1(params.index1()), mIndex2(params.index2()), mDistance1(params.distance1()), mDistance2(params.distance2()) {}
    JetMatchingParams(const JetMatchingParams &params) : mIndex1(params.index1()), mIndex2(params.index2()), mDistance1(params.distance1()), mDistance2(params.distance2()) {}

    /// access to the jet index
    int index1() const { return mIndex1; }
    int index2() const { return mIndex2; }

    /// access to the distance
    double distance1() const { return mDistance1; }
    double distance2() const { return mDistance2; }

    /// Set the indices and distance
    void setJetIndexDistance1(int index, double distance)
    {
        mIndex1 = index;
        mDistance1 = distance;
    }
    void setJetIndexDistance2(int index, double distance)
    {
        mIndex2 = index;
        mDistance2 = distance;
    }

protected:
    int mIndex1 = -1;        // jet1 index
    int mIndex2 = -1;        // jet2 index
    double mDistance1 = 999; // Distance to jet1
    double mDistance2 = 999; // Distance to jet2
};

class SW_JetArea : public fastjet::SelectorWorker
{
public:
    SW_JetArea(const double &minArea) : mMinArea(minArea) {}

    // the selector's description
    std::string description() const
    {
        std::ostringstream oss;
        oss << "Minimum area " << mMinArea;
        return oss.str();
    }

    bool pass(const fastjet::PseudoJet &p) const
    {
        return p.area() > mMinArea;
    }

private:
    double mMinArea; // the vertex number we're selecting
};

fastjet::Selector SelectorMinArea(const double &minArea)
{
    return fastjet::Selector(new SW_JetArea(minArea));
}

void AnalyzeJetsGrid(const char *inputdir = "./", const char *outputdir = "./")
{
    gSystem->Load("libpythia6_4_28.so");

    const Float_t Ecut = 2.0;
    const Float_t kPi0Mass = 0.135;
    const Float_t kPiPlusMass = 0.139;

    const Int_t nR = 3;
    const Float_t Rvals[nR] = {0.2, 0.4, 0.6}; // Cone radii
    // const Float_t Rvals[nR] = {0.2, 0.3, 0.4, 0.5, 0.6}; // Cone radii

    if (gSystem->AccessPathName(Form("%s/focalClusters.root", inputdir)))
    // if (gSystem->AccessPathName("focalClusters.root"))
    {
        cout << "AnalyzeJetsGrid() ERROR: focalClusters.root not found!" << endl;
        return;
    }
    TFile *clusterFile = new TFile(Form("%s/focalClusters.root", inputdir));
    // TFile *clusterFile = new TFile("focalClusters.root");
    cout << "Clusters from: " << clusterFile->GetName() << endl;

    TFile *fout = new TFile(Form("%s/analysisJets_0Mass.root", outputdir), "RECREATE"); // don't forget to change the name
    // TFile *fout = new TFile("analysisJets.root", "RECREATE");
    fout->cd();

    Float_t jetE, jetpT, jetpT_sub, jetPhi, jetEta, jetRap;
    Float_t jetE_match, jetpT_match, jetpT_sub_match, jetPhi_match, jetEta_match, jetRap_match;
    Float_t jet_distmatch;
    Int_t jetR;
    Int_t ievt;
    Int_t jetParts, jetParts_match, jetParts_truth;
    Float_t jetECALEnergyFrac, matchedjetECALEnergyFrac;

    TTree *results = new TTree("jetTree", "jets Info Tree");
    results->Branch("ievt", &ievt, "ievt/I");
    results->Branch("jetE", &jetE, "jetE/F");
    results->Branch("jetpT", &jetpT, "jetpT/F");
    results->Branch("jetpT_sub", &jetpT_sub, "jetpT_sub/F");
    results->Branch("jetPhi", &jetPhi, "jetPhi/F");
    results->Branch("jetEta", &jetEta, "jetEta/F");
    results->Branch("jetRap", &jetRap, "jetRap/F");
    results->Branch("jetParts", &jetParts, "jetParts/I");
    results->Branch("jetE_match", &jetE_match, "jetE_match/F");
    results->Branch("jetpT_match", &jetpT_match, "jetpT_match/F");
    results->Branch("jetPhi_match", &jetPhi_match, "jetPhi_match/F");
    results->Branch("jetEta_match", &jetEta_match, "jetEta_match/F");
    results->Branch("jetRap_match", &jetRap_match, "jetRap_match/F");
    results->Branch("jetParts_match", &jetParts_match, "jetParts_match/I");
    results->Branch("jet_distmatch", &jet_distmatch, "jet_distmatch/F");
    results->Branch("jetR", &jetR, "jetR/I");
    results->Branch("jetECALEnergyFrac", &jetECALEnergyFrac, "jetECALEnergyFrac/F");
    results->Branch("matchedjetECALEnergyFrac", &matchedjetECALEnergyFrac, "matchedjetECALEnergyFrac/F");

    Float_t truthjetE, truthjetpT, truthjetPhi, truthjetEta, truthjetRap;
    Float_t truthjetNeutralEnergyFrac;
    Int_t truthjetR;
    Int_t ievtTruth;

    TTree *results_truth = new TTree("truthjetTree", "MC jets Info Tree");
    results_truth->Branch("ievtTruth", &ievtTruth, "ievtTruth/I");
    results_truth->Branch("truthjetE", &truthjetE, "truthjetE/F");
    results_truth->Branch("truthjetpT", &truthjetpT, "truthjetpT/F");
    results_truth->Branch("truthjetPhi", &truthjetPhi, "truthjetPhi/F");
    results_truth->Branch("truthjetEta", &truthjetEta, "truthjetEta/F");
    results_truth->Branch("truthjetRap", &truthjetRap, "truthjetRap/F");
    results_truth->Branch("truthjetR", &truthjetR, "truthjetR/I");
    results_truth->Branch("jetParts_truth", &jetParts_truth, "jetParts_truth/I");
    results_truth->Branch("truthjetNeutralEnergyFrac", &truthjetNeutralEnergyFrac, "truthjetNeutralEnergyFrac/F");

    TH1D *nevt_h = new TH1D("nevt_h", "Number of Events", 1, 0, 1);
    TH1D *nEvPtHard_h = new TH1D("nEvPtHard_h", "pt hard of event", 1001, -0.5, 1000.5);
    TH1D *nLargeJets_h[nR];
    for (int iR = 0; iR < nR; iR++)
    {
        nLargeJets_h[iR] = new TH1D(Form("nLargeJets%d_h", int(Rvals[iR] * 10)), Form("Det and truth level jets with pT > 3*pthard of event, R=%0.1f", Rvals[iR]), 2, 0, 2);
    }

    TH1D *oneConstPDG_h = new TH1D("oneConstPDG_h", "PDG of 1-constituent jets, R=0.4", 12001, -6000.5, 6000.5);
    TH1D *twoConstPDG_h = new TH1D("twoConstPDG_h", "PDG of 2-constituent jets, R=0.4", 12001, -6000.5, 6000.5);
    TH1D *threeConstPDG_h = new TH1D("threeConstPDG_h", "PDG of 3-constituent jets, R=0.4", 12001, -6000.5, 6000.5);
    TH1D *totalConstPDG_h = new TH1D("totalConstPDG_h", "PDG of all jets' constituents, R=0.4", 12001, -6000.5, 6000.5);

    // Alice run loader
    AliRunLoader *fRunLoader = AliRunLoader::Open("galice.root");

    if (!fRunLoader->GetAliRun())
        fRunLoader->LoadgAlice();
    if (!fRunLoader->TreeE())
        fRunLoader->LoadHeader();
    if (!fRunLoader->TreeK())
        fRunLoader->LoadKinematics();

    TObjArray primTracks;

    std::vector<fastjet::JetDefinition> jet_def;
    std::vector<fastjet::JetDefinition> jet_defBkg;
    Double_t ghost_maxrap = 6.0;
    fastjet::GhostedAreaSpec area_spec(ghost_maxrap);
    fastjet::AreaDefinition area_def(fastjet::active_area, area_spec);

    for (Int_t iR = 0; iR < nR; iR++)
    {
        jet_def.push_back(fastjet::JetDefinition(fastjet::antikt_algorithm, Rvals[iR])); //, fastjet::pt_scheme)); //Etscheme, pTscheme, default is E scheme
        jet_defBkg.push_back(fastjet::JetDefinition(fastjet::kt_algorithm, Rvals[iR]));
    }

    Double_t evPtHard = 0;
    for (ievt = 0; ievt < fRunLoader->GetNumberOfEvents(); ievt++)
    {
        ievtTruth = ievt;
        // Get MC Stack
        Int_t ie = fRunLoader->GetEvent(ievt);

        if (ie != 0)
            continue;
        nevt_h->Fill(0.5); // count events

        AliStack *stack = fRunLoader->Stack();

        float Vertex[3] = {0x0};

        AliGenPythiaEventHeader *mcHeader = dynamic_cast<AliGenPythiaEventHeader *>(fRunLoader->GetHeader()->GenEventHeader());

        evPtHard = mcHeader->GetPtHard();
        nEvPtHard_h->Fill(evPtHard);

        std::vector<fastjet::PseudoJet> input_Particles;

        for (Int_t iTrk = 0; iTrk < stack->GetNtrack(); iTrk++)
        {
            TParticle *part = stack->Particle(iTrk);

            if (iTrk == 6)
            {
                Vertex[0] = part->Vx();
                Vertex[1] = part->Vy();
                Vertex[2] = part->Vz();
            }

            if (!stack->IsPhysicalPrimary(iTrk))
                continue;
            if (part->GetFirstDaughter() >= 0 && stack->IsPhysicalPrimary(part->GetFirstDaughter())) // if decay daughters are phys prim, only use those
                continue;

            if (part->Eta() > 5.8 || part->Eta() < 3.2)
                continue;

            input_Particles.push_back(fastjet::PseudoJet(part->Px(), part->Py(), part->Pz(), part->Energy()));
            input_Particles.back().set_user_index(iTrk);
        }

        // Get Clusters
        TTree *tClusters = 0;
        if (clusterFile->GetDirectory(Form("Event%i", ievt)))
        {
            clusterFile->GetDirectory(Form("Event%i", ievt))->GetObject("fTreeR", tClusters);
        }
        else
        {
            cout << "Cannot find event " << ievt << " in cluster file " << clusterFile->GetName() << endl;
            clusterFile->ls();
        }

        TBranch *bClusters;
        bClusters = tClusters->GetBranch("AliFOCALCluster"); // Branch for final ECAL clusters

        TClonesArray *clustersArray = 0;
        bClusters->SetAddress(&clustersArray);
        bClusters->GetEvent(0);

        std::vector<fastjet::PseudoJet> input_Clusters;

        for (int iClust = 0; iClust < clustersArray->GetEntries(); iClust++)
        {
            AliFOCALCluster *clust = (AliFOCALCluster *)clustersArray->At(iClust);

            TLorentzVector mom;
            Float_t vertex[3] = {0.};
            GetMomentum(mom, vertex, clust, 0.);
            input_Clusters.push_back(fastjet::PseudoJet(mom.Px(), mom.Py(), mom.Pz(), clust->E()));
            input_Clusters.back().set_user_index(1); // 1=isECAL
        }

        TBranch *AllClusters = tClusters->GetBranch("AliFOCALClusterItr");

        TClonesArray *FullclustersArray = 0;
        AllClusters->SetAddress(&FullclustersArray);
        AllClusters->GetEvent(0);

        // Selecting HCAL clusters
        for (int iClust = 0; iClust < FullclustersArray->GetEntries(); iClust++)
        {
            AliFOCALCluster *clust = (AliFOCALCluster *)FullclustersArray->At(iClust);
            if (clust->Segment() != 6)
            {
                continue;
            }

            if (clust->E() > 6500)
            {
                std::cout << "Rejecting fake cluster with large energy\n";
                continue;
            }

            TLorentzVector mom;
            Float_t vertex[3] = {0.};
            GetMomentum(mom, vertex, clust, 0. /*kPiPlusMass*/); // Don't forget to change it back
            input_Clusters.push_back(fastjet::PseudoJet(mom.Px(), mom.Py(), mom.Pz(), clust->E()));
            input_Clusters.back().set_user_index(0); // 0=isHCAL
        }

        for (Int_t iR = 0; iR < nR; iR++)
        {
            int isLargePt = 0; // If jet pthard = 3(2)*pthard of event or larger, this value is 1 and the event is skipped

            // These cuts are not used for now
            Double_t AreaCut = 0.6 * TMath::Pi() * TMath::Power(Rvals[iR], 2);
            fastjet::Selector selJet = SelectorMinArea(AreaCut) && fastjet::SelectorEMin(Ecut) && fastjet::SelectorEtaRange(3.5 + Rvals[iR], 5.5 - Rvals[iR]);

            fastjet::ClusterSequenceArea clust_seq(input_Clusters, jet_def[iR], area_def);
            fastjet::ClusterSequenceArea clust_seq_truth(input_Particles, jet_def[iR], area_def);
            double ptmin = 5.0; // GeV/c, might be too large in forward regions
            std::vector<fastjet::PseudoJet> inclusive_jets = sorted_by_pt(clust_seq.inclusive_jets(ptmin));
            std::vector<fastjet::PseudoJet> inclusive_jets_truth = sorted_by_pt(clust_seq_truth.inclusive_jets(ptmin));

            fastjet::ClusterSequenceArea clust_seqBkg(input_Clusters, jet_defBkg[iR], area_def);
            std::vector<fastjet::PseudoJet> inclusive_jetsBkg = sorted_by_pt(clust_seqBkg.inclusive_jets(0.));

            if (inclusive_jets.size() == 0)
                continue;

            Double_t Rho = 0.;

            // Perp Cone for underlying Event Subtraction
            Double_t PerpConePhi = inclusive_jets[0].phi() + TMath::Pi() / 2;
            PerpConePhi = (PerpConePhi > 2 * TMath::Pi()) ? PerpConePhi - 2 * TMath::Pi() : PerpConePhi; // fit to 0 < phi < 2pi

            Double_t PerpConeEta = inclusive_jets[0].eta();
            Double_t PerpConePt(0);

            for (unsigned int j = 0; j < input_Clusters.size(); j++)
            {

                Double_t deltaR(0);

                Double_t dPhi = TMath::Abs(input_Clusters[j].phi() - PerpConePhi);
                dPhi = (dPhi > TMath::Pi()) ? 2 * TMath::Pi() - dPhi : dPhi;

                Double_t dEta = TMath::Abs(input_Clusters[j].eta() - PerpConeEta);

                deltaR = TMath::Sqrt(TMath::Power(dEta, 2) + TMath::Power(dPhi, 2));

                if (deltaR <= Rvals[iR])
                    PerpConePt += input_Clusters[j].pt();
            }

            Double_t PerpConeRho = PerpConePt / (TMath::Pi() * TMath::Power(Rvals[iR], 2));
            Rho = PerpConeRho;

            // Additions 7.5.
            for (unsigned int i = 0; i < inclusive_jets.size(); i++)
            {
                Float_t bkgpT = inclusive_jets[i].perp() - (Rho * inclusive_jets[i].area());
                if (bkgpT >= 3.0 * evPtHard)
                {
                    isLargePt = 1;
                    nLargeJets_h[iR]->Fill(0.5);
                    break;
                }
            }
            // same thing for truth level jets, just to check
            for (unsigned int i = 0; i < inclusive_jets_truth.size(); i++)
            {
                Float_t trupT = inclusive_jets_truth[i].perp();
                if (trupT >= 3.0 * evPtHard)
                {
                    isLargePt = 1;
                    nLargeJets_h[iR]->Fill(1.5);
                    break;
                }
            }
            if (isLargePt == 1)
                continue;

            for (unsigned int i = 0; i < inclusive_jets_truth.size(); i++)
            {
                inclusive_jets_truth[i].set_user_index(i);
                inclusive_jets_truth[i].set_user_info(new JetMatchingParams());
                truthjetE = inclusive_jets_truth[i].E();
                truthjetpT = inclusive_jets_truth[i].perp();
                truthjetPhi = inclusive_jets_truth[i].phi();
                truthjetEta = inclusive_jets_truth[i].eta();
                truthjetRap = inclusive_jets_truth[i].rap();
                truthjetR = int(Rvals[iR] * 10);

                if (inclusive_jets_truth[i].has_constituents())
                {
                    jetParts_truth = inclusive_jets_truth[i].constituents().size();
                    float EnergyOfECALconstituents = 0.0;
                    float EnergyTotal = 0.0;
                    float truthjetECALEnergyFrac = -1.0;

                    std::vector<fastjet::PseudoJet> jetConstituents = inclusive_jets_truth[i].constituents();

                    for (int iH = 0; iH < jetConstituents.size(); iH++)
                    {

                        TParticle *particleA = ((TParticle *)(stack->Particle(jetConstituents[iH].user_index())));
                        EnergyTotal += particleA->Energy();
                        // std::cout << "Particle has E " << particleA->Energy() << "." << endl;
                        if (particleA->GetPDG()->Charge() == 0.0)
                        {
                            EnergyOfECALconstituents += (Float_t)particleA->Energy(); // cout << "Getting non-charged particle "<< iH << " for event " << ievt<< "." <<endl;
                        }
                    }
                    if (EnergyTotal == 0 || inclusive_jets_truth[i].eta() < 3.4 + Rvals[iR] || inclusive_jets_truth[i].eta() > 5.8 - Rvals[iR])
                    {
                        truthjetECALEnergyFrac = -1.0;
                        continue;
                    }
                    truthjetECALEnergyFrac = (Float_t)((Float_t)EnergyOfECALconstituents / (Float_t)EnergyTotal);
                    truthjetNeutralEnergyFrac = truthjetECALEnergyFrac;
                }

                results_truth->Fill();
            }

            for (unsigned int i = 0; i < inclusive_jets.size(); i++)
            {
                inclusive_jets[i].set_user_index(i);
                inclusive_jets[i].set_user_info(new JetMatchingParams());
            }

            DoJetMatching(inclusive_jets, inclusive_jets_truth);

            // adding ECAL fraction info to this code: Add a incl_jets_truth loop here, where pass thru all truth jets and assign ECAL fraction*100000 as user_index.

            for (unsigned int i = 0; i < inclusive_jets_truth.size(); i++)
            {
                if (inclusive_jets_truth[i].has_constituents())
                {

                    float EnergyOfECALconstituents = 0.0;
                    float EnergyTotal = 0.0;
                    float truthjetECALEnergyFrac = -1.0;

                    std::vector<fastjet::PseudoJet> jetConstituents = inclusive_jets_truth[i].constituents();

                    for (int iH = 0; iH < jetConstituents.size(); iH++)
                    {

                        TParticle *particleA = ((TParticle *)(stack->Particle(jetConstituents[iH].user_index())));
                        EnergyTotal += particleA->Energy();
                        // std::cout << "Particle has E " << particleA->Energy() << "." << endl;
                        if (particleA->GetPDG()->Charge() == 0.0)
                        {
                            EnergyOfECALconstituents += (Float_t)particleA->Energy(); // cout << "Getting non-charged particle "<< iH << " for event " << ievt<< "." <<endl;
                        }
                    }
                    if (EnergyTotal == 0 || inclusive_jets_truth[i].eta() < 3.4 + Rvals[iR] || inclusive_jets_truth[i].eta() > 5.8 - Rvals[iR])
                    {
                        truthjetECALEnergyFrac = -1.0;
                        continue;
                    }
                    truthjetECALEnergyFrac = (Float_t)((Float_t)EnergyOfECALconstituents / (Float_t)EnergyTotal);

                    Int_t ECALInt = (Int_t)(truthjetECALEnergyFrac * 100000);
                    inclusive_jets_truth[i].set_user_index(ECALInt);
                }
            }

            for (unsigned int i = 0; i < inclusive_jets.size(); i++)
            {
                jetE = inclusive_jets[i].E();
                jetpT = inclusive_jets[i].perp();
                jetpT_sub = inclusive_jets[i].perp() - (Rho * inclusive_jets[i].area());
                jetPhi = inclusive_jets[i].phi();
                jetEta = inclusive_jets[i].eta();
                jetRap = inclusive_jets[i].rap();
                jetR = int(Rvals[iR] * 10);

                if (inclusive_jets[i].has_constituents())
                {
                    jetParts = inclusive_jets[i].constituents().size();
                    // add det level ECAL fraction
                    std::vector<fastjet::PseudoJet> clustConstituents = inclusive_jets[i].constituents();
                    int IsHCAL; // 0 = is hcal
                    float EnergyOfECALconstituents = 0.0;
                    float EnergyTotal = 0.0;
                    for (int iH = 0; iH < clustConstituents.size(); iH++)
                    {
                        IsHCAL = clustConstituents[iH].user_index();
                        EnergyTotal += clustConstituents[iH].E();
                        if (IsHCAL == 1)
                        {
                            EnergyOfECALconstituents += clustConstituents[iH].E();
                        }
                    }
                    if (EnergyTotal > 0.0 && inclusive_jets[i].eta() >= 3.4 + Rvals[iR] && inclusive_jets[i].eta() <= 5.8 - Rvals[iR])
                    {
                        jetECALEnergyFrac = (Float_t)((Float_t)EnergyOfECALconstituents / (Float_t)EnergyTotal);
                    }
                    else
                    {
                        jetECALEnergyFrac = -1;
                    }
                }
                else
                {
                    jetParts = 0;
                }

                fastjet::PseudoJet *matchedJet;
                if (inclusive_jets[i].user_info<JetMatchingParams>().index1() >= 0)
                {
                    matchedJet = &(inclusive_jets_truth[inclusive_jets[i].user_info<JetMatchingParams>().index1()]);
                    // adding ECAL fraction part level matched: matchedjetECALEnergyFrac
                    matchedjetECALEnergyFrac = (Float_t)((Float_t)matchedJet->user_index() / 100000.0);
                }
                else
                {
                    // Assign an empty jets in case of no match
                    matchedJet = new fastjet::PseudoJet();
                    matchedjetECALEnergyFrac = -1.0;
                }

                jetE_match = matchedJet->E();
                jetpT_match = matchedJet->perp();
                jetPhi_match = matchedJet->phi();
                jetEta_match = matchedJet->eta();
                jetRap_match = matchedJet->rap();
                jet_distmatch = inclusive_jets[i].user_info<JetMatchingParams>().distance1();

                if (matchedJet->has_constituents())
                {
                    jetParts_match = matchedJet->constituents().size();

                    if (iR == 2)
                    {
                        std::vector<fastjet::PseudoJet> jetConstituents = matchedJet->constituents();

                        if (jetConstituents.size() == 1)
                        {
                            oneConstPDG_h->Fill(((TParticle *)(stack->Particle(jetConstituents[0].user_index())))->GetPdgCode());
                        }
                        if (jetConstituents.size() == 2)
                        {
                            twoConstPDG_h->Fill(((TParticle *)(stack->Particle(jetConstituents[0].user_index())))->GetPdgCode());
                            twoConstPDG_h->Fill(((TParticle *)(stack->Particle(jetConstituents[1].user_index())))->GetPdgCode());
                        }
                        if (jetConstituents.size() == 3)
                        {
                            threeConstPDG_h->Fill(((TParticle *)(stack->Particle(jetConstituents[0].user_index())))->GetPdgCode());
                            threeConstPDG_h->Fill(((TParticle *)(stack->Particle(jetConstituents[1].user_index())))->GetPdgCode());
                            threeConstPDG_h->Fill(((TParticle *)(stack->Particle(jetConstituents[2].user_index())))->GetPdgCode());
                        }
                        for (int j = 0; j < matchedJet->constituents().size(); j++)
                        {
                            totalConstPDG_h->Fill(((TParticle *)(stack->Particle(jetConstituents[j].user_index())))->GetPdgCode());
                        }
                    }
                }

                results->Fill();
            }
        }
    }

    fout->Write();
    fout->Close();
    clusterFile->Close();
    fRunLoader->Delete();
}

Float_t eta(TParticle *part)
{
    Double_t pt = sqrt(part->Px() * part->Px() + part->Py() * part->Py());
    if (pt == 0)
        return 9999;
    return -log(tan(TMath::ATan2(pt, part->Pz()) / 2));
}

Float_t phi(TParticle *part)
{
    return TMath::ATan2(part->Py(), -part->Px());
}

//_______________________________________________________________________
void GetMomentum(TLorentzVector &mom, const Float_t *vertex, AliFOCALCluster *foClust, Float_t mass)
{
    if (mass < 0)
        mass = 0;

    Double32_t energy = foClust->E();

    if (energy < mass)
        mass = 0;

    Double_t p = TMath::Sqrt(energy * energy - mass * mass);

    // std::cout<<"Total energy: "<<energy<<" HCAL Energy is: "<<foClust->GetSegmentEnergy(0)<<" "<<foClust->GetSegmentEnergy(1)<<" "<<foClust->GetSegmentEnergy(2)<<" "<<foClust->GetSegmentEnergy(3)<<" "<<foClust->GetSegmentEnergy(4)<<" "<<foClust->GetSegmentEnergy(5)<<" "<<foClust->GetSegmentEnergy(6)<<" "<< foClust->GetHCALEnergy()<<std::endl;
    Float_t pos[3] = {foClust->X(), foClust->Y(), foClust->Z()};

    if (vertex)
    { // calculate direction from vertex
        pos[0] -= vertex[0];
        pos[1] -= vertex[1];
        pos[2] -= vertex[2];
    }

    Double_t r = TMath::Sqrt(pos[0] * pos[0] + pos[1] * pos[1] + pos[2] * pos[2]);

    mom.SetPxPyPzE(p * pos[0] / r, p * pos[1] / r, p * pos[2] / r, energy);
}

void DoJetMatching(std::vector<fastjet::PseudoJet> &jetArray1, std::vector<fastjet::PseudoJet> &jetArray2)
{
    for (auto &jet1 : jetArray1)
    {
        JetMatchingParams matchingParams1(jet1.user_info<JetMatchingParams>());

        for (auto &jet2 : jetArray2)
        {
            JetMatchingParams matchingParams2(jet2.user_info<JetMatchingParams>());

            // This calculates the distance in the rap - phi
            double dist = jet1.delta_R(jet2);
            // In order to do calculation in eta-Phi uncomment this line
            // GetGeometricalMatchingLevel(jet1, jet2, dist);

            if (dist < matchingParams1.distance1())
            {
                matchingParams1.setJetIndexDistance2(matchingParams1.index1(), matchingParams1.distance1());
                matchingParams1.setJetIndexDistance1(jet2.user_index(), dist);
            }
            else if (dist < matchingParams1.distance2())
            {
                matchingParams1.setJetIndexDistance2(jet2.user_index(), dist);
            }

            if (dist < matchingParams2.distance1())
            {
                matchingParams2.setJetIndexDistance2(matchingParams2.index1(), matchingParams2.distance1());
                matchingParams2.setJetIndexDistance1(jet1.user_index(), dist);
            }
            else if (dist < matchingParams2.distance2())
            {
                matchingParams2.setJetIndexDistance2(jet1.user_index(), dist);
            }

            jet2.user_info_shared_ptr() = fastjet::SharedPtr<fastjet::PseudoJet::UserInfoBase>(new JetMatchingParams(matchingParams2));
        }

        jet1.user_info_shared_ptr() = fastjet::SharedPtr<fastjet::PseudoJet::UserInfoBase>(new JetMatchingParams(matchingParams1));
    }

    for (auto &jet1 : jetArray1)
    {
        fastjet::PseudoJet *jet2 = jet1.user_info<JetMatchingParams>().index1() >= 0 ? &(jetArray2[jet1.user_info<JetMatchingParams>().index1()]) : nullptr;
        fastjet::PseudoJet *jet3 = jet1.user_info<JetMatchingParams>().index2() >= 0 ? &(jetArray2[jet1.user_info<JetMatchingParams>().index2()]) : nullptr;

        if (jet2 && jet2->user_info<JetMatchingParams>().index1() != jet1.user_index())
        {
            if (!jet3 || jet3->user_info<JetMatchingParams>().index1() != jet1.user_index())
            {
                JetMatchingParams params;
                jet1.user_info_shared_ptr() = fastjet::SharedPtr<fastjet::PseudoJet::UserInfoBase>(new JetMatchingParams(params));
            }
            else
            {
                JetMatchingParams matchingParams(jet1.user_info<JetMatchingParams>());
                matchingParams.setJetIndexDistance1(jet1.user_info<JetMatchingParams>().index2(), jet1.user_info<JetMatchingParams>().distance2());
                matchingParams.setJetIndexDistance2(-1, 999);
                jet1.user_info_shared_ptr() = fastjet::SharedPtr<fastjet::PseudoJet::UserInfoBase>(new JetMatchingParams(matchingParams));
            }
        }
    }
}

void GetGeometricalMatchingLevel(fastjet::PseudoJet &jet1, fastjet::PseudoJet &jet2, Double_t &d)
{
    Double_t deta = jet2.eta() - jet1.eta();
    Double_t dphi = jet2.phi() - jet1.phi();
    dphi = TVector2::Phi_mpi_pi(dphi);
    d = TMath::Sqrt(deta * deta + dphi * dphi);
}
