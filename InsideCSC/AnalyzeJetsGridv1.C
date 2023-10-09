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
#include "AliFOCALCluster.h"

#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequenceArea.hh"
#ifndef __FJCORE__
#include "fastjet/GhostedAreaSpec.hh" // for area support
#endif								  // __FJCORE__

using std::cout;
using std::endl;

Float_t eta(TParticle *); // declaration, implementation at end of file
Float_t phi(TParticle *);

void GetMomentum(TLorentzVector &p, const Float_t *vertex, AliFOCALCluster *foClust, Float_t mass);
fastjet::PseudoJet GetMatchedJet(fastjet::PseudoJet jet1, std::vector<fastjet::PseudoJet> jetArray, Double_t &mindist, Float_t &ECALen);
void GetGeometricalMatchingLevel(fastjet::PseudoJet jet1, fastjet::PseudoJet jet2, Double_t &d);

void AnalyzeJetsGridv1(const char *inputdir = "./", const char *outputdir = "./")
{
	gSystem->Load("libpythia6_4_28.so");

	const Float_t Ecut = 2.0; //[GeV] cut on cluster energy for pi0 candidate pairs
	const Float_t kPi0Mass = 0.135;
	const Float_t kPiPlusMass = 0.139;

	const Int_t nR = 3;
	const Float_t Rvals[nR] = {0.2, 0.4, 0.6}; // Cone radii

	if (gSystem->AccessPathName(Form("%s/focalClusters.root", inputdir)))
	// if (gSystem->AccessPathName(Form("%s", inputdir)))
	{
		cout << "AnalyzeJetsGrid() ERROR: focalClusters.root not found!" << endl;
		return;
	}
	TFile *clusterFile = new TFile(Form("%s/focalClusters.root", inputdir));
	// TFile *clusterFile = new TFile(Form("%s", inputdir));
	cout << "Clusters from: " << clusterFile->GetName() << endl;

	TFile *fout = new TFile(Form("%s/analysisJets.root", outputdir), "RECREATE");
	fout->cd();

	Float_t jetE, jetpT, jetpT_sub, jetPhi, jetEta, jetRap;
	Float_t jetE_match, jetpT_match, jetpT_sub_match, jetPhi_match, jetEta_match, jetRap_match;
	Float_t jet_distmatch;
	Float_t jetHCALFrac, jetHCALEnergyFrac, jetECALEnergyFrac, matchedjetECALEnergyFrac;
	Int_t jetR;
	Int_t ievt;
	Int_t partPDG, jetParts, jetParts_match, jetParts_truth;

	TTree *results = new TTree("jetTree", "jets Info Tree");
	results->Branch("ievt", &ievt, "ievt/I");
	results->Branch("jetE", &jetE, "jetE/F");
	results->Branch("jetpT", &jetpT, "jetpT/F");
	results->Branch("jetpT_sub", &jetpT_sub, "jetpT_sub/F");
	results->Branch("jetPhi", &jetPhi, "jetPhi/F");
	results->Branch("jetEta", &jetEta, "jetEta/F");
	results->Branch("jetRap", &jetRap, "jetRap/F");
	results->Branch("jetParts", &jetParts, "jetParts/I");
	results->Branch("jetHCALFrac", &jetHCALFrac, "jetHCALFrac/F");
	results->Branch("jetHCALEnergyFrac", &jetHCALEnergyFrac, "jetHCALEnergyFrac/F");
	results->Branch("jetECALEnergyFrac", &jetECALEnergyFrac, "jetECALEnergyFrac/F");
	//results->Branch("truthjetNeutralEnergyFrac", &truthjetNeutralEnergyFrac, "truthjetNeutralEnergyFrac/F");
	//results->Branch("truthjetChargedEnergyFrac", &truthjetChargedEnergyFrac, "truthjetChargedEnergyFrac/F");
        results->Branch("matchedjetECALEnergyFrac", &matchedjetECALEnergyFrac, "matchedjetECALEnergyFrac/F");

	results->Branch("jetE_match", &jetE_match, "jetE_match/F");
	results->Branch("jetpT_match", &jetpT_match, "jetpT_match/F");
	results->Branch("jetPhi_match", &jetPhi_match, "jetPhi_match/F");
	results->Branch("jetEta_match", &jetEta_match, "jetEta_match/F");
	results->Branch("jetRap_match", &jetRap_match, "jetRap_match/F");
	results->Branch("jetParts_match", &jetParts_match, "jetParts_match/I");
	results->Branch("jet_distmatch", &jet_distmatch, "jet_distmatch/F");
	results->Branch("jetR", &jetR, "jetR/I");

	Float_t truthjetE, truthjetpT, truthjetPhi, truthjetEta, truthjetRap;
	Int_t truthjetR;
	Float_t truthjetHCALEnergyFrac, truthjetECALEnergyFrac;
	Float_t truthjetChargedEnergyFrac, truthjetNeutralEnergyFrac;
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
	results_truth->Branch("truthjetHCALEnergyFrac", &truthjetHCALEnergyFrac, "truthjetHCALEnergyFrac/F");
	results_truth->Branch("truthjetECALEnergyFrac", &truthjetECALEnergyFrac, "truthjetECALEnergyFrac/F");
	results_truth->Branch("truthjetNeutralEnergyFrac", &truthjetNeutralEnergyFrac, "truthjetNeutralEnergyFrac/F");
        results_truth->Branch("truthjetChargedEnergyFrac", &truthjetChargedEnergyFrac, "truthjetChargedEnergyFrac/F");

	Float_t partRap, partEta;

	TTree *in_PDG = new TTree("inPDGTree", "MC incoming particle PDG info");
	in_PDG->Branch("partPDG", &partPDG, "partPDG/I");
	in_PDG->Branch("partRap", &partRap, "partRap/F");
	in_PDG->Branch("partEta", &partEta, "partEta/F");

	TH1D *nevt_h = new TH1D("nevt_h", "Number of Events", 1, 0, 1);
/*
	TH1F *HCALClustersPerJet_h = new TH1F("HCALClustersPerJet_h", "Fraction of HCAL clusters out of all clusters in jet", 201, -0.05, 1.05);
	TH1F *HCALClusterEnergyPerJet_h = new TH1F("HCALClusterEnergyPerJet_h", "Fraction of HCAL cluster energy out of total cluster energy sum in jet", 201, -0.05, 1.05);
	TH2F *TotalEnergyDiff_h = new TH2F("TotalEnergyDiff_h", "Energy difference of particle and detector input, 3.5<eta<5.5 (particle level)", 400, -0.5, 1999.5, 400, -999.5, 999.5);
  TH2F *TotalEnergyDiffECAL_h = new TH2F("TotalEnergyDiffECAL_h", "Energy difference of em particle and detector input, 3.5<eta<5.5 (particle level)", 400, -0.5, 1999.5, 400, -999.5, 999.5);
  TH2F *TotalEnergyDiffHCAL_h = new TH2F("TotalEnergyDiffHCAL_h", "Energy difference of hadron particle and detector input, 3.5<eta<5.5 (particle level)", 400, -0.5, 1999.5, 400, -999.5, 999.5);
  TH2F *TotalDeltaEnergyDiff_h = new TH2F("TotalDeltaEnergyDiff_h", "Energy difference of particle and detector input, 3.5<eta<5.5 (particle level)", 400, -0.5, 1999.5, 1000, -99.5, 99.5);
	TH2F *TotalDeltaEnergyDiffECAL_h = new TH2F("TotalDeltaEnergyDiffECAL_h", "Energy difference of em particle and detector input, 3.5<eta<5.5 (particle level)", 400, -0.5, 1999.5, 1000, -99.5, 99.5);
	TH2F *TotalDeltaEnergyDiffHCAL_h = new TH2F("TotalDeltaEnergyDiffHCAL_h", "Energy difference of hadron particle and detector input, 3.5<eta<5.5 (particle level)", 400, -0.5, 1999.5, 1000, -99.5, 99.5);
	TH1F *HCALClustersPerEvent_h = new TH1F("HCALClustersPerEvent_h", "Fraction of HCAL clusters out of all clusters in event", 201, -0.05, 1.05);
	TH1F *HCALClusterEnergyPerEvent_h = new TH1F("HCALClusterEnergyPerEvent_h", "Fraction of HCAL cluster energy out of total cluster energy sum in event", 201, -0.05, 1.05);

	TH1F *TotalEnergyDiffECAL1d_h = new TH1F("TotalEnergyDiffECAL1d_h", "Photons and electrons, (E_det-E_part)/E_part, 3.4<eta<5.8 (particle level)", 400, -1, 1);
	TH1F *TotalEnergyDiffHCAL1d_h = new TH1F("TotalEnergyDiffHCAL1d_h", "Hadrons, (E_det-E_part)/E_part, 3.4<eta<5.8 (particle level)", 400, -1, 1);
	TH1F *TotalEnergyDiff1d_h = new TH1F("TotalEnergyDiff1d_h", "All particles, (E_det-E_part)/E_part, 3.4<eta<5.8 (particle level)", 400, -1, 1);

	TH1F *hdpT1D 	= new TH1F("hdpT1D","hdpT1D",1000,-10, 10);
	TH2F *hpT2Dpp 	= new TH2F("hpT2Dpp","hpT2Dpp",1000, 0, 100, 1000,-100,100);
	TH2F *hpT2Dcl 	= new TH2F("hpT2Dcl","hpT2Dcl", 1000, 0, 100, 1000, -100, 100);
	TH2F *hdpT2Dpp	= new TH2F("hdpT2Dpp","hdpT2Dpp",1000, 0, 100, 1000,-10,10);
	TH2F *hdpT2Dcl  = new TH2F("hdpT2Dcl","hdpT2Dcl",1000, 0, 100, 1000,-10,10);

	TH1F *hdpT1DHCAL= new TH1F("hdpT1DHCAL","hdpT1DHCAL",1000,-10, 10);
	TH1F *hdpT1DECAL= new TH1F("hdpT1DECAL","hdpT1DECAL",1000,-10, 10);
	TH2F *hpT2DppHCAL   = new TH2F("hpT2DppHCAL","hpT2DppHCAL",1000, 0, 100, 1000,-100,100);
	TH2F *hpT2DppECAL   = new TH2F("hpT2DppECAL","hpT2DppECAL",1000, 0, 100, 1000,-100,100);
	TH2F *hpT2DclHCAL   = new TH2F("hpT2DclHCAL","hpT2DclHCAL",1000, 0, 100, 1000,-100,100);
  TH2F *hpT2DclECAL   = new TH2F("hpT2DclECAL","hpT2DclECAL",1000, 0, 100, 1000,-100,100);
	TH2F *hdpT2DppHCAL   = new TH2F("hdpT2DppHCAL","hdpT2DppHCAL",1000, 0, 100, 1000,-100,100);
  TH2F *hdpT2DppECAL   = new TH2F("hdpT2DppECAL","hdpT2DppECAL",1000, 0, 100, 1000,-100,100);
  TH2F *hdpT2DclHCAL   = new TH2F("hdpT2DclHCAL","hdpT2DclHCAL",1000, 0, 100, 1000,-100,100);
  TH2F *hdpT2DclECAL   = new TH2F("hdpT2DclECAL","hdpT2DclECAL",1000, 0, 100, 1000,-100,100);
*/

TH2F *ECALfrac2D_part_h[nR], *HCALfrac2D_part_h[nR], *Chargedfrac2D_part_h[nR], *Neutralfrac2D_part_h[nR], *ECALfrac2D_det_h[nR], *HCALfrac2D_det_h[nR];
  //histograms for the energy fraction checks
  for(int iR=0; iR<nR; iR++){
		ECALfrac2D_part_h[iR] = new TH2F(Form("ECALfrac2D_part_%d_h",int(Rvals[iR]*10)),Form("Energy fraction of photon and electron energy of total jet energy vs. total jet energy, particle level, R=%0.1f", Rvals[iR]),2001, -0.5, 4000, 102, -0.05, 1.05);
		HCALfrac2D_part_h[iR] = new TH2F(Form("HCALfrac2D_part_%d_h",int(Rvals[iR]*10)),Form("Energy fraction of hadron energy of total jet energy vs. total jet energy, particle level, R=%0.1f",Rvals[iR]),2001, -0.5, 4000, 102, -0.05, 1.05);
		Chargedfrac2D_part_h[iR] = new TH2F(Form("Chargedfrac2D_part_%d_h",int(Rvals[iR]*10)),Form("Energy fraction of charged particle energy of total jet energy vs. total jet energy, particle level, R=%0.1f",Rvals[iR]),2001, -0.5, 4000, 102, -0.05, 1.05);
		Neutralfrac2D_part_h[iR] = new TH2F(Form("Neutralfrac2D_part_%d_h",int(Rvals[iR]*10)),Form("Energy fraction of neutral particle energy of total jet energy vs. total jet energy, particle level, R=%0.1f",Rvals[iR]),2001, -0.5, 4000, 102, -0.05, 1.05);
		ECALfrac2D_det_h[iR] = new TH2F(Form("ECALfrac2D_det_%d_h",int(Rvals[iR]*10)),Form("Energy fraction of ECAL cluster energy of total jet energy vs. total jet energy, detector level, R=%0.1f",Rvals[iR]),2001, -0.5, 4000, 102, -0.05, 1.05);
		HCALfrac2D_det_h[iR] = new TH2F(Form("HCALfrac2D_det_%d_h",int(Rvals[iR]*10)),Form("Energy fraction of HCAL energy of total jet energy vs. total jet energy, detector level, R=%0.1f", Rvals[iR]),2001, -0.5, 4000, 102, -0.05, 1.05);
  }
  //End of energy fraction histograms


	Float_t clusterEnergySum, particleEnergySum, clusterEnergySumHCAL, clusterEnergySumECAL, particleEnergySumHCAL, particleEnergySumECAL;
	Float_t clusterpTSum, particlepTSum, clusterpTSumHCAL, clusterpTSumECAL, particlepTSumHCAL, particlepTSumECAL;

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
		jet_def.push_back(fastjet::JetDefinition(fastjet::antikt_algorithm, Rvals[iR])); //, fastjet::pt_scheme));
		jet_defBkg.push_back(fastjet::JetDefinition(fastjet::kt_algorithm, Rvals[iR]));
	}

	for (ievt = 0; ievt <= fRunLoader->GetNumberOfEvents(); ievt++)
	{
		ievtTruth=ievt;
		particleEnergySum = 0.0;particlepTSum=0.0;
		clusterEnergySum = 0.0;clusterpTSum=0.0;
		clusterEnergySumECAL = 0.0;clusterpTSumECAL=0.0;
		particleEnergySumHCAL = 0.0;particlepTSumHCAL=0.0;
		particleEnergySumECAL = 0.0;particlepTSumECAL=0.0;
		// Get MC Stack
		Int_t ie = fRunLoader->GetEvent(ievt);

		if (ie != 0)
			continue;
		nevt_h->Fill(0.5); // count events

		AliStack *stack = fRunLoader->Stack();

		float Vertex[3] = {0x0};

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
			// check if removing muons neutrinos makes a difference :
			if (abs(part->GetPdgCode()) == 12 || abs(part->GetPdgCode()) == 13 || abs(part->GetPdgCode()) == 14 || abs(part->GetPdgCode()) == 16 || abs(part->GetPdgCode()) == 18) // not counting muons, neutrinos
				continue;
			if (!stack->IsPhysicalPrimary(iTrk))
				continue;
			if (part->GetFirstDaughter() >= 0 && stack->IsPhysicalPrimary(part->GetFirstDaughter())) // if decay daughters are phys prim, only use those
				continue;

			if (part->Eta() > 5.8 || part->Eta() < 3.4) // changed from 3.2//5.8 3.4
				continue;

			input_Particles.push_back(fastjet::PseudoJet(part->Px(), part->Py(), part->Pz(), part->Energy()));
			input_Particles.back().set_user_index(iTrk);
			partPDG = part->GetPdgCode();
			partEta = part->Eta();
			partRap = part->Y();
			in_PDG->Fill();

			bool bIsHadron = 1;
			if (abs(part->GetPdgCode()) == 11 || part->GetPdgCode() == 22)
				bIsHadron = 0;
			if (bIsHadron == 0)
			{
				particleEnergySumECAL += part->Energy();
				particlepTSumECAL += part->Pt();
			}
			if (bIsHadron == 1)
			{
				particleEnergySumHCAL += part->Energy();
				particlepTSumHCAL += part->Pt();
			}
			particleEnergySum += part->Energy();
			particlepTSum += part->Pt();
		}

		// Get Clusters
		int clustersHCALTotal = 0;
		int clustersTotal = 0;
		Float_t clustersHCALTotalEnergy = 0.0; Float_t clustersHCALTotalpT=0.0;

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
			
			Float_t trans = TMath::Sqrt(clust->X()*clust->X() + clust->Y()*clust->Y());
                        Float_t theta = TMath::Pi()/2.;
                        if (clust->Z() != 0) {
                          theta = TMath::ATan(trans/clust->Z());
                        }
                        double etaclust = 1e6;
                        if (theta > 1e-6)
                                  etaclust = -TMath::Log(TMath::Tan(theta/2.));

                        if(etaclust < 3.4 || etaclust > 5.8)
                                continue;

			clusterEnergySumECAL += clust->E(); clusterpTSumECAL += mom.Pt();
			clusterEnergySum += clust->E(); clusterpTSum += mom.Pt();
			clustersTotal++;
      input_Clusters.push_back(fastjet::PseudoJet(mom.Px(), mom.Py(), mom.Pz(), clust->E()));
      input_Clusters.back().set_user_index(1); //(iClust);
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
			GetMomentum(mom, vertex, clust, kPiPlusMass);

                        Float_t trans = TMath::Sqrt(clust->X()*clust->X() + clust->Y()*clust->Y());
                        Float_t theta = TMath::Pi()/2.;
                        if (clust->Z() != 0) {
                          theta = TMath::ATan(trans/clust->Z());
                        }
                        double etaclust = 1e6;
                        if (theta > 1e-6)
                                  etaclust = -TMath::Log(TMath::Tan(theta/2.));

                        if(etaclust < 3.4 || etaclust > 5.8)
                                continue;

			clusterEnergySum += clust->E(); clusterpTSum += mom.Pt();
			clustersHCALTotalEnergy += clust->E(); clustersHCALTotalpT += mom.Pt();
			clustersHCALTotal++;
			clustersTotal++;

			input_Clusters.push_back(fastjet::PseudoJet(mom.Px(), mom.Py(), mom.Pz(), clust->E()));
      input_Clusters.back().set_user_index(0); //(iClust);
		}
/*
		if(particlepTSum!=0)		hdpT1D->Fill((particlepTSum-clusterpTSum)/particlepTSum);
		if(particlepTSumHCAL!=0)	hdpT1DHCAL->Fill((particlepTSumHCAL-clustersHCALTotalpT)/particlepTSumHCAL);
		if(particlepTSumECAL!=0)	hdpT1DECAL->Fill((particlepTSumECAL-clusterpTSumECAL)/particlepTSumECAL);
		hpT2Dpp->Fill(particlepTSum,particlepTSum-clusterpTSum);
		hpT2DppHCAL->Fill(particlepTSumHCAL,particlepTSumHCAL-clustersHCALTotalpT);
		hpT2DppECAL->Fill(particlepTSumECAL,particlepTSumECAL-clusterpTSumECAL);
		hpT2Dcl->Fill(clusterpTSum,particlepTSum-clusterpTSum);
    hpT2DclHCAL->Fill(clustersHCALTotalpT,particlepTSumHCAL-clustersHCALTotalpT);
    hpT2DclECAL->Fill(clusterpTSumECAL,particlepTSumECAL-clusterpTSumECAL);

		hdpT2Dpp->Fill(particlepTSum,(particlepTSum-clusterpTSum)/particlepTSum);
    hdpT2DppHCAL->Fill(particlepTSumHCAL,(particlepTSumHCAL-clustersHCALTotalpT)/particlepTSumHCAL);
    hdpT2DppECAL->Fill(particlepTSumECAL,(particlepTSumECAL-clusterpTSumECAL)/particlepTSumECAL);
    hdpT2Dcl->Fill(clusterpTSum,(particlepTSum-clusterpTSum)/clusterpTSum);
    hdpT2DclHCAL->Fill(clustersHCALTotalpT,(particlepTSumHCAL-clustersHCALTotalpT)/clustersHCALTotalpT);
    hdpT2DclECAL->Fill(clusterpTSumECAL,(particlepTSumECAL-clusterpTSumECAL)/clusterpTSumECAL);



		HCALClustersPerEvent_h->Fill((Float_t)((Float_t)clustersHCALTotal / (Float_t)clustersTotal));
		HCALClusterEnergyPerEvent_h->Fill(clustersHCALTotalEnergy / clusterEnergySum);
		TotalEnergyDiff_h->Fill(particleEnergySum, clusterEnergySum - particleEnergySum);
		TotalEnergyDiffECAL_h->Fill(particleEnergySumECAL, clusterEnergySumECAL - particleEnergySumECAL);
		TotalEnergyDiffHCAL_h->Fill(particleEnergySumHCAL, clustersHCALTotalEnergy - particleEnergySumHCAL);
		if (particleEnergySumECAL != 0){
			TotalEnergyDiffECAL1d_h->Fill((clusterEnergySumECAL - particleEnergySumECAL) / particleEnergySumECAL);
			TotalDeltaEnergyDiffECAL_h->Fill(particleEnergySumECAL, (clusterEnergySumECAL - particleEnergySumECAL)/particleEnergySumECAL);
		}
		if (particleEnergySumHCAL != 0){
			TotalEnergyDiffHCAL1d_h->Fill((clustersHCALTotalEnergy - particleEnergySumHCAL) / particleEnergySumHCAL);
			TotalDeltaEnergyDiffHCAL_h->Fill(particleEnergySumHCAL, (clustersHCALTotalEnergy - particleEnergySumHCAL)/particleEnergySumHCAL);
		}
		if (particleEnergySum != 0){
			TotalEnergyDiff1d_h->Fill((clusterEnergySum - particleEnergySum) / particleEnergySum);
			TotalDeltaEnergyDiff_h->Fill(particleEnergySum, (clusterEnergySum - particleEnergySum)/particleEnergySum);
		}
*/
		for (Int_t iR = 0; iR < nR; iR++)
		{

			Double_t AreaCut = 0.6 * TMath::Pi() * TMath::Power(Rvals[iR], 2);

			fastjet::ClusterSequenceArea clust_seq(input_Clusters, jet_def[iR], area_def);
			fastjet::ClusterSequenceArea clust_seq_truth(input_Particles, jet_def[iR], area_def);
			double ptmin = 5.0; //0.
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

			// Truth jets loop
			for (unsigned int i = 0; i < inclusive_jets_truth.size(); i++)
			{
				//inclusive_jets_truth[i].set_user_index(i);
				truthjetE = inclusive_jets_truth[i].E();
				truthjetpT = inclusive_jets_truth[i].perp();
				truthjetPhi = inclusive_jets_truth[i].phi();
				truthjetEta = inclusive_jets_truth[i].eta();
				truthjetRap = inclusive_jets_truth[i].rap();
				truthjetR = int(Rvals[iR] * 10);

				if (!inclusive_jets_truth[i].has_constituents()) jetParts_truth=0;
				if (inclusive_jets_truth[i].has_constituents()){

					jetParts_truth = inclusive_jets_truth[i].constituents().size();

					float EnergyOfHCALconstituents = 0.0;float EnergyOfECALconstituents = 0.0;
					float EnergyOfNeutralconstituents = 0.0;float EnergyOfChargedconstituents = 0.0;
					float EnergyTotal = 0.0;

					std::vector<fastjet::PseudoJet> jetConstituents = inclusive_jets_truth[i].constituents();

					for (int iH = 0; iH < jetConstituents.size(); iH++){

						TParticle *particleA = ((TParticle *)(stack->Particle(jetConstituents[iH].user_index())));
						EnergyTotal += particleA->Energy();
						cout << "Particle has E "<< particleA->Energy() << "." <<endl; 
						if(abs(particleA->GetPdgCode()) != 11 && particleA->GetPdgCode() != 22){ EnergyOfHCALconstituents += (Float_t)particleA->Energy();//cout << "Getting hcal particle "<< iH << " for event " << ievt<< "." <<endl; 
						}
						if(abs(particleA->GetPdgCode()) == 11 || particleA->GetPdgCode() == 22){ EnergyOfECALconstituents += (Float_t)particleA->Energy();//cout << "Getting ecal particle "<< iH << " for event " << ievt<< "." <<endl; 
						}
						if(particleA->GetPDG()->Charge()==0.0) {EnergyOfNeutralconstituents += (Float_t)particleA->Energy();//cout << "Getting non-charged particle "<< iH << " for event " << ievt<< "." <<endl; 
						}
						if(particleA->GetPDG()->Charge()!=0.0) {EnergyOfChargedconstituents += (Float_t)particleA->Energy();//cout << "Getting charged particle "<< iH << " for event " << ievt<< "." <<endl; 
						}
					
					}

					//if(EnergyTotal<0.000005 || inclusive_jets_truth[i].E()<0.000005 || inclusive_jets_truth[i].eta()<3.4+Rvals[iR] || inclusive_jets_truth[i].eta()>5.8-Rvals[iR]){ 
					if(EnergyTotal==0 || inclusive_jets_truth[i].eta()<3.4+Rvals[iR] || inclusive_jets_truth[i].eta()>5.8-Rvals[iR]){ 
						truthjetHCALEnergyFrac = -1.0;truthjetChargedEnergyFrac = -1.0;
						truthjetECALEnergyFrac = -1.0;truthjetNeutralEnergyFrac = -1.0;
						continue;
					}
					truthjetHCALEnergyFrac = (Float_t)((Float_t)EnergyOfHCALconstituents/(Float_t)EnergyTotal);
					truthjetECALEnergyFrac = (Float_t)((Float_t)EnergyOfECALconstituents/(Float_t)EnergyTotal);
					truthjetNeutralEnergyFrac = (Float_t)((Float_t)EnergyOfNeutralconstituents/(Float_t)EnergyTotal);
					truthjetChargedEnergyFrac = (Float_t)((Float_t)EnergyOfChargedconstituents/(Float_t)EnergyTotal);
					
					Int_t ECALInt = (Int_t)(truthjetNeutralEnergyFrac*100000);
					inclusive_jets_truth[i].set_user_index(ECALInt);
//						cout << "Filling histos: "<< (Float_t)EnergyOfECALconstituents/(Float_t)inclusive_jets_truth[i].E() <<"?="<< (Float_t)EnergyOfChargedconstituents/(Float_t)inclusive_jets_truth[i].E() << " and " <<  EnergyOfHCALconstituents/(Float_t)inclusive_jets_truth[i].E() << ", sum is " <<EnergyTotal << "=?" <<(Float_t)inclusive_jets_truth[i].E() <<" ." <<endl; 
					ECALfrac2D_part_h[iR]->Fill((Float_t)inclusive_jets_truth[i].E(),truthjetECALEnergyFrac);
					HCALfrac2D_part_h[iR]->Fill((Float_t)inclusive_jets_truth[i].E(),truthjetHCALEnergyFrac);
					Chargedfrac2D_part_h[iR]->Fill((Float_t)inclusive_jets_truth[i].E(),truthjetChargedEnergyFrac);
					Neutralfrac2D_part_h[iR]->Fill((Float_t)inclusive_jets_truth[i].E(),truthjetNeutralEnergyFrac);

				}
				results_truth->Fill();
			}

			// Detector level jets loop + matching
			for (unsigned int i = 0; i < inclusive_jets.size(); i++)
			{
				inclusive_jets[i].set_user_index(i);
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

					std::vector<fastjet::PseudoJet> clustConstituents = inclusive_jets[i].constituents();
					int IsHCAL; // 0 = is hcal
					int nOfHCALconstituents = 0;
					float EnergyOfHCALconstituents = 0.0;float EnergyOfECALconstituents = 0.0;
					float EnergyTotal = 0.0;

					for (int iH = 0; iH < clustConstituents.size(); iH++)
					{
						IsHCAL = clustConstituents[iH].user_index();
						EnergyTotal += clustConstituents[iH].E();
						if (IsHCAL == 0){
							nOfHCALconstituents++;
							EnergyOfHCALconstituents += clustConstituents[iH].E();
						}
						if (IsHCAL == 1){
							EnergyOfECALconstituents += clustConstituents[iH].E();
						}
//						cout << "Getting "<< IsHCAL <<" particle "<< iH << " for event " << ievt<< "." <<endl; 
//						cout << "Particle has E "<< clustConstituents[iH].E() << "." <<endl; 

					}
//					HCALClustersPerJet_h->Fill((Float_t)((Float_t)nOfHCALconstituents / (Float_t)clustConstituents.size()));

//					if (EnergyTotal > 0)
//						HCALClusterEnergyPerJet_h->Fill((Float_t)EnergyOfHCALconstituents / (Float_t)EnergyTotal);
					// cout << "HCAL/ALL = " << (Float_t)((Float_t)nOfHCALconstituents/(Float_t)clustConstituents.size()) << endl;
					// cout << "HCAL E/ ALL E = " << EnergyOfHCALconstituents << "/" << EnergyTotal << endl;

					jetHCALFrac = (Float_t)nOfHCALconstituents / (Float_t)clustConstituents.size();
					if (EnergyTotal > 0.0 && inclusive_jets[i].eta()>3.4+Rvals[iR] && inclusive_jets[i].eta()<5.8-Rvals[iR])
					{
//						cout << "Filling histos: "<< (Float_t)EnergyOfECALconstituents/(Float_t)inclusive_jets[i].E() << " and " <<  (Float_t)EnergyOfHCALconstituents/(Float_t)inclusive_jets[i].E() << ", sum is " <<EnergyTotal << "=?" <<(Float_t)inclusive_jets[i].E() <<" ." <<endl; 
						jetHCALEnergyFrac = (Float_t)((Float_t)EnergyOfHCALconstituents / (Float_t)EnergyTotal);
						jetECALEnergyFrac = (Float_t)((Float_t)EnergyOfECALconstituents / (Float_t)EnergyTotal);
						HCALfrac2D_det_h[iR]->Fill((Float_t)inclusive_jets[i].E(), jetHCALEnergyFrac);
						ECALfrac2D_det_h[iR]->Fill((Float_t)inclusive_jets[i].E(), jetECALEnergyFrac);
					}
					else
					{
						jetHCALEnergyFrac = -1;jetECALEnergyFrac = -1;
					}
				}

				if (!inclusive_jets[i].has_constituents()) jetParts = 0;

				Double_t dist(-1.); Float_t ECALenergym = -1.0;
				fastjet::PseudoJet matchedJet = GetMatchedJet(inclusive_jets[i], inclusive_jets_truth, dist, ECALenergym);

				jetE_match = matchedJet.E();
				jetpT_match = matchedJet.perp();
				jetPhi_match = matchedJet.phi();
				jetEta_match = matchedJet.eta();
				jetRap_match = matchedJet.rap();
				// if(matchedJet.has_constituents())jetParts_match = matchedJet.constituents().size();
				jet_distmatch = dist;
				matchedjetECALEnergyFrac = ECALenergym;

				if (matchedJet.has_constituents())
				{
					jetParts_match = matchedJet.constituents().size();
				}
				else {jetParts_match=0;}

				results->Fill();
			}
		}
	}

				cout <<"writing" << endl;
	fout->Write();
				cout <<"closing" << endl;
	fout->Close();
				cout <<"closing cluster file" << endl;
	clusterFile->Close();
				cout <<"ALL DONE" << endl;
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

fastjet::PseudoJet GetMatchedJet(fastjet::PseudoJet jet1, std::vector<fastjet::PseudoJet> jetArray, Double_t &mindist, Float_t &ECALen)
{

	mindist = 999;
	Int_t matchedIndex = -1;

	for (unsigned int ijet = 0; ijet < jetArray.size(); ijet++)
	{
		Double_t dist(0.);
		GetGeometricalMatchingLevel(jet1, jetArray[ijet], dist);
		if (dist < mindist)
		{
			mindist = dist;
			matchedIndex = ijet;
		}
	}

	if (matchedIndex >= 0)
	{
		ECALen = (Float_t)(jetArray[matchedIndex].user_index()/100000.0);
		return jetArray[matchedIndex];
	}
	else
	{
		return fastjet::PseudoJet();
	}
}

void GetGeometricalMatchingLevel(fastjet::PseudoJet jet1, fastjet::PseudoJet jet2, Double_t &d)
{
	Double_t deta = jet2.eta() - jet1.eta();
	Double_t drap = jet2.rap() - jet1.rap();
	Double_t dphi = jet2.phi() - jet1.phi();
	dphi = TVector2::Phi_mpi_pi(dphi);
	d = TMath::Sqrt(drap * drap + dphi * dphi);
}
