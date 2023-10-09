#include "TFile.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TObjArray.h"
#include "TH1.h"
#include "TProfile.h"
#include "TH2.h"
#include "TH3.h"
#include "TParticle.h"
#include "TSystem.h"
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

//Int_t isFake;Int_t isFakeCut;
//Int_t isMatch;Int_t isMatchCut;
//Int_t isNoMatch;Int_t isNoMatchCut;

Float_t eta(TParticle *); // declaration, implementation at end of file
Float_t phi(TParticle *);
void GetMomentum(TLorentzVector &p, const Float_t *vertex, AliFOCALCluster *foClust, Float_t mass);

fastjet::PseudoJet GetMatchedJet(fastjet::PseudoJet jet1, std::vector<fastjet::PseudoJet> jetArray, Double_t &mindist);
void GetGeometricalMatchingLevel(fastjet::PseudoJet jet1, fastjet::PseudoJet jet2, Double_t &d);

Double_t SmearEMMomentum(Double_t mom, Double_t E);
Double_t SmearHadronMomentum(Double_t mom, Double_t E);
Double_t CalculateEnergy(Double_t dPx, Double_t dPy, Double_t dPz, bool isHadr);

void AnalyzeJets(Int_t startFolder, Int_t endFolder, Int_t startEvent, Int_t endEvent, const Char_t *simFolder, const Char_t *clustersFolder, const Char_t *dataSampleTag, const Char_t *outputdir)
{

	const Float_t Ecut = 2.0; //[GeV] cut on cluster energy for pi0 candidate pairs
	const Float_t kPi0Mass = 0.135;
	const Float_t kPiPlusMass = 0.139; //GeV

	const Int_t nR = 5;
	const Float_t Rvals[nR] = {0.2, 0.3, 0.4, 0.5, 0.6}; // Cone radii

	for (Int_t nfolder = startFolder; nfolder <= endFolder; nfolder++)
	{
		char filename[200];
		sprintf(filename, "%s/%03i/focalClusters.root", clustersFolder, nfolder);
		//
		//Due to Puhti things, we have to go through one file at a time, and the cluster files are each inside separate folders:
		//sprintf(filename, "%s/20230104pythia-gammajet-trig_ptmin2_%i-%i/clusters_%s_%i.root", clustersFolder,nfolder,nfolder, dataSampleTag, nfolder);
		//sprintf(filename, "%s/20230119_pythiamb-trig_ptmin-2_%i-%i/clusters_%s_%i.root", clustersFolder,nfolder,nfolder, dataSampleTag, nfolder);

		TFile *clusterFile = new TFile(filename);
		cout << "Clusters from: " << clusterFile->GetName() << endl;

		//TFile *fout = new TFile(Form("%s/analysis_%s_%d.root", outputdir, dataSampleTag, nfolder), "RECREATE");
		TFile *fout = new TFile(Form("%s/analysis_%d.root", outputdir, nfolder), "RECREATE");
		fout->cd();

		Float_t jetE, jetpT, jetpT_sub, jetPhi, jetEta, jetRap;
		Float_t jetE_match, jetpT_match, jetpT_sub_match, jetPhi_match, jetEta_match, jetRap_match;
		Float_t jet_distmatch;
		Int_t jetR;
		Int_t ievt;
		Int_t jetMatch;
		Int_t partPDG, jetParts, jetParts_match, jetParts_truth;
		//Float_t hcalX, hcalY, hcalZ;

		TTree *results = new TTree("jetTree", "jets Info Tree");
		results->Branch("nfolder", &nfolder, "nfolder/I");
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

		Float_t truthjetE, truthjetpT, truthjetPhi, truthjetEta, truthjetRap;
		Int_t truthjetR;

		TTree *results_truth = new TTree("truthjetTree", "MC jets Info Tree");
		results_truth->Branch("truthjetE", &truthjetE, "truthjetE/F");
		results_truth->Branch("truthjetpT", &truthjetpT, "truthjetpT/F");
		results_truth->Branch("truthjetPhi", &truthjetPhi, "truthjetPhi/F");
		results_truth->Branch("truthjetEta", &truthjetEta, "truthjetEta/F");
		results_truth->Branch("truthjetRap", &truthjetRap, "truthjetRap/F");
		results_truth->Branch("truthjetR", &truthjetR, "truthjetR/I");
		results_truth->Branch("jetParts_truth", &jetParts_truth, "jetParts_truth/I");

		TTree *in_PDG = new TTree("inPDGTree", "MC incoming particle PDG info");
		in_PDG->Branch("partPDG", &partPDG ,"partPDG/I");

		// cross section and trial information
		TProfile *xsec_h = new TProfile("xsec_h", "Cross section", 1, 0, 1);
		TH1D *trials_h = new TH1D("trials_h", "Trials", 1, 0, 1);
		TH1D *nevt_h = new TH1D("nevt_h", "Number of Events", 1, 0, 1);
		TH2F* hHitsXYHCAL = new TH2F("HitsXYHCAL", "X vs Y distribution of cluster hits in HCAL", 2200, -55.0, 55.0, 2200, -55.0, 55.0);
		TH1D *hHitsEnergyHCAL = new TH1D("hHitsEnergyHCAL", "Energy distr of HCAL clusters", 2200, -55.0, 55.0);

		AliRunLoader *fRunLoader = 0;
		cout << "FOLDER: " << nfolder << endl;

		TFile *fin_xsec = new TFile(Form("%s/%03d/pyxsec.root", simFolder, nfolder));

		if (fin_xsec->IsOpen())
		{
			TTree *xc_tree = (TTree *)fin_xsec->Get("Xsection");

			Double_t xsec;
			UInt_t ntrials;
			xc_tree->SetBranchAddress("xsection", &xsec);
			xc_tree->SetBranchAddress("ntrials", &ntrials);
			xc_tree->GetEvent(0);
			cout << "folder " << nfolder << " trials " << ntrials << " xsec " << xsec << endl;
			xsec_h->Fill(0.5, xsec, ntrials);
			trials_h->Fill(0.5, ntrials);
		}
		delete fin_xsec;

		sprintf(filename, "%s/%03i/%s", simFolder, nfolder, "galice.root");

		Long_t id = 0, size = 0, flags = 0, mt = 0;
		if (gSystem->GetPathInfo(filename, &id, &size, &flags, &mt) == 1)
		{
			cout << "ERROR: FOLDER: " << nfolder << endl;
			continue;
		}

		// Alice run loader
		fRunLoader = AliRunLoader::Open(filename);
		
		if (!fRunLoader)
		{
			cout << "ERROR: FOLDER: " << nfolder << endl;
			continue;
		}
		
		gSystem->Load("libAliPythia6.so");
		if (!fRunLoader->GetAliRun())
			fRunLoader->LoadgAlice();
		if (!fRunLoader->TreeE())
			fRunLoader->LoadHeader();
		if (!fRunLoader->TreeK())
			fRunLoader->LoadKinematics();
		vector<fastjet::JetDefinition> jet_def;
		vector<fastjet::JetDefinition> jet_defBkg;
		Double_t ghost_maxrap = 6.0;
		fastjet::GhostedAreaSpec area_spec(ghost_maxrap);
		fastjet::AreaDefinition area_def(fastjet::active_area, area_spec);

		for (Int_t iR = 0; iR < nR; iR++)
		{
			jet_def.push_back(fastjet::JetDefinition(fastjet::antikt_algorithm, Rvals[iR]));//, fastjet::pt_scheme)); //Default = E-scheme
			//jet_def.push_back(fastjet::JetDefinition(fastjet::cambridge_algorithm, Rvals[iR]));
			jet_defBkg.push_back(fastjet::JetDefinition(fastjet::kt_algorithm, Rvals[iR]));//, fastjet::pt_scheme));
		}

		Int_t maxEvent = TMath::Min(endEvent, fRunLoader->GetNumberOfEvents());
		for (ievt = startEvent; ievt <= maxEvent; ievt++)
		{
			Int_t ie = fRunLoader->GetEvent(ievt);
			cout << "Event: " << ievt << " folder " << nfolder << " event " << ievt << endl;

			if (ie != 0)
				continue;
			nevt_h->Fill(0.5); // count events

			// Get MC Stack
			AliStack *stack = fRunLoader->Stack();

			float Vertex[3] = {0x0};

			vector<fastjet::PseudoJet> input_Particles;
			//vector<fastjet::PseudoJet> input_Clusters;//smearing

			for (Int_t iTrk = 0; iTrk < stack->GetNtrack(); iTrk++)
			{
				TParticle *part = stack->Particle(iTrk);

				if (iTrk == 6)
				{
					Vertex[0] = part->Vx();
					Vertex[1] = part->Vy();
					Vertex[2] = part->Vz();
				}
				
				if(part->GetPdgCode() == 12 ||part->GetPdgCode() == 13 ||part->GetPdgCode() == 14 || part->GetPdgCode() == 16 || part->GetPdgCode() == 18)//not counting muons, neutrinos
					continue;
				if (!stack->IsPhysicalPrimary(iTrk))
					continue;
				if (part->GetFirstDaughter() >= 0 && stack->IsPhysicalPrimary(part->GetFirstDaughter())) // if decay daughters are phys prim, only use those
					continue;

				if (part->Eta() > 5.8 || part->Eta() < 3.2)
					continue;
				
				input_Particles.push_back(fastjet::PseudoJet(part->Px(), part->Py(), part->Pz(), part->Energy()));//the original input particles here
				//smearing
				Double_t dPxSmear, dPySmear, dPzSmear, dENew;
				bool bIsHadron = 1;
				if(part->GetPdgCode() == 11 || part->GetPdgCode() == 22) bIsHadron = 0;

	/*			dPxSmear = SmearEMMomentum(part->Px(), part->Energy());
				dPySmear = SmearEMMomentum(part->Py(), part->Energy());
				dPzSmear = SmearEMMomentum(part->Pz(), part->Energy());
				dENew    = CalculateEnergy(dPxSmear, dPySmear, dPzSmear, bIsHadron);
				
				cout << "Hadron?:" << bIsHadron << " , PDGcode " << part->GetPdgCode() << endl;
				cout << "px = " << part->Px() << " is now " << dPxSmear << endl;
				cout << "py = " << part->Py() << " is now " << dPySmear << endl;
				cout << "pz = " << part->Pz() << " is now " << dPzSmear << endl;
				cout << "E = " << part->Energy() << " is now " << dENew << endl;

				if(bIsHadron==0){ 
					//if(gRandom->Rndm()<0.05) continue;
					input_Clusters.push_back(fastjet::PseudoJet(dPxSmear, dPySmear, dPzSmear, dENew)); //smearing
				}
				if(bIsHadron==1){ 
					//if(gRandom->Rndm()<0.15) continue;
					input_Clusters.push_back(fastjet::PseudoJet(dPxSmear, dPySmear, dPzSmear, dENew)); //smearing
				}
	*/			
				partPDG=part->GetPdgCode();
				in_PDG->Fill();
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

			TBranch *bClusters = tClusters->GetBranch("AliFOCALCluster"); // Branch for final ECAL clusters

			TClonesArray *clustersArray = 0;
			bClusters->SetAddress(&clustersArray);
			bClusters->GetEvent(0);

			vector<fastjet::PseudoJet> input_Clusters;

			for (int iClust = 0; iClust < clustersArray->GetEntries(); iClust++)
			{
				AliFOCALCluster *clust = (AliFOCALCluster *)clustersArray->At(iClust);
				TLorentzVector mom;
				Float_t vertex[3] = {0.};
				GetMomentum(mom, vertex, clust, 0.);
				input_Clusters.push_back(fastjet::PseudoJet(mom.Px(), mom.Py(), mom.Pz(), clust->E()));
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

				hHitsXYHCAL->Fill(clust->X(),clust->Y());
				hHitsEnergyHCAL->Fill(clust->E());
				TLorentzVector mom;
				Float_t vertex[3] = {0.};
				GetMomentum(mom, vertex, clust, kPiPlusMass);
				input_Clusters.push_back(fastjet::PseudoJet(mom.Px(), mom.Py(), mom.Pz(), clust->E()));
			}
	
			for (Int_t iR = 0; iR < nR; iR++)
			{

				Double_t AreaCut = 0.6 * TMath::Pi() * TMath::Power(Rvals[iR], 2);

				fastjet::ClusterSequenceArea clust_seq(input_Clusters, jet_def[iR], area_def);
				fastjet::ClusterSequenceArea clust_seq_truth(input_Particles, jet_def[iR], area_def);
				double ptmin = 0.;
				vector<fastjet::PseudoJet> inclusive_jets = sorted_by_pt(clust_seq.inclusive_jets(ptmin));
				vector<fastjet::PseudoJet> inclusive_jets_truth = sorted_by_pt(clust_seq_truth.inclusive_jets(ptmin));

				fastjet::ClusterSequenceArea clust_seqBkg(input_Clusters, jet_defBkg[iR], area_def);
				vector<fastjet::PseudoJet> inclusive_jetsBkg = sorted_by_pt(clust_seqBkg.inclusive_jets(0.));

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
				
				for (unsigned int i = 0; i < inclusive_jets.size(); i++)
				{
					jetE = inclusive_jets[i].E();
					jetpT = inclusive_jets[i].perp();
					jetpT_sub = inclusive_jets[i].perp() - (Rho * inclusive_jets[i].area());
					jetPhi = inclusive_jets[i].phi();
					jetEta = inclusive_jets[i].eta();
					jetRap = inclusive_jets[i].rap();
					jetR = int(Rvals[iR] * 10);
					//jetParts = inclusive_jets[i].constituents().size();

					Double_t dist(-1.); 
					fastjet::PseudoJet matchedJet = GetMatchedJet(inclusive_jets[i], inclusive_jets_truth, dist);
					jetE_match = matchedJet.E();
					jetpT_match = matchedJet.perp();
					jetPhi_match = matchedJet.phi();
					jetEta_match = matchedJet.eta();
					jetRap_match = matchedJet.rap();
					jetParts_match = matchedJet.constituents().size();
					jet_distmatch = dist;

					results->Fill();
				}

				for (unsigned int i = 0; i < inclusive_jets_truth.size(); i++)
				{	
					truthjetE = inclusive_jets_truth[i].E();
					truthjetpT = inclusive_jets_truth[i].perp();
					truthjetPhi = inclusive_jets_truth[i].phi();
					truthjetEta = inclusive_jets_truth[i].eta();
					truthjetRap = inclusive_jets_truth[i].rap();
					truthjetR = int(Rvals[iR] * 10);
					
					jetParts_truth = inclusive_jets_truth[i].constituents().size();
					results_truth->Fill();
				}
			}
		}

		fout->Write();
		fout->Close();
		clusterFile->Close();
		fRunLoader->Delete();
	}
}


//We actually don't use these right now
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

fastjet::PseudoJet GetMatchedJet(fastjet::PseudoJet jet1, std::vector<fastjet::PseudoJet> jetArray, Double_t &mindist)
{
	mindist = 999; //should this be a more realistic distance value? Since now we are matching basically all jets to some jet, no matter how far away. 
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
/*
	//"removing" both jets if the matched, particle level jet has only one constituent. i.e. making the dist =999 but still returning the matched particle level jet
	if (jetArray[matchedIndex].constituents().size()==1){
		mindist = 999;
		return jetArray[matchedIndex];
	}
*/
	if (matchedIndex >= 0)
	{
		return jetArray[matchedIndex];
	}
	return fastjet::PseudoJet();
	
}

void GetGeometricalMatchingLevel(fastjet::PseudoJet jet1, fastjet::PseudoJet jet2, Double_t &d)
{
	Double_t deta = jet2.eta() - jet1.eta();
	//
	//HERE: changed from using eta to rapidity in calculating geometrical matching d
	Double_t drap = jet2.rap() - jet1.rap(); //should this be calculated or is it ok to use the fastjet rap()
	Double_t dphi = jet2.phi() - jet1.phi();
	dphi = TVector2::Phi_mpi_pi(dphi);
	d = TMath::Sqrt(deta * deta + dphi * dphi);
}

Double_t SmearEMMomentum(Double_t mom, Double_t E)
{
	Double_t dSmeared = 1./TMath::Sqrt(3)*TMath::Sqrt((0.27*0.27)/E + 0.01 * 0.01);//LOI
	//Double_t dSmeared = 1./TMath::Sqrt(3)*TMath::Sqrt((0.54*0.54)/E + 0.05 * 0.05);//Test beam
	//Double_t momSmeared = mom*gRandom->Gaus(1., dSmeared);
	//return momSmeared;
	Double_t dSigma = mom*dSmeared;
	Double_t pTrec  = gRandom->Gaus(mom, dSigma);
	return pTrec;
	//return mom;//no smearing
}

Double_t SmearHadronMomentum(Double_t mom, Double_t E)
{
	Double_t dSmeared = 1./TMath::Sqrt(3)*TMath::Sqrt((1.63*1.63)/E + 0.11 * 0.11); //LOI
	//Double_t dSmeared = 1./TMath::Sqrt(3)*TMath::Sqrt((2.50*2.50)/E + 0.12 * 0.12); //Test beam
        //Double_t momSmeared = mom*gRandom->Gaus(1., dSmeared);
	//return momSmeared;
	Double_t dSigma = mom*dSmeared;
        Double_t pTrec  = gRandom->Gaus(mom, dSigma);
        return pTrec;
	//return mom;//no smearing
}

Double_t CalculateEnergy(Double_t dPx, Double_t dPy, Double_t dPz, bool isHadr)
{
	Double_t dMass = 0.0;
	if (isHadr==1) dMass = 0.139;
	Double_t newEnergy = TMath::Sqrt(dMass * dMass + dPx*dPx + dPy*dPy + dPz*dPz);
	return newEnergy;
}


