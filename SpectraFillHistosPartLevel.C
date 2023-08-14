#include <TH1F.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TTree.h>

//Fills particle level jet pT E and const spectra
void SpectraFillHistosPartLevel(Int_t NormValue = 0) { //matched particle level jets

    const int nNorm = 8; //8 bins but we skip the first 5-10 for now
    const Float_t normalizations[nNorm] = {0.0610658,0.00716477,0.000557627,0.000107816,4.31694e-05,9.62255e-06,1.24904e-06,5.99517e-08};

    const int nBinspT = 100;
    const int nBinsE = 200;
    const double ptmin = 5;
    const double ptmax = 205;
    const double emin = 10;
    const double emax = 2010;
    const double cmin = -0.5;
    const double cmax = 49.5;
    const int nBinspTCon = (int)(cmax-cmin);

    // Arrays to hold histograms for different types and R values
    TH1F* hist_pT_R[3]; // hist_pT_R[0], hist_pT_R[1], hist_pT_R[2]
    TH1F* hist_E_R[3]; // hist_E_R[0], hist_E_R[1], hist_E_R[2]
    TH1F* hist_con_R[3]; // hist_con_R[0], hist_con_R[1], hist_con_R[2]

    TH1F* hist_pT_R_e[3][2]; // hist_pT_R_e[0][0], hist_pT_R_e[0][1], ..., hist_pT_R_e[2][2]
    TH1F* hist_E_R_e[3][2]; // hist_E_R_e[0][0], hist_E_R_e[0][1], ..., hist_E_R_e[2][2]
    TH1F* hist_con_R_e[3][2]; // hist_con_R_e[0][0], hist_con_R_e[0][1], ..., hist_con_R_e[2][2]

    //ADD TWO MORE HISTO ARRAYs to hold pT and E distribution for const = 1,2,3,5,10. 
    TH1F* hist_pTvscon_R_e[3][2][5];
    TH1F* hist_Evscon_R_e[3][2][5];
    TH1F* hist_pTvscon_R[3][5];
    TH1F* hist_Evscon_R[3][5];

    //Read in the data
    TFile* file;
    //if(NormValue==0)file = TFile::Open("/home/lmh/alice/fromPuhti/April2023/pTmin2_20230417_pythia8_JetJet_5-15GeV/Merged.root");
    if(NormValue==0)file = TFile::Open("/home/lmh/alice/fromPuhti/July2023/data/20230722_pythia8_JetJet_5-10GeV/Merged.root");

    if(NormValue==1)file = TFile::Open("/home/lmh/alice/fromPuhti/July2023/ptmin10/data/ptmin10_20230722_pythia8_JetJet_10-20GeV/Merged.root");
    if(NormValue==2)file = TFile::Open("/home/lmh/alice/fromPuhti/July2023/ptmin10/data/ptmin10_20230722_pythia8_JetJet_20-30GeV/Merged.root");
    if(NormValue==3)file = TFile::Open("/home/lmh/alice/fromPuhti/July2023/ptmin10/data/ptmin10_20230722_pythia8_JetJet_30-40GeV/Merged.root");
    if(NormValue==4)file = TFile::Open("/home/lmh/alice/fromPuhti/July2023/ptmin10/data/ptmin10_20230722_pythia8_JetJet_40-60GeV/Merged.root");
    if(NormValue==5)file = TFile::Open("/home/lmh/alice/fromPuhti/July2023/ptmin10/data/ptmin10_20230722_pythia8_JetJet_60-100GeV/Merged.root");
    if(NormValue==6)file = TFile::Open("/home/lmh/alice/fromPuhti/July2023/ptmin10/data/ptmin10_20230722_pythia8_JetJet_100-200GeV/Merged.root");
    if(NormValue==7)file = TFile::Open("/home/lmh/alice/fromPuhti/July2023/ptmin10/data/ptmin10_20230722_pythia8_JetJet_200-GeV/Merged.root");


    const Int_t nR = 3; 
    const Float_t Rvals[nR] = {0.2, 0.4, 0.6}; // Cone radii
    const Int_t nEtaBins  = 3;
    //const double EtaBinBorders[nEtaBins] = {3.4,5.0,5.5}; //
    const double EtaMin[2] = {3.4,3.9}; 
    const double EtaMax[2] = {5.0,5.5};
    const double EtaBinBorders[nR][nEtaBins] = {{3.6,4.5,5.3},{3.8,4.5,5.1},{4.0,4.5,4.9}}; //use this instead!!
    TString etaRange[2]     = {"3.4+R< #eta_{jet} < 4.5","4.5 < #eta_{jet} < 5.5-R"}; //consider making this different range for R=0.6, 0.4, 0.2. ...
    const Int_t nConst = 5;
    const Int_t ConstVals[nConst] = {1, 2, 3, 5, 10}; // constituent values to check

    // Create a ROOT file to save the histograms
    TFile *fout;
    //if(NormValue==0)fout = new TFile("JetJetOutput/NewBinning/TEST.root", "RECREATE");
    if(NormValue==0)fout = new TFile("JetJetOutput/July2023/data/ptmin10/match_20230722_pythia8_JetJet_5-10GeV_Merged_Output.root", "RECREATE");
    if(NormValue==1)fout = new TFile("JetJetOutput/July2023/data/ptmin10/match_20230722_pythia8_JetJet_10-20GeV_Merged_Output.root", "RECREATE");
    if(NormValue==2)fout = new TFile("JetJetOutput/July2023/data/ptmin10/match_20230722_pythia8_JetJet_20-30GeV_Merged_Output.root", "RECREATE");
    if(NormValue==3)fout = new TFile("JetJetOutput/July2023/data/ptmin10/match_20230722_pythia8_JetJet_30-40GeV_Merged_Output.root", "RECREATE");
    if(NormValue==4)fout = new TFile("JetJetOutput/July2023/data/ptmin10/match_20230722_pythia8_JetJet_40-60GeV_Merged_Output.root", "RECREATE");
    if(NormValue==5)fout = new TFile("JetJetOutput/July2023/data/ptmin10/match_20230722_pythia8_JetJet_60-100GeV_Merged_Output.root", "RECREATE");
    if(NormValue==6)fout = new TFile("JetJetOutput/July2023/data/ptmin10/match_20230722_pythia8_JetJet_100-200GeV_Merged_Output.root", "RECREATE");
    if(NormValue==7)fout = new TFile("JetJetOutput/July2023/data/ptmin10/match_20230722_pythia8_JetJet_200-GeV_Merged_Output.root", "RECREATE");

    TTree *tree = (TTree *)file->Get("jetTree");

    // Loop over different R values
    for (int R = 0; R < 3; ++R) {
        // Initialize hist_pT_R[3] histograms
        hist_pT_R[R] = new TH1F(Form("hist_pT_R%d", R), Form("pT Distribution for particle jets, R = %0.1f; pT [GeV/c]; Entries", Rvals[R]), nBinspT, ptmin, ptmax);
        tree->Draw(Form("jetpT_match >> hist_pT_R%d", R), Form("jetR == %d && jetEta_match>=%f && jetEta_match<%f", int(Rvals[R] * 10), EtaMin[0]+Rvals[R], EtaMax[1]-Rvals[R]));
        hist_pT_R[R]->Scale(normalizations[NormValue]/500.0);

        // Initialize hist_E_R[3] histograms
        hist_E_R[R] = new TH1F(Form("hist_E_R%d", R), Form("E Distribution for particle jets, R = %0.1f; E [GeV]; Entries", Rvals[R]), nBinsE, emin, emax);
        tree->Draw(Form("jetE_match >> hist_E_R%d", R), Form("jetR == %d && jetEta_match>=%f && jetEta_match<%f", int(Rvals[R] * 10), EtaMin[0]+Rvals[R], EtaMax[1]-Rvals[R]));
        hist_E_R[R]->Scale(normalizations[NormValue]/500.0);

        // Initialize hist_con_R[3] histograms
        hist_con_R[R] = new TH1F(Form("hist_con_R%d", R), Form("Constituent Distribution for particle jets, R = %0.1f; Constituent; Entries", Rvals[R]), nBinspTCon, cmin, cmax);
        tree->Draw(Form("jetParts_match>> hist_con_R%d", R), Form("jetR == %d && jetEta_match>=%f && jetEta_match<%f", int(Rvals[R] * 10), EtaMin[0]+Rvals[R], EtaMax[1]-Rvals[R]));
        hist_con_R[R]->Scale(normalizations[NormValue]/500.0);

        for (int e = 0; e < 2; ++e) {
            // Initialize hist_pT_R_e[3][3] histograms
            hist_pT_R_e[R][e] = new TH1F(Form("hist_pT_R%d_e%d", R, e), Form("pT Distribution for particle jets, R = %0.1f, %s; pT [GeV/c]; Entries", Rvals[R], etaRange[e].Data()), nBinspT, ptmin, ptmax);
            tree->Draw(Form("jetpT_match >> hist_pT_R%d_e%d", R, e), Form("jetR == %d && jetEta_match>=%f && jetEta_match<%f", int(Rvals[R] * 10), EtaBinBorders[R][e], EtaBinBorders[R][e+1]));
            hist_pT_R_e[R][e]->Scale(normalizations[NormValue]/500.0);

            // Initialize hist_E_R_e[3][3] histograms
            hist_E_R_e[R][e] = new TH1F(Form("hist_E_R%d_e%d", R, e), Form("E Distribution for particle jets, R = %0.1f, %s; E [GeV]; Entries", Rvals[R], etaRange[e].Data()), nBinsE, emin, emax);
            tree->Draw(Form("jetE_match >> hist_E_R%d_e%d", R, e), Form("jetR == %d && jetEta_match>=%f && jetEta_match<%f", int(Rvals[R] * 10), EtaBinBorders[R][e], EtaBinBorders[R][e+1]));
            hist_E_R_e[R][e]->Scale(normalizations[NormValue]/500.0);

            // Initialize hist_con_R_e[3][3] histograms
            hist_con_R_e[R][e] = new TH1F(Form("hist_con_R%d_e%d", R, e), Form("Constituent Distribution for particle jets, R = %0.1f, %s; Constituent; Entries", Rvals[R], etaRange[e].Data()), nBinspTCon, cmin, cmax);
            tree->Draw(Form("jetParts_match>> hist_con_R%d_e%d", R, e), Form("jetR == %d && jetEta_match>=%f && jetEta_match<%f", int(Rvals[R] * 10), EtaBinBorders[R][e], EtaBinBorders[R][e+1]));
            hist_con_R_e[R][e]->Scale(normalizations[NormValue]/500.0);

        //ADD LOOP over the constituent numbers chosen (5) and fill the final histos!
            for (int c = 0; c < 5; ++c) {
                hist_pTvscon_R_e[R][e][c] = new TH1F(Form("hist_pTvscon_R%d_e%d_c%d", R, e, c), Form("Constituent N = %d Distribution for particle jets, R = %0.1f, %s; pT [GeV/c]; Entries", ConstVals[c], Rvals[R], etaRange[e].Data()), nBinspT, ptmin, ptmax);
                tree->Draw(Form("jetpT_match >> hist_pTvscon_R%d_e%d_c%d", R, e, c), Form("jetR == %d && jetEta_match>=%f && jetEta_match<%f && jetParts_match== %d", int(Rvals[R] * 10), EtaBinBorders[R][e], EtaBinBorders[R][e+1], ConstVals[c]));
                hist_pTvscon_R_e[R][e][c]->Scale(normalizations[NormValue]/500.0);

                hist_Evscon_R_e[R][e][c] = new TH1F(Form("hist_Evscon_R%d_e%d_c%d", R, e, c), Form("Constituent N = %d Distribution for particle jets, R = %0.1f, %s; E [GeV]; Entries", ConstVals[c], Rvals[R], etaRange[e].Data()), nBinsE, emin, emax);
                tree->Draw(Form("jetE_match >> hist_Evscon_R%d_e%d_c%d", R, e, c), Form("jetR == %d && jetEta_match>=%f && jetEta_match<%f && jetParts_match== %d", int(Rvals[R] * 10), EtaBinBorders[R][e], EtaBinBorders[R][e+1], ConstVals[c]));
                hist_Evscon_R_e[R][e][c]->Scale(normalizations[NormValue]/500.0);
            }

        }
        for (int c = 0; c < 5; ++c) {
            hist_pTvscon_R[R][c] = new TH1F(Form("hist_pTvscon_R%d_c%d", R, c), Form("Constituent N = %d Distribution for particle jets, R = %0.1f; pT [GeV/c]; Entries", ConstVals[c], Rvals[R]), nBinspT, ptmin, ptmax);
            tree->Draw(Form("jetpT_match >> hist_pTvscon_R%d_c%d", R, c), Form("jetR == %d && jetEta_match>=%f && jetEta_match<%f && jetParts_match== %d", int(Rvals[R] * 10), EtaMin[0]+Rvals[R], EtaMax[1]-Rvals[R], ConstVals[c]));
            hist_pTvscon_R[R][c]->Scale(normalizations[NormValue]/500.0);

            hist_Evscon_R[R][c] = new TH1F(Form("hist_Evscon_R%d_c%d", R, c), Form("Constituent N = %d Distribution for particle jets, R = %0.1f; E [GeV]; Entries", ConstVals[c], Rvals[R]), nBinsE, emin, emax);
            tree->Draw(Form("jetE_match >> hist_Evscon_R%d_c%d", R, c), Form("jetR == %d && jetEta_match>=%f && jetEta_match<%f && jetParts_match== %d", int(Rvals[R] * 10), EtaMin[0]+Rvals[R], EtaMax[1]-Rvals[R], ConstVals[c]));
            hist_Evscon_R[R][c]->Scale(normalizations[NormValue]/500.0);
        }
    }

    // Write histograms to the ROOT file and close the file
    fout->cd();
    for (int R = 0; R < 3; ++R) {
        hist_pT_R[R]->Write();
        hist_E_R[R]->Write();
        hist_con_R[R]->Write();
        for (int e = 0; e < 2; ++e) {
            hist_pT_R_e[R][e]->Write();
            hist_E_R_e[R][e]->Write();
            hist_con_R_e[R][e]->Write();
            for (int c = 0; c < 5; ++c) {
                hist_Evscon_R_e[R][e][c]->Write();
                hist_pTvscon_R_e[R][e][c]->Write();
            }
        }
        for (int c = 0; c < 5; ++c) {
            hist_Evscon_R[R][c]->Write();
            hist_pTvscon_R[R][c]->Write();
        }
    }
    fout->Close();
}





    
