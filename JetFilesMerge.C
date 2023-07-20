/* Macro for normalizing some histograms from jet trees 
    2.21064
    0.0669805
    0.00182628
    0.000139462
    0.0000225822
*/
void JetFilesMerge()
{   
    const int nNorm = 5;
    const Float_t[nNorm] = {2.21064, 0.0669805, 0.00182628, 0.000139462, 0.0000225822};

    TFile *jetFile1 = TFile::Open("/home/lmh/alice/fromPuhti/April2023/20230417_pythia8_JetJet_0-5GeV/Merged.root");
    TFile *jetFile2 = TFile::Open("/home/lmh/alice/fromPuhti/April2023/20230417_pythia8_JetJet_5-15GeV/Merged.root");
    TFile *jetFile3 = TFile::Open("/home/lmh/alice/fromPuhti/April2023/20230417_pythia8_JetJet_15-30GeV/Merged.root");
    TFile *jetFile4 = TFile::Open("/home/lmh/alice/fromPuhti/April2023/20230417_pythia8_JetJet_30-50GeV/Merged.root");
    TFile *jetFile5 = TFile::Open("/home/lmh/alice/fromPuhti/April2023/20230417_pythia8_JetJet_50-1000GeV/Merged.root");

    TFile *fout1     = new TFile("JetJetOutput/MergedData/20230417_Merged_Output_0-5GeV.root");
    TFile *fout2     = new TFile("JetJetOutput/MergedData/20230417_Merged_Output_5-15GeV.root");
    TFile *fout3     = new TFile("JetJetOutput/MergedData/20230417_Merged_Output_15-30GeV.root");
    TFile *fout4     = new TFile("JetJetOutput/MergedData/20230417_Merged_Output_30-50GeV.root");
    TFile *fout5     = new TFile("JetJetOutput/MergedData/20230417_Merged_Output_50-1000GeV.root");

    TTree *jetTree = (TTree *)jetFile->Get("jetTree");
    TTree *TruthjetTree = (TTree *)jetFile->Get("truthjetTree");
    //TTree *PDGTree = (TTree *)jetFile->Get("inPDGTree");

    const Int_t nR = 3; //5
    const Float_t Rvals[nR] = {0.2, 0.4, 0.6};//{0.2, 0.3, 0.4, 0.5, 0.6}; // Cone radii
    int Rvalue = 1; // choose the index of the jet R you want to draw the main histos for !

    //histograms

    
    fout->Write();
    fout->Close();

}


