/* Macro for plotting some histograms from jet trees 
    2.21064
    0.0669805
    0.00182628
    0.000139462
    0.0000225822
*/

void JetPlottingMergedFractions2D(int NormValue=0)
{   
    const int nNorm = 8;
    const Float_t normalizations[nNorm] = {0.0611191, 0.00717001, 0.000558759, 0.000107936, 4.32163e-05, 9.57109e-06, 1.24606e-06, 6.01382e-08};//{0.0610658,0.00716477,0.000557627,0.000107816,4.31694e-05,9.62255e-06,1.24904e-06,5.99517e-08};

    //int NormValue = 0;

    TFile *jetFile;
    
    if(NormValue==0)jetFile = TFile::Open("/home/lmh/alice/fromPuhti/July2023/20230728JetJet/20230728_pythia8_JetJet_5-10GeV/Merged.root");
    if(NormValue==1)jetFile = TFile::Open("/home/lmh/alice/fromPuhti/July2023/20230728JetJet/20230728_pythia8_JetJet_10-20GeV/Merged.root");
    if(NormValue==2)jetFile = TFile::Open("/home/lmh/alice/fromPuhti/July2023/20230728JetJet/20230728_pythia8_JetJet_20-30GeV/Merged.root");
    if(NormValue==3)jetFile = TFile::Open("/home/lmh/alice/fromPuhti/July2023/20230728JetJet/20230728_pythia8_JetJet_30-40GeV/Merged.root");
    if(NormValue==4)jetFile = TFile::Open("/home/lmh/alice/fromPuhti/July2023/20230728JetJet/20230728_pythia8_JetJet_40-60GeV/Merged.root");
    if(NormValue==5)jetFile = TFile::Open("/home/lmh/alice/fromPuhti/July2023/20230728JetJet/20230728_pythia8_JetJet_60-100GeV/Merged.root");
    if(NormValue==6)jetFile = TFile::Open("/home/lmh/alice/fromPuhti/July2023/20230728JetJet/20230728_pythia8_JetJet_100-200GeV/Merged.root");
    if(NormValue==7)jetFile = TFile::Open("/home/lmh/alice/fromPuhti/July2023/20230728JetJet/20230728_pythia8_JetJet_200-GeV/Merged.root");

    const Int_t nR = 3; 
    const Float_t Rvals[nR] = {0.2, 0.4, 0.6}; // Cone radii
    //const double EtaMin[2] = {3.4,3.9}; 
    //const double EtaMax[2] = {5.0,5.5};
    

    const Float_t etaMin = 3.4; const Float_t etaMax = 5.5; //changed in June23, 3.4 to 5.5  

    //add eta range binning here, which depends on the value of R.   
    const Int_t nEtaBins  = 3;
    const double EtaBinBorders[nEtaBins][nR] = {{3.6, 4.5, 5.3}, {3.8, 4.5, 5.1}, {4.0, 4.5, 4.9}};
    TString etaRange[nEtaBins-1] = {"3.4+R < #eta_{jet} < 4.5", "4.5 < #eta_{jet} < 5.5-R"};
    const Int_t nEBins  = 15;//16;
    const Int_t npTBins  = 9;//16;
    const double JetEBorders[nEBins] =  {100.0, 200.0, 300.0, 400.0, 500.0, 600.0, 700.0, 800.0, 900.0, 1000.0, 1250.0, 1500.0, 1750.0, 2000.0, 3000.0};//{100.0, 400.0, 800.0, 1200.0, 1600.0, 2000.0}; 
    const double JetpTBorders[npTBins] = {5.0, 10.0, 20.0, 40.0, 60.0, 80.0, 100.0, 150.0, 200.0}; 
    const Int_t nConst = 5;
    const Int_t ConstVals[nConst] = {1, 2, 3, 5, 10}; // constituent values to check


    int constMin = 0; //min number of constituents in matchd jet 2

    TFile *fout;
    //if(NormValue==0)fout = new TFile("JetJetOutput/July2023/fractions/2D/TEST.root", "RECREATE");
    if(NormValue==0)fout = new TFile("Data20230728/2D/2D_20230722_pythia8_JetJet_5-10GeV_Merged_Output.root", "RECREATE");
    if(NormValue==1)fout = new TFile("Data20230728/2D/2D_20230722_pythia8_JetJet_10-20GeV_Merged_Output.root", "RECREATE");
    if(NormValue==2)fout = new TFile("Data20230728/2D/2D_20230722_pythia8_JetJet_20-30GeV_Merged_Output.root", "RECREATE");
    if(NormValue==3)fout = new TFile("Data20230728/2D/2D_20230722_pythia8_JetJet_30-40GeV_Merged_Output.root", "RECREATE");
    if(NormValue==4)fout = new TFile("Data20230728/2D/2D_20230722_pythia8_JetJet_40-60GeV_Merged_Output.root", "RECREATE");
    if(NormValue==5)fout = new TFile("Data20230728/2D/2D_20230722_pythia8_JetJet_60-100GeV_Merged_Output.root", "RECREATE");
    if(NormValue==6)fout = new TFile("Data20230728/2D/2D_20230722_pythia8_JetJet_100-200GeV_Merged_Output.root", "RECREATE");
    if(NormValue==7)fout = new TFile("Data20230728/2D/2D_20230722_pythia8_JetJet_200-GeV_Merged_Output.root", "RECREATE");

    TTree *jetTree = (TTree *)jetFile->Get("jetTree");

    TH2F *hRespMatrix[nR];
    TH2F *hRespMatrixE[nR][nEBins-1];
    TH2F *hRespMatrixpT[nR][npTBins-1];


    for(int Rvalue=0; Rvalue<nR; Rvalue++){
    //histograms
        hRespMatrix[Rvalue] = new TH2F(Form("hRespMatrix_R%d", Rvalue), Form("Detector ECAL fraction vs Particle neutral energy fraction, R=%0.1f;#it{E}_{ECAL}/#it{E}_{jet}^{det} (GeV);#it{E}_{neutral}/#it{E}_{jet}^{part} (GeV)", Rvals[Rvalue]), 41, -0.0125, 1.0125, 41, -0.0125, 1.0125); 

        jetTree->Draw(Form("matchedjetECALEnergyFrac:jetECALEnergyFrac>>hRespMatrix_R%d", Rvalue), Form("jetR==%d && jetEta>%f && jetEta<%f && jetEta_match>%f && jetEta_match<%f && jet_distmatch<%f && jetParts_match > %d&& jetParts > %d", int(Rvals[Rvalue] * 10), etaMin+Rvals[Rvalue], etaMax-Rvals[Rvalue], etaMin+Rvals[Rvalue], etaMax-Rvals[Rvalue], Rvals[Rvalue]*0.6, constMin, constMin), "goff");

        hRespMatrix[Rvalue]->Scale(normalizations[NormValue]);
        //hRespMatrix[Rvalue]->Scale(1./hRespMatrix[Rvalue]->Integral(), "width"); //only perform this normalization at final stage of plotting

        for (int iE = 0; iE < nEBins-1; ++iE) 
        {
            hRespMatrixE[Rvalue][iE] = new TH2F(Form("hRespMatrixE%d_R%d", iE, Rvalue), Form("Detector ECAL fraction vs Particle neutral energy fraction, E_{jet}^{det}: %d - %d GeV, R=%0.1f;#it{E}_{ECAL}/#it{E}_{jet}^{det} (GeV);#it{E}_{neutral}/#it{E}_{jet}^{part} (GeV)", int(JetEBorders[iE]), int(JetEBorders[iE+1]), Rvals[Rvalue]), 41, -0.0125, 1.0125, 41, -0.0125, 1.0125); 
            jetTree->Draw(Form("matchedjetECALEnergyFrac:jetECALEnergyFrac>>hRespMatrixE%d_R%d", iE, Rvalue), Form("jetR==%d && jetE>=%d  && jetE<%d && jetEta>%f && jetEta<%f && jetEta_match>%f && jetEta_match<%f && jet_distmatch<%f && jetParts_match > %d&& jetParts > %d", int(Rvals[Rvalue] * 10), int(JetEBorders[iE]), int(JetEBorders[iE+1]), etaMin+Rvals[Rvalue], etaMax-Rvals[Rvalue], etaMin+Rvals[Rvalue], etaMax-Rvals[Rvalue], Rvals[Rvalue]*0.6, constMin, constMin), "goff");

            hRespMatrixE[Rvalue][iE]->Scale(normalizations[NormValue]);
            //hRespMatrixE[Rvalue][iE]->Scale(1./hRespMatrixE[Rvalue][iE]->Integral(), "width"); //only perform this normalization at final stage of plotting
        }

        for (int ipT = 0; ipT < npTBins-1; ++ipT) 
        {
            hRespMatrixpT[Rvalue][ipT] = new TH2F(Form("hRespMatrixpT%d_R%d", ipT, Rvalue), Form("Detector ECAL fraction vs Particle neutral energy fraction, pT_{jet}^{det}: %d - %d GeV, R=%0.1f;#it{E}_{ECAL}/#it{E}_{jet}^{det} (GeV);#it{E}_{neutral}/#it{E}_{jet}^{part} (GeV)", int(JetpTBorders[ipT]), int(JetpTBorders[ipT+1]), Rvals[Rvalue]), 41, -0.0125, 1.0125, 41, -0.0125, 1.0125); 
            jetTree->Draw(Form("matchedjetECALEnergyFrac:jetECALEnergyFrac>>hRespMatrixpT%d_R%d", ipT, Rvalue), Form("jetR==%d && jetpT>=%d  && jetpT<%d && jetEta>%f && jetEta<%f && jetEta_match>%f && jetEta_match<%f && jet_distmatch<%f && jetParts_match > %d&& jetParts > %d", int(Rvals[Rvalue] * 10), int(JetpTBorders[ipT]), int(JetpTBorders[ipT+1]), etaMin+Rvals[Rvalue], etaMax-Rvals[Rvalue], etaMin+Rvals[Rvalue], etaMax-Rvals[Rvalue], Rvals[Rvalue]*0.6, constMin, constMin), "goff");

            hRespMatrixpT[Rvalue][ipT]->Scale(normalizations[NormValue]);
            //hRespMatrixpT[Rvalue][ipT]->Scale(1./hRespMatrixpT[Rvalue][ipT]->Integral(), "width"); //only perform this normalization at final stage of plotting
        }
    }

    fout->Write();
    fout->Close();
}