/* Macro for plotting some histograms from jet trees 
    2.21064
    0.0669805
    0.00182628
    0.000139462
    0.0000225822
*/

void JetPlottingMergedFractions(int NormValue=0,int Rvalue=0)
{   
    const int nNorm = 6;
    const Float_t normalizations[nNorm] = {2.21064, 0.0669805, 0.00182628, 0.000139462, 0.0000225822,1.0};

    //int NormValue = 0;

    TFile *jetFile;
    if(NormValue==0)jetFile = TFile::Open("/home/lmh/alice/fromPuhti/April2023/fractionsJuly_20230417_pythia8_JetJet_0-5GeV/Merged.root");
    if(NormValue==1)jetFile = TFile::Open("/home/lmh/alice/fromPuhti/April2023/fractionsJuly_20230417_pythia8_JetJet_5-15GeV/Merged.root");
    if(NormValue==2)jetFile = TFile::Open("/home/lmh/alice/fromPuhti/April2023/fractionsJuly_20230417_pythia8_JetJet_15-30GeV/Merged.root");
    if(NormValue==3)jetFile = TFile::Open("/home/lmh/alice/fromPuhti/April2023/fractionsJuly_20230417_pythia8_JetJet_30-50GeV/Merged.root");
    if(NormValue==4)jetFile = TFile::Open("/home/lmh/alice/fromPuhti/April2023/fractionsJuly_20230417_pythia8_JetJet_50-1000GeV/Merged.root");
    //TFile *jetFile = TFile::Open("JetJetOutput/MergedData/20230417_Merged/Merged.root");

    if(NormValue==5)jetFile = TFile::Open("/home/lmh/alice/jetAnalysis/NewPlotting/Fractions/Results/fractionsptmin2__20230127_PythiaMBTrig-2/Merged.root");

    const Int_t nR = 3; //5
    const Float_t Rvals[nR] = {0.2, 0.4, 0.6};//{0.2, 0.3, 0.4, 0.5, 0.6}; // Cone radii
    //int Rvalue = 0; // choose the index of the jet R you want to draw the main histos for !
    int constMin = 1; //min number of constituents in matchd jet 2

    TFile *fout;
    if(NormValue==0)fout = new TFile(Form("/home/lmh/alice/jetAnalysis/NewPlotting/Fractions/Results/100FIXEDfractionsJuly_20230417_0-5GeV_Merged_OutputR%d.root", int(Rvals[Rvalue] * 10)), "RECREATE");
    if(NormValue==1)fout = new TFile(Form("/home/lmh/alice/jetAnalysis/NewPlotting/Fractions/Results/100FIXEDfractionsJuly_20230417_5-15GeV_Merged_OutputR%d.root", int(Rvals[Rvalue] * 10)), "RECREATE");
    if(NormValue==2)fout = new TFile(Form("/home/lmh/alice/jetAnalysis/NewPlotting/Fractions/Results/100FIXEDfractionsJuly_20230417_15-30GeV_Merged_OutputR%d.root", int(Rvals[Rvalue] * 10)), "RECREATE");
    if(NormValue==3)fout = new TFile(Form("/home/lmh/alice/jetAnalysis/NewPlotting/Fractions/Results/100FIXEDfractionsJuly_20230417_30-50GeV_Merged_OutputR%d.root", int(Rvals[Rvalue] * 10)), "RECREATE");
    if(NormValue==4)fout = new TFile(Form("/home/lmh/alice/jetAnalysis/NewPlotting/Fractions/Results/100FIXEDfractionsJuly_20230417_50-1000GeV_Merged_OutputR%d.root", int(Rvals[Rvalue] * 10)), "RECREATE");

    if(NormValue==5)fout = new TFile(Form("/home/lmh/alice/jetAnalysis/NewPlotting/Fractions/Results/oneOfractionsptmin2__20230127_PythiaMBTrig-2_Merged_OutputR%d.root", int(Rvals[Rvalue] * 10)), "RECREATE");

    TTree *jetTree = (TTree *)jetFile->Get("jetTree");
    TTree *TruthjetTree = (TTree *)jetFile->Get("truthjetTree");
    //TTree *PDGTree = (TTree *)jetFile->Get("inPDGTree");
    const int nCol = 10;
    const int gcolors[nCol]={1,2,6,4,7,1,2,4,6,7};
    const int gmarkers[nCol]={4,8,25,21,8,21,25,4,8,21};

    const Float_t etaMin = 3.4; const Float_t etaMax = 5.5; //changed in June23, 3.4 to 5.5     

    const Int_t nEBins  = 6;//16;
    //const double JetEBorders[nEBins] = {0.0, 50.0, 100.0, 200.0, 300.0, 400.0, 500.0, 600.0, 700.0, 800.0, 900.0, 1000.0, 1250.0, 1500.0, 1750.0, 2000.0}; 
    const double JetEBorders[nEBins] = {100.0, 400.0, 800.0, 1200.0, 1600.0, 2000.0}; 

    //histograms

    TH1F *hDetHCALE[nEBins]; //HCALfracdet, drawn like delta E in bins, and then take the mean and stdev of the different bins. 
    TH1F *hDetECALE[nEBins]; //ECALfracdet

    TH1F *hPartNeutralE[nEBins]; //neutralfrac
    TH1F *hPartChargedE[nEBins]; //chargedfrac
    
    TH1F *hPartPhotonElectronE[nEBins]; //photon&electronfrac


    TH2F *hRespMatrix_DetHCALE = new TH2F("hRespMatrix_DetHCALE", Form("Detector HCAL fraction, R=%0.1f;#it{E}_{jet}^{det} (GeV);#it{E}_{HCAL}/#it{E}_{jet}^{det} (GeV)", Rvals[Rvalue]), 500, 0, 2000, 500, 0, 1);
    TH2F *hRespMatrix_DetECALE = new TH2F("hRespMatrix_DetECALE", Form("Detector ECAL fraction, R=%0.1f;#it{E}_{jet}^{det} (GeV);#it{E}_{ECAL}/#it{E}_{jet}^{det} (GeV)", Rvals[Rvalue]),  500, 0, 2000, 500, 0, 1);

    TH2F *hRespMatrix_PartNeutralE = new TH2F("hRespMatrix_PartNeutralE", Form("Particle neutral fraction, R=%0.1f;#it{E}_{jet}^{part} (GeV);#it{E}_{neutral}/#it{E}_{jet}^{part} (GeV)", Rvals[Rvalue]),  500, 0, 2000, 500, 0, 1);
    TH2F *hRespMatrix_PartChargedE = new TH2F("hRespMatrix_PartChargedE", Form("Particle charged fraction, R=%0.1f;#it{E}_{jet}^{part} (GeV);#it{E}_{charged}/#it{E}_{jet}^{part} (GeV)", Rvals[Rvalue]),  500, 0, 2000, 500, 0, 1);

    TH2F *hRespMatrix_PartPhotonElectronE = new TH2F("hRespMatrix_PartPhotonElectronE", Form("Particle photon and electron fraction, R=%0.1f;#it{E}_{jet}^{part} (GeV);#it{E}_{#gamma, e}/#it{E}_{jet}^{det} (GeV)", Rvals[Rvalue]),  500, 0, 2000, 500, 0, 1);

    //begin filling the histograms

    jetTree->Draw("jetHCALEnergyFrac:jetE>>hRespMatrix_DetHCALE", Form("jetR==%d && jetEta>%f && jetEta<%f && jetEta_match>%f && jetEta_match<%f && jet_distmatch<%f && jetParts_match > %d&& jetParts > %d", int(Rvals[Rvalue] * 10), etaMin+Rvals[Rvalue], etaMax-Rvals[Rvalue], etaMin+Rvals[Rvalue], etaMax-Rvals[Rvalue], Rvals[Rvalue]*0.6, constMin, constMin), "goff"); 

    jetTree->Draw("jetECALEnergyFrac:jetE>>hRespMatrix_DetECALE", Form("jetR==%d && jetEta>%f && jetEta<%f && jetEta_match>%f && jetEta_match<%f && jet_distmatch<%f && jetParts_match > %d&& jetParts > %d", int(Rvals[Rvalue] * 10), etaMin+Rvals[Rvalue], etaMax-Rvals[Rvalue], etaMin+Rvals[Rvalue], etaMax-Rvals[Rvalue], Rvals[Rvalue]*0.6, constMin, constMin), "goff");

    TruthjetTree->Draw("truthjetNeutralEnergyFrac:truthjetE>>hRespMatrix_PartNeutralE", Form("truthjetR==%d && truthjetEta>%f && truthjetEta<%f && jetParts_truth > %d", int(Rvals[Rvalue] * 10), etaMin+Rvals[Rvalue], etaMax-Rvals[Rvalue], constMin), "goff");

    TruthjetTree->Draw("truthjetChargedEnergyFrac:truthjetE>>hRespMatrix_PartChargedE", Form("truthjetR==%d && truthjetEta>%f && truthjetEta<%f && jetParts_truth > %d", int(Rvals[Rvalue] * 10), etaMin+Rvals[Rvalue], etaMax-Rvals[Rvalue], constMin), "goff");

    TruthjetTree->Draw("truthjetECALEnergyFrac:truthjetE>>hRespMatrix_PartPhotonElectronE", Form("truthjetR==%d && truthjetEta>%f && truthjetEta<%f && jetParts_truth > %d", int(Rvals[Rvalue] * 10), etaMin+Rvals[Rvalue], etaMax-Rvals[Rvalue], constMin), "goff");
    
    hRespMatrix_DetHCALE->Scale(normalizations[NormValue]);
    //hRespMatrix_E->Scale(1./hRespMatrix_E->GetEntries(), "width"); //only perform this normalization at final stage of plotting
    
    hRespMatrix_DetECALE->Scale(normalizations[NormValue]);
    hRespMatrix_PartNeutralE->Scale(normalizations[NormValue]);
    hRespMatrix_PartChargedE->Scale(normalizations[NormValue]);
    hRespMatrix_PartPhotonElectronE->Scale(normalizations[NormValue]);



    //energy loop
    for (int iE = 0; iE < nEBins-1; ++iE) 
    {
        // 1D histograms
        hDetHCALE[iE] = new TH1F(Form("hDetHCALE_%d", iE), Form("Jet-by-jet detector level HCAL energy fraction distribution, E_{jet}^{det}: %d - %d GeV", int(JetEBorders[iE]), int(JetEBorders[iE+1])), 50, 0.0, 1.0);
        hDetECALE[iE] = new TH1F(Form("hDetECALE_%d", iE), Form("Jet-by-jet detector level ECAL energy fraction distribution, E_{jet}^{det}: %d - %d GeV", int(JetEBorders[iE]), int(JetEBorders[iE+1])), 50, 0.0, 1.0);
        hPartNeutralE[iE] = new TH1F(Form("hPartNeutralE_%d", iE), Form("Jet-by-jet particle level neutral constituent energy fraction distribution, E_{jet}^{part}: %d - %d GeV", int(JetEBorders[iE]), int(JetEBorders[iE+1])), 50, 0.0, 1.0);
        hPartChargedE[iE] = new TH1F(Form("hPartChargedE_%d", iE), Form("Jet-by-jet particle level charged constituent energy fraction distribution, E_{jet}^{part}: %d - %d GeV", int(JetEBorders[iE]), int(JetEBorders[iE+1])), 50, 0.0, 1.0);     
        hPartPhotonElectronE[iE] = new TH1F(Form("hPartPhotonElectronE_%d", iE), Form("Jet-by-jet particle level #gamma,e constituent energy fraction distribution, matched jets only, E_{jet}^{part}: %d - %d GeV", int(JetEBorders[iE]), int(JetEBorders[iE+1])), 50, 0.0, 1.0);


        //Fill histos
        jetTree->Draw(Form("jetHCALEnergyFrac>>hDetHCALE_%d", iE), Form("jetR==%d && jetE>=%d  && jetE<%d  && jetEta>%f && jetEta<%f && jetEta_match>%f && jetEta_match<%f && jet_distmatch<%f && jetParts_match > %d && jetParts > %d", int(Rvals[Rvalue] * 10), int(JetEBorders[iE]), int(JetEBorders[iE+1]),etaMin+Rvals[Rvalue],etaMax-Rvals[Rvalue],etaMin+Rvals[Rvalue],etaMax-Rvals[Rvalue], Rvals[Rvalue]*0.6, constMin, constMin), "goff"); 
        jetTree->Draw(Form("jetECALEnergyFrac>>hDetECALE_%d", iE), Form("jetR==%d && jetE>=%d  && jetE<%d  && jetEta>%f && jetEta<%f && jetEta_match>%f && jetEta_match<%f && jet_distmatch<%f && jetParts_match > %d && jetParts > %d", int(Rvals[Rvalue] * 10), int(JetEBorders[iE]), int(JetEBorders[iE+1]),etaMin+Rvals[Rvalue],etaMax-Rvals[Rvalue],etaMin+Rvals[Rvalue],etaMax-Rvals[Rvalue], Rvals[Rvalue]*0.6, constMin, constMin), "goff"); 

        TruthjetTree->Draw(Form("truthjetNeutralEnergyFrac>>hPartNeutralE_%d", iE), Form("truthjetR==%d && truthjetE>=%d  && truthjetE<%d  && truthjetEta>%f && truthjetEta<%f && jetParts_truth > %d", int(Rvals[Rvalue] * 10), int(JetEBorders[iE]), int(JetEBorders[iE+1]),etaMin+Rvals[Rvalue],etaMax-Rvals[Rvalue], constMin), "goff"); 
        TruthjetTree->Draw(Form("truthjetChargedEnergyFrac>>hPartChargedE_%d", iE), Form("truthjetR==%d && truthjetE>=%d  && truthjetE<%d  && truthjetEta>%f && truthjetEta<%f && jetParts_truth > %d", int(Rvals[Rvalue] * 10), int(JetEBorders[iE]), int(JetEBorders[iE+1]),etaMin+Rvals[Rvalue],etaMax-Rvals[Rvalue], constMin), "goff"); 

        TruthjetTree->Draw(Form("truthjetECALEnergyFrac>>hPartPhotonElectronE_%d", iE), Form("truthjetR==%d && truthjetE>=%d  && truthjetE<%d  && truthjetEta>%f && truthjetEta<%f && jetParts_truth > %d", int(Rvals[Rvalue] * 10), int(JetEBorders[iE]), int(JetEBorders[iE+1]),etaMin+Rvals[Rvalue],etaMax-Rvals[Rvalue], constMin), "goff"); 


        if(hDetHCALE[iE]->GetEntries()==0) hDetHCALE[iE]->Fill(-2);
        if(hDetECALE[iE]->GetEntries()==0) hDetECALE[iE]->Fill(-2);
        if(hPartNeutralE[iE]->GetEntries()==0) hPartNeutralE[iE]->Fill(-2);
        if(hPartChargedE[iE]->GetEntries()==0) hPartChargedE[iE]->Fill(-2);
        if(hPartPhotonElectronE[iE]->GetEntries()==0) hPartPhotonElectronE[iE]->Fill(-2);

        hDetHCALE[iE]->Scale(normalizations[NormValue]);
        hDetECALE[iE]->Scale(normalizations[NormValue]);
        hPartNeutralE[iE]->Scale(normalizations[NormValue]);
        hPartChargedE[iE]->Scale(normalizations[NormValue]);
        hPartPhotonElectronE[iE]->Scale(normalizations[NormValue]);

        hDetHCALE[iE]->GetXaxis()->SetTitle("E_{HCAL}/(E_{jet}^{det}");
        hDetECALE[iE]->GetXaxis()->SetTitle("E_{ECAL}/(E_{jet}^{det}");
        hPartNeutralE[iE]->GetXaxis()->SetTitle("E_{neutral}/E_{jet}^{part}");
        hPartChargedE[iE]->GetXaxis()->SetTitle("E_{charged}/E_{jet}^{part}");
        hPartPhotonElectronE[iE]->GetXaxis()->SetTitle("E_{#gamma,e}/E_{jet}^{part}");

    }

    fout->Write();
    fout->Close();
}