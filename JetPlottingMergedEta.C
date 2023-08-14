/* Macro for plotting some histograms from jet trees 
    2.21064
    0.0669805
    0.00182628
    0.000139462
    0.0000225822
*/

double GetMedian( TH1D * h);

void JetPlottingMergedEta(int NormValue=0,int Rvalue=0)
{   
    const int nNorm = 8;
    const Float_t normalizations[nNorm] = {0.0610658,0.00716477,0.000557627,0.000107816,4.31694e-05,9.62255e-06,1.24904e-06,5.99517e-08};

    //int NormValue = 0;

    TFile *jetFile;
    if(NormValue==0)jetFile = TFile::Open("/home/lmh/alice/fromPuhti/July2023/ptmin10/data/ptmin10_20230722_pythia8_JetJet_5-10GeV/Merged.root");
    if(NormValue==1)jetFile = TFile::Open("/home/lmh/alice/fromPuhti/July2023/ptmin10/data/ptmin10_20230722_pythia8_JetJet_10-20GeV/Merged.root");
    if(NormValue==2)jetFile = TFile::Open("/home/lmh/alice/fromPuhti/July2023/ptmin10/data/ptmin10_20230722_pythia8_JetJet_20-30GeV/Merged.root");
    if(NormValue==3)jetFile = TFile::Open("/home/lmh/alice/fromPuhti/July2023/ptmin10/data/ptmin10_20230722_pythia8_JetJet_30-40GeV/Merged.root");
    if(NormValue==4)jetFile = TFile::Open("/home/lmh/alice/fromPuhti/July2023/ptmin10/data/ptmin10_20230722_pythia8_JetJet_40-60GeV/Merged.root");
    if(NormValue==5)jetFile = TFile::Open("/home/lmh/alice/fromPuhti/July2023/ptmin10/data/ptmin10_20230722_pythia8_JetJet_60-100GeV/Merged.root");
    if(NormValue==6)jetFile = TFile::Open("/home/lmh/alice/fromPuhti/July2023/ptmin10/data/ptmin10_20230722_pythia8_JetJet_100-200GeV/Merged.root");
    if(NormValue==7)jetFile = TFile::Open("/home/lmh/alice/fromPuhti/July2023/ptmin10/data/ptmin10_20230722_pythia8_JetJet_200-GeV/Merged.root");

    const Int_t nR = 3; //5
    const Float_t Rvals[nR] = {0.2, 0.4, 0.6};//{0.2, 0.3, 0.4, 0.5, 0.6}; // Cone radii
    //int Rvalue = 0; // choose the index of the jet R you want to draw the main histos for !
    int constMin = 0; //min number of constituents in matchd jet 1
    
    TFile *fout;
    if(NormValue==0)fout = new TFile(Form("Data20230728/JESEta/20230722_pythia8_JetJet_5-10GeV_Merged_OutputR%d.root", int(Rvals[Rvalue] * 10)), "RECREATE");
    if(NormValue==1)fout = new TFile(Form("Data20230728/JESEta/20230722_pythia8_JetJet_10-20GeV_Merged_OutputR%d.root", int(Rvals[Rvalue] * 10)), "RECREATE");
    if(NormValue==2)fout = new TFile(Form("Data20230728/JESEta/20230722_pythia8_JetJet_20-30GeV_Merged_OutputR%d.root", int(Rvals[Rvalue] * 10)), "RECREATE");
    if(NormValue==3)fout = new TFile(Form("Data20230728/JESEta/20230722_pythia8_JetJet_30-40GeV_Merged_OutputR%d.root", int(Rvals[Rvalue] * 10)), "RECREATE");
    if(NormValue==4)fout = new TFile(Form("Data20230728/JESEta/20230722_pythia8_JetJet_40-60GeV_Merged_OutputR%d.root", int(Rvals[Rvalue] * 10)), "RECREATE");
    if(NormValue==5)fout = new TFile(Form("Data20230728/JESEta/20230722_pythia8_JetJet_60-100GeV_Merged_OutputR%d.root", int(Rvals[Rvalue] * 10)), "RECREATE");
    if(NormValue==6)fout = new TFile(Form("Data20230728/JESEta/20230722_pythia8_JetJet_100-200GeV_Merged_OutputR%d.root", int(Rvals[Rvalue] * 10)), "RECREATE");
    if(NormValue==7)fout = new TFile(Form("Data20230728/JESEta/20230722_pythia8_JetJet_200-GeV_Merged_OutputR%d.root", int(Rvals[Rvalue] * 10)), "RECREATE");

    TTree *jetTree = (TTree *)jetFile->Get("jetTree");
    TTree *TruthjetTree = (TTree *)jetFile->Get("truthjetTree");
    //TTree *PDGTree = (TTree *)jetFile->Get("inPDGTree");
    const int nCol = 10;
    const int gcolors[nCol]={1,2,6,4,7,1,2,4,6,7};
    const int gmarkers[nCol]={4,8,25,21,8,21,25,4,8,21};

    const Float_t etaMin = 3.4; const Float_t etaMax = 5.5; //changed in June23, 5.5 

    const Int_t nPtBins = 7; 
    const Int_t nEBins  = 6;
    const double JetPtBorders[nPtBins] = {0.0, 20.0, 40.0, 60.0, 80.0, 100.0, 400.0}; 
    //const double JetPtBorders[nPtBins] = {0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0, 400.0}; 
    const double JetEBorders[nEBins] = {0.0, 400.0, 800.0, 1200.0, 1600.0, 2000.0}; 

    const Int_t nEtaBins  = 4;//10;
    //const double JetPtBorders[nPtBins] = {2.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0, 70.0};
    //const double JetPtBorders[nPtBins] = {2.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 80.0, 100.0, 400.0}; //setting new bins for largest pt binned pythia data
    //const double JetPtBorders[nPtBins] = {0.0, 20.0, 40.0, 60.0, 80.0, 100.0, 400.0}; 
    //const double JetPtBorders[nPtBins] = {0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0, 400.0}; 
    //const double JetEBorders[nEBins] = {0.0, 400.0, 800.0, 1200.0, 1600.0, 2000.0}; 
    //const double JetEBorders[nEBins] = {0.0, 200.0, 400.0, 600.0, 800.0, 1000.0, 1200.0, 1400.0, 1600.0, 1800.0, 2000.0}; 
    const double EtaBinBorders[nEtaBins] = {3.8, 4.2, 4.6, 5.1};//{3.6, 3.8, 4.0, 4.2, 4.4, 4.6, 4.8, 5.0, 5.2, 5.4}; //[3.8, 4.2, 4.6, 5.1]

    double medianspT[nEtaBins];
    double meanspT[nEtaBins];
    double mediansE[nEtaBins];
    double meansE[nEtaBins];
    double SDpT[nEtaBins];
    double SDE[nEtaBins];

    double EtamedianspT[nEtaBins][nPtBins];
    double EtameanspT[nEtaBins][nPtBins];
    double EtamediansE[nEtaBins][nEBins];
    double EtameansE[nEtaBins][nEBins];
    double EtaSDpT[nEtaBins][nPtBins];
    double EtaSDE[nEtaBins][nEBins];

    //histograms

    TH1D *hjetRatiopT_Eta[nEtaBins][nPtBins];
    TH1D *hjetRatioE_Eta[nEtaBins][nEBins];
    TH1D *hEtaMedianpT[nEtaBins], *hEtaMedianE[nEtaBins];
    TH1D *hEtaMeanpT[nEtaBins], *hEtaMeanE[nEtaBins];
    TH1D *hEtaSDpT[nEtaBins], *hEtaSDE[nEtaBins];

    TH2D *hRespMatrix_pT = new TH2D("hRespMatrix_pT", Form("#Delta p_{T}, R=%0.1f;#eta_{jet};#Delta p_{T}", Rvals[Rvalue]), 150, etaMin, etaMax, 50, -1.0, 1.0);
    TH2D *hRespMatrix_E = new TH2D("hRespMatrix_E", Form("#Delta E, R=%0.1f;#eta_{jet};#Delta E", Rvals[Rvalue]), 150, etaMin, etaMax, 50, -1.0, 1.0);

    TH2D *hRespMatrix_pT_many[nPtBins];
    TH2D *hRespMatrix_E_many[nEBins];

    jetTree->Draw("(jetpT-jetpT_match)/jetpT_match:jetEta>>hRespMatrix_pT", Form("jetR==%d && jetEta>%f && jetEta<%f && jetEta_match>%f && jetEta_match<%f && jet_distmatch<%f && jetParts_match > %d", int(Rvals[Rvalue] * 10), etaMin+Rvals[Rvalue], etaMax-Rvals[Rvalue], etaMin+Rvals[Rvalue], etaMax-Rvals[Rvalue], Rvals[Rvalue]*0.6, constMin), "goff");
    jetTree->Draw("(jetE-jetE_match)/jetE_match:jetEta>>hRespMatrix_E", Form("jetR==%d && jetEta>%f && jetEta<%f && jetEta_match>%f && jetEta_match<%f && jet_distmatch<%f && jetParts_match > %d", int(Rvals[Rvalue] * 10), etaMin+Rvals[Rvalue], etaMax-Rvals[Rvalue], etaMin+Rvals[Rvalue], etaMax-Rvals[Rvalue], Rvals[Rvalue]*0.6, constMin), "goff");

    hRespMatrix_E->Scale(normalizations[NormValue]/500.0);
    //hRespMatrix_E->Scale(1./hRespMatrix_E->GetEntries(), "width");
    
    hRespMatrix_pT->Scale(normalizations[NormValue]/500.0);
    //hRespMatrix_pT->Scale(1./hRespMatrix_pT->GetEntries(), "width");

    for (int iE = 0; iE < nEtaBins-1; ++iE) //eta loop
    {
        hEtaMedianpT[iE] = new TH1D(Form("hEtaMedianpT_%d", iE), Form("Mean and median of #Deltap_{T} distribution, eta %0.1f - %0.1f, R=%0.1f", EtaBinBorders[iE], EtaBinBorders[iE+1], Rvals[Rvalue]),nPtBins-1, JetPtBorders);
        hEtaMedianE[iE] = new TH1D(Form("hEtaMedianE_%d", iE), Form("Mean and median of #DeltaE distribution, #eta_{jet} %0.1f - %0.1f, R=%0.1f", EtaBinBorders[iE], EtaBinBorders[iE+1], Rvals[Rvalue]),nEBins-1, JetEBorders);
        hEtaMeanpT[iE] = new TH1D(Form("hEtaMeanpT_%d", iE), Form("Mean and median of #Deltap_{T} distribution, #eta_{jet} %0.1f - %0.1f, R=%0.1f", EtaBinBorders[iE], EtaBinBorders[iE+1], Rvals[Rvalue]),nPtBins-1, JetPtBorders); 
        hEtaMeanE[iE] = new TH1D(Form("hEtaMeanE_%d", iE), Form("Mean and median of #DeltaE distribution, #eta_{jet} %0.1f - %0.1f, R=%0.1f", EtaBinBorders[iE], EtaBinBorders[iE+1], Rvals[Rvalue]),nEBins-1, JetEBorders); 
        hEtaSDpT[iE] = new TH1D(Form("hEtaSDpT_%d", iE), Form("Standard deviation of #Deltap_{T} distribution, #eta_{jet} %0.1f - %0.1f, R=%0.1f", EtaBinBorders[iE], EtaBinBorders[iE+1], Rvals[Rvalue]),nPtBins-1, JetPtBorders);
        hEtaSDE[iE] = new TH1D(Form("hEtaSDE_%d", iE), Form("Standard deviation of #DeltaE distribution, #eta_{jet} %0.1f - %0.1f, R=%0.1f", EtaBinBorders[iE], EtaBinBorders[iE+1], Rvals[Rvalue]),nEBins-1, JetEBorders);  
    }

    for (int ipt = 0; ipt < nPtBins-1; ++ipt) {
        hRespMatrix_pT_many[ipt] = new TH2D(Form("hRespMatrix_pT%d", ipt), Form("#Delta p_{T}, p_{T}: %d - %d GeV/c, R=%0.1f;#eta_{jet};#Delta p_{T}",int(JetPtBorders[ipt]), int(JetPtBorders[ipt+1]), Rvals[Rvalue]), 150, etaMin, etaMax, 50, -1.0, 1.0);
        jetTree->Draw(Form("(jetpT-jetpT_match)/jetpT_match:jetEta>>hRespMatrix_pT%d", ipt), Form("jetR==%d && jetpT_match>=%d && jetpT_match<%d && jetEta>%f && jetEta<%f && jetEta_match>%f && jetEta_match<%f && jet_distmatch<%f && jetParts_match > %d", int(Rvals[Rvalue] * 10), int(JetPtBorders[ipt]), int(JetPtBorders[ipt+1]), etaMin+Rvals[Rvalue], etaMax-Rvals[Rvalue], etaMin+Rvals[Rvalue], etaMax-Rvals[Rvalue], Rvals[Rvalue]*0.6, constMin), "goff");
        hRespMatrix_pT_many[ipt]->Scale(normalizations[NormValue]/500.0);

    }
    for (int iE = 0; iE < nEBins-1; ++iE){
        hRespMatrix_E_many[iE] = new TH2D(Form("hRespMatrix_E%d", iE), Form("#Delta E, E: %d - %d GeV, R=%0.1f;#eta_{jet};#Delta E",  int(JetEBorders[iE]), int(JetEBorders[iE+1]), Rvals[Rvalue]), 150, etaMin, etaMax, 50, -1.0, 1.0);
        jetTree->Draw(Form("(jetE-jetE_match)/jetE_match:jetEta>>hRespMatrix_E%d", iE), Form("jetR==%d && jetE_match>=%d && jetE_match<%d && jetEta>%f && jetEta<%f && jetEta_match>%f && jetEta_match<%f && jet_distmatch<%f && jetParts_match > %d", int(Rvals[Rvalue] * 10), int(JetEBorders[iE]), int(JetEBorders[iE+1]), etaMin+Rvals[Rvalue], etaMax-Rvals[Rvalue], etaMin+Rvals[Rvalue], etaMax-Rvals[Rvalue], Rvals[Rvalue]*0.6, constMin), "goff");
        hRespMatrix_E_many[iE]->Scale(normalizations[NormValue]/500.0);
    }




    for (int iEta = 0; iEta < nEtaBins-1; ++iEta) //eta loop
    {
        //pt loop
        for (int ipt = 0; ipt < nPtBins-1; ++ipt) 
        {
            hjetRatiopT_Eta[iEta][ipt] = new TH1D(Form("hjetRatiopT_Eta_%d_%d", iEta, ipt), Form("Jet-by-jet #Deltap_{T} distribution, #eta: %0.1f - %0.1f, p_{T}: %d - %d GeV/c", EtaBinBorders[iEta], EtaBinBorders[iEta+1], int(JetPtBorders[ipt]), int(JetPtBorders[ipt+1])), 50, -1.0, 1.0);
            jetTree->Draw(Form("(jetpT-jetpT_match)/jetpT_match>>hjetRatiopT_Eta_%d_%d", iEta, ipt), Form("jetR==%d && jetpT_match>=%d && jetpT_match<%d && jetEta>%f && jetEta<%f && jetEta_match>%f && jetEta_match<%f && jet_distmatch<%f && jetParts_match > %d", int(Rvals[Rvalue] * 10), int(JetPtBorders[ipt]), int(JetPtBorders[ipt+1]),EtaBinBorders[iEta], EtaBinBorders[iEta+1],EtaBinBorders[iEta], EtaBinBorders[iEta+1], Rvals[Rvalue]*0.6, constMin), "goff");
            if(hjetRatiopT_Eta[iEta][ipt]->GetEntries()==0) hjetRatiopT_Eta[iEta][ipt]->Fill(-2);
            hjetRatiopT_Eta[iEta][ipt]->Scale(normalizations[NormValue]/500.0);

            EtamedianspT[iEta][ipt]=GetMedian(hjetRatiopT_Eta[iEta][ipt]);
            EtameanspT[iEta][ipt]=hjetRatiopT_Eta[iEta][ipt]->GetMean();
            EtaSDpT[iEta][ipt]=hjetRatiopT_Eta[iEta][ipt]->GetStdDev();

            hEtaMedianpT[iEta]->Fill(JetPtBorders[ipt+1]-((JetPtBorders[ipt+1]-JetPtBorders[ipt])/2.0),EtamedianspT[iEta][ipt]);
            hEtaMedianpT[iEta]->SetBinError(ipt+1,hjetRatiopT_Eta[iEta][ipt]->GetMeanError());

            hEtaMeanpT[iEta]->Fill(JetPtBorders[ipt+1]-((JetPtBorders[ipt+1]-JetPtBorders[ipt])/2.0),EtameanspT[iEta][ipt]);
            hEtaMeanpT[iEta]->SetBinError(ipt+1,hjetRatiopT_Eta[iEta][ipt]->GetMeanError());
            hEtaSDpT[iEta]->Fill(JetPtBorders[ipt+1]-((JetPtBorders[ipt+1]-JetPtBorders[ipt])/2.0),EtaSDpT[iEta][ipt]);
            hEtaSDpT[iEta]->SetBinError(ipt+1, hjetRatiopT_Eta[iEta][ipt]->GetStdDevError());

            //hjetRatiopT_Eta[iEta][ipt]->Scale(1./hjetRatiopT_Eta[iEta][ipt]->GetEntries(), "width");
            //cout << "In jet pT range " << JetPtBorders[ipt] << " - " << JetPtBorders[ipt+1] << "GeV/c, number of jets accepted = " << hjetRatiopT_Eta[iEta][ipt]->GetEntries() << endl;

            hjetRatiopT_Eta[iEta][ipt]->GetXaxis()->SetTitle("(p_{T}^{det}-p_{T}^{part})/p_{T}^{part}");
            hjetRatiopT_Eta[iEta][ipt]->GetYaxis()->SetTitle("probability/bin");
        }

        for (int iE = 0; iE < nEBins-1; ++iE)
        {
            hjetRatioE_Eta[iEta][iE] = new TH1D(Form("hjetRatioE_Eta_%d_%d", iEta, iE), Form("Jet-by-jet #Delta E distribution, #eta: %0.1f - %0.1f, E: %d - %d GeV", EtaBinBorders[iEta], EtaBinBorders[iEta+1], int(JetEBorders[iE]), int(JetEBorders[iE+1])), 50, -1.0, 1.0);
            jetTree->Draw(Form("(jetE-jetE_match)/jetE_match>>hjetRatioE_Eta_%d_%d", iEta, iE), Form("jetR==%d && jetE_match>=%d && jetE_match<%d && jetEta>%f && jetEta<%f && jetEta_match>%f && jetEta_match<%f && jet_distmatch<%f && jetParts_match > %d", int(Rvals[Rvalue] * 10), int(JetEBorders[iE]), int(JetEBorders[iE+1]),EtaBinBorders[iEta], EtaBinBorders[iEta+1],EtaBinBorders[iEta], EtaBinBorders[iEta+1], Rvals[Rvalue]*0.6, constMin), "goff"); 
            if(hjetRatioE_Eta[iEta][iE]->GetEntries()==0) hjetRatioE_Eta[iEta][iE]->Fill(-2);
            hjetRatioE_Eta[iEta][iE]->Scale(normalizations[NormValue]/500.0);

            EtamediansE[iEta][iE]=GetMedian(hjetRatioE_Eta[iEta][iE]);
            EtameansE[iEta][iE]=hjetRatioE_Eta[iEta][iE]->GetMean();
            EtaSDE[iEta][iE]=hjetRatioE_Eta[iEta][iE]->GetStdDev();

            hEtaMedianE[iEta]->Fill(JetEBorders[iE+1]-((JetEBorders[iE+1]-JetEBorders[iE])/2.0),EtamediansE[iEta][iE]);
            hEtaMedianE[iEta]->SetBinError(iE+1,hjetRatioE_Eta[iEta][iE]->GetMeanError());

            hEtaMeanE[iEta]->Fill(JetEBorders[iE+1]-((JetEBorders[iE+1]-JetEBorders[iE])/2.0),EtameansE[iEta][iE]);
            hEtaMeanE[iEta]->SetBinError(iE+1,hjetRatioE_Eta[iEta][iE]->GetMeanError());
            hEtaSDE[iEta]->Fill(JetEBorders[iE+1]-((JetEBorders[iE+1]-JetEBorders[iE])/2.0),EtaSDE[iEta][iE]);
            hEtaSDE[iEta]->SetBinError(iE+1, hjetRatioE_Eta[iEta][iE]->GetStdDevError());

            //hjetRatioE_Eta[iEta][iE]->Scale(1./hjetRatioE_Eta[iEta][iE]->GetEntries(), "width");
            //cout << "In jet pT range " << JetEBorders[iE] << " - " << JetEBorders[iE+1] << "GeV/c, number of jets accepted = " << hjetRatioE_Eta[iEta][iE]->GetEntries() << endl;

            hjetRatioE_Eta[iEta][iE]->GetXaxis()->SetTitle("(E^{det}-E^{part})/E^{part}");
            hjetRatioE_Eta[iEta][iE]->GetYaxis()->SetTitle("probability/bin");
        }
    }
    
    fout->Write();
    fout->Close();

}


double GetMedian( TH1D * h) { 

    Double_t x, q;
    q = 0.5; // 0.5 for "median"
    h->ComputeIntegral(); // just a precaution
    h->GetQuantiles(1, &x, &q);
    std::cout << "median = " << x << std::endl;
    return x;
}


