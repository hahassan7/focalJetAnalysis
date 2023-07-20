/* Macro for plotting some histograms from jet trees 
    2.21064
    0.0669805
    0.00182628
    0.000139462
    0.0000225822
*/

double GetMedian( TH1D * h);

void JetPlottingMergedFileEta(int Rvalue = 0)
{   
    //const int nNorm = 5;
    //const Float_t normalizations[nNorm] = {2.21064, 0.0669805, 0.00182628, 0.000139462, 0.0000225822};
    //int NormValue = 1;

    const Int_t nR = 3; //5
    const Float_t Rvals[nR] = {0.2, 0.4, 0.6};

    TFile *jetFile = TFile::Open(Form("JetJetOutput/FINALAN/oneETAMergedR%d.root", int(Rvals[Rvalue] * 10)));

    TTree *jetTree = (TTree *)jetFile->Get("jetTree");
    TTree *TruthjetTree = (TTree *)jetFile->Get("truthjetTree");
    //TTree *PDGTree = (TTree *)jetFile->Get("inPDGTree");

    const int nCol = 10;
    const int gcolors[nCol]={1,2,6,4,7,1,2,4,6,7};
    const int gmarkers[nCol]={4,8,25,21,8,21,25,4,8,21};

    const Float_t etaMin = 3.4; const Float_t etaMax = 5.5; 

    const Int_t nPtBins = 7; 
    const Int_t nEBins  = 6;
    const Int_t nEtaBins  = 10;
    //const double JetPtBorders[nPtBins] = {2.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0, 70.0};
    //const double JetPtBorders[nPtBins] = {2.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 80.0, 100.0, 400.0}; //setting new bins for largest pt binned pythia data
    const double JetPtBorders[nPtBins] = {0.0, 20.0, 40.0, 60.0, 80.0, 100.0, 400.0}; 
    //const double JetPtBorders[nPtBins] = {0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0, 400.0}; 
    const double JetEBorders[nEBins] = {0.0, 400.0, 800.0, 1200.0, 1600.0, 2000.0}; 
    //const double JetEBorders[nEBins] = {0.0, 200.0, 400.0, 600.0, 800.0, 1000.0, 1200.0, 1400.0, 1600.0, 1800.0, 2000.0}; 
    const double EtaBinBorders[nEtaBins] = {3.6, 3.8, 4.0, 4.2, 4.4, 4.6, 4.8, 5.0, 5.2, 5.4}; //

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

    TFile *fout = new TFile(Form("JetJetOutput/FINALAN/oneETAMerged20230417_0-1000GeV_OutputR%d.root", int(Rvals[Rvalue] * 10)), "RECREATE");
    //histograms

    TH1D *hjetRatiopT_Eta[nEtaBins][nPtBins];
    TH1D *hjetRatioE_Eta[nEtaBins][nEBins];
    TH1D *hEtaMedianpT[nEtaBins], *hEtaMedianE[nEtaBins];
    TH1D *hEtaMeanpT[nEtaBins], *hEtaMeanE[nEtaBins];
    TH1D *hEtaSDpT[nEtaBins], *hEtaSDE[nEtaBins];

    //TH2D *hRespMatrix_pT = new TH2D("hRespMatrix_pT", Form("#Delta {p}_{T}, R=%0.1f;#eta_{jet};#Delta{p}_{T}", Rvals[Rvalue]), 150, 0, 150, 150, 0, 150);
    //TH2D *hRespMatrix_E = new TH2D("hRespMatrix_E", Form("#Delta {E}, R=%0.1f;#eta_{jet};#Delta{E}", Rvals[Rvalue]), 2000, 0, 4000, 2000, 0, 4000);

    for (int iE = 0; iE < nEtaBins-1; ++iE)
    {
        hEtaMedianpT[iE] = new TH1D(Form("hEtaMedianpT_%d", iE), Form("Mean and median of #Delta p_{T} distribution, eta %0.1f - %0.1f, R=%0.1f", EtaBinBorders[iE], EtaBinBorders[iE+1], Rvals[Rvalue]),nPtBins-1, JetPtBorders);
        hEtaMedianE[iE] = new TH1D(Form("hEtaMedianE_%d", iE), Form("Mean and median of #Delta E distribution, #eta_{jet} %0.1f - %0.1f, R=%0.1f", EtaBinBorders[iE], EtaBinBorders[iE+1], Rvals[Rvalue]),nEBins-1, JetEBorders);
        hEtaMeanpT[iE] = new TH1D(Form("hEtaMeanpT_%d", iE), Form("Mean and median of #Deltap_{T} distribution, #eta_{jet} %0.1f - %0.1f, R=%0.1f", EtaBinBorders[iE], EtaBinBorders[iE+1], Rvals[Rvalue]),nPtBins-1, JetPtBorders); 
        hEtaMeanE[iE] = new TH1D(Form("hEtaMeanE_%d", iE), Form("Mean and median of #DeltaE distribution, #eta_{jet} %0.1f - %0.1f, R=%0.1f", EtaBinBorders[iE], EtaBinBorders[iE+1], Rvals[Rvalue]),nEBins-1, JetEBorders); 
        hEtaSDpT[iE] = new TH1D(Form("hEtaSDpT_%d", iE), Form("Standard deviation of #Deltap_{T} distribution, #eta_{jet} %0.1f - %0.1f, R=%0.1f", EtaBinBorders[iE], EtaBinBorders[iE+1], Rvals[Rvalue]),nPtBins-1, JetPtBorders);
        hEtaSDE[iE] = new TH1D(Form("hEtaSDE_%d", iE), Form("Standard deviation of #DeltaE distribution, #eta_{jet} %0.1f - %0.1f, R=%0.1f", EtaBinBorders[iE], EtaBinBorders[iE+1], Rvals[Rvalue]),nEBins-1, JetEBorders);  
    }

    TH2D *hRespMatrix_pT = (TH2D*)((TH2D*)jetFile->Get("hRespMatrix_pT"))->Clone();
    TH2D *hRespMatrix_E  = (TH2D*)((TH2D*)jetFile->Get("hRespMatrix_E"))->Clone();

    TH2D *hRespMatrix_pT_many[nPtBins];
    TH2D *hRespMatrix_E_many[nEBins];

    for (int ipt = 0; ipt < nPtBins-1; ++ipt) {
        hRespMatrix_pT_many[ipt] = (TH2D*)((TH2D*)jetFile->Get(Form("hRespMatrix_pT%d", ipt)))->Clone();

        //hRespMatrix_pT_many[ipt]->Scale(hRespMatrix_pT_many[ipt]->GetBinWidth(1));
        hRespMatrix_pT_many[ipt]->Scale(1./hRespMatrix_pT_many[ipt]->Integral());

    }
    for (int iE = 0; iE < nEBins-1; ++iE){
        hRespMatrix_E_many[iE] = (TH2D*)((TH2D*)jetFile->Get(Form("hRespMatrix_E%d", iE)))->Clone();
        //hRespMatrix_E_many[iE]->Scale(hRespMatrix_E_many[iE]->GetBinWidth(1));
        hRespMatrix_E_many[iE]->Scale(1./hRespMatrix_E_many[iE]->Integral());
    }

    for (int iEta = 0; iEta < nEtaBins-1; ++iEta) //eta loop
    {
        //pt loop
        for (int ipt = 0; ipt < nPtBins-1; ++ipt) 
        {
            hjetRatiopT_Eta[iEta][ipt] = (TH1D*)((TH1D*)jetFile->Get(Form("hjetRatiopT_Eta_%d_%d", iEta, ipt)))->Clone();

            EtamedianspT[iEta][ipt]=GetMedian(hjetRatiopT_Eta[iEta][ipt]);
            EtameanspT[iEta][ipt]=hjetRatiopT_Eta[iEta][ipt]->GetMean();
            EtaSDpT[iEta][ipt]=hjetRatiopT_Eta[iEta][ipt]->GetStdDev();

            hEtaMedianpT[iEta]->Fill(JetPtBorders[ipt+1]-((JetPtBorders[ipt+1]-JetPtBorders[ipt])/2.0),EtamedianspT[iEta][ipt]);
            hEtaMedianpT[iEta]->SetBinError(ipt+1,hjetRatiopT_Eta[iEta][ipt]->GetMeanError());

            hEtaMeanpT[iEta]->Fill(JetPtBorders[ipt+1]-((JetPtBorders[ipt+1]-JetPtBorders[ipt])/2.0),EtameanspT[iEta][ipt]);
            hEtaMeanpT[iEta]->SetBinError(ipt+1,hjetRatiopT_Eta[iEta][ipt]->GetMeanError());
            hEtaSDpT[iEta]->Fill(JetPtBorders[ipt+1]-((JetPtBorders[ipt+1]-JetPtBorders[ipt])/2.0),EtaSDpT[iEta][ipt]);
            hEtaSDpT[iEta]->SetBinError(ipt+1, hjetRatiopT_Eta[iEta][ipt]->GetStdDevError());

            //hjetRatiopT_Eta[iE][ipt]->Scale(1./hjetRatiopT_Eta[iE][ipt]->GetEntries());
            cout << "In jet pT range " << JetPtBorders[ipt] << " - " << JetPtBorders[ipt+1] << "GeV/c, number of jets accepted = " << hjetRatiopT_Eta[iEta][ipt]->GetEntries() << endl;

            hjetRatiopT_Eta[iEta][ipt]->Scale(1./hjetRatiopT_Eta[iEta][ipt]->GetEntries());
            hjetRatiopT_Eta[iEta][ipt]->GetXaxis()->SetTitle("(p_{T}^{det}-p_{T}^{part})/p_{T}^{part}");
            hjetRatiopT_Eta[iEta][ipt]->GetYaxis()->SetTitle("probability/bin");

            hEtaMedianpT[iEta]->GetXaxis()->SetTitle("p_{T}^{part} (GeV/c)");
            hEtaMedianpT[iEta]->GetYaxis()->SetTitle("Mean or median");

            hEtaMeanpT[iEta]->GetXaxis()->SetTitle("p_{T}^{part} (GeV/c)");
            hEtaMeanpT[iEta]->GetYaxis()->SetTitle("Mean or median");
         
            hEtaSDpT[iEta]->GetXaxis()->SetTitle("p_{T}^{part} (GeV/c)");
            hEtaSDpT[iEta]->GetYaxis()->SetTitle("Standard deviation");
        }

        //E loop
        for (int iE = 0; iE < nEBins-1; ++iE)
        {
            
            hjetRatioE_Eta[iEta][iE] = (TH1D*)((TH1D*)jetFile->Get(Form("hjetRatioE_Eta_%d_%d", iEta, iE)))->Clone();

            EtamediansE[iEta][iE]=GetMedian(hjetRatioE_Eta[iEta][iE]);
            EtameansE[iEta][iE]=hjetRatioE_Eta[iEta][iE]->GetMean();
            EtaSDE[iEta][iE]=hjetRatioE_Eta[iEta][iE]->GetStdDev();

            hEtaMedianE[iEta]->Fill(JetEBorders[iE+1]-((JetEBorders[iE+1]-JetEBorders[iE])/2.0),EtamediansE[iEta][iE]);
            hEtaMedianE[iEta]->SetBinError(iE+1,hjetRatioE_Eta[iEta][iE]->GetMeanError());

            hEtaMeanE[iEta]->Fill(JetEBorders[iE+1]-((JetEBorders[iE+1]-JetEBorders[iE])/2.0),EtameansE[iEta][iE]);
            hEtaMeanE[iEta]->SetBinError(iE+1,hjetRatioE_Eta[iEta][iE]->GetMeanError());
            hEtaSDE[iEta]->Fill(JetEBorders[iE+1]-((JetEBorders[iE+1]-JetEBorders[iE])/2.0),EtaSDE[iEta][iE]);
            hEtaSDE[iEta]->SetBinError(iE+1, hjetRatioE_Eta[iEta][iE]->GetStdDevError());

            hjetRatioE_Eta[iEta][iE]->Scale(1./hjetRatioE_Eta[iEta][iE]->GetEntries());

            hjetRatioE_Eta[iEta][iE]->GetXaxis()->SetTitle("(E^{det}-E^{part})/E^{part}");
            hjetRatioE_Eta[iEta][iE]->GetYaxis()->SetTitle("probability/bin");

            hEtaMedianE[iEta]->GetXaxis()->SetTitle("E^{part} (GeV)");
            hEtaMedianE[iEta]->GetYaxis()->SetTitle("Mean or median");

            hEtaMeanE[iEta]->GetXaxis()->SetTitle("E^{part} (GeV)");
            hEtaMeanE[iEta]->GetYaxis()->SetTitle("Mean or median");
         
            hEtaSDE[iEta]->GetXaxis()->SetTitle("E^{part} (GeV)");
            hEtaSDE[iEta]->GetYaxis()->SetTitle("Standard deviation");
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


