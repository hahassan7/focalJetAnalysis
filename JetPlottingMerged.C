/* Macro for plotting some histograms from jet trees 
    2.21064
    0.0669805
    0.00182628
    0.000139462
    0.0000225822
*/

double GetMedian( TH1D * h);

void JetPlottingMerged(int NormValue=0,int Rvalue=0)
{   
    const int nNorm = 8;
    const Float_t normalizations[nNorm] = {0.0611191, 0.00717001, 0.000558759, 0.000107936, 4.32163e-05, 9.57109e-06, 1.24606e-06, 6.01382e-08};//{0.0610658,0.00716477,0.000557627,0.000107816,4.31694e-05,9.62255e-06,1.24904e-06,5.99517e-08};

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
    int constMin = 0; //min number of constituents in matched jet

    TFile *fout;
    if(NormValue==0)fout = new TFile(Form("Data20230728/JES/20230728_pythia8_JetJet_5-10GeV_Merged_TotalOutputR%d.root", int(Rvals[Rvalue] * 10)), "RECREATE");
    if(NormValue==1)fout = new TFile(Form("Data20230728/JES/20230728_pythia8_JetJet_10-20GeV_Merged_TotalOutputR%d.root", int(Rvals[Rvalue] * 10)), "RECREATE");
    if(NormValue==2)fout = new TFile(Form("Data20230728/JES/20230728_pythia8_JetJet_20-30GeV_Merged_TotalOutputR%d.root", int(Rvals[Rvalue] * 10)), "RECREATE");
    if(NormValue==3)fout = new TFile(Form("Data20230728/JES/20230728_pythia8_JetJet_30-40GeV_Merged_TotalOutputR%d.root", int(Rvals[Rvalue] * 10)), "RECREATE");
    if(NormValue==4)fout = new TFile(Form("Data20230728/JES/20230728_pythia8_JetJet_40-60GeV_Merged_TotalOutputR%d.root", int(Rvals[Rvalue] * 10)), "RECREATE");
    if(NormValue==5)fout = new TFile(Form("Data20230728/JES/20230728_pythia8_JetJet_60-100GeV_Merged_TotalOutputR%d.root", int(Rvals[Rvalue] * 10)), "RECREATE");
    if(NormValue==6)fout = new TFile(Form("Data20230728/JES/20230728_pythia8_JetJet_100-200GeV_Merged_TotalOutputR%d.root", int(Rvals[Rvalue] * 10)), "RECREATE");
    if(NormValue==7)fout = new TFile(Form("Data20230728/JES/20230728_pythia8_JetJet_200-GeV_Merged_TotalOutputR%d.root", int(Rvals[Rvalue] * 10)), "RECREATE");

    TTree *jetTree = (TTree *)jetFile->Get("jetTree");
    TTree *TruthjetTree = (TTree *)jetFile->Get("truthjetTree");
    //TTree *PDGTree = (TTree *)jetFile->Get("inPDGTree");

    const int nCol = 10;
    const int gcolors[nCol]={1,2,6,4,7,1,2,4,6,7};
    const int gmarkers[nCol]={4,8,25,21,8,21,25,4,8,21};

    const Float_t etaMin = 3.4; const Float_t etaMax = 5.5; //changed in June23, 5.5 

    //Set the pT and energy bins for the projections, JES and JER. 
    const Int_t nPtBins = 14; 
    const Int_t nEBins  = 15;
    //const double JetPtBorders[nPtBins] = {5.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 80.0, 100.0, 400.0}; // nPtBins = 10
    const double JetPtBorders[nPtBins] = {5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0, 50.0, 60.0, 70.0, 80.0, 100.0, 150.0}; //nPtBins = 14
    const double JetEBorders[nEBins] = {100.0, 200.0, 300.0, 400.0, 500.0, 600.0, 700.0, 800.0, 900.0, 1000.0, 1250.0, 1500.0, 1750.0, 2000.0, 3000.0}; //const Int_t nEBins  = 15;
    double medianspT[nPtBins-1];
    double meanspT[nPtBins-1];
    double SDpT[nPtBins-1];
    double mediansE[nEBins-1];
    double meansE[nEBins-1];
    double SDE[nEBins-1];
    //eta binned
    const Int_t nEtaBins = 3;
    const Float_t EtaBinBorders[nEtaBins] = {4.0,4.5,4.9};//{3.8, 4.5, 5.1}; //CHANGED THESE {3.7, 4.5, 5.5}

    double EtamedianspT[nEtaBins-1][nPtBins-1];
    double EtameanspT[nEtaBins-1][nPtBins-1];
    double EtaSDpT[nEtaBins-1][nPtBins-1];

    double EtamediansE[nEtaBins-1][nEBins-1];
    double EtameansE[nEtaBins-1][nEBins-1];
    double EtaSDE[nEtaBins-1][nEBins-1];

    //histograms

    TH2D *hRespMatrix_pT_Eta[nEtaBins-1];
    TH2D *hRespMatrix_E_Eta[nEtaBins-1];
    TH1D *hjetRatiopT_Eta[nEtaBins-1][nPtBins-1];
    TH1D *hjetRatioE_Eta[nEtaBins-1][nEBins-1];
    TH1D *hEtaMedianpT[nEtaBins-1], *hEtaMedianE[nEtaBins-1];
    TH1D *hEtaMeanpT[nEtaBins-1], *hEtaMeanE[nEtaBins-1];
    TH1D *hEtaSDpT[nEtaBins-1], *hEtaSDE[nEtaBins-1];

    //Not-eta-binned histograms
    TH1F *hjetpT[nR], *hjetE[nR];
    TH1F *hTruthjetpT[nR], *hTruthjetE[nR];

    TH1D *hjetRatiopT[nPtBins-1], *hjetRatioE[nEBins-1];
    TH1D *hMedianpT, *hMedianE;
    TH1D *hMeanpT, *hMeanE;
    TH1D *hSDpT, *hSDE;

    hMedianpT = new TH1D("hMedianpT", "Mean and median of #Deltap_{T} distribution",nPtBins-1, JetPtBorders);
    hMedianE = new TH1D("hMedianE", "Mean and median of #DeltaE distribution",nEBins-1, JetEBorders);
    hMeanpT = new TH1D("hMeanpT", "Mean and median of #Deltap_{T} distribution",nPtBins-1, JetPtBorders); 
    hMeanE = new TH1D("hMeanE", "Mean and median of #DeltaE distribution",nEBins-1, JetEBorders); 
    hSDpT = new TH1D("hSDpT", "Standard deviation of #Deltap_{T} distribution",nPtBins-1, JetPtBorders);
    hSDE = new TH1D("hSDE", "Standard deviation of #DeltaE distribution",nEBins-1, JetEBorders);


    //begin filling the histograms

    for (int iR = 0; iR < nR; iR++)
    {
        hjetpT[iR] = new TH1F(Form("hjetpT_%d", iR), Form("Jet pT spectrum R=%0.1f, jet #eta %0.1f - %0.1f", Rvals[iR], etaMin+Rvals[iR], etaMax-Rvals[iR]), 150, 0, 150);
        hjetE[iR] = new TH1F(Form("hjetE_%d", iR), Form("Jet energy spectrum R=%0.1f, jet #eta %0.1f - %0.1f", Rvals[iR], etaMin+Rvals[iR], etaMax-Rvals[iR]), 2000, 0, 4000);

        hTruthjetpT[iR] = new TH1F(Form("hTruthjetpT_%d", iR), Form("Truth jet pT spectrum R=%0.1f", Rvals[iR]), 150, 0, 150);
        hTruthjetE[iR] = new TH1F(Form("hTruthjetE_%d", iR), Form("Truth jet energy spectrum R=%0.1f", Rvals[iR]), 2000, 0, 4000);

        jetTree->Draw(Form("jetpT>>hjetpT_%d", iR), Form("jetR==%d && jetEta>=%f && jetEta<%f", int(Rvals[iR] * 10), etaMin+Rvals[iR], etaMax-Rvals[iR]), "goff");
        hjetpT[iR]->Scale(normalizations[NormValue]/500.0);

        jetTree->Draw(Form("jetE>>hjetE_%d", iR), Form("jetR==%d && jetEta>=%f && jetEta<%f", int(Rvals[iR] * 10), etaMin+Rvals[iR], etaMax-Rvals[iR]), "goff");
        hjetE[iR]->Scale(normalizations[NormValue]/500.0);

        TruthjetTree->Draw(Form("truthjetpT>>hTruthjetpT_%d", iR), Form("truthjetR==%d && truthjetEta>=%f && truthjetEta<%f", int(Rvals[iR] * 10), etaMin+Rvals[iR], etaMax-Rvals[iR]), "goff");
        hTruthjetpT[iR]->Scale(normalizations[NormValue]/500.0);

        TruthjetTree->Draw(Form("truthjetE>>hTruthjetE_%d", iR), Form("truthjetR==%d && truthjetEta>=%f && truthjetEta<%f", int(Rvals[iR] * 10), etaMin+Rvals[iR], etaMax-Rvals[iR]), "goff");
        hTruthjetE[iR]->Scale(normalizations[NormValue]/500.0);
    }

    TH2D *hRespMatrix_pT = new TH2D("hRespMatrix_pT", Form("Detector response matrix, R=%0.1f;#it{p}_{T}^{det} (GeV/c);#it{p}_{T}^{part} (GeV/c)", Rvals[Rvalue]), 150, 0, 150, 150, 0, 150);
    TH2D *hRespMatrix_E = new TH2D("hRespMatrix_E", Form("Detector response matrix, R=%0.1f;#it{E}_{jet}^{det} (GeV);#it{E}_{jet}^{part} (GeV)", Rvals[Rvalue]), 2000, 0, 4000, 2000, 0, 4000);

    jetTree->Draw("jetE_match:jetE>>hRespMatrix_E", Form("jetR==%d && jetEta>=%f && jetEta<%f && jetEta_match>%f && jetEta_match<%f && jet_distmatch<%f && jetParts_match > %d", int(Rvals[Rvalue] * 10), etaMin+Rvals[Rvalue], etaMax-Rvals[Rvalue], etaMin+Rvals[Rvalue], etaMax-Rvals[Rvalue], Rvals[Rvalue]*0.6, constMin), "goff");
    jetTree->Draw("jetpT_match:jetpT>>hRespMatrix_pT", Form("jetR==%d && jetEta>=%f && jetEta<%f && jetEta_match>%f && jetEta_match<%f && jet_distmatch<%f && jetParts_match > %d", int(Rvals[Rvalue] * 10), etaMin+Rvals[Rvalue], etaMax-Rvals[Rvalue], etaMin+Rvals[Rvalue], etaMax-Rvals[Rvalue], Rvals[Rvalue]*0.6, constMin), "goff");
    
    hRespMatrix_pT->Scale(normalizations[NormValue]/500.0);
    hRespMatrix_E->Scale(normalizations[NormValue]/500.0);

    for (int iE = 0; iE < nEtaBins-1; ++iE) //eta loop
    {
        hRespMatrix_pT_Eta[iE] = new TH2D(Form("hRespMatrix_pT_Eta_%d", iE), Form("Detector response matrix, eta %0.1f - %0.1f, R=%0.1f;#it{p}_{T}^{det} (GeV/c);#it{p}_{T}^{part} (GeV/c)", EtaBinBorders[iE], EtaBinBorders[iE+1], Rvals[Rvalue]), 150, 0, 150, 150, 0, 150);
        jetTree->Draw(Form("jetpT_match:jetpT>>hRespMatrix_pT_Eta_%d", iE), Form("jetR==%d && jetEta>=%f && jetEta<%f && jetEta_match>%f && jetEta_match<%f && jet_distmatch<%f && jetParts_match > %d", int(Rvals[Rvalue] * 10), EtaBinBorders[iE], EtaBinBorders[iE+1], EtaBinBorders[iE], EtaBinBorders[iE+1], Rvals[Rvalue]*0.6, constMin), "goff");   
        hRespMatrix_pT_Eta[iE]->Scale(normalizations[NormValue]/500.0);

        hRespMatrix_E_Eta[iE] = new TH2D(Form("hRespMatrix_E_Eta_%d", iE), Form("Detector response matrix, eta %0.1f - %0.1f, R=%0.1f;#it{E}_{jet}^{det} (GeV);#it{E}_{jet}^{part} (GeV)", EtaBinBorders[iE], EtaBinBorders[iE+1], Rvals[Rvalue]), 2000, 0, 4000, 2000, 0, 4000);
        jetTree->Draw(Form("jetE_match:jetE>>hRespMatrix_E_Eta_%d", iE), Form("jetR==%d && jetEta>=%f && jetEta<%f && jetEta_match>%f && jetEta_match<%f && jet_distmatch<%f && jetParts_match > %d", int(Rvals[Rvalue] * 10), EtaBinBorders[iE], EtaBinBorders[iE+1], EtaBinBorders[iE], EtaBinBorders[iE+1], Rvals[Rvalue]*0.6, constMin), "goff"); 
        hRespMatrix_E_Eta[iE]->Scale(normalizations[NormValue]/500.0);

        hEtaMedianpT[iE] = new TH1D(Form("hEtaMedianpT_%d", iE), Form("Mean and median of #Deltap_{T} distribution, eta %0.1f - %0.1f, R=%0.1f", EtaBinBorders[iE], EtaBinBorders[iE+1], Rvals[Rvalue]),nPtBins-1, JetPtBorders);
        hEtaMedianE[iE] = new TH1D(Form("hEtaMedianE_%d", iE), Form("Mean and median of #DeltaE distribution, eta %0.1f - %0.1f, R=%0.1f", EtaBinBorders[iE], EtaBinBorders[iE+1], Rvals[Rvalue]),nEBins-1, JetEBorders);
        hEtaMeanpT[iE] = new TH1D(Form("hEtaMeanpT_%d", iE), Form("Mean and median of #Deltap_{T} distribution, eta %0.1f - %0.1f, R=%0.1f", EtaBinBorders[iE], EtaBinBorders[iE+1], Rvals[Rvalue]),nPtBins-1, JetPtBorders); 
        hEtaMeanE[iE] = new TH1D(Form("hEtaMeanE_%d", iE), Form("Mean and median of #DeltaE distribution, eta %0.1f - %0.1f, R=%0.1f", EtaBinBorders[iE], EtaBinBorders[iE+1], Rvals[Rvalue]),nEBins-1, JetEBorders); 
        hEtaSDpT[iE] = new TH1D(Form("hEtaSDpT_%d", iE), Form("Standard deviation of #Deltap_{T} distribution, eta %0.1f - %0.1f, R=%0.1f", EtaBinBorders[iE], EtaBinBorders[iE+1], Rvals[Rvalue]),nPtBins-1, JetPtBorders);
        hEtaSDE[iE] = new TH1D(Form("hEtaSDE_%d", iE), Form("Standard deviation of #DeltaE distribution, eta %0.1f - %0.1f, R=%0.1f", EtaBinBorders[iE], EtaBinBorders[iE+1], Rvals[Rvalue]),nEBins-1, JetEBorders);  
    }


    //pt loop
    for (int ipt = 0; ipt < nPtBins-1; ++ipt) 
    {
        hjetRatiopT[ipt] = new TH1D(Form("hjetRatiopT_%d", ipt), Form("Jet-by-jet #Deltap_{T} distribution, p_{T}: %d - %d GeV/c", int(JetPtBorders[ipt]), int(JetPtBorders[ipt+1])), 50, -1.0, 1.0);
        
        //Fill deltapt histos, no eta cut
        //Add here filling from 
        jetTree->Draw(Form("(jetpT-jetpT_match)/jetpT_match>>hjetRatiopT_%d", ipt), Form("jetR==%d && jetpT_match>=%d  && jetpT_match<%d  && jetEta>3.4+%f && jetEta<5.5-%f && jetEta_match>3.4+%f && jetEta_match<5.5-%f && jet_distmatch<%f && jetParts_match > %d", int(Rvals[Rvalue] * 10), int(JetPtBorders[ipt]), int(JetPtBorders[ipt+1]),Rvals[Rvalue],Rvals[Rvalue],Rvals[Rvalue],Rvals[Rvalue], Rvals[Rvalue]*0.6, constMin), "goff");
        //if(hjetRatiopT[ipt]->GetEntries()==0) hjetRatiopT[ipt]->Fill(-2);
        hjetRatiopT[ipt]->Scale(normalizations[NormValue]/500.0);

        medianspT[ipt]=GetMedian(hjetRatiopT[ipt]);
        meanspT[ipt]=hjetRatiopT[ipt]->GetMean();
        SDpT[ipt]=hjetRatiopT[ipt]->GetStdDev();

        hMedianpT->Fill(JetPtBorders[ipt+1]-((JetPtBorders[ipt+1]-JetPtBorders[ipt])/2.0),medianspT[ipt]);
        hMedianpT->SetBinError(ipt+1,hjetRatiopT[ipt]->GetMeanError());

        hMeanpT->Fill(JetPtBorders[ipt+1]-((JetPtBorders[ipt+1]-JetPtBorders[ipt])/2.0),meanspT[ipt]);
        hMeanpT->SetBinError(ipt+1,hjetRatiopT[ipt]->GetMeanError());
        hSDpT->Fill(JetPtBorders[ipt+1]-((JetPtBorders[ipt+1]-JetPtBorders[ipt])/2.0),SDpT[ipt]);
        hSDpT->SetBinError(ipt+1, hjetRatiopT[ipt]->GetStdDevError());

        cout << "In jet pT range " << JetPtBorders[ipt] << " - " << JetPtBorders[ipt+1] << "GeV/c, number of jets accepted = " << hjetRatiopT[ipt]->GetEntries() << endl;

        hjetRatiopT[ipt]->GetXaxis()->SetTitle("(p_{T}^{det}-p_{T}^{part})/p_{T}^{part}");
        hjetRatiopT[ipt]->GetYaxis()->SetTitle("probability");

        //same for eta cut
        for (int iE = 0; iE < nEtaBins-1; ++iE)
        {
            hjetRatiopT_Eta[iE][ipt] = new TH1D(Form("hjetRatiopT_Eta_%d_%d", iE, ipt), Form("Jet-by-jet #Deltap_{T} distribution, #eta: %0.1f - %0.1f (det jet), p_{T}: %d - %d GeV/c", EtaBinBorders[iE], EtaBinBorders[iE+1], int(JetPtBorders[ipt]), int(JetPtBorders[ipt+1])), 50, -1.0, 1.0);
            jetTree->Draw(Form("(jetpT-jetpT_match)/jetpT_match>>hjetRatiopT_Eta_%d_%d", iE, ipt), Form("jetR==%d && jetpT_match>=%d && jetpT_match<%d && jetEta>=%f && jetEta<%f && jetEta_match>%f && jetEta_match<%f && jet_distmatch<%f && jetParts_match > %d", int(Rvals[Rvalue] * 10), int(JetPtBorders[ipt]), int(JetPtBorders[ipt+1]),EtaBinBorders[iE], EtaBinBorders[iE+1],EtaBinBorders[iE], EtaBinBorders[iE+1], Rvals[Rvalue]*0.6, constMin), "goff"); 
            //if(hjetRatiopT_Eta[iE][ipt]->GetEntries()==0) hjetRatiopT_Eta[iE][ipt]->Fill(-2);
            hjetRatiopT_Eta[iE][ipt]->Scale(normalizations[NormValue]/500.0);

            EtamedianspT[iE][ipt]=GetMedian(hjetRatiopT_Eta[iE][ipt]);
            EtameanspT[iE][ipt]=hjetRatiopT_Eta[iE][ipt]->GetMean();
            EtaSDpT[iE][ipt]=hjetRatiopT_Eta[iE][ipt]->GetStdDev();

            hEtaMedianpT[iE]->Fill(JetPtBorders[ipt+1]-((JetPtBorders[ipt+1]-JetPtBorders[ipt])/2.0),EtamedianspT[iE][ipt]);
            hEtaMedianpT[iE]->SetBinError(ipt+1,hjetRatiopT_Eta[iE][ipt]->GetMeanError());

            hEtaMeanpT[iE]->Fill(JetPtBorders[ipt+1]-((JetPtBorders[ipt+1]-JetPtBorders[ipt])/2.0),EtameanspT[iE][ipt]);
            hEtaMeanpT[iE]->SetBinError(ipt+1,hjetRatiopT_Eta[iE][ipt]->GetMeanError());
            hEtaSDpT[iE]->Fill(JetPtBorders[ipt+1]-((JetPtBorders[ipt+1]-JetPtBorders[ipt])/2.0),EtaSDpT[iE][ipt]);
            hEtaSDpT[iE]->SetBinError(ipt+1, hjetRatiopT_Eta[iE][ipt]->GetStdDevError());

            //cout << "In jet pT range " << JetPtBorders[ipt] << " - " << JetPtBorders[ipt+1] << "GeV/c, number of jets accepted = " << hjetRatiopT_Eta[iE][ipt]->GetEntries() << endl;

            hjetRatiopT_Eta[iE][ipt]->GetXaxis()->SetTitle("(p_{T}^{det}-p_{T}^{part})/p_{T}^{part}");
            hjetRatiopT_Eta[iE][ipt]->GetYaxis()->SetTitle("probability");
        }       
    }

    //energy loop
    for (int ipt = 0; ipt < nEBins-1; ++ipt) 
    {
        hjetRatioE[ipt] = new TH1D(Form("hjetRatioE_%d", ipt), Form("Jet-by-jet #Delta E distribution, E: %d - %d GeV", int(JetEBorders[ipt]), int(JetEBorders[ipt+1])), 50, -1.0, 1.0);

        //Fill deltaE histos, no eta cut
        jetTree->Draw(Form("(jetE-jetE_match)/jetE_match>>hjetRatioE_%d", ipt), Form("jetR==%d && jetE_match>=%d  && jetE_match<%d  && jetEta>3.4+%f && jetEta<5.5-%f && jetEta_match>3.4+%f && jetEta_match<5.5-%f && jet_distmatch<%f && jetParts_match > %d", int(Rvals[Rvalue] * 10), int(JetEBorders[ipt]), int(JetEBorders[ipt+1]),Rvals[Rvalue],Rvals[Rvalue],Rvals[Rvalue],Rvals[Rvalue], Rvals[Rvalue]*0.6, constMin), "goff"); 
        //if(hjetRatioE[ipt]->GetEntries()==0) hjetRatioE[ipt]->Fill(-2);
        hjetRatioE[ipt]->Scale(normalizations[NormValue]/500.0);

        mediansE[ipt]=GetMedian(hjetRatioE[ipt]);
        meansE[ipt]=hjetRatioE[ipt]->GetMean();
        SDE[ipt]=hjetRatioE[ipt]->GetStdDev();

        hMedianE->Fill(JetEBorders[ipt+1]-((JetEBorders[ipt+1]-JetEBorders[ipt])/2),mediansE[ipt]);
        hMedianE->SetBinError(ipt+1,hjetRatioE[ipt]->GetMeanError());
        hMeanE->Fill(JetEBorders[ipt+1]-((JetEBorders[ipt+1]-JetEBorders[ipt])/2),meansE[ipt]);
        hMeanE->SetBinError(ipt+1,hjetRatioE[ipt]->GetMeanError());
        hSDE->Fill(JetEBorders[ipt+1]-((JetEBorders[ipt+1]-JetEBorders[ipt])/2),SDE[ipt]);
        hSDE->SetBinError(ipt+1, hjetRatioE[ipt]->GetStdDevError());

        cout << "In jet E range " << JetEBorders[ipt] << " - " << JetEBorders[ipt+1] << "GeV, number of jets accepted = " << hjetRatioE[ipt]->GetEntries() << endl;


        hjetRatioE[ipt]->GetXaxis()->SetTitle("(E^{det}-E^{part})/E^{part}");
        hjetRatioE[ipt]->GetYaxis()->SetTitle("probability");

        //same for eta cut
        for (int iE = 0; iE < nEtaBins-1; ++iE)
        {
            hjetRatioE_Eta[iE][ipt] = new TH1D(Form("hjetRatioE_Eta_%d_%d", iE, ipt), Form("Jet-by-jet #Delta E distribution, #eta: %0.1f - %0.1f (det jet), E: %d - %d GeV", EtaBinBorders[iE], EtaBinBorders[iE+1], int(JetEBorders[ipt]), int(JetEBorders[ipt+1])), 50, -1.0, 1.0);
            jetTree->Draw(Form("(jetE-jetE_match)/jetE_match>>hjetRatioE_Eta_%d_%d", iE, ipt), Form("jetR==%d && jetE_match>=%d && jetE_match<%d && jetEta>=%f && jetEta<%f && jetEta_match>%f && jetEta_match<%f && jet_distmatch<%f && jetParts_match > %d", int(Rvals[Rvalue] * 10), int(JetEBorders[ipt]), int(JetEBorders[ipt+1]),EtaBinBorders[iE], EtaBinBorders[iE+1],EtaBinBorders[iE], EtaBinBorders[iE+1], Rvals[Rvalue]*0.6, constMin), "goff"); 
            //if(hjetRatioE_Eta[iE][ipt]->GetEntries()==0) hjetRatioE_Eta[iE][ipt]->Fill(-2);
            hjetRatioE_Eta[iE][ipt]->Scale(normalizations[NormValue]/500.0);

            EtamediansE[iE][ipt]=GetMedian(hjetRatioE_Eta[iE][ipt]);
            EtameansE[iE][ipt]=hjetRatioE_Eta[iE][ipt]->GetMean();
            EtaSDE[iE][ipt]=hjetRatioE_Eta[iE][ipt]->GetStdDev();

            hEtaMedianE[iE]->Fill(JetEBorders[ipt+1]-((JetEBorders[ipt+1]-JetEBorders[ipt])/2.0),EtamediansE[iE][ipt]);
            hEtaMedianE[iE]->SetBinError(ipt+1,hjetRatioE_Eta[iE][ipt]->GetMeanError());

            hEtaMeanE[iE]->Fill(JetEBorders[ipt+1]-((JetEBorders[ipt+1]-JetEBorders[ipt])/2.0),EtameansE[iE][ipt]);
            hEtaMeanE[iE]->SetBinError(ipt+1,hjetRatioE_Eta[iE][ipt]->GetMeanError());
            hEtaSDE[iE]->Fill(JetEBorders[ipt+1]-((JetEBorders[ipt+1]-JetEBorders[ipt])/2.0),EtaSDE[iE][ipt]);
            hEtaSDE[iE]->SetBinError(ipt+1, hjetRatioE_Eta[iE][ipt]->GetStdDevError());

            //cout << "In jet pT range " << JetEBorders[ipt] << " - " << JetEBorders[ipt+1] << "GeV/c, number of jets accepted = " << hjetRatioE_Eta[iE][ipt]->GetEntries() << endl;

            hjetRatioE_Eta[iE][ipt]->GetXaxis()->SetTitle("(E^{det}-E^{part})/E^{part}");
            hjetRatioE_Eta[iE][ipt]->GetYaxis()->SetTitle("probability");
        }       
    }

    hMedianpT->GetXaxis()->SetTitle("p_{T}^{part} (GeV/c)");
    hMedianpT->GetYaxis()->SetTitle("Mean or median");

    hMeanpT->GetXaxis()->SetTitle("p_{T}^{part} (GeV/c)");
    hMeanpT->GetYaxis()->SetTitle("Mean or median");
 
    hSDpT->GetXaxis()->SetTitle("p_{T}^{part} (GeV/c)");
    hSDpT->GetYaxis()->SetTitle("Standard deviation");

    hMedianE->GetXaxis()->SetTitle("E^{part} (GeV)");
    hMedianE->GetYaxis()->SetTitle("Mean or median");

    hMeanE->GetXaxis()->SetTitle("E^{part} (GeV)");
    hMeanE->GetYaxis()->SetTitle("Mean or median");

    hSDE->GetXaxis()->SetTitle("E^{part} (GeV)");
    hSDE->GetYaxis()->SetTitle("Standard deviation");

    
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


