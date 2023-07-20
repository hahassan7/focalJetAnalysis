/* Macro for plotting some histograms from jet trees 
    THIS IS NOT USED SINCE YOU CAN TAKE THE MERGED R FILE DIRECTLY IN PlotDetResFractions.C
*/

double GetMedian( TH1D * h);

void JetPlottingMergedFileFractions(int Rvalue = 0)
{  
    const Int_t nR = 3; //5
    const Float_t Rvals[nR] = {0.2, 0.4, 0.6};//{0.2, 0.3, 0.4, 0.5, 0.6}; // Cone radii

    TFile *jetFile = TFile::Open(Form("/home/lmh/alice/jetAnalysis/NewPlotting/Fractions/Results/OneMergedR%d.root", int(Rvals[Rvalue] * 10)));


    const int nCol = 10;
    const int gcolors[nCol]={1,2,6,4,7,1,2,4,6,7};
    const int gmarkers[nCol]={4,8,25,21,8,21,25,4,8,21};

    const Float_t etaMin = 3.4; const Float_t etaMax = 5.5; 

    const Int_t nEBins  = 6;//16;
    //const double JetEBorders[nEBins] = {0.0, 50.0, 100.0, 200.0, 300.0, 400.0, 500.0, 600.0, 700.0, 800.0, 900.0, 1000.0, 1250.0, 1500.0, 1750.0, 2000.0}; 
    const double JetEBorders[nEBins] = {0.0, 400.0, 800.0, 1200.0, 1600.0, 2000.0}; 

    double sMeanDetHCALE[nEBins], sSDEDetHCALE[nEBins], sMeanDetECALE[nEBins], sSDEDetECALE[nEBins];

    double sMeanPartNeutralE[nEBins], sSDEPartNeutralE[nEBins], sMeanPartChargedE[nEBins], sSDEPartChargedE[nEBins];

    double sMeanPartPhotonElectronE[nEBins], sSDEPartPhotonElectronE[nEBins];

    TFile *fout = new TFile(Form("/home/lmh/alice/jetAnalysis/NewPlotting/Fractions/Results/OneMerged20230417_0-1000GeV_OutputR%d.root", int(Rvals[Rvalue] * 10)), "RECREATE");
    //histograms

    TH1F *hDetHCALE[nEBins]; //HCALfracdet, drawn like delta E in bins, and then take the mean and stdev of the different bins. 
    TH1F *hDetECALE[nEBins]; //ECALfracdet

    TH1F *hPartNeutralE[nEBins]; //neutralfrac
    TH1F *hPartChargedE[nEBins]; //chargedfrac
    
    TH1F *hPartPhotonElectronE[nEBins]; //photon&electronfrac
    
//////////////////////////////////////////////////////////////////////////////////////////
    TH1F *hMeanDetHCALE, *hMeanDetECALE, *hSDEDetHCALE, *hSDEDetECALE;                  //
    TH1F *hMeanPartNeutralE, *hMeanPartChargedE, *hSDEPartNeutralE, *hSDEPartChargedE;  //
    TH1F *hMeanPartPhotonElectronE, *hSDEPartPhotonElectronE;                           //

    hMeanDetHCALE = new TH1F("hMeanDetHCALE", "Mean of detector HCAL energy fraction distribution",nEBins-1, JetEBorders); 
    hMeanDetECALE = new TH1F("hMeanDetECALE", "Mean of detector HCAL energy fraction distribution",nEBins-1, JetEBorders); 
    hSDEDetHCALE = new TH1F("hSDEDetHCALE", "Standard deviation of detector HCAL energy fraction distribution",nEBins-1, JetEBorders);
    hSDEDetECALE = new TH1F("hSDEDetECALE", "Standard deviation of detector ECAL energy fraction distribution",nEBins-1, JetEBorders);

    hMeanPartNeutralE = new TH1F("hMeanPartNeutralE", "Mean of particle neutral energy fraction distribution",nEBins-1, JetEBorders); 
    hMeanPartChargedE = new TH1F("hMeanPartChargedE", "Mean of particle charged energy fraction distribution",nEBins-1, JetEBorders); 
    hSDEPartNeutralE = new TH1F("hSDEPartNeutralE", "Standard deviation of particle neutral energy fraction distribution",nEBins-1, JetEBorders);
    hSDEPartChargedE = new TH1F("hSDEPartChargedE", "Standard deviation of particle charged energy fraction distribution",nEBins-1, JetEBorders);

    hMeanPartPhotonElectronE = new TH1F("hMeanPartPhotonElectronE", "Mean of particle photon and electron energy fraction distribution",nEBins-1, JetEBorders); 
    hSDEPartPhotonElectronE = new TH1F("hSDEPartPhotonElectronE", "Standard deviation of particle photon and electron energy fraction distribution",nEBins-1, JetEBorders);

//////////////////////////////////////////////////////////////////////////////////////////

    TH2F *hRespMatrix_DetHCALE = (TH2F*)((TH2F*)jetFile->Get("hRespMatrix_DetHCALE"))->Clone();
    TH2F *hRespMatrix_DetECALE = (TH2F*)((TH2F*)jetFile->Get("hRespMatrix_DetECALE"))->Clone();

    TH2F *hRespMatrix_PartNeutralE = (TH2F*)((TH2F*)jetFile->Get("hRespMatrix_PartNeutralE"))->Clone();
    TH2F *hRespMatrix_PartChargedE = (TH2F*)((TH2F*)jetFile->Get("hRespMatrix_PartChargedE"))->Clone();

    TH2F *hRespMatrix_PartPhotonElectronE = (TH2F*)((TH2F*)jetFile->Get("hRespMatrix_PartPhotonElectronE"))->Clone();

    //E loop
    for (int iE = 0; iE < nEBins-1; ++iE) 
    {

        hDetHCALE[iE] = (TH1F*)((TH1F*)jetFile->Get(Form("hDetHCALE_%d", iE)))->Clone();
        hDetECALE[iE] = (TH1F*)((TH1F*)jetFile->Get(Form("hDetECALE_%d", iE)))->Clone();
        hPartNeutralE[iE] = (TH1F*)((TH1F*)jetFile->Get(Form("hPartNeutralE_%d", iE)))->Clone();
        hPartChargedE[iE] = (TH1F*)((TH1F*)jetFile->Get(Form("hPartChargedE_%d", iE)))->Clone();
        hPartPhotonElectronE[iE] = (TH1F*)((TH1F*)jetFile->Get(Form("hPartPhotonElectronE_%d", iE)))->Clone();
        
        //////////////
        sMeanDetHCALE[iE]=hDetHCALE[iE]->GetMean();
        sSDEDetHCALE[iE]=hDetHCALE[iE]->GetStdDev();
        sMeanDetECALE[iE]=hDetECALE[iE]->GetMean();
        sSDEDetECALE[iE]=hDetECALE[iE]->GetStdDev();

        sMeanPartNeutralE[iE]=hPartNeutralE[iE]->GetMean();
        sSDEPartNeutralE[iE]=hPartNeutralE[iE]->GetStdDev();
        sMeanPartChargedE[iE]=hPartChargedE[iE]->GetMean();
        sSDEPartChargedE[iE]=hPartChargedE[iE]->GetStdDev();

        sMeanPartPhotonElectronE[iE]=hPartPhotonElectronE[iE]->GetMean();
        sSDEPartPhotonElectronE[iE]=hPartPhotonElectronE[iE]->GetStdDev();
        //////////////

        hMeanDetHCALE->Fill(JetEBorders[iE+1]-((JetEBorders[iE+1]-JetEBorders[iE])/2),sMeanDetHCALE[iE]);
        hMeanDetHCALE->SetBinError(iE+1,hDetHCALE[iE]->GetMeanError());
        hSDEDetHCALE->Fill(JetEBorders[iE+1]-((JetEBorders[iE+1]-JetEBorders[iE])/2),sSDEDetHCALE[iE]);
        hSDEDetHCALE->SetBinError(iE+1, hDetHCALE[iE]->GetStdDevError());

        hMeanDetECALE->Fill(JetEBorders[iE+1]-((JetEBorders[iE+1]-JetEBorders[iE])/2),sMeanDetECALE[iE]);
        hMeanDetECALE->SetBinError(iE+1,hDetECALE[iE]->GetMeanError());
        hSDEDetECALE->Fill(JetEBorders[iE+1]-((JetEBorders[iE+1]-JetEBorders[iE])/2),sSDEDetECALE[iE]);
        hSDEDetECALE->SetBinError(iE+1, hDetECALE[iE]->GetStdDevError());

        hMeanPartNeutralE->Fill(JetEBorders[iE+1]-((JetEBorders[iE+1]-JetEBorders[iE])/2),sMeanPartNeutralE[iE]);
        hMeanPartNeutralE->SetBinError(iE+1,hPartNeutralE[iE]->GetMeanError());
        hSDEPartNeutralE->Fill(JetEBorders[iE+1]-((JetEBorders[iE+1]-JetEBorders[iE])/2),sSDEPartNeutralE[iE]);
        hSDEPartNeutralE->SetBinError(iE+1, hPartNeutralE[iE]->GetStdDevError());

        hMeanPartChargedE->Fill(JetEBorders[iE+1]-((JetEBorders[iE+1]-JetEBorders[iE])/2),sMeanPartChargedE[iE]);
        hMeanPartChargedE->SetBinError(iE+1,hPartChargedE[iE]->GetMeanError());
        hSDEPartChargedE->Fill(JetEBorders[iE+1]-((JetEBorders[iE+1]-JetEBorders[iE])/2),sSDEPartChargedE[iE]);
        hSDEPartChargedE->SetBinError(iE+1, hPartChargedE[iE]->GetStdDevError());

        hMeanPartPhotonElectronE->Fill(JetEBorders[iE+1]-((JetEBorders[iE+1]-JetEBorders[iE])/2),sMeanPartPhotonElectronE[iE]);
        hMeanPartPhotonElectronE->SetBinError(iE+1,hPartPhotonElectronE[iE]->GetMeanError());
        hSDEPartPhotonElectronE->Fill(JetEBorders[iE+1]-((JetEBorders[iE+1]-JetEBorders[iE])/2),sSDEPartPhotonElectronE[iE]);
        hSDEPartPhotonElectronE->SetBinError(iE+1, hPartPhotonElectronE[iE]->GetStdDevError());


        /////////////

        hDetHCALE[iE]->Scale(1/hDetHCALE[iE]->GetEntries());
        hDetECALE[iE]->Scale(1/hDetECALE[iE]->GetEntries());
        hPartNeutralE[iE]->Scale(1/hPartNeutralE[iE]->GetEntries());
        hPartChargedE[iE]->Scale(1/hPartChargedE[iE]->GetEntries());
        hPartPhotonElectronE[iE]->Scale(1/hPartPhotonElectronE[iE]->GetEntries());

        hDetHCALE[iE]->GetXaxis()->SetTitle("E_{HCAL}/(E_{jet}^{det}");
        hDetECALE[iE]->GetXaxis()->SetTitle("E_{ECAL}/(E_{jet}^{det}");
        hPartNeutralE[iE]->GetXaxis()->SetTitle("E_{neutral}/E_{jet}^{part}");
        hPartChargedE[iE]->GetXaxis()->SetTitle("E_{charged}/E_{jet}^{part}");
        hPartPhotonElectronE[iE]->GetXaxis()->SetTitle("E_{#gamma,e}/E_{jet}^{part}");

    }

    hMeanDetHCALE->GetYaxis()->SetTitle("Mean"); 
    hMeanDetHCALE->GetXaxis()->SetTitle("E^{det} (GeV)");
    hMeanDetECALE->GetYaxis()->SetTitle("Mean"); 
    hMeanDetECALE->GetXaxis()->SetTitle("E^{det} (GeV)");
    hSDEDetHCALE->GetYaxis()->SetTitle("Standard deviation"); 
    hSDEDetHCALE->GetXaxis()->SetTitle("E^{part} (GeV)");
    hSDEDetECALE->GetYaxis()->SetTitle("Standard deviation"); 
    hSDEDetECALE->GetXaxis()->SetTitle("E^{part} (GeV)");

    hMeanPartNeutralE->GetYaxis()->SetTitle("Mean"); 
    hMeanPartNeutralE->GetXaxis()->SetTitle("E^{part} (GeV)");
    hMeanPartChargedE->GetYaxis()->SetTitle("Mean");  
    hMeanPartChargedE->GetXaxis()->SetTitle("E^{part} (GeV)");
    hSDEPartNeutralE->GetYaxis()->SetTitle("Standard deviation"); 
    hSDEPartNeutralE->GetXaxis()->SetTitle("E^{part} (GeV)");
    hSDEPartChargedE->GetYaxis()->SetTitle("Standard deviation"); 
    hSDEPartChargedE->GetXaxis()->SetTitle("E^{part} (GeV)");

    hMeanPartPhotonElectronE->GetYaxis()->SetTitle("Mean"); 
    hMeanPartPhotonElectronE->GetXaxis()->SetTitle("E^{part} (GeV)");
    hSDEPartPhotonElectronE->GetYaxis()->SetTitle("Standard deviation"); 
    hSDEPartPhotonElectronE->GetXaxis()->SetTitle("E^{part} (GeV)");

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


