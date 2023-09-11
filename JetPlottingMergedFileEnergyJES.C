/* Macro for plotting some histograms from jet trees
    2.21064
    0.0669805
    0.00182628
    0.000139462
    0.0000225822
*/
#include "iostream"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"

double GetMedian(TH1D *h);

void JetPlottingMergedFileEnergyJES(int Rvalue = 2)
{
    const Int_t nR = 3;                        // 5
    const Float_t Rvals[nR] = {0.2, 0.4, 0.6}; // Cone radii

    TFile *jetFile = TFile::Open(Form("Data20230816/JES/EnMergedR%d_0Mass_20PercPart.root", int(Rvals[Rvalue] * 10)));

    TTree *jetTree = (TTree *)jetFile->Get("jetTree");
    TTree *TruthjetTree = (TTree *)jetFile->Get("truthjetTree");
    // TTree *PDGTree = (TTree *)jetFile->Get("inPDGTree");

    const int nCol = 10;
    const int gcolors[nCol] = {1, 2, 6, 4, 7, 1, 2, 4, 6, 7};
    const int gmarkers[nCol] = {4, 8, 25, 21, 8, 21, 25, 4, 8, 21};

    const Float_t etaMin = 3.4;
    const Float_t etaMax = 5.5;

    const Int_t nEBins = 15;
    const double JetEBorders[nEBins] = {100.0, 200.0, 300.0, 400.0, 500.0, 600.0, 700.0, 800.0, 900.0, 1000.0, 1250.0, 1500.0, 1750.0, 2000.0, 3000.0}; // const Int_t nEBins  = 15;

    // eta binned
    const Int_t nEtaBins = 3; // 11;
    const Float_t EtaBinBorders[nEtaBins] = {3.8, 4.5, 5.1};

    double mediansE[nEBins - 1];
    double meansE[nEBins - 1];
    double SDE[nEBins - 1];

    TH1D *hjetRatioE[nEBins - 1];
    TH1D *hMedianE;
    TH1D *hMeanE;
    TH1D *hSDE;

    const int NEFbins = 2;

    double mediansENEF[nEBins - 1][NEFbins];
    double meansENEF[nEBins - 1][NEFbins];
    double SDENEF[nEBins - 1][NEFbins];

    TH1D *hjetRatioENEF[nEBins - 1][NEFbins];
    TH1D *hMedianENEF[NEFbins];
    TH1D *hMeanENEF[NEFbins];
    TH1D *hSDENEF[NEFbins];

    TFile *fout = new TFile(Form("Data20230816/JES/EneMerged_OutputR%d_0Mass_20PercPart.root", int(Rvals[Rvalue] * 10)), "RECREATE");

    hMedianE = new TH1D("hMedianE", "Mean and median of #DeltaE distribution;E^{part} (GeV);Mean or median", nEBins - 1, JetEBorders);
    hMeanE = new TH1D("hMeanE", "Mean and median of #DeltaE distribution;E^{part} (GeV);Mean or median", nEBins - 1, JetEBorders);
    hSDE = new TH1D("hSDE", "Standard deviation of #DeltaE distribution;E^{part} (GeV);Standard deviation", nEBins - 1, JetEBorders);
    for (int nef = 0; nef < NEFbins; nef++)
    {
        hMedianENEF[nef] = new TH1D(Form("hMedianENEF%d", nef), Form("Mean and median of #DeltaE distribution. NEF bin%d;E^{part} (GeV);Mean or median", nef), nEBins - 1, JetEBorders);
        hMeanENEF[nef] = new TH1D(Form("hMeanENEF%d", nef), Form("Mean and median of #DeltaE distribution. NEF bin%d;E^{part} (GeV);Mean or median", nef), nEBins - 1, JetEBorders);
        hSDENEF[nef] = new TH1D(Form("hSDENEF%d", nef), Form("Standard deviation of #DeltaE distribution. NEF bin%d;E^{part} (GeV);Standard deviation", nef), nEBins - 1, JetEBorders);
    }

    // E loop
    for (int ipt = 0; ipt < nEBins - 1; ++ipt)
    {
        hjetRatioE[ipt] = (TH1D *)((TH1D *)jetFile->Get(Form("hjetRatioE_%d", ipt)))->Clone();

        hjetRatioE[ipt]->Scale(1. / hjetRatioE[ipt]->Integral(), "");

        mediansE[ipt] = GetMedian(hjetRatioE[ipt]);
        meansE[ipt] = hjetRatioE[ipt]->GetMean();
        SDE[ipt] = hjetRatioE[ipt]->GetStdDev();

        hMedianE->Fill(JetEBorders[ipt + 1] - ((JetEBorders[ipt + 1] - JetEBorders[ipt]) / 2), mediansE[ipt]);
        hMedianE->SetBinError(ipt + 1, hjetRatioE[ipt]->GetMeanError());
        hMeanE->Fill(JetEBorders[ipt + 1] - ((JetEBorders[ipt + 1] - JetEBorders[ipt]) / 2), meansE[ipt]);
        hMeanE->SetBinError(ipt + 1, hjetRatioE[ipt]->GetMeanError());
        hSDE->Fill(JetEBorders[ipt + 1] - ((JetEBorders[ipt + 1] - JetEBorders[ipt]) / 2), SDE[ipt]);
        hSDE->SetBinError(ipt + 1, hjetRatioE[ipt]->GetStdDevError());

        hjetRatioENEF[ipt][0] = (TH1D *)((TH1D *)jetFile->Get(Form("hjetRatioENEF_%d_0", ipt)))->Clone();
        hjetRatioENEF[ipt][1] = (TH1D *)((TH1D *)jetFile->Get(Form("hjetRatioENEF_%d_1", ipt)))->Clone();

        for (int nef = 0; nef < NEFbins; nef++)
        {
            hjetRatioENEF[ipt][nef]->Scale(1. / hjetRatioENEF[ipt][nef]->Integral(), "");

            mediansENEF[ipt][nef] = GetMedian(hjetRatioENEF[ipt][nef]);
            meansENEF[ipt][nef] = hjetRatioENEF[ipt][nef]->GetMean();
            SDENEF[ipt][nef] = hjetRatioENEF[ipt][nef]->GetStdDev();

            hMedianENEF[nef]->Fill(JetEBorders[ipt + 1] - ((JetEBorders[ipt + 1] - JetEBorders[ipt]) / 2), mediansENEF[ipt][nef]);
            hMedianENEF[nef]->SetBinError(ipt + 1, hjetRatioENEF[ipt][nef]->GetMeanError());
            hMeanENEF[nef]->Fill(JetEBorders[ipt + 1] - ((JetEBorders[ipt + 1] - JetEBorders[ipt]) / 2), meansENEF[ipt][nef]);
            hMeanENEF[nef]->SetBinError(ipt + 1, hjetRatioENEF[ipt][nef]->GetMeanError());
            hSDENEF[nef]->Fill(JetEBorders[ipt + 1] - ((JetEBorders[ipt + 1] - JetEBorders[ipt]) / 2), SDENEF[ipt][nef]);
            hSDENEF[nef]->SetBinError(ipt + 1, hjetRatioENEF[ipt][nef]->GetStdDevError());

            std::cout << "In jet E range " << JetEBorders[ipt] << " - " << JetEBorders[ipt + 1] << "GeV, number of jets accepted = " << hjetRatioENEF[ipt][nef]->GetEntries() << std::endl;

        }
    }

    fout->Write();
    fout->Close();
}

double GetMedian(TH1D *h)
{

    Double_t x, q;
    q = 0.5;              // 0.5 for "median"
    h->ComputeIntegral(); // just a precaution
    h->GetQuantiles(1, &x, &q);
    std::cout << "median = " << x << std::endl;
    return x;
}
