#include <TROOT.h>
#include <TStyle.h>
#include <TColor.h>
#include <TH1.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TLine.h>
#include <TFile.h>
#include <TF1.h>
#include <TGraph.h>
#include <TLegend.h>
#include <TSpline.h>
#include <TLatex.h>
#include <TArrow.h>
#include <TGaxis.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include "paperPlotsHeader.h"

Double_t textSize = 0.05; // 5;
const int SettingLogZ = 1;

void NormalizeToProb(TH2F *hHisto);

void ChangeTitleOfCanvas(TCanvas *cCanvas, TString sTitle);

void JetPlottingMergedFractions2DPlot()
{
  TGaxis::SetMaxDigits(3);

    // setting NormValue = 8 if merged file in use

    const double limitYconst = 0.2;
    const double limitYpT = 0.2;
    const double limitYEta = 0.2;
    const double limitYconstmin = 0.0;
    const double limitYpTmin = 0.0;
    const double limitYEtamin = 0.0;

    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    // Open the input file
    TFile *fin;
    fin = new TFile("Data20230816/2D/2D_Merged_0Mass.root", "READ");

    const int nR = 3;
    const int nEBins = 15;                                                                                                                            // 16;
    const int npTBins = 9;                                                                                                                            // 16;
    const double JetEBorders[nEBins] = {100.0, 200.0, 300.0, 400.0, 500.0, 600.0, 700.0, 800.0, 900.0, 1000.0, 1250.0, 1500.0, 1750.0, 2000.0, 3000.0}; //{100.0, 400.0, 800.0, 1200.0, 1600.0, 2000.0};
    const double JetpTBorders[npTBins] = {5.0, 10.0, 20.0, 40.0, 60.0, 80.0, 100.0, 150.0, 200.0};

    TH2F *hRespMatrix[nR];
    TH2F *hRespMatrixE[nR][nEBins - 1];
    TH2F *hRespMatrixpT[nR][npTBins - 1];

    int exampleBins[3] = {0, 2, 3};

    for (int Rvalue = 2; Rvalue < nR; Rvalue++)
    {
        // histograms
        hRespMatrix[Rvalue] = (TH2F *)fin->Get(Form("hRespMatrix_R%d", Rvalue));
        NormalizeToProb(hRespMatrix[Rvalue]);

        for (int iE = 0; iE < nEBins - 1; ++iE)
        {
            hRespMatrixE[Rvalue][iE] = (TH2F *)fin->Get(Form("hRespMatrixE%d_R%d", iE, Rvalue));
            NormalizeToProb(hRespMatrixE[Rvalue][iE]);
        }

        for (int ipT = 0; ipT < npTBins - 1; ++ipT)
        {
            hRespMatrixpT[Rvalue][ipT] = (TH2F *)fin->Get(Form("hRespMatrixpT%d_R%d", ipT, Rvalue));
            NormalizeToProb(hRespMatrixpT[Rvalue][ipT]);
        }
    }

    // Style projection
    TCanvas *cCorr = new TCanvas("cCorr", "cCorr", 2 * 900, 2 * 700);
    DrawPaperCanvasSettings(cCorr, 0.105, 0.12, 0.09, 0.11);
    // cCorr->SetLogz(1);
    gPad->SetTopMargin(0.11);
    // gPad->SetLeftMargin(0.1);
    // gPad->SetRightMargin(0.05);
    // gPad->SetBottomMargin(0.1);

    for (int e = 0; e < 3; e++)
    {
        for (int r = 2; r < nR; r++)
        {
            cCorr->Clear();
            SetStyleHistoTH2ForGraphs(hRespMatrixpT[r][exampleBins[e]], "#it{E}_{ECAL}/#it{E}_{jet}^{det} (GeV)", "#it{E}_{neutral}/#it{E}_{jet}^{part} (GeV)", 0.85 * textSize, textSize, 0.85 * textSize, textSize, 0.95, 0.95, 510, 510, 42, 62);
            hRespMatrixpT[r][exampleBins[e]]->GetZaxis()->SetLabelSize(0.85 * textSize);
            // hRespMatrixpT[r][exampleBins[e]]->GetZaxis()->SetRangeUser(0.00001, 0.5);
            // hRespMatrixpT[r][exampleBins[e]]->GetXaxis()->SetRangeUser(JESMinX, JESMaxX);
            // hRespMatrixpT[r][exampleBins[e]]->GetYaxis()->SetRangeUser(JESMinX, JESMaxX);
            hRespMatrixpT[r][exampleBins[e]]->Draw("colz");

            drawLatexAdd("ALICE simulation, pp #sqrt{#it{s}} = 14 TeV", 0.1, 1.01 - 1 * 0.9 * textSize, 0.85 * textSize, kFALSE, kFALSE, kFALSE);
            drawLatexAdd("FoCal upgrade", 0.1, 1.01 - 2 * 0.9 * textSize, 0.85 * textSize, kFALSE, kFALSE, kFALSE);
            drawLatexAdd("jets, anti-#it{k}_{T}, #it{R} = 0.6", 0.9, 1.01 - 1 * 0.9 * textSize, 0.85 * textSize, kFALSE, kFALSE, kTRUE);
            drawLatexAdd(Form("%.0f< #it{p}_{T,det} < %.0f GeV/#it{c}", JetpTBorders[exampleBins[e]], JetpTBorders[exampleBins[e] + 1]), 0.9, 1.01 - 2 * 0.9 * textSize, 0.85 * textSize, kFALSE, kFALSE, kTRUE);
            cCorr->Update();
            cCorr->SaveAs(Form("figs/2D/ResponseMatrix_%d_6.pdf", e));
        }
    }
}

void NormalizeToProb(TH2F *hHisto)
{
    Double_t factor = 1.;
    hHisto->Rebin2D();
    if (hHisto != NULL)
        hHisto->Scale(factor / hHisto->Integral()); //, "width");
    // Add the following line to scale to picobarns
    // if (hHisto!=NULL) hHisto->Scale(factor/scalingFactor, "width");
}

void ChangeTitleOfCanvas(TCanvas *cCanvas, TString sTitle)
{
    cCanvas->cd();
    TLatex *tex = new TLatex(0.31, 0.92, Form("%s", sTitle.Data()));
    tex->SetNDC();
    tex->SetTextSize(0.038);
    tex->Draw();
    if (SettingLogZ == 1)
        cCanvas->SetLogz();
}