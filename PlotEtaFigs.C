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

Double_t textSize = 0.045;

void PlotEtaFigs(
                TString inputfileR02 = "JetJetOutput/FINALAN/oneETAMerged20230417_0-1000GeV_OutputR6.root", 
                TString inputfileR04 = "JetJetOutput/FINALAN/oneETAMerged20230417_0-1000GeV_OutputR4.root")
{
    StyleSettingsPaper();
    TGaxis::SetMaxDigits(4);

    Color_t colors[4]     = {kBlack,  kRed+1, kBlue+2, kGray+1};
    Color_t colorsMC[4]   = {kGray+1, kRed+3, kBlue+3, kBlack+1};
    Style_t marker[4]     = {20, 21, 33, 20};
    Style_t markerMC[4]   = {24, 25, 27, 24};
    Style_t style[4]      =  {1, 5, 7, 1};
    Size_t  markerS[4]    = { 2, 2, 2.4, 2};
    
    TString radiusOut[2]    = {"R06", "R04"};
    TString radiusLabel[2]  = {"#it{R} = 0.6", "#it{R} = 0.4"};
    TString etaOut[10]       = {"3.6", "3.8", "4.0", "4.2", "4.4", "4.6", "4.8", "5.0", "5.2", "5.4"};
    TString etaRange[9]     = {"3.6< #eta_{jet} < 3.8", "3.8 < #eta_{jet} < 4.0", "4.0 < #eta_{jet} < 4.2","4.2< #eta_{jet} < 4.4","4.4< #eta_{jet} < 4.6","4.6< #eta_{jet} < 4.8","4.8< #eta_{jet} < 5.0","5.0< #eta_{jet} < 5.2","5.2< #eta_{jet} < 5.4"}; //{3.8, 4.6, 5.4}
    Double_t rangeJES[2][2] = { {-0.7, -0.05}, {-0.7, 0.1} };

    const Int_t nEtaBinsN  = 10;
    const double EtaBinBorders[nEtaBinsN] = {3.6, 3.8, 4.0, 4.2, 4.4, 4.6, 4.8, 5.0, 5.2, 5.4};

    const Int_t nEtaBins  = 10;
    Int_t maxNEbins        = 5;
    const Int_t nNEBins = 6;

    Int_t exampleBins[3]    = {1,3,6};
    //Double_t binningE[nNEBins]  = {0.0, 200.0, 400.0, 600.0, 800.0, 1000.0, 1200.0, 1400.0, 1600.0, 1800.0, 2000.0}; 
    Double_t binningE[nNEBins]  = {0.0, 400.0, 800.0, 1200.0, 1600.0, 2000.0};  
    
    TH2D* histResponseMat_E[2][nNEBins];
    TH1D* histMean_E[2][nNEBins];
    TH1D* histMedian_E[2][nNEBins];
    TH1D* histSigma_E[2][nNEBins];
    TH1D* histDeltaE_bins[2][nEtaBinsN][nNEBins];

    TFile* fileR02 = new TFile(inputfileR02.Data());
    for (Int_t i = 0; i < nEtaBins-1; i++){
        histResponseMat_E[0][i]  = (TH2D*)fileR02->Get(Form("hRespMatrix_E%d",i));
        histMean_E[0][i]         = (TH1D*)fileR02->Get(Form("hEtaMeanE_%d",i));
        histMedian_E[0][i]       = (TH1D*)fileR02->Get(Form("hEtaMedianE_%d",i));
        histSigma_E[0][i]        = (TH1D*)fileR02->Get(Form("hEtaSDE_%d",i));
        for (Int_t p = 0; p < nNEBins-1; p++){
          histDeltaE_bins[0][i][p]   = (TH1D*)fileR02->Get(Form("hjetRatioE_Eta_%d_%d", i, p));
          histDeltaE_bins[0][i][p]->Scale(histDeltaE_bins[0][i][p]->GetBinWidth(1));
          histDeltaE_bins[0][i][p]->Scale(1./histDeltaE_bins[0][i][p]->Integral()); 
        }
      cout << "R = 0.2 \t" << histMean_E[0][i] << "\t"<<histSigma_E[0][i] << "\t"<<histResponseMat_E[0][i] << endl;
    }
    
    TFile* fileR04 = new TFile(inputfileR04.Data());
    for (Int_t i = 0; i < nEtaBins-1; i++){
       histResponseMat_E[1][i]  = (TH2D*)fileR04->Get(Form("hRespMatrix_E%d",i));
        histMean_E[1][i]         = (TH1D*)fileR04->Get(Form("hEtaMeanE_%d",i));
        histMedian_E[1][i]       = (TH1D*)fileR04->Get(Form("hEtaMedianE_%d",i));
        histSigma_E[1][i]        = (TH1D*)fileR04->Get(Form("hEtaSDE_%d",i));
        for (Int_t p = 0; p < nNEBins-1; p++){
          histDeltaE_bins[1][i][p]   = (TH1D*)fileR04->Get(Form("hjetRatioE_Eta_%d_%d", i, p));
          histDeltaE_bins[1][i][p]->Scale(histDeltaE_bins[1][i][p]->GetBinWidth(1));
          histDeltaE_bins[1][i][p]->Scale(1./histDeltaE_bins[1][i][p]->Integral()); 
        }  
      cout << "R = 0.4 \t" << histMean_E[1][i] << "\t"<<histSigma_E[1][i]  << "\t"<<histResponseMat_E[1][i] << endl;
    }

    TCanvas * cJES = new TCanvas("cJES", "cJES", 2*600,2*400);
    DrawPaperCanvasSettings(cJES, 0.081, 0.01, 0.013, 0.099 );  
    cJES->cd();
    
    TH2D* jesframe = new TH2D("jesframe", "", 5000, 100, 2002, 2000, -0.7, 0.2);
    SetStyleHistoTH2ForGraphs(jesframe, "E_{part} (GeV)","JES", 0.85*textSize,textSize, 0.85*textSize,textSize, 0.92, 0.95, 510, 505, 42, 62);


    TCanvas * cJER = new TCanvas("cJER", "cJER", 2*600,2*400);
    DrawPaperCanvasSettings(cJER, 0.081, 0.01, 0.013, 0.097 );  
    cJER->cd();
    
    TH2D* jerframe = new TH2D("jerframe", "", 5000, 100, 2002, 2000, 0.0001, 0.4);
    SetStyleHistoTH2ForGraphs(jerframe, "E_{part} (GeV)","JER", 0.85*textSize,textSize, 0.85*textSize,textSize, 0.92, 0.95, 510, 505, 42, 62);



    for (Int_t r = 0; r < 2; r++){
        cJES->Clear();
        jesframe->GetYaxis()->SetRangeUser(rangeJES[r][0], rangeJES[r][1]);
        jesframe->Draw("axis");

        TLegend* leg2 = GetAndSetLegend2(0.1, 0.93-3*1.0*textSize, 0.3, 0.93, textSize, 1, "", 42, 0.35);     
        for (Int_t e = 0; e < 3; e++){
            DrawSetMarker(histMean_E[r][2*e+1], marker[e], markerS[e], colors[e], colors[e]);
            histMean_E[r][2*e+1]->Draw("same");
            leg2->AddEntry(histMean_E[r][2*e+1], etaRange[2*e+1].Data(), "p");
        } 
        for (Int_t e = 3; e < 4; e++){
            DrawSetMarker(histMean_E[r][2*e], marker[e], markerS[e], colors[e], colors[e]);
            histMean_E[r][2*e]->Draw("same");
            leg2->AddEntry(histMean_E[r][2*e], etaRange[2*e].Data(), "p");
        } 
        leg2->Draw("same");
        drawLatexAdd("ALICE simulation, pp #sqrt{#it{s}} = 14 TeV",0.95,0.965-1*1.1*textSize, textSize,kFALSE, kFALSE, kTRUE);  
        drawLatexAdd("FoCal upgrade",0.95,0.965-2*1.1*textSize, textSize,kFALSE, kFALSE, kTRUE);  
        drawLatexAdd(Form("jets, anti-#it{k}_{T}, %s",radiusLabel[r].Data()),0.95,0.965-3*1.1*textSize, textSize,kFALSE, kFALSE, kTRUE);  
        jesframe->Draw("axis,same");
        
        cJES->cd();
        cJES->Update();
        cJES->SaveAs(Form("EtaJES_%s.png", radiusOut[r].Data()));

        cJER->Clear();
        jerframe->Draw("axis");
        TLegend* leg3 = GetAndSetLegend2(0.74, 0.94-6*1.1*textSize, 0.96, 0.94-3*1.1*textSize, textSize, 1, "", 42, 0.25);     
        for (Int_t e = 0; e < 3; e++){
            DrawSetMarker(histSigma_E[r][2*e+1], marker[e], markerS[e], colors[e], colors[e]);
            histSigma_E[r][2*e+1]->Draw("same");
            leg3->AddEntry(histSigma_E[r][2*e+1], etaRange[2*e+1].Data(), "p");
        } 
        for (Int_t e = 3; e < 4; e++){
            DrawSetMarker(histSigma_E[r][2*e], marker[e], markerS[e], colors[e], colors[e]);
            histSigma_E[r][2*e]->Draw("same");
            leg3->AddEntry(histSigma_E[r][2*e], etaRange[2*e].Data(), "p");
        } 
        leg3->Draw("same");
        drawLatexAdd("ALICE simulation, pp #sqrt{#it{s}} = 14 TeV",0.95,0.965-1*1.1*textSize, textSize,kFALSE, kFALSE, kTRUE);  
        drawLatexAdd("FoCal upgrade",0.95,0.965-2*1.1*textSize, textSize,kFALSE, kFALSE, kTRUE);  
        drawLatexAdd(Form("jets, anti-#it{k}_{T}, %s",radiusLabel[r].Data()),0.95,0.965-3*1.1*textSize, textSize,kFALSE, kFALSE, kTRUE);  
        jerframe->Draw("axis,same");
        
        cJER->cd();
        cJER->Update();
        cJER->SaveAs(Form("EtaJER_%s.png", radiusOut[r].Data()));
    }



    TCanvas *c1 = new TCanvas("c1", "c1: #Delta E for single energy bin, varying eta range, R=0.4", 1000, 1600);
    c1->Divide(3, 3);

    for (Int_t ibin = 1; ibin < 10; ibin++)
    {
        c1->cd(ibin);
        histDeltaE_bins[1][ibin-1][0]->SetTitle(Form("#Delta E , %4.0f < E < %4.0f GeV, %s", binningE[0], binningE[1], etaRange[ibin-1].Data()));
        //hDetHCALE[ibin-1]->SetLineColor(kBlack);
        histDeltaE_bins[1][ibin-1][0]->Draw();
    }

    c1->SaveAs("DeltaE_R04.png");

}