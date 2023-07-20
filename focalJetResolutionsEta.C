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

void focalJetResolutionsEta( 
                          TString inputfileR02 = "JetJetOutput/FINALAN/oneETAMerged20230417_0-1000GeV_OutputR6.root", 
                          TString inputfileR04 = "JetJetOutput/FINALAN/oneETAMerged20230417_0-1000GeV_OutputR4.root"
                        ){
    StyleSettingsPaper();
    TGaxis::SetMaxDigits(4);

    Color_t colors[3]     = {kBlack,  kRed+1, kBlue+2};
    Color_t colorsMC[3]   = {kGray+1, kRed+3, kBlue+3};
    Style_t marker[3]     = {20, 21, 33};
    Style_t markerMC[3]   = {24, 25, 27};
    Style_t style[3]      =  {1, 5, 7};
    Size_t  markerS[3]    = { 2, 2, 2.4};
    
    TString radiusOut[2]    = {"R06", "R04"};
    TString radiusLabel[2]  = {"#it{R} = 0.6", "#it{R} = 0.4"};
    TString etaOut[3]       = {"Full", "38to45", "45to51"};
    TString etaRange[9]     = {"3.6< #eta_{jet} < 3.8", "3.8 < #eta_{jet} < 4.0", "4.0 < #eta_{jet} < 4.2","4.2< #eta_{jet} < 4.4","4.4< #eta_{jet} < 4.6","4.6< #eta_{jet} < 4.8","4.8< #eta_{jet} < 5.0","5.0< #eta_{jet} < 5.2","5.2< #eta_{jet} < 5.4"}; //{3.8, 4.6, 5.4}
    Double_t rangeJES[2][2] = { {-0.7, -0.05}, {-0.7, 0.1} };

    Int_t nEtaBins  = 10;
    Int_t maxNEbins        = 5;
    Int_t exampleBins[3]    = {1,3,6};
    //Double_t binningE[11]  = {0.0, 200.0, 400.0, 600.0, 800.0, 1000.0, 1200.0, 1400.0, 1600.0, 1800.0, 2000.0}; 
    Double_t binningE[6]  = {0.0, 400.0, 800.0, 1200.0, 1600.0, 2000.0};  
    
    TH2D* histResponseMat_E[2][11];
    TH1D* histMean_E[2][10];
    TH1D* histMedian_E[2][10];
    TH1D* histSigma_E[2][10];
    TH1D* histDeltaE_bins[2][11][10];

    TFile* fileR02 = new TFile(inputfileR02.Data());
    for (Int_t i = 0; i < maxNEbins; i++){
        histResponseMat_E[0][i]  = (TH2D*)fileR02->Get(Form("hRespMatrix_E%d",i));
        histMean_E[0][i]         = (TH1D*)fileR02->Get(Form("hEtaMeanE_%d",i));
        histMedian_E[0][i]       = (TH1D*)fileR02->Get(Form("hEtaMedianE_%d",i));
        histSigma_E[0][i]        = (TH1D*)fileR02->Get(Form("hEtaSDE_%d",i));
        for (Int_t p = 0; p < nEtaBins-1; p++){
          histDeltaE_bins[0][i][p]   = (TH1D*)fileR02->Get(Form("hjetRatioE_Eta_%d_%d", p, i));
          histDeltaE_bins[0][i][p]->Scale(histDeltaE_bins[0][i][p]->GetBinWidth(1));
          histDeltaE_bins[0][i][p]->Scale(1./histDeltaE_bins[0][i][p]->Integral()); 
        }
      cout << "R = 0.2 \t" << histMean_E[0][i] << "\t"<<histSigma_E[0][i] << "\t"<<histResponseMat_E[0][i] << endl;
    }
    
    TFile* fileR04 = new TFile(inputfileR04.Data());
    for (Int_t i = 0; i < maxNEbins; i++){
       histResponseMat_E[1][i]  = (TH2D*)fileR04->Get(Form("hRespMatrix_E%d",i));
        histMean_E[1][i]         = (TH1D*)fileR04->Get(Form("hEtaMeanE_%d",i));
        histMedian_E[1][i]       = (TH1D*)fileR04->Get(Form("hEtaMedianE_%d",i));
        histSigma_E[1][i]        = (TH1D*)fileR04->Get(Form("hEtaSDE_%d",i));
        for (Int_t p = 0; p < nEtaBins-1; p++){
          histDeltaE_bins[1][i][p]   = (TH1D*)fileR04->Get(Form("hjetRatioE_Eta_%d_%d", p, i));
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
            DrawSetMarker(histMean_E[r][e], marker[e], markerS[e], colors[e], colors[e]);
            histMean_E[r][e]->Draw("same");
            leg2->AddEntry(histMean_E[r][e], etaRange[e].Data(), "p");
        } 
        leg2->Draw("same");
        drawLatexAdd("ALICE simulation, pp #sqrt{#it{s}} = 14 TeV",0.95,0.965-1*1.1*textSize, textSize,kFALSE, kFALSE, kTRUE);  
        drawLatexAdd("FoCal upgrade",0.95,0.965-2*1.1*textSize, textSize,kFALSE, kFALSE, kTRUE);  
        drawLatexAdd(Form("jets, anti-#it{k}_{T}, %s",radiusLabel[r].Data()),0.95,0.965-3*1.1*textSize, textSize,kFALSE, kFALSE, kTRUE);  
        jesframe->Draw("axis,same");
        
        cJES->cd();
        cJES->Update();
        cJES->SaveAs(Form("JES_%s.pdf", radiusOut[r].Data()));

        cJER->Clear();
        jerframe->Draw("axis");
        TLegend* leg3 = GetAndSetLegend2(0.74, 0.94-6*1.1*textSize, 0.96, 0.94-3*1.1*textSize, textSize, 1, "", 42, 0.25);     
        for (Int_t e = 0; e < 3; e++){
            DrawSetMarker(histSigma_E[r][e], marker[e], markerS[e], colors[e], colors[e]);
            histSigma_E[r][e]->Draw("same");
            leg3->AddEntry(histSigma_E[r][e], etaRange[e].Data(), "p");
        } 
        leg3->Draw("same");
        drawLatexAdd("ALICE simulation, pp #sqrt{#it{s}} = 14 TeV",0.95,0.965-1*1.1*textSize, textSize,kFALSE, kFALSE, kTRUE);  
        drawLatexAdd("FoCal upgrade",0.95,0.965-2*1.1*textSize, textSize,kFALSE, kFALSE, kTRUE);  
        drawLatexAdd(Form("jets, anti-#it{k}_{T}, %s",radiusLabel[r].Data()),0.95,0.965-3*1.1*textSize, textSize,kFALSE, kFALSE, kTRUE);  
        jerframe->Draw("axis,same");
        
        cJER->cd();
        cJER->Update();
        cJER->SaveAs(Form("JER_%s.pdf", radiusOut[r].Data()));
    }

    for (Int_t e = 0; e < 3; e++){
        cJES->Clear();
        jesframe->GetYaxis()->SetRangeUser(rangeJES[0][0], rangeJES[1][1]);
        jesframe->Draw("axis");

        TLegend* leg2 = GetAndSetLegend2(0.1, 0.93-2*1.0*textSize, 0.3, 0.93, textSize, 1, "", 42, 0.35);     
        for (Int_t r = 0; r < 2; r++){
            if (r == 1) DrawSetMarker(histMean_E[r][e], marker[e], markerS[e], colors[e], colors[e]);
            else        DrawSetMarker(histMean_E[r][e], markerMC[e], markerS[e]*1.2, colorsMC[e], colorsMC[e]);
            histMean_E[r][e]->Draw("same");
            leg2->AddEntry(histMean_E[r][e], radiusLabel[r].Data(), "p");
        } 
        leg2->Draw("same");
        drawLatexAdd("ALICE simulation, pp #sqrt{#it{s}} = 14 TeV",0.95,0.965-1*1.1*textSize, textSize,kFALSE, kFALSE, kTRUE);  
        drawLatexAdd("FoCal upgrade",0.95,0.965-2*1.1*textSize, textSize,kFALSE, kFALSE, kTRUE);  
        drawLatexAdd(Form("jets, anti-#it{k}_{T}, %s",etaRange[e].Data()),0.95,0.965-3*1.1*textSize, textSize,kFALSE, kFALSE, kTRUE);  
        jesframe->Draw("axis,same");
        
        cJES->cd();
        cJES->Update();
        cJES->SaveAs(Form("JES_%s.pdf", etaOut[e].Data()));

        cJER->Clear();
        jerframe->Draw("axis");
        TLegend* leg3 = GetAndSetLegend2(0.82, 0.94-5*1.1*textSize, 0.96, 0.94-3*1.1*textSize, textSize, 1, "", 42, 0.3);     
        for (Int_t r = 0; r < 2; r++){
            if (r == 1) DrawSetMarker(histSigma_E[r][e], marker[e], markerS[e], colors[e], colors[e]);
            else        DrawSetMarker(histSigma_E[r][e], markerMC[e], markerS[e]*1.2, colorsMC[e], colorsMC[e]);
            histSigma_E[r][e]->Draw("same");
            leg3->AddEntry(histSigma_E[r][e], radiusLabel[r].Data(), "p");
        } 
        leg3->Draw("same");
        drawLatexAdd("ALICE simulation, pp #sqrt{#it{s}} = 14 TeV",0.95,0.965-1*1.1*textSize, textSize,kFALSE, kFALSE, kTRUE);  
        drawLatexAdd("FoCal upgrade",0.95,0.965-2*1.1*textSize, textSize,kFALSE, kFALSE, kTRUE);  
        drawLatexAdd(Form("jets, anti-#it{k}_{T}, %s",etaRange[e].Data()),0.95,0.965-3*1.1*textSize, textSize,kFALSE, kFALSE, kTRUE);  
        jerframe->Draw("axis,same");
        
        cJER->cd();
        cJER->Update();
        cJER->SaveAs(Form("JER_%s.pdf", etaOut[e].Data()));
    }
    
        // Style projection
    TCanvas * cCorr = new TCanvas("cCorr", "cCorr", 2*900,2*700);
    DrawPaperCanvasSettings(cCorr, 0.095, 0.11, 0.02, 0.1 );  
    cCorr->SetLogz(1);

    for (Int_t e = 0; e < 3; e++){
      for (Int_t r = 0; r < 2; r++){
        cCorr->Clear();
         SetStyleHistoTH2ForGraphs(histResponseMat_E[r][e], "E_{det} (GeV)", "E_{part} (GeV)", 0.85*textSize,textSize, 0.85*textSize,textSize, 0.95, 0.9, 510, 510, 42, 62);
         histResponseMat_E[r][e]->GetZaxis()->SetLabelSize( 0.85*textSize);
         //histResponseMat_E[r][e]->Rebin2D(6,6);
         //histResponseMat_E[r][e]->GetZaxis()->SetRangeUser(0.001, 10);
         histResponseMat_E[r][e]->GetXaxis()->SetRangeUser(0, 1300);
         histResponseMat_E[r][e]->GetYaxis()->SetRangeUser(0, 1300);
         histResponseMat_E[r][e]->Draw("colz");

        drawLatexAdd("ALICE simulation, pp #sqrt{#it{s}} = 14 TeV",0.14,0.945-1*0.9*textSize, 0.85*textSize,kFALSE, kFALSE, kFALSE);  
        drawLatexAdd("FoCal upgrade",0.14,0.945-2*0.9*textSize, 0.85*textSize,kFALSE, kFALSE, kFALSE);  
        drawLatexAdd(Form("jets, anti-#it{k}_{T}, %s",radiusLabel[r].Data()),0.85,0.14+0.9*textSize, 0.85*textSize,kFALSE, kFALSE, kTRUE);  
        drawLatexAdd(Form("%s", etaRange[e].Data()),0.85,0.14, 0.85*textSize,kFALSE, kFALSE, kTRUE);  
        cCorr->Update();
        cCorr->SaveAs(Form("ResponseMatrix_%s_%s.pdf",etaOut[e].Data(), radiusOut[r].Data()));
      }
    }
    
    
    Double_t maxY   = 0.5;
    Int_t yDivs     = 505;
    Double_t pixelX = 1400*2;
    Double_t pixelY = 500*2;
    Double_t arrayBoundariesX1_4[4];
    Double_t arrayBoundariesY1_4[2];
    Double_t relativeMarginsX[3];
    Double_t relativeMarginsY[3];
    Int_t textSizeLabelsPixel = (Int_t)(textSize*pixelY);
    ReturnCorrectValuesForCanvasScaling(pixelX,pixelY, 3, 1, 0.048, 0.005, 0.01,0.1,arrayBoundariesX1_4,arrayBoundariesY1_4,relativeMarginsX,relativeMarginsY);
    Double_t textSizeLabels[3]  = {0.045, 0.045, 0.045};
    Double_t textSizeFac[3]     = {0.045, 0.045, 0.045};
    Double_t margin             = relativeMarginsX[0]*pixelX;
    
    TCanvas * cSlicesJES = new TCanvas("JESslice", "JESslice", pixelX,pixelY);
    DrawPaperCanvasSettings(cSlicesJES, 0, 0, 0, 0 );  

    cout <<  arrayBoundariesX1_4[0] << "\t"<<  arrayBoundariesY1_4[0] << "\t"<<  arrayBoundariesX1_4[1] << "\t"<<  arrayBoundariesY1_4[1] << endl;
    TPad* padRight = new TPad("padRight", "padRight", arrayBoundariesX1_4[0], arrayBoundariesY1_4[1], arrayBoundariesX1_4[1], arrayBoundariesY1_4[0]);
    DrawPaperPadSettings(padRight, relativeMarginsX[0], relativeMarginsX[1], relativeMarginsY[0], relativeMarginsY[2] );  
    ReturnCorrectValuesTextSize(   padRight, textSizeLabels[0], textSizeFac[0], textSizeLabelsPixel, 0.08);

    cout <<  arrayBoundariesX1_4[1] << "\t"<<  arrayBoundariesY1_4[0] << "\t"<<  arrayBoundariesX1_4[2] << "\t"<<  arrayBoundariesY1_4[1] << endl;
    TPad* padMiddle = new TPad("padMiddle", "padMiddle", arrayBoundariesX1_4[1], arrayBoundariesY1_4[1], arrayBoundariesX1_4[2], arrayBoundariesY1_4[0]);
    DrawPaperPadSettings(padMiddle, relativeMarginsX[1], relativeMarginsX[1], relativeMarginsY[0], relativeMarginsY[2] );  
    ReturnCorrectValuesTextSize(   padMiddle, textSizeLabels[1], textSizeFac[1], textSizeLabelsPixel, 0.08);

    cout <<  arrayBoundariesX1_4[2] << "\t"<<  arrayBoundariesY1_4[0] << "\t"<<  arrayBoundariesX1_4[3] << "\t"<<  arrayBoundariesY1_4[1] << endl;
    TPad* padLeft = new TPad("padLeft", "padLeft", arrayBoundariesX1_4[2], arrayBoundariesY1_4[1], arrayBoundariesX1_4[3], arrayBoundariesY1_4[0]);
    DrawPaperPadSettings(padLeft, relativeMarginsX[1], relativeMarginsX[2], relativeMarginsY[0], relativeMarginsY[2] );  
    ReturnCorrectValuesTextSize(   padLeft, textSizeLabels[2], textSizeFac[2], textSizeLabelsPixel, 0.08);

    
    for (Int_t r = 0; r < 2; r++){
      padRight->Draw();
      padMiddle->Draw();
      padLeft->Draw();

      
      padRight->cd();
      TH2D* esframeR = new TH2D("effframe", "", 5000,  -1.01, 0.5, 2000, 0., maxY);
      SetStyleHistoTH2ForGraphs(esframeR, "#Delta E = (E_{part}-E_{part})/E_{part} ","1/#it{k} d#it{N}/d#Delta E", 0.85*textSizeLabels[0],textSizeLabels[0], 0.85*textSizeLabels[0],textSizeLabels[0], 0.955, 1.45, 510, yDivs, 42, 62);
      esframeR->GetYaxis()->SetLabelOffset(0.01);
      esframeR->GetYaxis()->SetTickLength(0.035);
      esframeR->Draw("axis");

        for (Int_t e = 0; e < 3; e++){
            DrawSetMarker(histDeltaE_bins[r][e][exampleBins[0]], marker[e], markerS[e]*1.5, colors[e], colors[e]);
            histDeltaE_bins[r][e][exampleBins[0]]->Draw("same");
        } 
        
        drawLatexAdd("ALICE simulation",0.17,0.965-1*1.1*textSizeLabels[0], textSizeLabels[0],kFALSE, kFALSE, kFALSE);  
        drawLatexAdd("FoCal upgrade",0.17,0.965-2*1.1*textSizeLabels[0], textSizeLabels[0],kFALSE, kFALSE, kFALSE);  
        drawLatexAdd("pp #sqrt{#it{s}} = 13 TeV",0.17,0.965-3*1.1*textSizeLabels[0], textSizeLabels[0],kFALSE, kFALSE, kFALSE);  
        drawLatexAdd(Form("%.0f< #it{E}_{part}<%.0f GeV", binningE[exampleBins[0]], binningE[exampleBins[0]+1]),0.96,0.97-1*1.1*textSizeLabels[0], textSizeLabels[0],kFALSE, kFALSE, kTRUE);  

        DrawLines(0,0,0,0.05,2,kGray+2,7);
        
      esframeR->Draw("axis,same");
      padMiddle->cd();
        TH2D* esframeM = new TH2D("effframe", "", 5000, -1.01, 0.5, 2000, 0., maxY);
        SetStyleHistoTH2ForGraphs(esframeM, "#Delta E = (E_{part}-E_{part})/E_{part} ","1/#it{k} d#it{N}/d#Delta E", 0.85*textSizeLabels[1],textSizeLabels[1], 0.85*textSizeLabels[1],textSizeLabels[1], 0.84, 0.95,  510, yDivs, 42, 62);
        esframeM->GetYaxis()->SetTickLength(0.04);
        esframeM->GetXaxis()->SetLabelOffset(-0.001);
        esframeM->Draw("axis");

        
          TLegend* leg2 = GetAndSetLegend2(0.03, 0.945-3*1.*textSizeLabels[1], 0.4, 0.945, textSizeLabels[1], 1, "", 42, 0.35);     
          for (Int_t e = 0; e < 3; e++){
            DrawSetMarker(histDeltaE_bins[r][e][exampleBins[1]], marker[e], markerS[e]*1.5, colors[e], colors[e]);
            histDeltaE_bins[r][e][exampleBins[1]]->Draw("same");
            leg2->AddEntry(histDeltaE_bins[r][e][exampleBins[1]], etaRange[e].Data(), "p");
          } 
          leg2->Draw("same");
          drawLatexAdd(Form("%.0f< #it{E}_{part}<%.0f GeV", binningE[exampleBins[1]], binningE[exampleBins[1]+1]),0.95,0.97-1*1.1*textSizeLabels[1], textSizeLabels[1],kFALSE, kFALSE, kTRUE);  
          
          DrawLines(0,0,0,0.05,2,kGray+2,7);
          esframeM->Draw("axis,same");
          
      padLeft->cd();
        TH2D* esframeL = new TH2D("esframeL", "", 5000,  -1.01, 0.5, 2000, 0., maxY);
        SetStyleHistoTH2ForGraphs(esframeL, "#Delta E = (E_{part}-E_{part})/E_{part} ","#it{E}", 0.85*textSizeLabels[1],textSizeLabels[1], 0.85*textSizeLabels[1],textSizeLabels[1], 0.83, 0.95,  510, yDivs, 42, 62);
        esframeL->GetYaxis()->SetTickLength(0.038);
        esframeL->GetXaxis()->SetLabelOffset(-0.001);
        esframeL->Draw("axis");
      
        for (Int_t e = 0; e < 3; e++){
            DrawSetMarker(histDeltaE_bins[r][e][exampleBins[2]], marker[e], markerS[e]*1.5, colors[e], colors[e]);
            histDeltaE_bins[r][e][exampleBins[2]]->Draw("same");
        } 
        drawLatexAdd(Form("jets, anti-#it{k}_{T}, %s", radiusLabel[r].Data()),0.03,0.965-1*1.1*textSizeLabels[2], textSizeLabels[2],kFALSE, kFALSE, kFALSE);  
        drawLatexAdd(Form("%.0f< #it{E}_{part}<%.0f GeV", binningE[exampleBins[2]], binningE[exampleBins[2]+1]),0.93,0.965-1*1.1*textSizeLabels[2], textSizeLabels[2],kFALSE, kFALSE, kTRUE);  

        DrawLines(0,0,0,0.05,2,kGray+2,7);
        esframeL->Draw("axis,same");
        
      cSlicesJES->cd();
      cSlicesJES->Update();
      cSlicesJES->SaveAs(Form("JetEscaleProj_%s.pdf",radiusOut[r].Data()));
    }

    for (Int_t e = 0; e < 3; e++){
      padRight->Draw();
      padMiddle->Draw();
      padLeft->Draw();

      
      padRight->cd();
      TH2D* esframeR = new TH2D("effframe", "", 5000,  -1.01, 0.5, 2000, 0., maxY);
      SetStyleHistoTH2ForGraphs(esframeR, "#Delta E = (E_{part}-E_{part})/E_{part} ","1/#it{k} d#it{N}/d#Delta E", 0.85*textSizeLabels[0],textSizeLabels[0], 0.85*textSizeLabels[0],textSizeLabels[0], 0.955, 1.45, 510, yDivs, 42, 62);
      esframeR->GetYaxis()->SetLabelOffset(0.01);
      esframeR->GetYaxis()->SetTickLength(0.035);
      esframeR->Draw("axis");

        for (Int_t r = 0; r < 2; r++){
            if (r == 1)   DrawSetMarker(histDeltaE_bins[r][e][exampleBins[0]], marker[e], markerS[e]*1.5, colors[e], colors[e]);
            else          DrawSetMarker(histDeltaE_bins[r][e][exampleBins[0]], markerMC[e], markerS[e]*1.9, colorsMC[e], colorsMC[e]);
            histDeltaE_bins[r][e][exampleBins[0]]->Draw("same");
        } 
        
        drawLatexAdd("ALICE simulation",0.17,0.965-1*1.1*textSizeLabels[0], textSizeLabels[0],kFALSE, kFALSE, kFALSE);  
        drawLatexAdd("FoCal upgrade",0.17,0.965-2*1.1*textSizeLabels[0], textSizeLabels[0],kFALSE, kFALSE, kFALSE);  
        drawLatexAdd("pp #sqrt{#it{s}} = 13 TeV",0.17,0.965-3*1.1*textSizeLabels[0], textSizeLabels[0],kFALSE, kFALSE, kFALSE);  
        drawLatexAdd(Form("%.0f< E<%.0f GeV", binningE[exampleBins[0]], binningE[exampleBins[0]+1]),0.96,0.97-1*1.1*textSizeLabels[0], textSizeLabels[0],kFALSE, kFALSE, kTRUE);  

        DrawLines(0,0,0,0.05,2,kGray+2,7);
        
      esframeR->Draw("axis,same");
      padMiddle->cd();
        TH2D* esframeM = new TH2D("effframe", "", 5000, -1.01, 0.5, 2000, 0., maxY);
        SetStyleHistoTH2ForGraphs(esframeM, "#Delta E = (E_{part}-E_{part})/E_{part} ","1/#it{k} d#it{N}/d#Delta E", 0.85*textSizeLabels[1],textSizeLabels[1], 0.85*textSizeLabels[1],textSizeLabels[1], 0.84, 0.95,  510, yDivs, 42, 62);
        esframeM->GetYaxis()->SetTickLength(0.04);
        esframeM->GetXaxis()->SetLabelOffset(-0.001);
        esframeM->Draw("axis");

        
          TLegend* leg2 = GetAndSetLegend2(0.03, 0.945-2*1.*textSizeLabels[1], 0.4, 0.945, textSizeLabels[1], 1, "", 42, 0.35);     
          for (Int_t r = 0; r < 2; r++){
            if (r == 1)   DrawSetMarker(histDeltaE_bins[r][e][exampleBins[1]], marker[e], markerS[e]*1.5, colors[e], colors[e]);
            else          DrawSetMarker(histDeltaE_bins[r][e][exampleBins[1]], markerMC[e], markerS[e]*1.9, colorsMC[e], colorsMC[e]);
            histDeltaE_bins[r][e][exampleBins[1]]->Draw("same");
            leg2->AddEntry(histDeltaE_bins[r][e][exampleBins[1]], radiusLabel[r].Data(), "p");
          } 
          leg2->Draw("same");
          drawLatexAdd(Form("%.0f< E<%.0f GeV", binningE[exampleBins[1]], binningE[exampleBins[1]+1]),0.95,0.97-1*1.1*textSizeLabels[1], textSizeLabels[1],kFALSE, kFALSE, kTRUE);  
          
          DrawLines(0,0,0,0.05,2,kGray+2,7);
          esframeM->Draw("axis,same");
          
      padLeft->cd();
        TH2D* esframeL = new TH2D("esframeL", "", 5000,  -1.01, 0.5, 2000, 0., maxY);
        SetStyleHistoTH2ForGraphs(esframeL, "#Delta E = (E_{part}-E_{part})/E_{part} ","#it{E}", 0.85*textSizeLabels[1],textSizeLabels[1], 0.85*textSizeLabels[1],textSizeLabels[1], 0.83, 0.95,  510, yDivs, 42, 62);
        esframeL->GetYaxis()->SetTickLength(0.038);
        esframeL->GetXaxis()->SetLabelOffset(-0.001);
        esframeL->Draw("axis");
      
        for (Int_t r = 0; r < 2; r++){
            if (r == 1)   DrawSetMarker(histDeltaE_bins[r][e][exampleBins[2]], marker[e], markerS[e]*1.5, colors[e], colors[e]);
            else          DrawSetMarker(histDeltaE_bins[r][e][exampleBins[2]], markerMC[e], markerS[e]*1.9, colorsMC[e], colorsMC[e]);

            histDeltaE_bins[r][e][exampleBins[2]]->Draw("same");
        } 
        drawLatexAdd("jets, anti-#it{k}_{T}",0.03,0.965-1*1.1*textSizeLabels[2], textSizeLabels[2],kFALSE, kFALSE, kFALSE);  
        drawLatexAdd(Form("%s", etaRange[e].Data()),0.03,0.965-2*1.1*textSizeLabels[2], textSizeLabels[2],kFALSE, kFALSE, kFALSE);  
        drawLatexAdd(Form("%.0f< E<%.0f GeV", binningE[exampleBins[2]], binningE[exampleBins[2]+1]),0.93,0.965-1*1.1*textSizeLabels[2], textSizeLabels[2],kFALSE, kFALSE, kTRUE);  

        DrawLines(0,0,0,0.05,2,kGray+2,7);
        esframeL->Draw("axis,same");
        
      cSlicesJES->cd();
      cSlicesJES->Update();
      cSlicesJES->SaveAs(Form("JetEscaleProj_%s.pdf",etaOut[e].Data()));
    }
    
    
}
