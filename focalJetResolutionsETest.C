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

void focalJetResolutionsETest( 
                          TString inputfileR02 = "JetJetOutput/FINALAN/Merged20230417_0-1000GeV_OutputR6.root", 
                          TString inputfileR04 = "JetJetOutput/FINALAN/Merged20230417_0-1000GeV_OutputR4.root"
                        ){
    StyleSettingsPaper();
    TGaxis::SetMaxDigits(4);

    Color_t colors[3]     = {kBlack,  kRed+1, kBlue+2};
    Color_t colorsMC[3]   = {kGray+1, kRed+3, kBlue+3};
    Style_t marker[3]     = {20, 21, 33};
    Style_t markerMC[3]   = {24, 25, 27};
    Style_t style[3]      =  {1, 5, 7};
    Size_t  markerS[3]    = { 2, 2, 2};
    
    TString radiusOut[2]    = {"R06", "R04"};
    TString radiusLabel[2]  = {"#it{R} = 0.6", "#it{R} = 0.4"};
    TString etaOut[3]       = {"Full", "38to45", "45to51"};
    TString etaRange[3]     = {"3.4+R < #eta_{jet} < 5.5-R", "3.8 < #eta_{jet} < 4.5", "4.5 < #eta_{jet} < 5.1"}; //{3.8, 4.6, 5.4}
    Double_t rangeJES[2][2] = { {-0.7, -0.05}, {-0.7, 0.1} };
    Int_t maxNPtbins        = 15;
    Int_t exampleBins[3]    = {1,2,3};
    Double_t binningPt[16]  = {0.0, 50.0, 100.0, 200.0, 300.0, 400.0, 500.0, 600.0, 700.0, 800.0, 900.0, 1000.0, 1250.0, 1500.0, 1750.0, 2000.0}; 
     
    TH2D* histResponseMat_pT[2][3];
    TH1D* histMean_pT[2][3];
    TH1D* histMedian_pT[2][3];
    TH1D* histSigma_pT[2][3];
    TH1D* histDeltaPt_bins[2][3][9];
    TFile* fileR02 = new TFile(inputfileR02.Data());
    for (Int_t i = 0; i < 3; i++){
      if (i == 0){
        histResponseMat_pT[0][i]  = (TH2D*)fileR02->Get("hRespMatrix_E");
        histMean_pT[0][i]         = (TH1D*)fileR02->Get("hMeanE");
        histMedian_pT[0][i]       = (TH1D*)fileR02->Get("hMedianE");
        histSigma_pT[0][i]        = (TH1D*)fileR02->Get("hSDE");
        for (Int_t p = 0; p< maxNPtbins; p++){
          histDeltaPt_bins[0][i][p]   = (TH1D*)fileR02->Get(Form("hjetRatioE_%d", p));
          histDeltaPt_bins[0][i][p]->Scale(histDeltaPt_bins[0][i][p]->GetBinWidth(1)); 
          histDeltaPt_bins[0][i][p]->Scale(1./histDeltaPt_bins[0][i][p]->Integral()); 
        }
      } else {
        histResponseMat_pT[0][i]  = (TH2D*)fileR02->Get(Form("hRespMatrix_pT_Eta_%i",i-1));
        histMean_pT[0][i]         = (TH1D*)fileR02->Get(Form("hEtaMeanpT_%i",i-1));
        histMedian_pT[0][i]       = (TH1D*)fileR02->Get(Form("hEtaMedianpT_%i",i-1));
        histSigma_pT[0][i]        = (TH1D*)fileR02->Get(Form("hEtaSDpT_%i",i-1));
        for (Int_t p = 0; p< maxNPtbins; p++){
          histDeltaPt_bins[0][i][p]   = (TH1D*)fileR02->Get(Form("hjetRatioE_Eta_%i_%i",i-1, p));
          histDeltaPt_bins[0][i][p]->Scale(histDeltaPt_bins[0][i][p]->GetBinWidth(1));
          histDeltaPt_bins[0][i][p]->Scale(1./histDeltaPt_bins[0][i][p]->Integral()); 
        }
      }
      cout << "R = 0.2 \t" << histMean_pT[0][i] << "\t"<<histSigma_pT[0][i] << "\t"<<histResponseMat_pT[0][i] << endl;
    }
    
    TFile* fileR04 = new TFile(inputfileR04.Data());
    for (Int_t i = 0; i < 3; i++){
      if (i == 0){
        histResponseMat_pT[1][i]  = (TH2D*)fileR04->Get("hRespMatrix_pT_Eta_");
        histMean_pT[1][i]         = (TH1D*)fileR04->Get("hMedianE");
        histMedian_pT[1][i]       = (TH1D*)fileR04->Get("hMedianpT");
        histSigma_pT[1][i]        = (TH1D*)fileR04->Get("hSDpT");
        for (Int_t p = 0; p< maxNPtbins; p++){
          histDeltaPt_bins[1][i][p]   = (TH1D*)fileR04->Get(Form("hjetRatioE_%d", p));
          histDeltaPt_bins[1][i][p]->Scale(histDeltaPt_bins[1][i][p]->GetBinWidth(1));
          histDeltaPt_bins[1][i][p]->Scale(1./histDeltaPt_bins[1][i][p]->Integral()); 

        }
      } else {
        histResponseMat_pT[1][i]  = (TH2D*)fileR04->Get(Form("hRespMatrix_pT_Eta_%i",i-1));
        histMean_pT[1][i]         = (TH1D*)fileR04->Get(Form("hEtaMeanpT_%i",i-1));
        histMedian_pT[1][i]       = (TH1D*)fileR04->Get(Form("hEtaMedianpT_%i",i-1));
        histSigma_pT[1][i]        = (TH1D*)fileR04->Get(Form("hEtaSDpT_%i",i-1));
        for (Int_t p = 0; p< maxNPtbins; p++){
          histDeltaPt_bins[1][i][p]   = (TH1D*)fileR04->Get(Form("hjetRatioE_Eta_%i_%i",i-1, p));
          histDeltaPt_bins[1][i][p]->Scale(histDeltaPt_bins[1][i][p]->GetBinWidth(1));
          histDeltaPt_bins[1][i][p]->Scale(1./histDeltaPt_bins[1][i][p]->Integral()); 

        }        
      }
      cout << "R = 0.4 \t" << histMean_pT[1][i] << "\t"<<histSigma_pT[1][i]  << "\t"<<histResponseMat_pT[1][i] << endl;
    }
        
    Double_t maxY   = 0.25;//0.151;
    Int_t yDivs     = 505;
    Double_t pixelX = 500*2;
    Double_t pixelY = 500*2;
    Double_t textSizeLabels[3]  = {0.025, 0.025, 0.045};
    
    TCanvas * cSlicesJES = new TCanvas("JESslice", "JESslice", pixelX,pixelY);
       
    cSlicesJES->cd();


    Int_t r = 1; Int_t e=0;
    TLegend* leg2 = GetAndSetLegend2(0.63, 0.53, 0.84, 0.75, textSizeLabels[1], 1, "", 42, 0.35);     
  
    for (Int_t p = 0; p <3; p++){
      DrawSetMarker(histDeltaPt_bins[r][e][exampleBins[p]], marker[p], markerS[p]*1.2, colors[p], colors[p]);
      histDeltaPt_bins[r][e][exampleBins[p]]->Rebin(2);
      histDeltaPt_bins[r][e][exampleBins[p]]->GetYaxis()->SetRangeUser(0,maxY);
      histDeltaPt_bins[r][e][exampleBins[p]]->GetXaxis()->SetRangeUser(-1,0.8);
      //histDeltaPt_bins[r][e][exampleBins[p]]->Fit("gaus");
      histDeltaPt_bins[r][e][exampleBins[p]]->Draw("same");
      histDeltaPt_bins[r][e][exampleBins[p]]->SetTitle("");
      leg2->AddEntry(histDeltaPt_bins[r][e][exampleBins[p]], Form("%.0f< #it{E}_{part}<%.0f GeV", binningPt[exampleBins[p]], binningPt[exampleBins[p]+1]), "p");
    }
 
    leg2->Draw("same");  


    drawLatexAdd("ALICE simulation",0.7,0.945-1*1.1*textSizeLabels[0], textSizeLabels[0],kFALSE, kFALSE, kFALSE);  
    drawLatexAdd("FoCal upgrade",0.7,0.945-2*1.1*textSizeLabels[0], textSizeLabels[0],kFALSE, kFALSE, kFALSE);  
    drawLatexAdd("pp #sqrt{#it{s}} = 14 TeV",0.7,0.945-3*1.1*textSizeLabels[0], textSizeLabels[0],kFALSE, kFALSE, kFALSE);  
    drawLatexAdd("#it{R} = 0.4",0.7,0.945-4*1.1*textSizeLabels[0], textSizeLabels[0],kFALSE, kFALSE, kFALSE); 
    //drawLatexAdd(Form("%.0f< #it{p}_{T,part}<%.0f GeV/#it{c}", binningPt[exampleBins[0]], binningPt[exampleBins[0]+1]),0.96,0.97-1*1.1*textSizeLabels[0], textSizeLabels[0],kFALSE, kFALSE, kTRUE);  

    //DrawLines(0,0,0,0.05,2,kGray+2,7);
      
    //esframeR->Draw("axis,same");
        
    cSlicesJES->Update();
    cSlicesJES->SaveAs(Form("EJetEscaleProj_%s.pdf",radiusOut[r].Data()));
        
}