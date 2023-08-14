#include <TH1F.h>
#include <TCanvas.h>
#include <TFile.h>
#include "paperPlotsHeader.h"

void NormalizeToProb(TH1F* hHisto);

void ChangeTitleOfCanvas(TCanvas* cCanvas, TString sTitle);

const double scalingFactor = 1e-9;

const int SettingLogY = 1;
Double_t textSize = 0.03;//5;


//Draws particle level jet pT E and const spectra
void SpectraBoth() {

    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    //gStyle->SetHistLineWidth(20);



    const double limitYconst = 1e10;//0.2;
    const double limitYpT    = 2e10;//0.5;
    const double limitYE     = 1e11;//0.011;
    const double limitYconstmin = 0.000001;
    const double limitYpTmin    = 0.0100;
    const double limitYEmin     = 10000;
    const double XpTmin     = 5; //change to 10 unless lowest bin, then 5
    const double XpTmax     = 200;
    const double XpTminFit     = 10; //change to 10 unless lowest bin, then 5
    const double XpTmaxFit     = 100;
    const double XEmin     = 10; 
    const double XEmax     = 2000;
    //->GetYaxis()->SetRangeUser

    // Arrays to hold histograms for different types and R values
    TH1F* hist_pT_R[3]; // hist_pT_R[0], hist_pT_R[1], hist_pT_R[2]
    TH1F* hist_E_R[3]; // hist_E_R[0], hist_E_R[1], hist_E_R[2]

    TH1F* hist_pT_R_e[3][2]; // hist_pT_R_e[0][0], hist_pT_R_e[0][1], ..., hist_pT_R_e[2][2]
    TH1F* hist_E_R_e[3][2]; // hist_E_R_e[0][0], hist_E_R_e[0][1], ..., hist_E_R_e[2][2]
    // Arrays to hold histograms for different types and R values
    TH1F* Dethist_pT_R[3]; // hist_pT_R[0], hist_pT_R[1], hist_pT_R[2]
    TH1F* Dethist_E_R[3]; // hist_E_R[0], hist_E_R[1], hist_E_R[2]

    TH1F* Dethist_pT_R_e[3][2]; // hist_pT_R_e[0][0], hist_pT_R_e[0][1], ..., hist_pT_R_e[2][2]
    TH1F* Dethist_E_R_e[3][2]; // hist_E_R_e[0][0], hist_E_R_e[0][1], ..., hist_E_R_e[2][2]


    Color_t colors[5]     = {kBlack,  kRed+1, kViolet, kCyan-1, kGray+1};
    Color_t colorsMC[3]   = {kGray+1, kRed+3, kBlue+3};
    Style_t marker[5]     = {20, 21, 33, 20, 21};
    Style_t markerMC[5]   = {24, 25, 27, 24, 25};
    Style_t style[5]      =  {1, 5, 7, 1, 5};
    Size_t  markerS[5]    = { 2, 2, 2.4, 2, 2};

    // Read histograms from the ROOT file
    TFile* file = TFile::Open("Data20230728/Spectra/ScaledMerged.root", "READ");
    TFile* fileDet = TFile::Open("Data20230728/Spectra/ScaledDetMerged.root", "READ");


    const Int_t nR = 3; 
    const Float_t Rvals[nR] = {0.2, 0.4, 0.6}; // Cone radii
    const Int_t nEtaBins  = 3;
    //const double EtaBinBorders[nEtaBins] = {3.4,5.0,5.5}; //
    //const double EtaMin[2] = {3.4,3.9}; 
    //const double EtaMax[2] = {5.0,5.5};
    TString etaRange[nEtaBins-1]     =  {"4.0< #eta_{jet} < 4.5","4.5 < #eta_{jet} < 4.9"}; //consider making this different range for R=0.6, 0.4, 0.2. ...
    const Int_t nConst = 5;
    const Int_t ConstVals[nConst] = {1, 2, 3, 5, 10}; // constituent values to check


    for (int R = 0; R < 3; ++R) {
        // Read hist_pT_R[3] histograms from the file using Clone method
        hist_pT_R[R] = (TH1F*)(file->Get(Form("hist_pT_R%d", R))->Clone());
        hist_pT_R[R]->SetLineColor(colors[R]);
        NormalizeToProb(hist_pT_R[R]);
        hist_pT_R[R]->SetTitle(Form("pT Distribution for R = %d; pT [GeV/c]; d#sigma/dp_{T}dy (pb(GeV/c)^{-1})", R));
        hist_pT_R[R]->GetYaxis()->SetRangeUser(limitYpTmin,limitYpT);
        hist_pT_R[R]->GetXaxis()->SetRangeUser(XpTmin,XpTmax);
        //normalize all the histograms to the NofEntries (check if 1/Integral is enough etc. or nentries or what not)

        Dethist_pT_R[R] = (TH1F*)(fileDet->Get(Form("hist_pT_R%d", R))->Clone());
        Dethist_pT_R[R]->SetLineColor(colors[R]);
        NormalizeToProb(Dethist_pT_R[R]);
        Dethist_pT_R[R]->SetTitle(Form("pT Distribution for R = %d; pT [GeV/c]; d#sigma/dp_{T}dy (pb(GeV/c)^{-1})", R));
        Dethist_pT_R[R]->GetYaxis()->SetRangeUser(limitYpTmin,limitYpT);
        Dethist_pT_R[R]->GetXaxis()->SetRangeUser(XpTmin,XpTmax);

        // Read hist_E_R[3] histograms from the file using Clone method
        hist_E_R[R] = (TH1F*)(file->Get(Form("hist_E_R%d", R))->Clone());
        hist_E_R[R]->SetLineColor(colors[R]);
        NormalizeToProb(hist_E_R[R]);
        hist_E_R[R]->SetTitle(Form("E Distribution for R = %d; E [GeV]; d#sigma/dEdy (pb(GeV)^{-1})", R));
        hist_E_R[R]->GetYaxis()->SetRangeUser(limitYEmin,limitYE);
        hist_E_R[R]->GetXaxis()->SetRangeUser(XEmin,XEmax);


        Dethist_E_R[R] = (TH1F*)(fileDet->Get(Form("hist_E_R%d", R))->Clone());
        Dethist_E_R[R]->SetLineColor(colors[R]);
        NormalizeToProb(Dethist_E_R[R]);
        Dethist_E_R[R]->SetTitle(Form("E Distribution for R = %d; E [GeV]; d#sigma/dEdy (pb(GeV)^{-1})", R));
        Dethist_E_R[R]->GetYaxis()->SetRangeUser(limitYEmin,limitYE);
        Dethist_E_R[R]->GetXaxis()->SetRangeUser(XEmin,XEmax);

        for (int e = 0; e < 2; ++e) {
            // Read hist_pT_R_e[3][3] histograms from the file using Clone method
            hist_pT_R_e[R][e] = (TH1F*)(file->Get(Form("hist_pT_R%d_e%d", R, e))->Clone());
            hist_pT_R_e[R][e]->SetLineColor(colors[e]);
            hist_pT_R_e[R][e]->SetLineStyle(style[0]);
            NormalizeToProb(hist_pT_R_e[R][e]);
            hist_pT_R_e[R][e]->SetTitle(Form("pT Distribution for R = %d, e = %d; p_{T} (GeV/c); d#sigma/dp_{T}dy (pb(GeV/c)^{-1})", R, e));
            hist_pT_R_e[R][e]->GetYaxis()->SetRangeUser(limitYpTmin,limitYpT);
            hist_pT_R_e[R][e]->GetXaxis()->SetRangeUser(XpTmin,XpTmax);


            Dethist_pT_R_e[R][e] = (TH1F*)(fileDet->Get(Form("Dethist_pT_R%d_e%d", R, e))->Clone());
            Dethist_pT_R_e[R][e]->SetLineColor(colors[e]);
            Dethist_pT_R_e[R][e]->SetLineStyle(style[1]);
            NormalizeToProb(Dethist_pT_R_e[R][e]);
            //Dethist_pT_R_e[R][e]->SetTitle(Form("pT Distribution for R = %d, e = %d; pT [GeV/c]; d#sigma/dp_{T}dy", R, e));
            Dethist_pT_R_e[R][e]->GetYaxis()->SetRangeUser(limitYpTmin,limitYpT);
            Dethist_pT_R_e[R][e]->GetXaxis()->SetRangeUser(XpTmin,XpTmax);

            // Read hist_E_R_e[3][3] histograms from the file using Clone method
            hist_E_R_e[R][e] = (TH1F*)(file->Get(Form("hist_E_R%d_e%d", R, e))->Clone());
            hist_E_R_e[R][e]->SetLineColor(colors[e]);
            hist_E_R_e[R][e]->SetLineStyle(style[0]);
            NormalizeToProb(hist_E_R_e[R][e]);
            hist_E_R_e[R][e]->SetTitle(Form("E Distribution for R = %d, e = %d; E [GeV]; d#sigma/dEdy (pb(GeV)^{-1})", R, e));
            hist_E_R_e[R][e]->GetYaxis()->SetRangeUser(limitYEmin,limitYE);
            hist_E_R_e[R][e]->GetXaxis()->SetRangeUser(XEmin,XEmax);


            Dethist_E_R_e[R][e] = (TH1F*)(fileDet->Get(Form("hist_E_R%d_e%d", R, e))->Clone());
            Dethist_E_R_e[R][e]->SetLineColor(colors[e]);
            Dethist_E_R_e[R][e]->SetLineStyle(style[1]);
            NormalizeToProb(Dethist_E_R_e[R][e]);
            //Dethist_E_R_e[R][e]->SetTitle(Form("E Distribution for R = %d, e = %d; E [GeV]; d#sigma/dEdy", R, e));
            Dethist_E_R_e[R][e]->GetYaxis()->SetRangeUser(limitYEmin,limitYE);
            Dethist_E_R_e[R][e]->GetXaxis()->SetRangeUser(XEmin,XEmax);
        }
    }

/*
    // Draw histograms for hist_pT_R[3]
    TCanvas* canvas_pT_R = new TCanvas("canvas_pT_R", "Histogram Canvas (pT_R)", 900, 600);
    canvas_pT_R->SetWindowSize(900, 600);
    canvas_pT_R->SetCanvasSize(900, 600);

    canvas_pT_R->cd();


    hist_pT_R[0]->Draw("");

    TLegend* legend_pT_R = new TLegend(0.58, 0.7, 0.75, 0.85);
    legend_pT_R->AddEntry(hist_pT_R[0], Form("R = %.1f", Rvals[0]), "l");
    for (int R = 1; R < 3; ++R) {
        hist_pT_R[R]->Draw(" SAME");
        legend_pT_R->AddEntry(hist_pT_R[R], Form("R = %.1f", Rvals[R]), "l");

    }

    legend_pT_R->SetTextSize(0.035);
    legend_pT_R->SetBorderSize(0);
    //legend_pT_R->SetFillColor(0);

    legend_pT_R->Draw(); // Add legend to the canvas
    ChangeTitleOfCanvas(canvas_pT_R, "particle level jet pT");
    canvas_pT_R->SaveAs("figs/Spectra/TruthpartLevelcanvas_pT_R.pdf");*/

   

    // Draw the hist_pT_R_e histograms for each R value in separate canvases and one plot each
    for (int R = 0; R < 3; ++R) {
        TCanvas* canvas_pT_R_e = new TCanvas(Form("canvas_pT_R%d_e", R), Form("Histogram Canvas (pT_R = %d, e)", R), 900, 600);
        canvas_pT_R_e->SetWindowSize(900, 600);
        canvas_pT_R_e->SetCanvasSize(900, 600);

        TF1  *f1 = new TF1("f1","[0]*(1./TMath::Power(x,[1]))");
        TF1  *f2 = new TF1("f2","[0]*(1./TMath::Power(x,[1]))");
        TF1  *f3 = new TF1("f3","[0]*(1./TMath::Power(x,[1]))");
        TF1  *f4 = new TF1("f4","[0]*(1./TMath::Power(x,[1]))");

        canvas_pT_R_e->cd();
        hist_pT_R_e[R][0]->Draw("hist");
        hist_pT_R_e[R][0]->Fit("f1","","",XpTminFit,XpTmaxFit);

        TLegend* legend_pT_R_e = new TLegend(0.62, 0.72, 0.82, 0.87);
        TLegend* Detlegend_pT_R_e = new TLegend(0.62, 0.54, 0.82, 0.69);

        legend_pT_R_e->AddEntry((TObject*)0,"Particle level jets", "");
        legend_pT_R_e->AddEntry(hist_pT_R_e[R][0], Form("%s : n = %0.1f", etaRange[0].Data(), f1->GetParameter(1)), "l");

        for (int e = 1; e < 2; ++e) {
            hist_pT_R_e[R][e]->Draw("hist SAME");
        }
        hist_pT_R_e[R][1]->Fit("f2","","",XpTminFit,70);
        legend_pT_R_e->AddEntry(hist_pT_R_e[R][1], Form("%s : n = %0.1f", etaRange[1].Data(), f2->GetParameter(1)), "l");

        legend_pT_R_e->SetTextSize(0.032);
        legend_pT_R_e->SetBorderSize(0);
        //legend_pT_R_e->SetFillColor(0);
        legend_pT_R_e->Draw(); // Add legend to the canvas



        canvas_pT_R_e->cd();
        Detlegend_pT_R_e->AddEntry((TObject*)0,"Detector level jets", "");

        for (int e = 0; e < 2; ++e) {
            Dethist_pT_R_e[R][e]->Draw("hist SAME");
        }
        
        Dethist_pT_R_e[R][0]->Fit("f3","","",XpTminFit,XpTmaxFit);
        Dethist_pT_R_e[R][1]->Fit("f4","","",XpTminFit,70);

        Detlegend_pT_R_e->AddEntry(Dethist_pT_R_e[R][0], Form("%s : n = %0.1f", etaRange[0].Data(), f3->GetParameter(1)), "l");    
        Detlegend_pT_R_e->AddEntry(Dethist_pT_R_e[R][1], Form("%s : n = %0.1f", etaRange[1].Data(), f4->GetParameter(1)), "l");    

        Detlegend_pT_R_e->SetTextSize(0.032);
        Detlegend_pT_R_e->SetBorderSize(0);
        //legend_pT_R_e->SetFillColor(0);
        Detlegend_pT_R_e->Draw(); // Add legend to the canvas



        ChangeTitleOfCanvas(canvas_pT_R_e, Form("particle level jet pT, R = %0.1f", Rvals[R]));
        canvas_pT_R_e->SaveAs(Form("figs/Spectra/BothpartLevelcanvas_pT_R%d_e.pdf", R));
    }



    
}



void NormalizeToProb(TH1F* hHisto){
    Double_t factor = 1.;

    //if (hHisto!=NULL) hHisto->Scale(factor/hHisto->Integral(), "");
    //if (hHisto!=NULL) hHisto->Scale(factor/8, "");
    //Add the following line to scale to picobarns
    if (hHisto!=NULL) hHisto->Scale(factor/scalingFactor, "width");//hHisto->SetLineWidth(1);
    if (hHisto!=NULL) hHisto->Rebin(2);
}

void ChangeTitleOfCanvas(TCanvas* cCanvas, TString sTitle){
    cCanvas->cd();
    TLatex *   tex = new TLatex(0.31,0.92,Form("%s", sTitle.Data()));
    tex->SetNDC();
    tex->SetTextSize(0.038);
    //tex->Draw();

        drawLatexAdd("ALICE simulation, pp #sqrt{#it{s}} = 14 TeV",0.45,0.365-1*1.1*textSize, textSize,kFALSE, kFALSE, kTRUE);  
        drawLatexAdd("FoCal upgrade",0.45,0.365-2*1.1*textSize, textSize,kFALSE, kFALSE, kTRUE);  
        drawLatexAdd("jets, anti-#it{k}_{T}, R=0.6",0.45,0.365-3*1.1*textSize, textSize,kFALSE, kFALSE, kTRUE);  
        //drawLatexAdd("R=0.6",0.95,0.965-4.2*1.1*textSize, textSize,kFALSE, kFALSE, kTRUE);  

    if(SettingLogY==1) cCanvas->SetLogy();
}
