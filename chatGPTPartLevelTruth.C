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
void chatGPTPartLevelTruth() {

    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gStyle->SetHistLineWidth(5);


    const double limitYconst = 1e10;//0.2;
    const double limitYpT    = 1e12;//0.5;
    const double limitYE     = 1e11;//0.011;
    const double limitYconstmin = 0.000001;
    const double limitYpTmin    = 100;
    const double limitYEmin     = 10000;
    const double XpTmin     = 10; //change to 10 unless lowest bin, then 5
    const double XpTmax     = 200;
    const double XEmin     = 10; 
    const double XEmax     = 2000;
    //->GetYaxis()->SetRangeUser

    // Arrays to hold histograms for different types and R values
    TH1F* hist_pT_R[3]; // hist_pT_R[0], hist_pT_R[1], hist_pT_R[2]
    TH1F* hist_E_R[3]; // hist_E_R[0], hist_E_R[1], hist_E_R[2]
    TH1F* hist_con_R[3]; // hist_con_R[0], hist_con_R[1], hist_con_R[2]

    TH1F* hist_pT_R_e[3][2]; // hist_pT_R_e[0][0], hist_pT_R_e[0][1], ..., hist_pT_R_e[2][2]
    TH1F* hist_E_R_e[3][2]; // hist_E_R_e[0][0], hist_E_R_e[0][1], ..., hist_E_R_e[2][2]
    TH1F* hist_con_R_e[3][2]; // hist_con_R_e[0][0], hist_con_R_e[0][1], ..., hist_con_R_e[2][2]

    //ADD TWO MORE HISTO ARRAYs to hold pT and E distribution for const = 1,2,3,5,10. 
    TH1F* hist_pTvscon_R_e[3][2][5];
    TH1F* hist_Evscon_R_e[3][2][5];
    TH1F* hist_pTvscon_R[3][5];
    TH1F* hist_Evscon_R[3][5];

    Color_t colors[5]     = {kBlack,  kRed+1, kViolet, kCyan-1, kGray+1};
    Color_t colorsMC[3]   = {kGray+1, kRed+3, kBlue+3};
    Style_t marker[5]     = {20, 21, 33, 20, 21};
    Style_t markerMC[5]   = {24, 25, 27, 24, 25};
    Style_t style[5]      =  {1, 5, 7, 1, 5};
    Size_t  markerS[5]    = { 2, 2, 2.4, 2, 2};

    // Read histograms from the ROOT file
    TFile* file = TFile::Open("Data20230728/Spectra/Merged.root", "READ");
    //TFile* file = TFile::Open("JetJetOutput/July2023/data/truth_20230722_pythia8_JetJet_5-10GeV_Merged_Output.root", "READ");
    //TFile* file = TFile::Open("JetJetOutput/July2023/data/truth_20230722_pythia8_JetJet_10-20GeV_Merged_Output.root", "READ");
    //TFile* file = TFile::Open("JetJetOutput/July2023/data/truth_20230722_pythia8_JetJet_20-30GeV_Merged_Output.root", "READ");
    //TFile* file = TFile::Open("JetJetOutput/July2023/data/truth_20230722_pythia8_JetJet_30-40GeV_Merged_Output.root", "READ");
    //TFile* file = TFile::Open("JetJetOutput/July2023/data/truth_20230722_pythia8_JetJet_40-60GeV_Merged_Output.root", "READ");
    //TFile* file = TFile::Open("JetJetOutput/July2023/data/truth_20230722_pythia8_JetJet_60-100GeV_Merged_Output.root", "READ");
    //TFile* file = TFile::Open("JetJetOutput/July2023/data/truth_20230722_pythia8_JetJet_100-200GeV_Merged_Output.root", "READ");
    //TFile* file = TFile::Open("JetJetOutput/July2023/data/truth_20230722_pythia8_JetJet_200-GeV_Merged_Output.root", "READ");


    const Int_t nR = 3; 
    const Float_t Rvals[nR] = {0.2, 0.4, 0.6}; // Cone radii
    const Int_t nEtaBins  = 3;
    //const double EtaBinBorders[nEtaBins] = {3.4,5.0,5.5}; //
    //const double EtaMin[2] = {3.4,3.9}; 
    //const double EtaMax[2] = {5.0,5.5};
    TString etaRange[nEtaBins-1]     =  {"3.4+R< #eta_{jet} < 4.5","4.5 < #eta_{jet} < 5.5-R"}; //consider making this different range for R=0.6, 0.4, 0.2. ...
    const Int_t nConst = 5;
    const Int_t ConstVals[nConst] = {1, 2, 3, 5, 10}; // constituent values to check


    for (int R = 0; R < 3; ++R) {
        // Read hist_pT_R[3] histograms from the file using Clone method
        hist_pT_R[R] = (TH1F*)(file->Get(Form("hist_pT_R%d", R))->Clone());
        hist_pT_R[R]->SetLineColor(colors[R]);
        NormalizeToProb(hist_pT_R[R]);
        hist_pT_R[R]->SetTitle(Form("pT Distribution for R = %d; pT [GeV/c]; d#sigma/dp_{T}dy", R));
        hist_pT_R[R]->GetYaxis()->SetRangeUser(limitYpTmin,limitYpT);
        hist_pT_R[R]->GetXaxis()->SetRangeUser(XpTmin,XpTmax);
        //normalize all the histograms to the NofEntries (check if 1/Integral is enough etc. or nentries or what not)

        // Read hist_E_R[3] histograms from the file using Clone method
        hist_E_R[R] = (TH1F*)(file->Get(Form("hist_E_R%d", R))->Clone());
        hist_E_R[R]->SetLineColor(colors[R]);
        NormalizeToProb(hist_E_R[R]);
        hist_E_R[R]->SetTitle(Form("E Distribution for R = %d; E [GeV]; d#sigma/dEdy", R));
        hist_E_R[R]->GetYaxis()->SetRangeUser(limitYEmin,limitYE);
        hist_E_R[R]->GetXaxis()->SetRangeUser(XEmin,XEmax);

        // Read hist_con_R[3] histograms from the file using Clone method
        hist_con_R[R] = (TH1F*)(file->Get(Form("hist_con_R%d", R))->Clone());
        hist_con_R[R]->SetLineColor(colors[R]);
        NormalizeToProb(hist_con_R[R]);
        hist_con_R[R]->SetTitle(Form("Constituent Distribution for R = %d; Constituent N; d#sigma/dconstdy", R));
        hist_con_R[R]->GetYaxis()->SetRangeUser(limitYconstmin,limitYconst);

        for (int e = 0; e < 2; ++e) {
            // Read hist_pT_R_e[3][3] histograms from the file using Clone method
            hist_pT_R_e[R][e] = (TH1F*)(file->Get(Form("hist_pT_R%d_e%d", R, e))->Clone());
            hist_pT_R_e[R][e]->SetLineColor(colors[e]);
            NormalizeToProb(hist_pT_R_e[R][e]);
            hist_pT_R_e[R][e]->SetTitle(Form("pT Distribution for R = %d, e = %d; pT [GeV/c]; d#sigma/dp_{T}dy", R, e));
            hist_pT_R_e[R][e]->GetYaxis()->SetRangeUser(limitYpTmin,limitYpT);
            hist_pT_R_e[R][e]->GetXaxis()->SetRangeUser(XpTmin,XpTmax);

            // Read hist_E_R_e[3][3] histograms from the file using Clone method
            hist_E_R_e[R][e] = (TH1F*)(file->Get(Form("hist_E_R%d_e%d", R, e))->Clone());
            hist_E_R_e[R][e]->SetLineColor(colors[e]);
            NormalizeToProb(hist_E_R_e[R][e]);
            hist_E_R_e[R][e]->SetTitle(Form("E Distribution for R = %d, e = %d; E [GeV]; d#sigma/dEdy", R, e));
            hist_E_R_e[R][e]->GetYaxis()->SetRangeUser(limitYEmin,limitYE);
            hist_E_R_e[R][e]->GetXaxis()->SetRangeUser(XEmin,XEmax);

            // Read hist_con_R_e[3][3] histograms from the file using Clone method
            hist_con_R_e[R][e] = (TH1F*)(file->Get(Form("hist_con_R%d_e%d", R, e))->Clone());
            hist_con_R_e[R][e]->SetLineColor(colors[e]);
            NormalizeToProb(hist_con_R_e[R][e]);
            hist_con_R_e[R][e]->SetTitle(Form("Constituent Distribution for R = %d, e = %d; Constituent N; d#sigma/dconstdy", R, e));
            hist_con_R_e[R][e]->GetYaxis()->SetRangeUser(limitYconstmin,limitYconst);

            for (int c = 0; c < 5; ++c) {
                hist_pTvscon_R_e[R][e][c] = (TH1F*)(file->Get(Form("hist_pTvscon_R%d_e%d_c%d", R, e, c))->Clone());
                hist_pTvscon_R_e[R][e][c]->SetLineColor(colors[c]);
                NormalizeToProb(hist_pTvscon_R_e[R][e][c]);
                hist_pTvscon_R_e[R][e][c]->SetTitle(Form("Constituent N = %d pT Distribution for particle jets, R = %0.1f, %s; pT [GeV/c]; d#sigma/dp_{T}dy", ConstVals[c], Rvals[R], etaRange[e].Data()));
                hist_pTvscon_R_e[R][e][c]->GetYaxis()->SetRangeUser(limitYpTmin,limitYpT);
                hist_pTvscon_R_e[R][e][c]->GetXaxis()->SetRangeUser(XpTmin,XpTmax);

                hist_Evscon_R_e[R][e][c] = (TH1F*)(file->Get(Form("hist_Evscon_R%d_e%d_c%d", R, e, c))->Clone());
                hist_Evscon_R_e[R][e][c]->SetLineColor(colors[c]);
                NormalizeToProb(hist_Evscon_R_e[R][e][c]);
                hist_Evscon_R_e[R][e][c]->SetTitle(Form("Constituent N = %d E Distribution for particle jets, R = %0.1f, %s; E [GeV]; d#sigma/dEdy", ConstVals[c], Rvals[R], etaRange[e].Data()));
                hist_Evscon_R_e[R][e][c]->GetYaxis()->SetRangeUser(limitYEmin,limitYE);
                hist_Evscon_R_e[R][e][c]->GetXaxis()->SetRangeUser(XEmin,XEmax);

            }
        }
        for (int c = 0; c < 5; ++c) {
            hist_pTvscon_R[R][c] = (TH1F*)(file->Get(Form("hist_pTvscon_R%d_c%d", R, c))->Clone());
            hist_pTvscon_R[R][c]->SetLineColor(colors[c]);
            NormalizeToProb(hist_pTvscon_R[R][c]);
            hist_pTvscon_R[R][c]->SetTitle(Form("Constituent N = %d pT Distribution for particle jets, R = %0.1f; pT [GeV/c]; d#sigma/dp_{T}dy", ConstVals[c], Rvals[R]));
            hist_pTvscon_R[R][c]->GetYaxis()->SetRangeUser(XpTmin,limitYpT);

            hist_Evscon_R[R][c] = (TH1F*)(file->Get(Form("hist_Evscon_R%d_c%d", R, c))->Clone());
            hist_Evscon_R[R][c]->SetLineColor(colors[c]);
            NormalizeToProb(hist_Evscon_R[R][c]);
            hist_Evscon_R[R][c]->SetTitle(Form("Constituent N = %d E Distribution for particle jets, R = %0.1f; E [GeV]; d#sigma/dEdy", ConstVals[c], Rvals[R]));
            hist_Evscon_R[R][c]->GetYaxis()->SetRangeUser(XpTmin,limitYpT);

        }

    }

    // Close the ROOT file
    //file->Close();

    // Draw histograms for hist_pT_R[3]
    TCanvas* canvas_pT_R = new TCanvas("canvas_pT_R", "Histogram Canvas (pT_R)", 900, 600);
    canvas_pT_R->SetWindowSize(900, 600);
    canvas_pT_R->SetCanvasSize(900, 600);

    canvas_pT_R->cd();

    /*TF1  *f1 = new TF1("f1","[0]*(1./TMath::Power(x,[1]))"); //[0]*TMath::Power(x-[1],[2]) ,25,XpTmax
    TF1  *f2 = new TF1("f2","[0]*(1./TMath::Power(x,[1]))");
    TF1  *f3 = new TF1("f3","[0]*(1./TMath::Power(x,[1]))");

    //hist_pT_R[0]->Fit("f1","","",25,XpTmax);
    hist_pT_R[1]->Fit("f2","","",20,120);
    //hist_pT_R[2]->Fit("f3","","",25,XpTmax);*/

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
    canvas_pT_R->SaveAs("figs/Spectra/TruthpartLevelcanvas_pT_R.png");

    // Draw histograms for hist_E_R[3]
    TCanvas* canvas_E_R = new TCanvas("canvas_E_R", "Histogram Canvas (E_R)", 900, 600);
    canvas_E_R->SetWindowSize(900, 600);
    canvas_E_R->SetCanvasSize(900, 600);

    canvas_E_R->cd();
    hist_E_R[0]->Draw("");
    TLegend* legend_E_R = new TLegend(0.58, 0.7, 0.75, 0.85);
    legend_E_R->AddEntry(hist_E_R[0], Form("R = %.1f", Rvals[0]), "l");

    for (int R = 1; R < 3; ++R) {
        hist_E_R[R]->Draw(" SAME");
        legend_E_R->AddEntry(hist_E_R[R], Form("R = %.1f", Rvals[R]), "l");
    }

    legend_E_R->SetTextSize(0.035);
    legend_E_R->SetBorderSize(0);
    //legend_E_R->SetFillColor(0);
    legend_E_R->Draw(); // Add legend to the canvas
    ChangeTitleOfCanvas(canvas_E_R, "particle level jet E");
    canvas_E_R->SaveAs("figs/Spectra/TruthpartLevelcanvas_E_R.png");

    /*// Draw histograms for hist_con_R[3]
    TCanvas* canvas_con_R = new TCanvas("canvas_con_R", "Histogram Canvas (con_R)", 900, 600);
    canvas_con_R->SetWindowSize(900, 600);
    canvas_con_R->SetCanvasSize(900, 600);

    canvas_con_R->cd();
    hist_con_R[0]->Draw("");
    TLegend* legend_con_R = new TLegend(0.7, 0.7, 0.9, 0.9);
    legend_con_R->AddEntry(hist_con_R[0], Form("R = %.1f", Rvals[0]), "l");

    for (int R = 1; R < 3; ++R) {
        hist_con_R[R]->Draw(" SAME");
        legend_con_R->AddEntry(hist_con_R[R], Form("R = %.1f", Rvals[R]), "l");
    }

    legend_con_R->Draw(); // Add legend to the canvas
    ChangeTitleOfCanvas(canvas_con_R, "particle level jet constituent N");
    canvas_con_R->SaveAs("figs/Spectra/TruthpartLevelcanvas_con_R.png");*/

    // Draw the hist_pT_R_e histograms for each R value in separate canvases and one plot each
    for (int R = 0; R < 3; ++R) {
        TCanvas* canvas_pT_R_e = new TCanvas(Form("canvas_pT_R%d_e", R), Form("Histogram Canvas (pT_R = %d, e)", R), 900, 600);
        canvas_pT_R_e->SetWindowSize(900, 600);
        canvas_pT_R_e->SetCanvasSize(900, 600);

        canvas_pT_R_e->cd();
        hist_pT_R_e[R][0]->Draw("");
        TLegend* legend_pT_R_e = new TLegend(0.58, 0.7, 0.75, 0.85);
        legend_pT_R_e->AddEntry(hist_pT_R_e[R][0], Form("%s", etaRange[0].Data()), "l");

        for (int e = 1; e < 2; ++e) {
            hist_pT_R_e[R][e]->Draw(" SAME");
            legend_pT_R_e->AddEntry(hist_pT_R_e[R][e], Form("%s", etaRange[e].Data()), "l");
        }

        legend_pT_R_e->SetTextSize(0.035);
        legend_pT_R_e->SetBorderSize(0);
        //legend_pT_R_e->SetFillColor(0);
        legend_pT_R_e->Draw(); // Add legend to the canvas
        ChangeTitleOfCanvas(canvas_pT_R_e, Form("particle level jet pT, R = %0.1f", Rvals[R]));
        canvas_pT_R_e->SaveAs(Form("figs/Spectra/TruthpartLevelcanvas_pT_R%d_e.png", R));
    }

    // Draw the hist_E_R_e histograms for each R value in separate canvases and one plot each
    for (int R = 0; R < 3; ++R) {
        TCanvas* canvas_E_R_e = new TCanvas(Form("canvas_E_R%d_e", R), Form("Histogram Canvas (E_R = %d, e)", R), 900, 600);
        canvas_E_R_e->SetWindowSize(900, 600);
        canvas_E_R_e->SetCanvasSize(900, 600);

        canvas_E_R_e->cd();
        hist_E_R_e[R][0]->Draw("");
        TLegend* legend_E_R_e = new TLegend(0.58, 0.7, 0.75, 0.85);
        legend_E_R_e->AddEntry(hist_E_R_e[R][0], Form("%s", etaRange[0].Data()), "l");

        for (int e = 1; e < 2; ++e) {
            hist_E_R_e[R][e]->Draw(" SAME");
            legend_E_R_e->AddEntry(hist_E_R_e[R][e], Form("%s", etaRange[e].Data()), "l");
        }

        legend_E_R_e->SetTextSize(0.035);
        legend_E_R_e->SetBorderSize(0);
        //legend_E_R_e->SetFillColor(0);
        legend_E_R_e->Draw(); // Add legend to the canvas
        ChangeTitleOfCanvas(canvas_E_R_e, Form("particle level jet E, R = %0.1f", Rvals[R]));
        canvas_E_R_e->SaveAs(Form("figs/Spectra/TruthpartLevelcanvas_E_R%d_e.png", R));
    }

    /*// Draw the hist_con_R_e histograms for each R value in separate canvases and one plot each
    for (int R = 0; R < 3; ++R) {
        TCanvas* canvas_con_R_e = new TCanvas(Form("canvas_con_R%d_e", R), Form("Histogram Canvas (con_R = %d, e)", R), 900, 600);
        canvas_con_R_e->SetWindowSize(900, 600);
        canvas_con_R_e->SetCanvasSize(900, 600);

        canvas_con_R_e->cd();
        hist_con_R_e[R][0]->Draw("");
        TLegend* legend_con_R_e = new TLegend(0.7, 0.7, 0.9, 0.9);
        legend_con_R_e->AddEntry(hist_con_R_e[R][0], Form("%s", etaRange[0].Data()), "l");

        for (int e = 1; e < 2; ++e) {
            hist_con_R_e[R][e]->Draw(" SAME");
            legend_con_R_e->AddEntry(hist_con_R_e[R][e], Form("%s", etaRange[e].Data()), "l");
        }

        legend_con_R_e->Draw(); // Add legend to the canvas
        ChangeTitleOfCanvas(canvas_con_R_e, Form("particle level jet constituents, R = %0.1f", Rvals[R]));
        canvas_con_R_e->SaveAs(Form("figs/Spectra/TruthpartLevelcanvas_con_R%d_e.png", R));
    }*/

   /* // Draw histograms for hist_pTvscon_R_e[3][2][5]
    for (int R = 0; R < 3; ++R) {
        for (int e = 0; e < 2; ++e) {
            TCanvas* canvas_pTvscon_R_e = new TCanvas(Form("canvas_pTvscon_R%d_e%d.png", R, e), Form("Histogram Canvas (pT vs con_R = %d, e)", R), 900, 600);
            canvas_pTvscon_R_e->SetWindowSize(900, 600);
            canvas_pTvscon_R_e->SetCanvasSize(900, 600);

            canvas_pTvscon_R_e->cd();
            hist_pTvscon_R_e[R][e][0]->Draw("");
            TLegend* legend_pTvscon_R_e = new TLegend(0.7, 0.7, 0.9, 0.9);
            legend_pTvscon_R_e->AddEntry(hist_pTvscon_R_e[R][e][0], Form("Constituent N = %d", ConstVals[0]), "l");
            for (int c = 1; c < 5; ++c) {
                hist_pTvscon_R_e[R][e][c]->Draw(" SAME");
                legend_pTvscon_R_e->AddEntry(hist_pTvscon_R_e[R][e][c], Form("Constituent N = %d", ConstVals[c]), "l");
            }
            legend_pTvscon_R_e->Draw(); // Add legend to the canvas
            ChangeTitleOfCanvas(canvas_pTvscon_R_e, Form("particle level jet pT, R = %0.1f, %s", Rvals[R], etaRange[e].Data()));
            canvas_pTvscon_R_e->SaveAs(Form("figs/Spectra/TruthpartLevelcanvas_pTvscon_R%d_e%d.png", R, e));
        }
    }

    // Draw histograms for hist_Evscon_R_e[3][2][5]
    for (int R = 0; R < 3; ++R) {
        for (int e = 0; e < 2; ++e) {
            TCanvas* canvas_Evscon_R_e = new TCanvas(Form("canvas_Evscon_R%d_e%d", R, e), Form("Histogram Canvas (E vs con_R = %d, e)", R), 900, 600);
            canvas_Evscon_R_e->SetWindowSize(900, 600);
            canvas_Evscon_R_e->SetCanvasSize(900, 600);

            canvas_Evscon_R_e->cd();
            hist_Evscon_R_e[R][e][0]->Draw("");
            TLegend* legend_Evscon_R_e = new TLegend(0.7, 0.7, 0.9, 0.9);
            legend_Evscon_R_e->AddEntry(hist_Evscon_R_e[R][e][0], Form("Constituent N = %d", ConstVals[0]), "l");
            for (int c = 1; c < 5; ++c) {
                hist_Evscon_R_e[R][e][c]->Draw(" SAME");
                legend_Evscon_R_e->AddEntry(hist_Evscon_R_e[R][e][c], Form("Constituent N = %d", ConstVals[c]), "l");
            }
            legend_Evscon_R_e->Draw(); // Add legend to the canvas
            ChangeTitleOfCanvas(canvas_Evscon_R_e, Form("particle level jet E, R = %0.1f, %s", Rvals[R], etaRange[e].Data()));
            canvas_Evscon_R_e->SaveAs(Form("figs/Spectra/TruthpartLevelcanvas_Evscon_R%d_e%d.png", R, e));
        }
    }

    // Draw histograms for hist_pTvscon_R[3][5]
    for (int R = 0; R < 3; ++R) {
        TCanvas* canvas_pTvscon_R = new TCanvas(Form("canvas_pTvscon_R%d", R), Form("Histogram Canvas (pT vs con_R = %d)", R), 900, 600);
        canvas_pTvscon_R->SetWindowSize(900, 600);
        canvas_pTvscon_R->SetCanvasSize(900, 600);

        canvas_pTvscon_R->cd();
        hist_pTvscon_R[R][0]->Draw("");
        TLegend* legend_pTvscon_R = new TLegend(0.7, 0.7, 0.9, 0.9);
        legend_pTvscon_R->AddEntry(hist_pTvscon_R[R][0], Form("Constituent N = %d", ConstVals[0]), "l");

        for (int c = 1; c < 5; ++c) {
            hist_pTvscon_R[R][c]->Draw(" SAME");
            legend_pTvscon_R->AddEntry(hist_pTvscon_R[R][c], Form("Constituent N = %d", ConstVals[c]), "l");
        }

        legend_pTvscon_R->Draw(); // Add legend to the canvas
        ChangeTitleOfCanvas(canvas_pTvscon_R, Form("particle level jet pT, R = %0.1f", Rvals[R]));
        canvas_pTvscon_R->SaveAs(Form("figs/Spectra/TruthpartLevelcanvas_pTvscon_R%d.png", R));
    }

    // Draw histograms for hist_Evscon_R[3][5]
    for (int R = 0; R < 3; ++R) {
        TCanvas* canvas_Evscon_R = new TCanvas(Form("canvas_Evscon_R%d", R), Form("Histogram Canvas (E vs con_R = %d)", R), 900, 600);
        canvas_Evscon_R->SetWindowSize(900, 600);
        canvas_Evscon_R->SetCanvasSize(900, 600);

        canvas_Evscon_R->cd();
        hist_Evscon_R[R][0]->Draw("");
        TLegend* legend_Evscon_R = new TLegend(0.7, 0.7, 0.9, 0.9);
        legend_Evscon_R->AddEntry(hist_Evscon_R[R][0], Form("Constituent N = %d", ConstVals[0]), "l");

        for (int c = 1; c < 5; ++c) {
            hist_Evscon_R[R][c]->Draw(" SAME");
            legend_Evscon_R->AddEntry(hist_Evscon_R[R][c], Form("Constituent N = %d", ConstVals[c]), "l");
        }

        legend_Evscon_R->Draw(); // Add legend to the canvas
        ChangeTitleOfCanvas(canvas_Evscon_R, Form("particle level jet E, R = %0.1f", Rvals[R]));
        canvas_Evscon_R->SaveAs(Form("figs/Spectra/TruthpartLevelcanvas_Evscon_R%d.png", R));
    }
*/

//file->Close();
}



void NormalizeToProb(TH1F* hHisto){
    Double_t factor = 1.;
    //if (hHisto!=NULL) hHisto->Scale(factor/hHisto->Integral(), "width");
    //Add the following line to scale to picobarns
    if (hHisto!=NULL) hHisto->Scale(factor/scalingFactor, "width");
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