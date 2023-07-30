#include <iostream>
#include <TFile.h>
#include <TH1F.h>
#include <TCanvas.h>
using namespace std;

const Int_t nR = 3;
const Float_t Rvals[nR] = {0.2, 0.4, 0.6}; // Cone radii

const Int_t nEBins = 6;
const Int_t nEtaBins = 3;
const double JetEBorders[nEBins] = {100.0, 400.0, 800.0, 1200.0, 1600.0, 2000.0};
const double EtaBinBorders[nEtaBins][nR] = {{3.6, 4.5, 5.3}, {3.8, 4.5, 5.1}, {4.0, 4.5, 4.9}};
TString etaRange[nEtaBins-1] = {"3.4+R < #eta_{jet} < 4.5", "4.5 < #eta_{jet} < 5.5-R"};

const int nCol = 10;
const int gcolors[nCol] = {1, 2, 6, 4, 7, 1, 2, 4, 6, 7};
const int gmarkers[nCol] = {4, 8, 25, 21, 8, 21, 25, 4, 8, 21};
const Int_t nConst = 5;
const Int_t ConstVals[nConst] = {1, 2, 3, 5, 10}; // constituent values to check

void PlotHistogramsMultiPanel(TH1F* histos[nEBins - 1][nR], const TString& title, const TString& xTitle, const TString& outputFileName);
void PlotHistogramsMultiPad(TH1F* histos[nEBins-1][nR][nEtaBins-1], TString canvasTitle, TString xAxisTitle, TString outputFileName);
void PlotHistogramsMultiPadConst(TH1F* histos[nEBins - 1][nR][nConst], TString title, TString xTitle, TString outputFileName);
void NormalizeToProb(TH1F* hHisto);

void JetPlottingMergedFractionsEightPlotEnergyBins(int NormValue=8) {


    const double limitYconst = 0.2;
    const double limitYpT    = 0.2;
    const double limitYEta     = 0.2;
    const double limitYconstmin = 0.0;
    const double limitYpTmin    = 0.0;
    const double limitYEtamin     = 0.0;


    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    // Open the input file
    TFile *fin;
    if(NormValue==0)fin = new TFile("JetJetOutput/July2023/ptmin10/fractions/20230722_pythia8_JetJet_5-10GeV_Merged_Output.root", "READ");
    if(NormValue==1)fin = new TFile("JetJetOutput/July2023/ptmin10/fractions/20230722_pythia8_JetJet_10-20GeV_Merged_Output.root", "READ");
    if(NormValue==2)fin = new TFile("JetJetOutput/July2023/ptmin10/fractions/20230722_pythia8_JetJet_20-30GeV_Merged_Output.root", "READ");
    if(NormValue==3)fin = new TFile("JetJetOutput/July2023/ptmin10/fractions/20230722_pythia8_JetJet_30-40GeV_Merged_Output.root", "READ");
    if(NormValue==4)fin = new TFile("JetJetOutput/July2023/ptmin10/fractions/20230722_pythia8_JetJet_40-60GeV_Merged_Output.root", "READ");
    if(NormValue==5)fin = new TFile("JetJetOutput/July2023/ptmin10/fractions/20230722_pythia8_JetJet_60-100GeV_Merged_Output.root", "READ");
    if(NormValue==6)fin = new TFile("JetJetOutput/July2023/ptmin10/fractions/20230722_pythia8_JetJet_100-200GeV_Merged_Output.root", "READ");
    if(NormValue==7)fin = new TFile("JetJetOutput/July2023/ptmin10/fractions/20230722_pythia8_JetJet_200-GeV_Merged_Output.root", "READ");
    if(NormValue==8)fin = new TFile("JetJetOutput/July2023/ptmin10/fractions/Merged.root", "READ");

    // Read histograms
    TH1F* hDetECALE[nEBins - 1][nR];
    TH1F* hPartNeutralMatchE[nEBins - 1][nR];

	TH1F* hDetECALEEta[nEBins - 1][nR][nEtaBins - 1];
	TH1F* hPartNeutralMatchEEta[nEBins - 1][nR][nEtaBins - 1];

	TH1F *hPartNeutralMatchEConst[nEBins-1][nR][nConst]; 

    for (int Rvalue = 0; Rvalue < nR; Rvalue++) {
        for (int eb = 0; eb < nEBins - 1; eb++) {
            hDetECALE[eb][Rvalue] = (TH1F*)fin->Get(Form("hDetECALE_%d_R%d", eb, Rvalue));NormalizeToProb(hDetECALE[eb][Rvalue]);
            hDetECALE[eb][Rvalue]->GetYaxis()->SetRangeUser(limitYpTmin,limitYpT);
            hPartNeutralMatchE[eb][Rvalue] = (TH1F*)fin->Get(Form("hPartNeutralMatchE_%d_R%d", eb, Rvalue));NormalizeToProb(hPartNeutralMatchE[eb][Rvalue]);
            hPartNeutralMatchE[eb][Rvalue]->GetYaxis()->SetRangeUser(limitYpTmin,limitYpT);

            for (int ieta = 0; ieta < nEtaBins - 1; ++ieta) {
		        hDetECALEEta[eb][Rvalue][ieta] = (TH1F*)fin->Get(Form("hDetECALE_%dEta_%d_R%d", eb, ieta, Rvalue));NormalizeToProb(hDetECALEEta[eb][Rvalue][ieta]);
		        hDetECALEEta[eb][Rvalue][ieta]->GetYaxis()->SetRangeUser(limitYEtamin,limitYEta);

		        hPartNeutralMatchEEta[eb][Rvalue][ieta] = (TH1F*)fin->Get(Form("hPartNeutralMatchE_%dEta_%d_R%d", eb, ieta, Rvalue));NormalizeToProb(hPartNeutralMatchEEta[eb][Rvalue][ieta]);
		        hPartNeutralMatchEEta[eb][Rvalue][ieta]->GetYaxis()->SetRangeUser(limitYEtamin,limitYEta);
            }

            for (int iconst = 0; iconst < nConst; ++iconst) {
		        hPartNeutralMatchEConst[eb][Rvalue][iconst] = (TH1F*)fin->Get(Form("hPartNeutralMatchE_%dConst_%d_R%d", eb, iconst, Rvalue));NormalizeToProb(hPartNeutralMatchEConst[eb][Rvalue][iconst]);
		        hPartNeutralMatchEConst[eb][Rvalue][iconst]->GetYaxis()->SetRangeUser(limitYconstmin,limitYconst);
            }
        }
    }
    PlotHistogramsMultiPad(hDetECALEEta, "ECAL Fraction vs. Eta for Different R Values", "E_{jet}/E_{jet}^{det}", "fractionFigs/ptmin10/fractionsEnergy/EtaDetECALMultiPanel");
    PlotHistogramsMultiPad(hPartNeutralMatchEEta, "Neutral Fraction vs. Eta for Different R Values", "E_{neutral}/E_{jet}^{part}", "fractionFigs/ptmin10/fractionsEnergy/EtaPartNeutralMatchMultiPanel");

    PlotHistogramsMultiPadConst(hPartNeutralMatchEConst, "Neutral Fraction vs. Constituent Bins for Different R Values", "E_{neutral}/E_{jet}^{part}", "fractionFigs/ptmin10/fractionsEnergy/ConstPartNeutralMatchMultiPanel");

    // Plot histograms using multi-panel
    PlotHistogramsMultiPanel(hDetECALE, "Matched detector level jet ECAL energy fraction distribution", "E_{ECAL}/E_{jet}^{det}", "fractionFigs/ptmin10/fractionsEnergy/DetECALMultiPanel.png");
    PlotHistogramsMultiPanel(hPartNeutralMatchE, "Matched particle level jet neutral energy fraction distribution", "E_{neutral}/E_{jet}^{part}", "fractionFigs/ptmin10/fractionsEnergy/PartNeutralMatchMultiPanel.png");

}

void PlotHistogramsMultiPanel(TH1F* histos[nEBins - 1][nR], const TString& title, const TString& xTitle, const TString& outputFileName) {
    TCanvas* canvas = new TCanvas("canvas", title, 1200, 800);
    canvas->Divide(3, 2);

    for (int eb = 0; eb < nEBins - 1; eb++) {
        canvas->cd(eb + 1);
        // Create and draw the legend for this pad
        TLegend* legend = new TLegend(0.6, 0.75, 0.9, 0.9);
        legend->SetTextSize(0.04);
        legend->SetBorderSize(0);

        for (int r = 0; r < nR; r++) {
            gPad->SetLeftMargin(0.15);
            gPad->SetRightMargin(0.05);
            gPad->SetBottomMargin(0.15);

            if (histos[eb][r] != nullptr) {
                histos[eb][r]->SetLineColor(gcolors[r]);
                histos[eb][r]->SetMarkerStyle(gmarkers[r]);
                histos[eb][r]->GetXaxis()->SetTitle(xTitle);
                histos[eb][r]->Draw("same");

                if (r == 0) {
                    histos[eb][r]->GetYaxis()->SetTitle("Normalized");
                    histos[eb][r]->GetYaxis()->SetTitleOffset(1.3);
                }
            }

            if (histos[eb][r] != nullptr) {
                legend->AddEntry(histos[eb][r], Form("R = %0.1f", Rvals[r]), "lp");
            }
	        
	        legend->Draw();
    	}
        // Add title with the corresponding energy bin
        TLatex* titleLatex = new TLatex(0.15, 0.92, Form("E_{jet} = %0.0f - %0.0f", JetEBorders[eb], JetEBorders[eb + 1]));
        titleLatex->SetNDC();
        titleLatex->SetTextSize(0.04);
        titleLatex->Draw();

        // Update the pad
        canvas->Update();

    }

    canvas->SaveAs(outputFileName);
    delete canvas;
}

void PlotHistogramsMultiPad(TH1F* histos[nEBins-1][nR][nEtaBins-1], TString canvasTitle, TString xAxisTitle, TString outputFileName) {
    // Loop over R values
    for (int Rvalue = 0; Rvalue < nR; Rvalue++) {
        // Create a canvas for each R value
        TCanvas* c = new TCanvas(Form("c_R%d", Rvalue), canvasTitle, 1200, 800);
        c->Divide(3, 2); // Three columns and two rows for nEBins-1

        // Loop over energy bins
        for (int eb = 0; eb < nEBins - 1; ++eb) {
            c->cd(eb + 1);

            // Create a legend for each energy bin
            TLegend* legend = new TLegend(0.7, 0.7, 0.9, 0.9);
            legend->SetTextSize(0.035);
            legend->SetBorderSize(0);
            legend->SetFillColor(0);

            // Plot histograms for both eta bins in the same pad
            histos[eb][Rvalue][0]->SetLineColor(gcolors[0]);
            histos[eb][Rvalue][0]->SetMarkerStyle(gmarkers[0]);
            histos[eb][Rvalue][0]->SetMarkerColor(gcolors[0]);
            histos[eb][Rvalue][0]->GetXaxis()->SetTitle(xAxisTitle);
            histos[eb][Rvalue][0]->Draw("EP");
            legend->AddEntry(histos[eb][Rvalue][0], Form("%s",etaRange[0].Data()), "lep");

            histos[eb][Rvalue][1]->SetLineColor(gcolors[1]);
            histos[eb][Rvalue][1]->SetMarkerStyle(gmarkers[1]);
            histos[eb][Rvalue][1]->SetMarkerColor(gcolors[1]);
            histos[eb][Rvalue][1]->Draw("EP SAME");
            legend->AddEntry(histos[eb][Rvalue][1], Form("%s",etaRange[1].Data()), "lep");

            // Draw the legend
	            legend->Draw();

	        // Add title with the corresponding energy bin
	        TLatex* titleLatex = new TLatex(0.15, 0.92, Form("E_{jet} = %0.0f - %0.0f", JetEBorders[eb], JetEBorders[eb + 1]));
	        titleLatex->SetNDC();
	        titleLatex->SetTextSize(0.04);
	        titleLatex->Draw();

            // Update the pad
            c->Update();
        }

        // Save the canvas as an image file
        TString outputFile = outputFileName + Form("_R%d.png", Rvalue);
        c->Print(outputFile);
        delete c;
    }
}

void PlotHistogramsMultiPadConst(TH1F* histos[nEBins - 1][nR][nConst], TString title, TString xTitle, TString outputFileName) {
    const int nEBins = 6; // Number of energy bins

    // Loop over all R values
    for (int Rvalue = 0; Rvalue < nR; Rvalue++) {
        // Create a canvas for each R value
        TCanvas* cMultiPad = new TCanvas(Form("cMultiPad_R%d", Rvalue), title, 1200, 800);
        cMultiPad->Divide(2, 3); // Divide canvas into 2x3 pads (for 6 energy bins)

        // Loop over all energy bins
        for (int eb = 0; eb < nEBins - 1; eb++) {
            cMultiPad->cd(eb + 1); // Select the pad for this energy bin

            // Create a legend for this pad
            TLegend* legend = new TLegend(0.6, 0.7, 0.85, 0.85);
            legend->SetBorderSize(0);
            legend->SetTextSize(0.03);

            // Loop over all constituent bins
            for (int iconst = 0; iconst < nConst; iconst++) {
                // Set line color and marker style using gcolors and gmarkers arrays
                histos[eb][Rvalue][iconst]->SetLineColor(gcolors[iconst]);
                histos[eb][Rvalue][iconst]->SetMarkerStyle(gmarkers[iconst]);

                // Draw the histograms on the same pad
                if (iconst == 0) {
                    histos[eb][Rvalue][iconst]->Draw();
                } else {
                    histos[eb][Rvalue][iconst]->Draw("same");
                }

                // Add legend entries
                legend->AddEntry(histos[eb][Rvalue][iconst], Form("Const %d", ConstVals[iconst]), "lp");
            }

            // Draw the legend in this pad
            legend->Draw();

            // Add title with energy bin info
            TLatex* titleLatex = new TLatex(0.15, 0.9, Form("E_{jet} = %0.0f - %0.0f", JetEBorders[eb], JetEBorders[eb + 1]));
            titleLatex->SetNDC();
            titleLatex->SetTextSize(0.04);
            titleLatex->Draw();

            // Set X and Y axis titles for the pad
            histos[eb][Rvalue][0]->GetXaxis()->SetTitle(xTitle);
            histos[eb][Rvalue][0]->GetYaxis()->SetTitle("Normalized");

            // Update the pad
            cMultiPad->Update();
        }

        // Save the canvas as a PNG image with a unique name based on the R value
        TString outputFileName_R = Form("%s_R%d.png", outputFileName.Data(), Rvalue);
        cMultiPad->SaveAs(outputFileName_R);

        // Delete the canvas to free up memory
        delete cMultiPad;
    }
}


void NormalizeToProb(TH1F* hHisto){
    Double_t factor = 1.;
    if (hHisto!=NULL) hHisto->Scale(factor/hHisto->Integral());//, "width");
}
