#include "TFile.h"
#include "TH1.h"

void PlotDetResolution()
{
    //TFile *fileSpaghetti = TFile::Open("AnalysisClusters_BoxPiPlus_Spaghetti.root");
    //TFile *fileSandwich = TFile::Open("AnalysisClusters_BoxPiPlus_Sandwich.root");
    TFile *fileSandwich = TFile::Open("/home/lmh/alice/fromPuhti/v19/20230127_PythiaMBTrig-2_v4_eta35-55/Merged.root");
    // TFile *fileSpaghetti = TFile::Open("AnalysisClusters_BoxElectrons_Spaghetti.root");
    // TFile *fileSandwich = TFile::Open("AnalysisClusters_BoxElectrons_Sandwich.root");

    TH2F *hTotalEnergyDiff = (TH2F*)fileSandwich->Get("TotalEnergyDiffHCAL_h");
    TH2F *hTotalDeltaEnergyDiff = (TH2F*)fileSandwich->Get("TotalDeltaEnergyDiffHCAL_h");

    hTotalEnergyDiff->SetTitle("HCAL (ClustE - TruthE)/TruthE for Sandwich HCAL;Truth Energy (GeV);(ClustE - TruthE)/TruthE");

    TH1D *hMeanDeltaE_Sandwich = new TH1D("hMeanDeltaE_Sandwich", "Mean #Delta E versus Truth E;Truth E (GeV);<#Delta E>", 20, 0, 1000);
    hMeanDeltaE_Sandwich->SetLineWidth(2);


    TH1D *hResolution_Sandwich = new TH1D("hResolution_Sandwich", "Resolution;Truth E (GeV);#sigma(E)/E", 20, 0, 1000);
    hResolution_Sandwich->GetYaxis()->SetRangeUser(0, 0.4);
    // hResolution_Sandwich->GetYaxis()->SetRangeUser(0, 0.12);

    TCanvas *c1 = new TCanvas("c1", "c1: response slice", 1200, 1200);
    c1->Divide(4, 5);

    for (Int_t ibin = 1; ibin <= hResolution_Sandwich->GetNbinsX(); ibin++)
    {
        c1->cd(ibin);
        // gPad->SetLogy();

        auto eresp_Sandwich = (TH1F *)hTotalDeltaEnergyDiff->ProjectionY(Form("eresp_Sandwich_%d", ibin),1+(ibin-1)*10,ibin*10);
        eresp_Sandwich->Rebin();
        eresp_Sandwich->Scale(1./eresp_Sandwich->Integral());
        eresp_Sandwich->SetTitle(Form("(ClustE - TruthE)/TruthE for %d < E < %d GeV", (ibin - 1) * 50, (ibin)*50));
        eresp_Sandwich->SetLineColor(kBlue);
        //eresp_Sandwich->Draw();
        eresp_Sandwich->GetXaxis()->SetRangeUser(-3., 6.);

        auto f2 = new TF1("f2", "crystalball", -2., 5.);
        f2->SetLineColor(kRed);
        //f2->SetParameters(eresp_Sandwich->GetMaximum(), eresp_Sandwich->GetMean(), 50, 0.6, 1.5);
        f2->SetParameters(eresp_Sandwich->GetMaximum(), 0, 0.2, -0.6, 2.);

        eresp_Sandwich->Fit("f2", "", "same", -2., 5.);


        hMeanDeltaE_Sandwich->SetBinContent(ibin, eresp_Sandwich->GetFunction("f2")->GetParameter(1));
        hMeanDeltaE_Sandwich->SetBinError(ibin, eresp_Sandwich->GetFunction("f2")->GetParError(1));

        hResolution_Sandwich->SetBinContent(ibin, eresp_Sandwich->GetFunction("f2")->GetParameter(2));
        hResolution_Sandwich->SetBinError(ibin, eresp_Sandwich->GetFunction("f2")->GetParError(2));
        
    }

    TLine *UnityLine = new TLine(0, 0, 2000, 0);
    UnityLine->SetLineWidth(2);
    UnityLine->SetLineColor(kBlack);
    UnityLine->SetLineStyle(2);

    TCanvas *c2 = new TCanvas("c2", "DeltaE vs TruthE", 800, 600);

    hTotalDeltaEnergyDiff->Draw("colz");
    hTotalDeltaEnergyDiff->GetYaxis()->SetRangeUser(-3., 6.);
    hMeanDeltaE_Sandwich->SetLineColor(kRed);
    hMeanDeltaE_Sandwich->SetLineWidth(2);
    hMeanDeltaE_Sandwich->Draw("E1,same");
    UnityLine->Draw("same");

}