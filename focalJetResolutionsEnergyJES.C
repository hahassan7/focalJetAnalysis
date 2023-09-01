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

double textSize = 0.05; // 5;

// Macro to plot the energy binned JES and JER as well as response matrices and projections
// Base macro provided to me by F. Bock

void focalJetResolutionsEnergyJES(TString inputfileR02 = "Data20230816/JES/EneMerged_OutputR6_0Mass.root",
                                  TString inputfileR04 = "Data20230816/JES/EneMerged_OutputR6_0Mass.root")
{
  StyleSettingsPaper();
  TGaxis::SetMaxDigits(4);

  int RMax = 1;       // how many R values you want to draw
  int JESMinX = 200;  // xmin for the jes and jer plots, also used for other E min
  int JESMaxX = 3050; // xmin for the jes and jer plots, also used for other E max

  Color_t colors[3] = {kBlue + 2, kRed + 1, kBlack};   // kBlack,  kRed+1, kBlue+2
  Color_t colorsMC[3] = {kBlack, kRed + 1, kBlue + 2}; //{kBlue+2,  kRed+1, kBlack}; //{kGray+1, kRed+3, kBlue+3};
  Style_t marker[3] = {20, 21, 33};
  Style_t markerMC[3] = {33, 20, 21}; // 24, 25, 27
  Style_t style[3] = {1, 5, 7};
  Size_t markerS[3] = {2, 2, 2}; // 2.4

  TString radiusOut[2] = {"R06", "R04"};
  TString radiusLabel[2] = {"#it{R} = 0.6", "#it{R} = 0.4"};
  TString etaOut[3] = {"Full", "40to45", "45to49"};
  TString etaRange[3] = {"4.0 < #it{#eta}_{jet} < 4.9", "4.0 < #it{#eta}_{jet} < 4.5", "4.5 < #it{#eta}_{jet} < 4.9"}; //{3.8, 4.5, 5.1}  = {"3.4+R < #it{#eta}_{jet} < 5.5-R", "4.0 < #it{#eta}_{jet} < 4.5", "4.5 < #it{#eta}_{jet} < 4.9"};
  double rangeJES[2][2] = {{-0.35, 0.1}, {-0.35, 0.1}};
  const int maxNEbins = 14;
  int exampleBins[3] = {7, 9, 11};
  // double binningE[16]  = {0.0, 50.0, 100.0, 200.0, 300.0, 400.0, 500.0, 600.0, 700.0, 800.0, 900.0, 1000.0, 1250.0, 1500.0, 1750.0, 2000.0};
  double binningE[15] = {100.0, 200.0, 300.0, 400.0, 500.0, 600.0, 700.0, 800.0, 900.0, 1000.0, 1250.0, 1500.0, 1750.0, 2000.0, 3000.0}; // const int nEBins  = 15;

  TH2D *histResponseMat_E[2][3];
  TH1D *histMean_E[2][3];
  TH1D *histMedian_E[2][3];
  TH1D *histSigma_E[2][3];
  TH1D *histDeltaE_bins[2][3][maxNEbins];

  TH1D *hMeanENEF[2];
  TH1D *hSDENEF[2];

  TFile *fileR02 = new TFile(inputfileR02.Data());
  hMeanENEF[0] = (TH1D *)fileR02->Get("hMeanENEF0");
  hMeanENEF[1] = (TH1D *)fileR02->Get("hMeanENEF1");
  hSDENEF[0] = (TH1D *)fileR02->Get("hSDENEF0");
  hSDENEF[1] = (TH1D *)fileR02->Get("hSDENEF1");
  for (int i = 0; i < 3; i++)
  {
    if (i == 0)
    {

      histMean_E[0][i] = (TH1D *)fileR02->Get("hMeanE");
      histSigma_E[0][i] = (TH1D *)fileR02->Get("hSDE");
    }
    std::cout << "R = 0.2 \t" << histMean_E[0][i] << "\t" << histSigma_E[0][i] << "\t" << histResponseMat_E[0][i] << std::endl;
  }

  TCanvas *cJES = new TCanvas("cJES", "cJES", 2 * 600, 2 * 400);
  DrawPaperCanvasSettings(cJES, 0.081, 0.011, 0.013, 0.099);
  cJES->cd();

  TH2D *jesframe = new TH2D("jesframe", "", 5000, JESMinX, JESMaxX, 2000, -0.7, 0.2);
  SetStyleHistoTH2ForGraphs(jesframe, "#it{E}_{part} (GeV)", "JES", 0.85 * textSize, textSize, 0.85 * textSize, textSize, 0.92, 0.8, 510, 505, 42, 62);
  jesframe->GetXaxis()->SetLabelOffset(0.004);
  jesframe->GetYaxis()->SetLabelOffset(0.006);

  TCanvas *cJER = new TCanvas("cJER", "cJER", 2 * 600, 2 * 400);
  DrawPaperCanvasSettings(cJER, 0.081, 0.011, 0.013, 0.097);
  cJER->cd();

  TH2D *jerframe = new TH2D("jerframe", "", 5000, JESMinX, JESMaxX, 2000, 0.075, 0.22);
  SetStyleHistoTH2ForGraphs(jerframe, "#it{E}_{part} (GeV)", "JER", 0.85 * textSize, textSize, 0.85 * textSize, textSize, 0.92, 0.8, 510, 505, 42, 62);
  jerframe->GetXaxis()->SetLabelOffset(0.004);
  jerframe->GetYaxis()->SetLabelOffset(0.006);

  for (int e = 0; e < 1; e++)
  {
    cJES->Clear();
    jesframe->GetYaxis()->SetRangeUser(rangeJES[0][0], rangeJES[1][1]);
    jesframe->Draw("axis");

    TLegend *leg2 = GetAndSetLegend2(0.15, 0.9 - 3 * 1.0 * textSize, 0.25, 0.93, textSize, 1, "", 42, 0.35);
    for (int r = 0; r < RMax; r++)
    {
      if (r == 1)
        DrawSetMarker(histMean_E[r][e], marker[e], markerS[e], colors[e], colors[e]);
      else
      {
        DrawSetMarker(histMean_E[r][e], markerMC[e], markerS[e] * 1.2, colorsMC[e], colorsMC[e]);
        DrawSetMarker(hMeanENEF[0], marker[e], markerS[e] * 1.2, colors[e], colors[e]);
        DrawSetMarker(hMeanENEF[1], marker[e + 1], markerS[e] * 1.2, colors[e + 1], colors[e + 1]);
      }
      histMean_E[r][e]->Draw("same");
      hMeanENEF[0]->Draw("same");
      hMeanENEF[1]->Draw("same");
      leg2->AddEntry(histMean_E[r][e], "mean", "p"); // leg2->AddEntry(histMean_E[r][e], radiusLabel[r].Data(), "p");
      leg2->AddEntry(hMeanENEF[0], "mean, det NEF < 1/3", "p");
      leg2->AddEntry(hMeanENEF[1], "mean, det NEF > 2/3", "p");
    }
    leg2->Draw("same");
    drawLatexAdd("ALICE simulation, pp #sqrt{#it{s}} = 14 TeV", 0.95, 0.965 - 1 * 1.1 * textSize, textSize, kFALSE, kFALSE, kTRUE);
    drawLatexAdd("FoCal upgrade", 0.95, 0.965 - 2 * 1.1 * textSize, textSize, kFALSE, kFALSE, kTRUE);
    drawLatexAdd(Form("jets, anti-#it{k}_{T}, #it{R} = 0.6, %s", etaRange[e].Data()), 0.95, 0.965 - 3 * 1.1 * textSize, textSize, kFALSE, kFALSE, kTRUE);
    // drawLatexAdd("#it{R} = 0.6", 0.95, 0.965 - 4.2 * 1.1 * textSize, textSize, kFALSE, kFALSE, kTRUE);
    drawLatexAdd("#Delta#it{E} = (#it{E}_{det} #font[122]{-} #it{E}_{part})/#it{E}_{part}", 0.95, 0.965 - 4.5 * 1.1 * textSize, textSize, kFALSE, kFALSE, kTRUE);

    jesframe->Draw("axis,same");

    cJES->cd();
    cJES->Update();
    cJES->SaveAs(Form("figs/Energy/EnJES_%s.pdf", etaOut[e].Data()));

    cJER->Clear();
    jerframe->Draw("axis");
    TLegend *leg3 = GetAndSetLegend2(0.3, 0.75 - 4 * 1.1 * textSize, 0.5, 0.75 - 1 * 1.1 * textSize, textSize, 1, "", 42, 0.3);
    for (int r = 0; r < RMax; r++)
    {
      if (r == 1)
        DrawSetMarker(histSigma_E[r][e], marker[e], markerS[e], colors[e], colors[e]);
      else
      {
        DrawSetMarker(histSigma_E[r][e], markerMC[e], markerS[e] * 1.2, colorsMC[e], colorsMC[e]);
        DrawSetMarker(hSDENEF[0], marker[e], markerS[e] * 1.2, colors[e], colors[e]);
        DrawSetMarker(hSDENEF[1], marker[e + 1], markerS[e] * 1.2, colors[e + 1], colors[e + 1]);
      }
      histSigma_E[r][e]->Draw("same");
      hSDENEF[0]->Draw("same");
      hSDENEF[1]->Draw("same");
      // leg3->AddEntry((TObject *)0, radiusLabel[r].Data(), "");
      leg3->AddEntry(histSigma_E[r][e], "standard deviation", "p");
      leg3->AddEntry(hSDENEF[0], "standard deviation, det NEF < 1/3 ", "p");
      leg3->AddEntry(hSDENEF[1], "standard deviation, det NEF > 2/3 ", "p");
    }
    leg3->Draw("same");
    drawLatexAdd("ALICE simulation, pp #sqrt{#it{s}} = 14 TeV", 0.95, 0.965 - 1 * 1.1 * textSize, textSize, kFALSE, kFALSE, kTRUE);
    drawLatexAdd("FoCal upgrade", 0.95, 0.965 - 2 * 1.1 * textSize, textSize, kFALSE, kFALSE, kTRUE);
    drawLatexAdd(Form("jets, anti-#it{k}_{T}, #it{R} = 0.6, %s", etaRange[e].Data()), 0.95, 0.965 - 3 * 1.1 * textSize, textSize, kFALSE, kFALSE, kTRUE);
    jerframe->Draw("axis,same");

    cJER->cd();
    cJER->Update();
    cJER->SaveAs(Form("figs/Energy/EnJER_%s.pdf", etaOut[e].Data()));
  }
}
