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

void focalJetResolutionsEnergyJES_30PercentNEF(TString inputfilePart = "Data20230816/JES/EneMerged_OutputR6_0Mass_20PercPart.root",
                                               TString inputfileDet = "Data20230816/JES/EneMerged_OutputR6_0Mass_20PercDet.root")
{
  StyleSettingsPaper();
  TGaxis::SetMaxDigits(4);

  int RMax = 1;       // how many R values you want to draw
  int JESMinX = 200;  // xmin for the jes and jer plots, also used for other E min
  int JESMaxX = 3050; // xmin for the jes and jer plots, also used for other E max

  Color_t colors[3] = {kBlue + 2, kRed + 1, kBlack};     // kBlack,  kRed+1, kBlue+2
  Color_t colorsMC[3] = {kBlack, kGreen + 2, kBlue + 2}; //{kBlue+2,  kRed+1, kBlack}; //{kGray+1, kRed+3, kBlue+3};
  Style_t marker[3] = {20, 24, 33};
  Style_t markerMC[3] = {21, 25, 33}; // 24, 25, 27
  Style_t style[3] = {1, 5, 7};
  Size_t markerS[3] = {2, 2, 2}; // 2.4

  TString radiusOut[2] = {"R06", "R04"};
  TString radiusLabel[2] = {"#it{R} = 0.6", "#it{R} = 0.4"};
  TString etaOut[3] = {"Full", "40to45", "45to49"};
  TString etaRange[3] = {"4.0 < #it{#eta}_{jet} < 4.9", "4.0 < #it{#eta}_{jet} < 4.5", "4.5 < #it{#eta}_{jet} < 4.9"}; //{3.8, 4.5, 5.1}  = {"3.4+R < #it{#eta}_{jet} < 5.5-R", "4.0 < #it{#eta}_{jet} < 4.5", "4.5 < #it{#eta}_{jet} < 4.9"};
  // double rangeJES[2][2] = {{-0.35, 0.1}, {-0.35, 0.1}};
  double rangeJES[2][2] = {{-0.35, -0.05}, {-0.35, -0.05}};

  TFile *fileDet = TFile::Open(inputfileDet.Data());
  TFile *filePart = TFile::Open(inputfilePart.Data());

  TH1D *hJES_30Perc_Lowest_Det = (TH1D *)fileDet->Get("hMeanENEF0");
  TH1D *hJES_30Perc_Highest_Det = (TH1D *)fileDet->Get("hMeanENEF1");
  TH1D *hJES_30Perc_Lowest_Part = (TH1D *)filePart->Get("hMeanENEF0");
  TH1D *hJES_30Perc_Highest_Part = (TH1D *)filePart->Get("hMeanENEF1");

  TCanvas *cJES = new TCanvas("cJES", "cJES", 2 * 600, 2 * 400);
  DrawPaperCanvasSettings(cJES, 0.081, 0.011, 0.013, 0.099);
  cJES->cd();

  TH2D *jesframe = new TH2D("jesframe", "", 5000, JESMinX, JESMaxX, 2000, -0.7, 0.2);
  SetStyleHistoTH2ForGraphs(jesframe, "#it{E}_{part} (GeV)", "JES", 0.85 * textSize, textSize, 0.85 * textSize, textSize, 0.92, 0.85, 510, 505, 42, 62);
  jesframe->GetXaxis()->SetLabelOffset(0.004);
  jesframe->GetYaxis()->SetLabelOffset(0.006);

  cJES->Clear();
  jesframe->GetYaxis()->SetRangeUser(rangeJES[0][0], rangeJES[1][1]);
  jesframe->Draw("axis");

  TLegend *leg2 = GetAndSetLegend2(0.15, 0.8 - 3 * 1.0 * textSize, 0.3, 0.93, textSize, 1, "", 42, 0.35);

  DrawSetMarker(hJES_30Perc_Highest_Det, marker[0], markerS[0] * 1.2, colors[1], colors[1]);
  DrawSetMarker(hJES_30Perc_Lowest_Det, marker[1], markerS[1] * 1.2, colors[1], colors[1]);
  DrawSetMarker(hJES_30Perc_Highest_Part, markerMC[0], markerS[0] * 1.2, colorsMC[0], colorsMC[0]);
  DrawSetMarker(hJES_30Perc_Lowest_Part, markerMC[1], markerS[0] * 1.2, colorsMC[0], colorsMC[0]);

  hJES_30Perc_Highest_Det->Draw("same");
  hJES_30Perc_Lowest_Det->Draw("same");
  hJES_30Perc_Highest_Part->Draw("same");
  hJES_30Perc_Lowest_Part->Draw("same");

  leg2->AddEntry(hJES_30Perc_Highest_Det, "Highest 20\% NEF_{det}", "p"); // leg2->AddEntry(histMean_E[r][0], radiusLabel[r].Data(), "p");
  leg2->AddEntry(hJES_30Perc_Lowest_Det, "Lowest 20\% NEF_{det}", "p");
  leg2->AddEntry(hJES_30Perc_Highest_Part, "Highest 20\% NEF_{part}", "p");
  leg2->AddEntry(hJES_30Perc_Lowest_Part, "Lowest 20\% NEF_{part}", "p");

  leg2->Draw("same");

  drawLatexAdd("ALICE simulation, pp #sqrt{#it{s}} = 14 TeV", 0.95, 0.965 - 1 * 1.1 * textSize, textSize, kFALSE, kFALSE, kTRUE);
  drawLatexAdd("FoCal upgrade", 0.95, 0.965 - 2 * 1.1 * textSize, textSize, kFALSE, kFALSE, kTRUE);
  drawLatexAdd(Form("jets, anti-#it{k}_{T}, #it{R} = 0.6, %s", etaRange[0].Data()), 0.95, 0.965 - 3 * 1.1 * textSize, textSize, kFALSE, kFALSE, kTRUE);
  drawLatexAdd("#Delta#it{E} = (#it{E}_{det} #font[122]{-} #it{E}_{part})/#it{E}_{part}", 0.95, 0.965 - 4.5 * 1.1 * textSize, textSize, kFALSE, kFALSE, kTRUE);

  jesframe->Draw("axis,same");

  cJES->cd();
  cJES->Update();
  cJES->SaveAs(Form("figs/Energy/EnJES_%s_20Perc.pdf", etaOut[0].Data()));
}
