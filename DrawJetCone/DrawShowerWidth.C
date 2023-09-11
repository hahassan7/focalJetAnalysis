#include "TFile.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TEllipse.h"
#include "TLatex.h"
#include "TSystem.h"
#include "AliFOCALGeometry.h"
#include "TExec.h"
#include "/Users/hadi/focalJetAnalysis/paperPlotsHeader.h"

TH2D *ConvertColRowToXY(TH2D *myhist, AliFOCALGeometry *geom, std::string name);

float textSize = 0.05;
float textSizeLeg = 0.04;

void DrawShowerWidth()
{

    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);

    TFile *myfile = TFile::Open("ShowerProfile_PiPlus_500GeV.root");

    TH2D *myhist49 = (TH2D *)gFile->Get("hColRowHCAL_energy_Total49");
    TH2D *myhist45 = (TH2D *)gFile->Get("hColRowHCAL_energy_Total45");
    TH2D *myhist40 = (TH2D *)gFile->Get("hColRowHCAL_energy_Total40");

    AliFOCALGeometry *gFocalGeometry = new AliFOCALGeometry();
    gFocalGeometry->Init(Form("%s", gSystem->ExpandPathName("$FOCAL/geometryFiles/geometry.txt")));

    TH2D *myhist49_XY = ConvertColRowToXY(myhist49, gFocalGeometry, Form("%s_XY", myhist49->GetName()));
    SetStyleHistoTH2ForGraphs(myhist49_XY, "X (cm)", "Y (cm)", 0.9 * textSize, textSize, 0.9 * textSize, textSize, 0.92, 0.95, 510, 510, 42, 62);
    myhist49_XY->SetTitle("Shower profile at #it{#eta} = 4.9;#it{X} (cm);#it{Y} (cm)");
    myhist49_XY->GetZaxis()->SetLabelOffset(0.);
    TH2D *myhist45_XY = ConvertColRowToXY(myhist45, gFocalGeometry, Form("%s_XY", myhist45->GetName()));
    SetStyleHistoTH2ForGraphs(myhist45_XY, "X (cm)", "Y (cm)", 0.9 * textSize, textSize, 0.9 * textSize, textSize, 0.92, 0.95, 510, 510, 42, 62);
    myhist45_XY->SetTitle("Shower profile at #it{#eta} = 4.5;#it{X} (cm);#it{Y} (cm)");
    myhist45_XY->GetZaxis()->SetLabelOffset(0.);
    TH2D *myhist40_XY = ConvertColRowToXY(myhist40, gFocalGeometry, Form("%s_XY", myhist40->GetName()));
    SetStyleHistoTH2ForGraphs(myhist40_XY, "X (cm)", "Y (cm)", 0.9 * textSize, textSize, 0.9 * textSize, textSize, 0.92, 0.95, 510, 510, 42, 62);
    myhist40_XY->SetTitle("Shower profile at #it{#eta} = 4.0;#it{X} (cm);#it{Y} (cm)");
    myhist40_XY->GetZaxis()->SetLabelOffset(0.);

    TFile *JetConesFile = TFile::Open("JetConesEta.root");
    TH2D *jetConeEta49 = (TH2D *)JetConesFile->Get("JetCone_eta49");
    TH2D *jetConeEta45 = (TH2D *)JetConesFile->Get("JetCone_eta45");
    TH2D *jetConeEta40 = (TH2D *)JetConesFile->Get("JetCone_eta40");

    TExec *ex1 = new TExec("ex1", "gStyle->SetPalette(57);");
    TExec *ex2 = new TExec("ex2", "gStyle->SetPalette(52);");
    TExec *ex3 = new TExec("ex3", "gStyle->SetPalette(51);");

    TCanvas *c1 = new TCanvas("c1", "Shower Shape #it{#eta} = 4.9", 800, 854);
    c1->cd();
    gPad->SetTickx();
    gPad->SetTicky();
    gPad->SetLogz();
    gPad->SetRightMargin(0.105);
    gPad->SetTopMargin(0.145);
    // myhist49->Draw("colz");
    myhist49_XY->Draw("colz");
    ex1->Draw();
    myhist49_XY->SetMinimum(1e-2);
    myhist49_XY->Draw("col,same");
    ex2->Draw();
    jetConeEta49->Draw("col,same");
    ex1->Draw();

    // TEllipse *eX0_49 = new TEllipse(23.5, 20.5, 6);
    TEllipse *eX0_49 = new TEllipse(11.25, 1.25, 15.32);
    eX0_49->SetFillStyle(0);
    eX0_49->SetLineColor(kRed);
    eX0_49->SetLineWidth(3);
    eX0_49->Draw();

    // TLatex *TexX0_49 = new TLatex(-6.0, 20, "#it{r} = 17.6cm");
    TLatex *TexX0_49 = new TLatex(-8.0, 20, "#it{r} = 15.32cm");
    TexX0_49->SetTextColor(kRed);
    // TexX0_49->SetTextSizePixels(30);
    TexX0_49->SetTextSize(textSize);
    TexX0_49->Draw();

    // TEllipse *eLong = new TEllipse(23.5, 20.5, 3.42, 2.8, -90, 90);
    TEllipse *eLong = new TEllipse(10.425, 1.25, 19.53 - 10.425 - 1.25, 7.17, -90, 90);
    eLong->SetFillStyle(0);
    eLong->SetLineWidth(3);
    // TEllipse *eShort = new TEllipse(23.5, 20.5, 1.88, 2.8, 90, 270);
    TEllipse *eShort = new TEllipse(10.425, 1.25, 10.425 - 5.88, 7.17, 90, 270);
    eShort->SetLineWidth(3);
    eShort->SetFillStyle(0);

    // eLong->Draw();
    // eShort->Draw();

    TLatex *TexEta0 = new TLatex(-55, 50.5, "Jet #it{R} = 0.6, #it{#eta} = 4.9, #it{#varphi} = 0.");
    // TLatex *TexEta0 = new TLatex(-55, 50.5, "Jet \\mathit{R}=0.6, \\mathit{\\eta}=4.9, \\mathscr{\\phi}=0.");
    TexEta0->SetTextColor(kBlack);
    TexEta0->SetTextSize(textSize);
    TexEta0->Draw();

    drawLatexAdd("ALICE simulation, Single #pi^{+}", 0.12, 1.01 - 1 * 1.1 * textSizeLeg, textSizeLeg, kFALSE, kFALSE, kFALSE);
    drawLatexAdd("FoCal upgrade", 0.12, 1.01 - 2 * 1.1 * textSizeLeg, textSizeLeg, kFALSE, kFALSE, kFALSE);
    drawLatexAdd("#it{E}_{#pi^{+}} = 500 GeV", 0.12, 1.01 - 3 * 1.1 * textSizeLeg, textSizeLeg, kFALSE, kFALSE, kFALSE);

    TCanvas *c2 = new TCanvas("c2", "Shower Shape #it{#eta} = 4.5", 800, 854);
    c2->cd();
    gPad->SetTickx();
    gPad->SetTicky();
    gPad->SetLogz();
    gPad->SetRightMargin(0.105);
    gPad->SetTopMargin(0.145);

    // myhist45->Draw("colz");
    ex1->Draw();
    myhist45_XY->Draw("colz");
    ex1->Draw();
    myhist45_XY->SetMinimum(1e-2);
    myhist45_XY->Draw("col,same");
    ex2->Draw();
    jetConeEta45->Draw("col,same");
    ex1->Draw();

    // TEllipse *eX0_45 = new TEllipse(25.5, 20.5, 6);
    TEllipse *eX0_45 = new TEllipse(16.25, 1.25, 15.32);
    eX0_45->SetFillStyle(0);
    eX0_45->SetLineColor(kRed);
    eX0_45->SetLineWidth(3);
    eX0_45->Draw();

    // TLatex *TexX0_45 = new TLatex(-1, 20, "#it{r} = 17.6cm");
    TLatex *TexX0_45 = new TLatex(-3, 20, "#it{r} = 15.32cm");
    TexX0_45->SetTextColor(kRed);
    TexX0_45->SetTextSize(textSize);
    TexX0_45->Draw();

    // TEllipse *eLong1 = new TEllipse(25.5, 20.5, 5.11, 4.1, -90, 90);
    TEllipse *eLong1 = new TEllipse(15.55, 1.25, 29.2 - 15.55 - 1.25, 10.7, -90, 90);
    eLong1->SetFillStyle(0);
    eLong1->SetLineWidth(3);
    eLong1->SetLineColor(kMagenta);
    // TEllipse *eShort1 = new TEllipse(25.5, 20.5, 2.8, 4.1, 90, 270);
    TEllipse *eShort1 = new TEllipse(15.55, 1.25, 15.55 - 8.77, 10.7, 90, 270);
    eShort1->SetLineWidth(3);
    eShort1->SetFillStyle(0);
    eShort1->SetLineColor(kMagenta);

    // eLong1->Draw();
    // eShort1->Draw();

    TLatex *TexEta1 = new TLatex(-55, 50.5, "Jet #it{R} = 0.6, #it{#eta} = 4.5, #it{#varphi} = 0.");
    TexEta1->SetTextColor(kBlack);
    TexEta1->SetTextSize(textSize);
    TexEta1->Draw();

    drawLatexAdd("ALICE simulation, Single #pi^{+}", 0.12, 1.01 - 1 * 1.1 * textSizeLeg, textSizeLeg, kFALSE, kFALSE, kFALSE);
    drawLatexAdd("FoCal upgrade", 0.12, 1.01 - 2 * 1.1 * textSizeLeg, textSizeLeg, kFALSE, kFALSE, kFALSE);
    drawLatexAdd("#it{E}_{#pi^{+}} = 500 GeV", 0.12, 1.01 - 3 * 1.1 * textSizeLeg, textSizeLeg, kFALSE, kFALSE, kFALSE);

    TCanvas *c3 = new TCanvas("c3", "Shower Shape #it{#eta} = 4.0", 800, 854);
    c3->cd();
    gPad->SetTickx();
    gPad->SetTicky();
    gPad->SetLogz();
    gPad->SetRightMargin(0.105);
    gPad->SetTopMargin(0.145);

    // myhist40->Draw("colz");
    ex1->Draw();
    myhist40_XY->Draw("colz");
    ex1->Draw();
    myhist40_XY->SetMinimum(1e-2);
    myhist40_XY->Draw("col,same");
    ex2->Draw();
    jetConeEta40->Draw("col,same");
    ex1->Draw();

    // TEllipse *eX0_40 = new TEllipse(29.5, 20.5, 6);
    TEllipse *eX0_40 = new TEllipse(26.25, 1.25, 15.32);
    eX0_40->SetFillStyle(0);
    eX0_40->SetLineColor(kRed);
    eX0_40->SetLineWidth(3);
    eX0_40->Draw();

    // TLatex *TexX0_40 = new TLatex(9, 20, "#it{r} = 17.6cm");
    TLatex *TexX0_40 = new TLatex(7, 20, "#it{r} = 15.32cm");
    TexX0_40->SetTextColor(kRed);
    // TexX0_40->SetTextSizePixels(30);
    TexX0_40->SetTextSize(textSize);
    TexX0_40->Draw();

    // TEllipse *eLong2 = new TEllipse(29.5, 20.5, 8.4, 6.8, -90, 90);
    TEllipse *eLong2 = new TEllipse(25.64, 1.25, 48.05 - 25.64 - 1.25, 17.6, -90, 90);
    eLong2->SetFillStyle(0);
    eLong2->SetLineWidth(3);
    eLong2->SetLineColor(kGreen);
    // TEllipse *eShort2 = new TEllipse(29.5, 20.5, 4.6, 6.8, 90, 270);
    TEllipse *eShort2 = new TEllipse(25.64, 1.25, 25.64 - 14.5, 17.6, 90, 270);
    eShort2->SetLineWidth(3);
    eShort2->SetFillStyle(0);
    eShort2->SetLineColor(kGreen);

    // eLong2->Draw();
    // eShort2->Draw();

    TLatex *TexEta2 = new TLatex(-55, 50.5, "Jet #it{R} = 0.6, #it{#eta} = 4.0, #it{#varphi} = 0.");
    TexEta2->SetTextColor(kBlack);
    TexEta2->SetTextSize(textSize);
    TexEta2->Draw();

    drawLatexAdd("ALICE simulation, Single #pi^{+}", 0.12, 1.01 - 1 * 1.1 * textSizeLeg, textSizeLeg, kFALSE, kFALSE, kFALSE);
    drawLatexAdd("FoCal upgrade", 0.12, 1.01 - 2 * 1.1 * textSizeLeg, textSizeLeg, kFALSE, kFALSE, kFALSE);
    drawLatexAdd("#it{E}_{#pi^{+}} = 500 GeV", 0.12, 1.01 - 3 * 1.1 * textSizeLeg, textSizeLeg, kFALSE, kFALSE, kFALSE);
}

TH2D *ConvertColRowToXY(TH2D *myhist, AliFOCALGeometry *geom, std::string name)
{

    TH2D *convertedHist = new TH2D(Form("%s", name.c_str()), "", 48, -60, 60, 48, -60, 60);

    for (int ibinX = 1; ibinX <= myhist->GetNbinsX(); ibinX++)
    {
        for (int ibinY = 1; ibinY <= myhist->GetNbinsY(); ibinY++)
        {

            double Energy = myhist->GetBinContent(ibinX, ibinY);

            float Xpos(0.), Ypos(0.), Zpos(0.);
            geom->GetXYZFromColRowSeg(ibinX - 1, ibinY - 1, 6, Xpos, Ypos, Zpos);

            convertedHist->Fill(Xpos, Ypos, Energy);
        }
    }

    return convertedHist;
}
