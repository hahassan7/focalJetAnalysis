#include "TFile.h"
#include "TH1.h"

void PlotDetResolutionFractions(int Rvalue = 0)
{
    const Int_t nR = 3; //5
    const Float_t Rvals[nR] = {0.2, 0.4, 0.6};//{0.2, 0.3, 0.4, 0.5, 0.6}; // Cone radii

    const Int_t nEBins  = 6;//16;
    //const double JetEBorders[nEBins] = {0.0, 50.0, 100.0, 200.0, 300.0, 400.0, 500.0, 600.0, 700.0, 800.0, 900.0, 1000.0, 1250.0, 1500.0, 1750.0, 2000.0}; 
    const double JetEBorders[nEBins] = {0.0, 400.0, 800.0, 1200.0, 1600.0, 2000.0}; 

    //TFile *jetFile = TFile::Open(Form("/home/lmh/alice/jetAnalysis/NewPlotting/Fractions/Results/Ofractionsptmin2__20230127_PythiaMBTrig-2_Merged_OutputR%d.root", int(Rvals[Rvalue] * 10)));
    TFile *jetFile = TFile::Open(Form("/home/lmh/alice/jetAnalysis/NewPlotting/Fractions/Results/OneMergedR%d.root", int(Rvals[Rvalue] * 10)));
    //TFile *fout = new TFile(Form("/home/lmh/alice/jetAnalysis/NewPlotting/Fractions/Results/FractionsAllMerged_OutputR%d.root", int(Rvals[Rvalue] * 10)), "RECREATE");

    TH1F *hDetHCALE[nEBins]; //HCALfracdet
    TH1F *hDetECALE[nEBins]; //ECALfracdet

    TH1F *hPartNeutralE[nEBins]; //neutralfrac
    TH1F *hPartChargedE[nEBins]; //chargedfrac
    
    TH1F *hPartPhotonElectronE[nEBins]; //photon&electronfrac

    //E loop
    for (int iE = 0; iE < nEBins-1; ++iE) 
    {

        hDetHCALE[iE] = (TH1F*)((TH1F*)jetFile->Get(Form("hDetHCALE_%d", iE)))->Clone();
        hDetECALE[iE] = (TH1F*)((TH1F*)jetFile->Get(Form("hDetECALE_%d", iE)))->Clone();
        hPartNeutralE[iE] = (TH1F*)((TH1F*)jetFile->Get(Form("hPartNeutralE_%d", iE)))->Clone();
        hPartChargedE[iE] = (TH1F*)((TH1F*)jetFile->Get(Form("hPartChargedE_%d", iE)))->Clone();
        hPartPhotonElectronE[iE] = (TH1F*)((TH1F*)jetFile->Get(Form("hPartPhotonElectronE_%d", iE)))->Clone();
    }


    TCanvas *c1 = new TCanvas("c1", "c1: the detector HCAL jet energy fractions", 1200, 1200);
    c1->Divide(2, 3);

    for (Int_t ibin = 1; ibin < nEBins; ibin++)
    {
        c1->cd(ibin);
        if(hDetHCALE[ibin-1]->GetEntries()==0) hDetHCALE[ibin-1]->Fill(-1);
        hDetHCALE[ibin-1]->Rebin();
        hDetHCALE[ibin-1]->Scale(1./hDetHCALE[ibin-1]->Integral());
        hDetHCALE[ibin-1]->SetTitle(Form("Jet energy fraction of HCAL energy, detector level jets, for %4.0f < E < %4.0f GeV", JetEBorders[ibin - 1], JetEBorders[ibin]));
        //hDetHCALE[ibin-1]->SetLineColor(kBlack);
        hDetHCALE[ibin-1]->Draw();

        //eresp_Spaghetti->GetYaxis()->SetRangeUser(0, std::max(eresp_Spaghetti->GetMaximum(), eresp_Sandwich->GetMaximum()));
        // eresp_Spaghetti->GetYaxis()->SetRangeUser(1e-5, std::max(eresp_Spaghetti->GetMaximum(), eresp_Sandwich->GetMaximum()));
    }
    c1->SaveAs(Form("MBTrighDetHCALE_%d.pdf", int(Rvals[Rvalue] * 10)));


    TCanvas *c2 = new TCanvas("c2", "c2: the detector ECAL jet energy fractions", 1200, 1200);
    c2->Divide(2, 3);

    for (Int_t ibin = 1; ibin < nEBins; ibin++)
    {
        c2->cd(ibin);
        if(hDetECALE[ibin-1]->GetEntries()==0) hDetECALE[ibin-1]->Fill(-1);
        hDetECALE[ibin-1]->Rebin();
        hDetECALE[ibin-1]->Scale(1./hDetECALE[ibin-1]->Integral());
        hDetECALE[ibin-1]->SetTitle(Form("Jet energy fraction of ËCAL energy, detector level jets, for %4.0f < E < %4.0f GeV", JetEBorders[ibin - 1], JetEBorders[ibin]));
        hDetECALE[ibin-1]->Draw();
    }
    c2->SaveAs(Form("MBTrighDetECALE_%d.pdf", int(Rvals[Rvalue] * 10)));



    TCanvas *c3 = new TCanvas("c3", "c3: the particle level jet neutral energy fractions", 1200, 1200);
    c3->Divide(2, 3);

    for (Int_t ibin = 1; ibin < nEBins; ibin++)
    {
        c3->cd(ibin);
        if(hPartNeutralE[ibin-1]->GetEntries()==0) hPartNeutralE[ibin-1]->Fill(-1);
        hPartNeutralE[ibin-1]->Rebin();
        hPartNeutralE[ibin-1]->Scale(1./hPartNeutralE[ibin-1]->Integral());
        hPartNeutralE[ibin-1]->SetTitle(Form("Jet energy fraction of neutral energy, particle level jets, for %4.0f < E < %4.0f GeV", JetEBorders[ibin - 1], JetEBorders[ibin]));
        hPartNeutralE[ibin-1]->Draw();
    }
    c3->SaveAs(Form("MBTrighPartNeutralE_%d.pdf", int(Rvals[Rvalue] * 10)));


    TCanvas *c4 = new TCanvas("c4", "c4: the particle level jet charged energy fractions", 1200, 1200);
    c4->Divide(2, 3);

    for (Int_t ibin = 1; ibin < nEBins; ibin++)
    {
        c4->cd(ibin);
        if(hPartChargedE[ibin-1]->GetEntries()==0) hPartChargedE[ibin-1]->Fill(-1);
        hPartChargedE[ibin-1]->Rebin();
        hPartChargedE[ibin-1]->Scale(1./hPartChargedE[ibin-1]->Integral());
        hPartChargedE[ibin-1]->SetTitle(Form("Jet energy fraction of charged energy, particle level jets, for %4.0f < E < %4.0f GeV", JetEBorders[ibin - 1], JetEBorders[ibin]));
        hPartChargedE[ibin-1]->Draw();
    }
    c4->SaveAs(Form("MBTrighPartChargedE_%d.pdf", int(Rvals[Rvalue] * 10)));


    TCanvas *c5 = new TCanvas("c5", "c5: the particle level jet photon and electron energy fractions", 1200, 1200);
    c5->Divide(2, 3);

    for (Int_t ibin = 1; ibin < nEBins; ibin++)
    {
        c5->cd(ibin);
        if(hPartPhotonElectronE[ibin-1]->GetEntries()==0) hPartPhotonElectronE[ibin-1]->Fill(-1);
        hPartPhotonElectronE[ibin-1]->Rebin();
        hPartPhotonElectronE[ibin-1]->Scale(1./hPartPhotonElectronE[ibin-1]->Integral());
        hPartPhotonElectronE[ibin-1]->SetTitle(Form("Jet energy fraction of photon and electron energy, particle level jets, for %4.0f < E < %4.0f GeV", JetEBorders[ibin - 1], JetEBorders[ibin]));
        hPartPhotonElectronE[ibin-1]->Draw();
    }
    c5->SaveAs(Form("MBTrighPartPhotonElectronE_%d.pdf", int(Rvals[Rvalue] * 10)));


//gStyle->SetOptStat(1111111);
// Set stat options
gStyle->SetStatY(0.9);                
// Set y-position (fraction of pad size)
gStyle->SetStatX(0.9);                
// Set x-position (fraction of pad size)
gStyle->SetStatW(0.2);                
// Set width of stat-box (fraction of pad size)
//gStyle->SetStatH(0.2);                
// Set height of stat-box (fraction of pad size)


    TCanvas *ca = new TCanvas("ca", "ca: Detector HCAL fraction", 1800, 1200);
    TH2F *hRespMatrix_DetHCALE = (TH2F*)((TH2F*)jetFile->Get("hRespMatrix_DetHCALE"))->Clone();
    hRespMatrix_DetHCALE->Rebin2D(6);
    hRespMatrix_DetHCALE->RebinY(4);
    hRespMatrix_DetHCALE->Scale(1./hRespMatrix_DetHCALE->Integral());
    hRespMatrix_DetHCALE->Draw("colz");
    ca->SetLogz();
    hRespMatrix_DetHCALE->GetXaxis()->SetRangeUser(0,1000);
    ca->SaveAs(Form("hRespMatrix_DetHCALE%d.pdf", int(Rvals[Rvalue] * 10)));

    TCanvas *cb = new TCanvas("cb", "cb: Detector ECAL fraction", 1800, 1200);
    TH2F *hRespMatrix_DetECALE =(TH2F*)((TH2F*)jetFile->Get("hRespMatrix_DetECALE"))->Clone();
    hRespMatrix_DetECALE->Rebin2D(6);
    hRespMatrix_DetECALE->RebinY(4);
    hRespMatrix_DetECALE->Scale(1./hRespMatrix_DetECALE->Integral());
    hRespMatrix_DetECALE->Draw("colz");
    cb->SetLogz();
    hRespMatrix_DetECALE->GetXaxis()->SetRangeUser(0,1000);
    cb->SaveAs(Form("hRespMatrix_DetECALE%d.pdf", int(Rvals[Rvalue] * 10)));


    TCanvas *cc = new TCanvas("cc", "cc: Particle neutral fraction", 1800, 1200);
    TH2F *hRespMatrix_PartNeutralE =  (TH2F*)((TH2F*)jetFile->Get("hRespMatrix_PartNeutralE"))->Clone();
    hRespMatrix_PartNeutralE->Rebin2D(6);
    hRespMatrix_PartNeutralE->RebinY(4);
    hRespMatrix_PartNeutralE->Scale(1./hRespMatrix_PartNeutralE->Integral());
    hRespMatrix_PartNeutralE->Draw("colz");
    cc->SetLogz();
    hRespMatrix_PartNeutralE->GetXaxis()->SetRangeUser(0,1000);
    cc->SaveAs(Form("hRespMatrix_PartNeutralE%d.pdf", int(Rvals[Rvalue] * 10)));



    TCanvas *cd = new TCanvas("cd", "cd:Particle charged fraction Particle charged fraction", 1800, 1200);
    TH2F *hRespMatrix_PartChargedE = (TH2F*)((TH2F*)jetFile->Get("hRespMatrix_PartChargedE"))->Clone();
    hRespMatrix_PartChargedE->Rebin2D(6);
    hRespMatrix_PartChargedE->RebinY(4);
    hRespMatrix_PartChargedE->Scale(1./hRespMatrix_PartChargedE->Integral());
    hRespMatrix_PartChargedE->Draw("colz");
    cd->SetLogz();
    hRespMatrix_PartChargedE->GetXaxis()->SetRangeUser(0,1000);
    cd->SaveAs(Form("hRespMatrix_PartChargedE%d.pdf", int(Rvals[Rvalue] * 10)));



    TCanvas *ce = new TCanvas("ce", "ce: Particle photon and electron fraction", 1800, 1200);
    TH2F *hRespMatrix_PartPhotonElectronE = (TH2F*)((TH2F*)jetFile->Get("hRespMatrix_PartPhotonElectronE"))->Clone();
    hRespMatrix_PartPhotonElectronE->Rebin2D(6);
    hRespMatrix_PartPhotonElectronE->RebinY(4);
    hRespMatrix_PartPhotonElectronE->Scale(1./hRespMatrix_PartPhotonElectronE->Integral());
    hRespMatrix_PartPhotonElectronE->Draw("colz");
    ce->SetLogz();
    hRespMatrix_PartPhotonElectronE->GetXaxis()->SetRangeUser(0,1000);
    ce->SaveAs(Form("hRespMatrix_PartPhotonElectronE%d.pdf", int(Rvals[Rvalue] * 10)));





}