#include "TH2D.h"
#include "TMath.h"
#include "iostream"
#include "TVector2.h"
#include "TGraph.h"
#include "TEllipse.h"
#include "TStyle.h"
#include "TExec.h"
#include "TFile.h"
#include "TCanvas.h"

TH2D *GetXYFromPhiEta(double eta0, double phi0, double radius, double distace, std::string name);

void DrawCircle()
{
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);

    bool writeCircles = false;

    double eta(0.), phi(0.), distance(700);

    double phi0 = 0.;
    double eta0 = 4.9;
    double eta1 = 4.5;
    double eta2 = 4.0;
    double radius = 0.6;

    double xPos0 = ((2 * distance) / (TMath::Exp(eta0))) * TMath::Cos(phi0); // Eta=4.9, Phi 0, 10.425 Col=23 Col+3.42 Col-1.88  // Eta=4.5, Phi=0, 15.55 Col=25 Col+5.11 Col-2.8 // Eta=4.0, Phi=0, 25.64 Col=29 Col+8.4 Col-4.6
    double yPos0 = ((2 * distance) / (TMath::Exp(eta0))) * TMath::Sin(phi0); // Eta=4.9, Phi 0, 0. Row=20 Row+-2.8               // Eta=4.5, Phi=0, 0. Row=20 Row+-4.1            // Eta=4.0, Phi=0, 0.   Row=20 Row-+6.8

    double xPos1 = ((2 * distance) / (TMath::Exp(eta1))) * TMath::Cos(phi0);
    double yPos1 = ((2 * distance) / (TMath::Exp(eta1))) * TMath::Sin(phi0);
    
    double xPos2 = ((2 * distance) / (TMath::Exp(eta2))) * TMath::Cos(phi0);
    double yPos2 = ((2 * distance) / (TMath::Exp(eta2))) * TMath::Sin(phi0);

    std::cout << "orig xPos: " << xPos0 << " orig yPos: " << yPos0 << std::endl;

    TEllipse *eLong1 = new TEllipse(xPos0, yPos0 + 1.25, 19.53 - xPos0, 7.17, -90, 90);
    TEllipse *eShort1 = new TEllipse(xPos0, yPos0 + 1.25, xPos0 - 5.88, 7.17, 90, 270);
    TEllipse *eLong2 = new TEllipse(xPos1, yPos1 + 1.25, 29.2 - xPos1, 10.7, -90, 90);
    TEllipse *eShort2 = new TEllipse(xPos1, yPos1 + 1.25, xPos1 - 8.77, 10.7, 90, 270);
    TEllipse *eLong3 = new TEllipse(xPos2, yPos2 + 1.25, 48.05 - xPos2, 17.6, -90, 90);
    TEllipse *eShort3 = new TEllipse(xPos2, yPos2 + 1.25, xPos2 - 14.5, 17.6, 90, 270);
    eLong1->SetFillStyle(0);
    eLong1->SetLineColor(kBlack);
    eLong1->SetLineWidth(2);
    eShort1->SetFillStyle(0);
    eShort1->SetLineColor(kBlack);
    eShort1->SetLineWidth(2);
    eLong2->SetFillStyle(0);
    eLong2->SetLineColor(kRed);
    eLong2->SetLineWidth(2);
    eShort2->SetFillStyle(0);
    eShort2->SetLineColor(kRed);
    eShort2->SetLineWidth(2);
    eLong3->SetFillStyle(0);
    eLong3->SetLineColor(kGreen);
    eLong3->SetLineWidth(2);
    eShort3->SetFillStyle(0);
    eShort3->SetLineColor(kGreen);
    eShort3->SetLineWidth(2);

    // gStyle->SetPalette(kInvertedDarkBodyRadiator);
    // myhist->SetFillColorAlpha(kBlack, 1.);

    TH2D *myhist_49 = GetXYFromPhiEta(eta0, phi0, radius, distance, "JetCone_eta49");
    TH2D *myhist_45 = GetXYFromPhiEta(eta1, phi0, radius, distance, "JetCone_eta45");
    TH2D *myhist_40 = GetXYFromPhiEta(eta2, phi0, radius, distance, "JetCone_eta40");

    TExec *ex1 = new TExec("ex1", "gStyle->SetPalette(0);");
    TExec *ex2 = new TExec("ex2", "gStyle->SetPalette(1);");
    TExec *ex3 = new TExec("ex3", "gStyle->SetPalette(2);");

    TCanvas* c1 = new TCanvas("c1", "Jet Cones", 700, 700);
    c1->cd();

    myhist_49->Draw("colz");
    ex1->Draw();
    myhist_45->Draw("col,same");
    ex2->Draw();
    myhist_40->Draw("col,same");
    ex3->Draw();

    eLong1->Draw();
    eShort1->Draw();

    eLong2->Draw();
    eShort2->Draw();

    eLong3->Draw();
    eShort3->Draw();

    if (writeCircles)
    {
        TFile *outputFile = new TFile("JetConesEta.root", "RECREATE");
        outputFile->cd();
        myhist_49->Write();
        myhist_45->Write();
        myhist_40->Write();
    }
}

TH2D *GetXYFromPhiEta(double eta0, double phi0, double radius, double distance, std::string name)
{
    TH2D *myhist = new TH2D(Form("%s", name.c_str()), "Circle from eta-phi to x-y;x (cm);y (cm)", 240, -47.5, 47.5, 240, -50, 50);

    double xPos(0.), yPos(0.);
    double phi(0.), eta(0.);

    int number1 = 0;
    int number2 = 24;
    for (phi = 0; phi < TMath::TwoPi(); phi += 0.000001)
    {
        eta = TMath::Sqrt(radius * radius - (phi0 - TVector2::Phi_mpi_pi(phi)) * (phi0 - TVector2::Phi_mpi_pi(phi))) + eta0;

        xPos = ((2 * distance) / (TMath::Exp(eta))) * TMath::Cos(phi);
        yPos = ((2 * distance) / (TMath::Exp(eta))) * TMath::Sin(phi);

        if (!std::isnan(eta))
        {
            // std::cout << "xPos: " << xPos << " yPos: " << yPos << std::endl;
            number1++;
        }

        myhist->Fill(xPos, yPos + 1.25);

        eta = -1 * TMath::Sqrt(radius * radius - (phi0 - TVector2::Phi_mpi_pi(phi)) * (phi0 - TVector2::Phi_mpi_pi(phi))) + eta0;

        xPos = ((2 * distance) / (TMath::Exp(eta))) * TMath::Cos(phi);
        yPos = ((2 * distance) / (TMath::Exp(eta))) * TMath::Sin(phi);
        // graph->SetPoint(number + 31, xPos, yPos);

        if (!std::isnan(eta))
        {
            // std::cout << "xPos: " << xPos << " yPos: " << yPos << std::endl;
            number2--;
        }

        myhist->Fill(xPos, yPos + 1.25);
    }

    for (int ibinX = 1; ibinX <= myhist->GetNbinsX(); ibinX++)
    {
        for (int ibinY = 1; ibinY <= myhist->GetNbinsY(); ibinY++)
        {
            if (myhist->GetBinContent(ibinX, ibinY) > 0)
            {
                myhist->SetBinContent(ibinX, ibinY, 1.1e-2);
            }
        }
    }

    return myhist;
}