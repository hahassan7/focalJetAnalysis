#include <TH1F.h>
#include <TCanvas.h>
#include <TFile.h>

void NormalizeToProb(TH1F* hHisto);

void ChangeTitleOfCanvas(TCanvas* cCanvas, TString sTitle);

//Draw fraction histos from v1 code (no 2D fraction plot)
void JetPlottingMergedFractionsEightPlot(int NormValue=8)
{   
    const double limitYconst = 0.3;
    const double limitYpT    = 0.15;
    const double limitYEta     = 0.13;
    const double limitYconstmin = 0.0;
    const double limitYpTmin    = 0.0;
    const double limitYEtamin     = 0.0;

    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
//setting NormValue = 8 if merged file in use

    if(NormValue!=8) return;
    TFile *fin;
    /*if(NormValue==0)fin = new TFile("JetJetOutput/July2023/ptmin10/fractions/20230722_pythia8_JetJet_5-10GeV_Merged_Output.root", "READ");
    if(NormValue==1)fin = new TFile("JetJetOutput/July2023/ptmin10/fractions/20230722_pythia8_JetJet_10-20GeV_Merged_Output.root", "READ");
    if(NormValue==2)fin = new TFile("JetJetOutput/July2023/ptmin10/fractions/20230722_pythia8_JetJet_20-30GeV_Merged_Output.root", "READ");
    if(NormValue==3)fin = new TFile("JetJetOutput/July2023/ptmin10/fractions/20230722_pythia8_JetJet_30-40GeV_Merged_Output.root", "READ");
    if(NormValue==4)fin = new TFile("JetJetOutput/July2023/ptmin10/fractions/20230722_pythia8_JetJet_40-60GeV_Merged_Output.root", "READ");
    if(NormValue==5)fin = new TFile("JetJetOutput/July2023/ptmin10/fractions/20230722_pythia8_JetJet_60-100GeV_Merged_Output.root", "READ");
    if(NormValue==6)fin = new TFile("JetJetOutput/July2023/ptmin10/fractions/20230722_pythia8_JetJet_100-200GeV_Merged_Output.root", "READ");
    if(NormValue==7)fin = new TFile("JetJetOutput/July2023/ptmin10/fractions/20230722_pythia8_JetJet_200-GeV_Merged_Output.root", "READ");*/
    if(NormValue==8)fin = new TFile("Data20230728/FRAC/Merged.root", "READ");

    //Add normalization to probability density function: 1/Int, "width"

    // ... Constants and binning definitions ...
    const Int_t nR = 3;
    const Float_t Rvals[nR] = {0.2, 0.4, 0.6}; // Cone radii
    const double EtaMin[2] = {3.4, 3.9}; 
    const double EtaMax[2] = {5.0, 5.5};
    const Float_t etaMin = 3.4; 
    const Float_t etaMax = 5.5;

    // Eta range binning
    const Int_t nEtaBins = 3;
    const double EtaBinBorders[nR][nEtaBins] = {{3.6, 4.5, 5.3}, {3.8, 4.5, 5.1}, {4.0, 4.5, 4.9}};
    TString etaRange[nEtaBins-1] = {"3.4+R < #eta_{jet} < 4.5", "4.5 < #eta_{jet} < 5.5-R"};

    // Energy and pT binning
    const Int_t nEBins = 6;
    const Int_t npTBins = 8;
    const double JetEBorders[nEBins] = {100.0, 400.0, 800.0, 1200.0, 1600.0, 2000.0}; 
    const double JetpTBorders[npTBins] = {5.0, 10.0, 20.0, 40.0, 60.0, 80.0, 100.0, 150.0};//, 200.0}; 

    const Int_t nConst = 5;
    const Int_t ConstVals[nConst] = {1, 2, 3, 5, 10}; // Constituent values to check



    // Histograms
    TH1F *hDetECAL[nR]; // ECALfracdet matched
    TH1F *hPartNeutralMatch[nR]; // neutralfrac matched
    TH1F *hDetECALAll[nR]; // ECALfracdet all not just matched
    TH1F *hPartNeutral[nR]; // neutralfrac all not just matched

    TH1F *hDetECALEta[nR][nEtaBins-1]; // ECALfracdet matched
    TH1F *hPartNeutralMatchEta[nR][nEtaBins-1]; // neutralfrac matched
    TH1F *hDetECALAllEta[nR][nEtaBins-1]; // ECALfracdet all not just matched
    TH1F *hPartNeutralEta[nR][nEtaBins-1]; // neutralfrac all not just matched

    TH1F *hDetECALEtaConst[nR][nEtaBins-1][nConst]; // ECALfracdet matched
    TH1F *hPartNeutralMatchEtaConst[nR][nEtaBins-1][nConst]; // neutralfrac matched
    TH1F *hDetECALAllEtaConst[nR][nEtaBins-1][nConst]; // ECALfracdet all not just matched
    TH1F *hPartNeutralEtaConst[nR][nEtaBins-1][nConst]; // neutralfrac all not just matched

    TH1F *hDetECALConst[nR][nConst]; // ECALfracdet matched
    TH1F *hPartNeutralMatchConst[nR][nConst]; // neutralfrac matched
    TH1F *hDetECALAllConst[nR][nConst]; // ECALfracdet all not just matched
    TH1F *hPartNeutralConst[nR][nConst]; // neutralfrac all not just matched


    // ... Open the input file and read histograms ...
    for (int Rvalue = 0; Rvalue < nR; Rvalue++) {
        // Read histograms
        hDetECAL[Rvalue] = (TH1F*)(fin->Get(Form("hDetECAL_R%d", Rvalue))->Clone());NormalizeToProb(hDetECAL[Rvalue]);
        hDetECAL[Rvalue]->GetYaxis()->SetRangeUser(limitYpTmin,limitYpT);
        hPartNeutralMatch[Rvalue] = (TH1F*)(fin ->Get(Form("hPartNeutralMatch_R%d", Rvalue))->Clone());NormalizeToProb(hPartNeutralMatch[Rvalue]);
        hPartNeutralMatch[Rvalue]->GetYaxis()->SetRangeUser(limitYpTmin,limitYpT);
        hDetECALAll[Rvalue] = (TH1F*)(fin ->Get(Form("hDetECALAll_R%d", Rvalue))->Clone());NormalizeToProb(hDetECALAll[Rvalue]);
        hDetECALAll[Rvalue]->GetYaxis()->SetRangeUser(limitYpTmin,limitYpT);
        hPartNeutral[Rvalue] = (TH1F*)(fin ->Get(Form("hPartNeutral_R%d", Rvalue))->Clone());NormalizeToProb(hPartNeutral[Rvalue]);
        hPartNeutral[Rvalue]->GetYaxis()->SetRangeUser(limitYpTmin,limitYpT);

        for (int ieta = 0; ieta < nEtaBins-1; ++ieta) {
            hDetECALEta[Rvalue][ieta] = (TH1F*)(fin ->Get(Form("hDetECALEta_%d_R%d", ieta, Rvalue))->Clone());NormalizeToProb(hDetECALEta[Rvalue][ieta]);
            hDetECALEta[Rvalue][ieta]->GetYaxis()->SetRangeUser(limitYEtamin,limitYEta);
            hPartNeutralMatchEta[Rvalue][ieta] = (TH1F*)(fin ->Get(Form("hPartNeutralMatchEta_%d_R%d", ieta, Rvalue))->Clone());NormalizeToProb(hPartNeutralMatchEta[Rvalue][ieta]);
            hPartNeutralMatchEta[Rvalue][ieta]->GetYaxis()->SetRangeUser(limitYEtamin,limitYEta);
            hDetECALAllEta[Rvalue][ieta] = (TH1F*)(fin ->Get(Form("hDetECALAllEta_%d_R%d", ieta, Rvalue))->Clone());NormalizeToProb(hDetECALAllEta[Rvalue][ieta]);
            hDetECALAllEta[Rvalue][ieta]->GetYaxis()->SetRangeUser(limitYEtamin,limitYEta);
            hPartNeutralEta[Rvalue][ieta] = (TH1F*)(fin ->Get(Form("hPartNeutralEta_%d_R%d", ieta, Rvalue))->Clone());NormalizeToProb(hPartNeutralEta[Rvalue][ieta]);
            hPartNeutralEta[Rvalue][ieta]->GetYaxis()->SetRangeUser(limitYEtamin,limitYEta);

            for (int iconst = 0; iconst < nConst; ++iconst) {
                hDetECALEtaConst[Rvalue][ieta][iconst] = (TH1F*)(fin ->Get(Form("hDetECALEta_%d_Const_%d_R%d", ieta, iconst, Rvalue))->Clone());NormalizeToProb(hDetECALEtaConst[Rvalue][ieta][iconst]);
            	hDetECALEtaConst[Rvalue][ieta][iconst]->GetYaxis()->SetRangeUser(limitYconstmin,limitYconst);
                hPartNeutralMatchEtaConst[Rvalue][ieta][iconst] = (TH1F*)(fin ->Get(Form("hPartNeutralMatchEta_%d_Const_%d_R%d", ieta, iconst, Rvalue))->Clone());NormalizeToProb(hPartNeutralMatchEtaConst[Rvalue][ieta][iconst]);
            	hPartNeutralMatchEtaConst[Rvalue][ieta][iconst]->GetYaxis()->SetRangeUser(limitYconstmin,limitYconst);
                hDetECALAllEtaConst[Rvalue][ieta][iconst] = (TH1F*)(fin ->Get(Form("hDetECALAllEta_%d_Const_%d_R%d", ieta, iconst, Rvalue))->Clone());NormalizeToProb(hDetECALAllEtaConst[Rvalue][ieta][iconst]);
            	hDetECALAllEtaConst[Rvalue][ieta][iconst]->GetYaxis()->SetRangeUser(limitYconstmin,limitYconst);
                hPartNeutralEtaConst[Rvalue][ieta][iconst] = (TH1F*)(fin ->Get(Form("hPartNeutralEta_%d_Const_%d_R%d", ieta, iconst, Rvalue))->Clone());NormalizeToProb(hPartNeutralEtaConst[Rvalue][ieta][iconst]);
            	hPartNeutralEtaConst[Rvalue][ieta][iconst]->GetYaxis()->SetRangeUser(limitYconstmin,limitYconst);
            }
        }
        for (int iconst = 0; iconst < nConst; ++iconst) {
            hDetECALConst[Rvalue][iconst] = (TH1F*)(fin ->Get(Form("hDetECALConst_%d_R%d", iconst, Rvalue))->Clone());NormalizeToProb(hDetECALConst[Rvalue][iconst]);
            hDetECALConst[Rvalue][iconst]->GetYaxis()->SetRangeUser(limitYconstmin,limitYconst);
            hPartNeutralMatchConst[Rvalue][iconst] = (TH1F*)(fin ->Get(Form("hPartNeutralMatchConst_%d_R%d", iconst, Rvalue))->Clone());NormalizeToProb(hPartNeutralMatchConst[Rvalue][iconst]);
            hPartNeutralMatchConst[Rvalue][iconst]->GetYaxis()->SetRangeUser(limitYconstmin,limitYconst);
            hDetECALAllConst[Rvalue][iconst] = (TH1F*)(fin ->Get(Form("hDetECALAllConst_%d_R%d", iconst, Rvalue))->Clone());NormalizeToProb(hDetECALAllConst[Rvalue][iconst]);
            hDetECALAllConst[Rvalue][iconst]->GetYaxis()->SetRangeUser(limitYconstmin,limitYconst);
            hPartNeutralConst[Rvalue][iconst] = (TH1F*)(fin ->Get(Form("hPartNeutralConst_%d_R%d", iconst, Rvalue))->Clone());NormalizeToProb(hPartNeutralConst[Rvalue][iconst]);
            hPartNeutralConst[Rvalue][iconst]->GetYaxis()->SetRangeUser(limitYconstmin,limitYconst);
        }
    }

    const int nCol = 10;
    const int gcolors[nCol] = {1, 2, 4, 6, 4, 7, 1, 2, 4, 6};
    const int gmarkers[nCol] = {4, 8, 21, 21, 8, 21, 25, 4, 8, 21};

    // Create separate canvases for each histogram array
    TCanvas *cDetECAL = new TCanvas("cDetECAL", "Matched detector level jet ECAL energy fraction distribution", 800, 600);
    TCanvas *cPartNeutralMatch = new TCanvas("cPartNeutralMatch", "Matched particle level jet neutral energy fraction distribution", 800, 600);
    TCanvas *cDetECALAll = new TCanvas("cDetECALAll", "All detector level jet ECAL energy fraction distribution", 800, 600);
    TCanvas *cPartNeutral = new TCanvas("cPartNeutral", "All particle level jet neutral energy fraction distribution", 800, 600);


    	// Add legends for each canvas
	    TLegend *legendDetECAL = new TLegend(0.7, 0.7, 0.9, 0.9);
	    TLegend *legendPartNeutralMatch = new TLegend(0.7, 0.7, 0.9, 0.9);
	    TLegend *legendDetECALAll = new TLegend(0.7, 0.7, 0.9, 0.9);
	    TLegend *legendPartNeutral = new TLegend(0.7, 0.7, 0.9, 0.9);

    // Loop over R values and draw histograms on the canvases
    for (int Rvalue = 0; Rvalue < nR; Rvalue++) {
        // Set line color and marker style based on Rvalue
        hDetECAL[Rvalue]->SetLineColor(gcolors[Rvalue]);
        hDetECAL[Rvalue]->SetMarkerStyle(gmarkers[Rvalue]);
        hPartNeutralMatch[Rvalue]->SetLineColor(gcolors[Rvalue]);
        hPartNeutralMatch[Rvalue]->SetMarkerStyle(gmarkers[Rvalue]);
        hDetECALAll[Rvalue]->SetLineColor(gcolors[Rvalue]);
        hDetECALAll[Rvalue]->SetMarkerStyle(gmarkers[Rvalue]);
        hPartNeutral[Rvalue]->SetLineColor(gcolors[Rvalue]);
        hPartNeutral[Rvalue]->SetMarkerStyle(gmarkers[Rvalue]);

        // Set x-axis titles
        hDetECAL[Rvalue]->GetXaxis()->SetTitle("E_{ECAL}/E_{jet}^{det}");
        hPartNeutralMatch[Rvalue]->GetXaxis()->SetTitle("E_{neutral}/E_{part}");
        hDetECALAll[Rvalue]->GetXaxis()->SetTitle("E_{ECAL}/E_{jet}^{det}");
        hPartNeutral[Rvalue]->GetXaxis()->SetTitle("E_{neutral}/E_{part}");

        // Draw histograms on the corresponding canvas
        cDetECAL->cd();
        if (Rvalue == 0) hDetECAL[Rvalue]->Draw();
        else hDetECAL[Rvalue]->Draw("same");
	    legendDetECAL->AddEntry(hDetECAL[Rvalue], Form("R = %0.1f", Rvals[Rvalue]), "pl");
	    legendDetECAL->Draw();

        cPartNeutralMatch->cd();
        if (Rvalue == 0) hPartNeutralMatch[Rvalue]->Draw();
        else hPartNeutralMatch[Rvalue]->Draw("same");
	    legendPartNeutralMatch->AddEntry(hPartNeutralMatch[Rvalue], Form("R = %0.1f", Rvals[Rvalue]), "pl");
	    legendPartNeutralMatch->Draw();

        cDetECALAll->cd();
        if (Rvalue == 0) hDetECALAll[Rvalue]->Draw();
        else hDetECALAll[Rvalue]->Draw("same");
	    legendDetECALAll->AddEntry(hDetECALAll[Rvalue], Form("R = %0.1f", Rvals[Rvalue]), "pl");
	    legendDetECALAll->Draw();

        cPartNeutral->cd();
        if (Rvalue == 0) hPartNeutral[Rvalue]->Draw();
        else hPartNeutral[Rvalue]->Draw("same");
	    legendPartNeutral->AddEntry(hPartNeutral[Rvalue], Form("R = %0.1f", Rvals[Rvalue]), "pl");
    	legendPartNeutral->Draw();


	    TCanvas *cDetECALEta = new TCanvas("cDetECALEta", "Matched detector level jet ECAL energy fraction distribution, eta", 800, 600);
	    TCanvas *cPartNeutralMatchEta = new TCanvas("cPartNeutralMatchEta", "Matched particle level jet neutral energy fraction distribution, eta", 800, 600);
        TCanvas *cPartNeutralEta = new TCanvas("cPartNeutralEta", "Total particle level jet neutral energy fraction distribution, eta", 800, 600);
		TLegend *legendDetECALEta = new TLegend(0.7, 0.7, 0.9, 0.9);
		TLegend *legendPartNeutralMatchEta = new TLegend(0.7, 0.7, 0.9, 0.9);
        TLegend *legendPartNeutralEta = new TLegend(0.65, 0.7, 0.85, 0.85);legendPartNeutralEta->SetBorderSize(0);
    	for (int ieta = 0; ieta < nEtaBins-1; ++ieta) {

            // Set line color and marker style based on Rvalue
            hDetECALEta[Rvalue][ieta]->SetLineColor(gcolors[ieta]);
            hDetECALEta[Rvalue][ieta]->SetMarkerStyle(gmarkers[ieta]);
            hPartNeutralMatchEta[Rvalue][ieta]->SetLineColor(gcolors[ieta]);
            hPartNeutralMatchEta[Rvalue][ieta]->SetMarkerStyle(gmarkers[ieta]);
            hPartNeutralEta[Rvalue][ieta]->SetLineColor(gcolors[ieta]);
            hPartNeutralEta[Rvalue][ieta]->SetMarkerStyle(gmarkers[ieta]);

            // Set x-axis titles
            hDetECALEta[Rvalue][ieta]->GetXaxis()->SetTitle("E_{ECAL}/E_{jet}^{det}");
            hPartNeutralMatchEta[Rvalue][ieta]->GetXaxis()->SetTitle("E_{neutral}/E_{part}");
            hPartNeutralEta[Rvalue][ieta]->GetXaxis()->SetTitle("E_{neutral}/E_{part}");

            // Draw histograms on the corresponding canvas
            cDetECALEta->cd();
            if (ieta == 0) hDetECALEta[Rvalue][ieta]->Draw();
            else hDetECALEta[Rvalue][ieta]->Draw("same");
		    legendDetECALEta->AddEntry(hDetECALEta[Rvalue][ieta], Form("%s", etaRange[ieta].Data()), "pl");
		    legendDetECALEta->Draw();

            cPartNeutralMatchEta->cd();
            if (ieta == 0) hPartNeutralMatchEta[Rvalue][ieta]->Draw();
            else hPartNeutralMatchEta[Rvalue][ieta]->Draw("same");
		    legendPartNeutralMatchEta->AddEntry(hPartNeutralMatchEta[Rvalue][ieta], Form("%s", etaRange[ieta].Data()), "pl");
		    legendPartNeutralMatchEta->Draw();


            cPartNeutralEta->cd();
            if (ieta == 0) hPartNeutralEta[Rvalue][ieta]->Draw();
            else hPartNeutralEta[Rvalue][ieta]->Draw("same");
            legendPartNeutralEta->AddEntry(hPartNeutralEta[Rvalue][ieta], Form("%s", etaRange[ieta].Data()), "pl");
            legendPartNeutralEta->Draw();


		    TCanvas *cDetECALEtaConst = new TCanvas("cDetECALEtaConst", "Matched detector level jet ECAL energy fraction distribution, eta, const", 800, 600);
		    TCanvas *cPartNeutralMatchEtaConst = new TCanvas("cPartNeutralMatchEtaConst", "Matched particle level jet neutral energy fraction distribution, eta, const", 800, 600);
			TLegend *legendDetECALEtaConst = new TLegend(0.7, 0.7, 0.9, 0.9);
			TLegend *legendPartNeutralMatchEtaConst = new TLegend(0.7, 0.7, 0.9, 0.9);
            // Loop over iconst and draw histograms on corresponding canvases
            for (int iconst = 0; iconst < nConst; ++iconst) {

                // Set line color and marker style based on Rvalue
                hDetECALEtaConst[Rvalue][ieta][iconst]->SetLineColor(gcolors[iconst]);
                hDetECALEtaConst[Rvalue][ieta][iconst]->SetMarkerStyle(gmarkers[iconst]);
                hPartNeutralMatchEtaConst[Rvalue][ieta][iconst]->SetLineColor(gcolors[iconst]);
                hPartNeutralMatchEtaConst[Rvalue][ieta][iconst]->SetMarkerStyle(gmarkers[iconst]);

                // Set x-axis titles
                hDetECALEtaConst[Rvalue][ieta][iconst]->GetXaxis()->SetTitle("E_{ECAL}/E_{jet}^{det}");
                hPartNeutralMatchEtaConst[Rvalue][ieta][iconst]->GetXaxis()->SetTitle("E_{neutral}/E_{part}");

                // Draw histograms on the corresponding canvas
	            cDetECALEtaConst->cd();
	            if (ieta == 0) hDetECALEtaConst[Rvalue][ieta][iconst]->Draw();
	            else hDetECALEtaConst[Rvalue][ieta][iconst]->Draw("same");
			    legendDetECALEtaConst->AddEntry(hDetECALEtaConst[Rvalue][ieta][iconst], Form("Constituent N = %d", ConstVals[iconst]), "pl");
			    legendDetECALEtaConst->Draw();

	            cPartNeutralMatchEtaConst->cd();
	            if (ieta == 0) hPartNeutralMatchEtaConst[Rvalue][ieta][iconst]->Draw();
	            else hPartNeutralMatchEtaConst[Rvalue][ieta][iconst]->Draw("same");
			    legendPartNeutralMatchEtaConst->AddEntry(hPartNeutralMatchEtaConst[Rvalue][ieta][iconst], Form("Constituent N = %d", ConstVals[iconst]), "pl");
			    legendPartNeutralMatchEtaConst->Draw();
            }
            ChangeTitleOfCanvas(cDetECALEtaConst, Form("Matched detector level jet ECAL energy fraction, %s, R=%0.1f", etaRange[ieta].Data(), Rvals[Rvalue]));
		    cDetECALEtaConst->SaveAs(Form("figs/FRAC/cDetECALEta_%dConst_R%d.png",ieta, Rvalue));
            ChangeTitleOfCanvas(cPartNeutralMatchEtaConst, Form("Matched particle level jet neutral energy fraction, %s, R=%0.1f", etaRange[ieta].Data(), Rvals[Rvalue]));
		    cPartNeutralMatchEtaConst->SaveAs(Form("figs/FRAC/cPartNeutralMatchEta_%dConst_R%d.png", ieta, Rvalue));

        }
        ChangeTitleOfCanvas(cDetECALEta, Form("Matched detector level jet ECAL energy fraction, R=%0.1f", Rvals[Rvalue]));
	    cDetECALEta->SaveAs(Form("figs/FRAC/cDetECALEta_R%d.png", Rvalue));
        ChangeTitleOfCanvas(cPartNeutralMatchEta, Form("Matched particle level jet neutral energy fraction, R=%0.1f", Rvals[Rvalue]));
	    cPartNeutralMatchEta->SaveAs(Form("figs/FRAC/cPartNeutralMatchEta_R%d.png", Rvalue));
        ChangeTitleOfCanvas(cPartNeutralEta, Form("Total particle level jet neutral energy fraction, R=%0.1f", Rvals[Rvalue]));
        cPartNeutralEta->SaveAs(Form("figs/FRAC/cPartNeutralEta_R%d.png", Rvalue));
/*
	    TCanvas *cDetECALConst = new TCanvas("cDetECALConst", "Matched detector level jet ECAL energy fraction distribution, const", 800, 600);
	    TCanvas *cPartNeutralMatchConst = new TCanvas("cPartNeutralMatchConst", "Matched particle level jet neutral energy fraction distribution, const", 800, 600);
		TLegend *legendDetECALConst = new TLegend(0.7, 0.7, 0.9, 0.9);
		TLegend *legendPartNeutralMatchConst = new TLegend(0.7, 0.7, 0.9, 0.9);
        // Loop over iconst and draw histograms on corresponding canvases
        for (int iconst = 0; iconst < nConst; ++iconst) {

            // Set line color and marker style based on Rvalue
            hDetECALConst[Rvalue][iconst]->SetLineColor(gcolors[iconst]);
            hDetECALConst[Rvalue][iconst]->SetMarkerStyle(gmarkers[iconst]);
            hPartNeutralMatchConst[Rvalue][iconst]->SetLineColor(gcolors[iconst]);
            hPartNeutralMatchConst[Rvalue][iconst]->SetMarkerStyle(gmarkers[iconst]);

            // Set x-axis titles
            hDetECALConst[Rvalue][iconst]->GetXaxis()->SetTitle("E_{ECAL}/E_{jet}^{det}");
            hPartNeutralMatchConst[Rvalue][iconst]->GetXaxis()->SetTitle("E_{neutral}/E_{part}");

            // Draw histograms on the corresponding canvas
            cDetECALConst->cd();
            if (iconst == 0) hDetECALConst[Rvalue][iconst]->Draw();
            else hDetECALConst[Rvalue][iconst]->Draw("same");
		    legendDetECALConst->AddEntry(hDetECALConst[Rvalue][iconst], Form("Constituent N = %d", ConstVals[iconst]), "pl");
		    legendDetECALConst->Draw();

            cPartNeutralMatchConst->cd();
            if (iconst == 0) hPartNeutralMatchConst[Rvalue][iconst]->Draw();
            else hPartNeutralMatchConst[Rvalue][iconst]->Draw("same");
		    legendPartNeutralMatchConst->AddEntry(hPartNeutralMatchConst[Rvalue][iconst], Form("Constituent N = %d", ConstVals[iconst]), "pl");
		    legendPartNeutralMatchConst->Draw();

        } 
        ChangeTitleOfCanvas(cDetECALConst, Form("Matched detector level jet ECAL energy fraction, R=%0.1f", Rvals[Rvalue]));
	    cDetECALConst->SaveAs(Form("figs/FRAC/cDetECALConst_R%d.png", Rvalue));
        ChangeTitleOfCanvas(cPartNeutralMatchConst, Form("Matched particle level jet neutral energy fraction, R=%0.1f", Rvals[Rvalue]));
	    cPartNeutralMatchConst->SaveAs(Form("figs/FRAC/cPartNeutralMatchConst_R%d.png", Rvalue));

*///skip constituent plots for now
    
    }

    // Save canvases as PNG images
    ChangeTitleOfCanvas(cDetECAL, "Matched detector level jet ECAL energy fraction");
    cDetECAL->SaveAs("figs/FRAC/cDetECAL.png");
    ChangeTitleOfCanvas(cPartNeutralMatch, "Matched particle level jet neutral energy fraction");
    cPartNeutralMatch->SaveAs("figs/FRAC/cPartNeutralMatch.png");
    ChangeTitleOfCanvas(cDetECALAll, "Total detector level jet ECAL energy fraction");
    cDetECALAll->SaveAs("figs/FRAC/cDetECALAll.png");
    ChangeTitleOfCanvas(cPartNeutral, "Total particle level jet neutral energy fraction");
    cPartNeutral->SaveAs("figs/FRAC/cPartNeutral.png");


}


void NormalizeToProb(TH1F* hHisto){
    Double_t factor = 1.;
    if (hHisto!=NULL) hHisto->Scale(factor/hHisto->Integral()); //, "width"
}

void ChangeTitleOfCanvas(TCanvas* cCanvas, TString sTitle){
    cCanvas->cd();
    TLatex *   tex = new TLatex(0.25,0.92,Form("%s", sTitle.Data()));
    tex->SetNDC();
    tex->SetTextSize(0.03);
    tex->Draw();
    //if(SettingLogY==1) cCanvas->SetLogy();
}