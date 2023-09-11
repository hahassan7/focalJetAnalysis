/* Macro for plotting some histograms from jet trees
    2.21064
    0.0669805
    0.00182628
    0.000139462
    0.0000225822
*/
#include "iostream"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "TSystem.h"

double GetMedian(TH1D *h);

void JetPlottingMergedEnergyJES(int Rvalue = 2)
{
    const int nR = 3;
    const float Rvals[nR] = {0.2, 0.4, 0.6}; // Cone radii
    int constMin = 0;                        // min number of constituents in matched jet

    const int nNorm = 8;
    const float normalizations[nNorm] = {0.0611191, 0.00717001, 0.000558759, 0.000107936, 4.32163e-05, 9.57109e-06, 1.24606e-06, 6.01382e-08}; //{0.0610658,0.00716477,0.000557627,0.000107816,4.31694e-05,9.62255e-06,1.24904e-06,5.99517e-08};

    int pTHardRange[] = {5, 10, 20, 30, 40, 60, 100, 200, 1000};

    const int nCol = 10;
    const int gcolors[nCol] = {1, 2, 6, 4, 7, 1, 2, 4, 6, 7};
    const int gmarkers[nCol] = {4, 8, 25, 21, 8, 21, 25, 4, 8, 21};

    const float etaMin = 3.4;
    const float etaMax = 5.5; // changed in June23, 5.5

    // Set the pT and energy bins for the projections, JES and JER.
    const int nEBins = 15;
    // const double JetPtBorders[nPtBins] = {5.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 80.0, 100.0, 400.0}; // nPtBins = 10
    const double JetEBorders[nEBins] = {100.0, 200.0, 300.0, 400.0, 500.0, 600.0, 700.0, 800.0, 900.0, 1000.0, 1250.0, 1500.0, 1750.0, 2000.0, 3000.0}; // const int nEBins  = 15;

    for (int iNorm = 0; iNorm < nNorm; iNorm++)
    {
        TFile *jetFile;
        jetFile = TFile::Open(Form("/Users/hadi/focalJetAnalysis/jetOutput/20230728_pythia8_JetJet_%d-%dGeV/Merged_0Mass.root", pTHardRange[iNorm], pTHardRange[iNorm + 1]));

        TFile *fout;
        fout = new TFile(Form("Data20230816/JES/En20230816_pythia8_JetJet_%d-%dGeV_Merged_TotalOutputR%d_0Mass.root", pTHardRange[iNorm], pTHardRange[iNorm + 1], int(Rvals[Rvalue] * 10)), "RECREATE");

        TTree *jetTree = (TTree *)jetFile->Get("jetTree");
        TTree *TruthjetTree = (TTree *)jetFile->Get("truthjetTree");
        // TTree *PDGTree = (TTree *)jetFile->Get("inPDGTree");

        double mediansE[nEBins - 1];
        double meansE[nEBins - 1];
        double SDE[nEBins - 1];

        TH1D *hjetRatioE[nEBins - 1];
        TH1D *hMedianE;
        TH1D *hMeanE;
        TH1D *hSDE;

        const int NEFbins = 2;

        double mediansENEF[nEBins - 1][NEFbins];
        double meansENEF[nEBins - 1][NEFbins];
        double SDENEF[nEBins - 1][NEFbins];

        TH1D *hjetRatioENEF[nEBins - 1][NEFbins];
        TH1D *hMedianENEF[NEFbins];
        TH1D *hMeanENEF[NEFbins];
        TH1D *hSDENEF[NEFbins];

        hMedianE = new TH1D("hMedianE", "Mean and median of #DeltaE distribution;E^{part} (GeV);Mean or median", nEBins - 1, JetEBorders);
        hMeanE = new TH1D("hMeanE", "Mean and median of #DeltaE distribution;E^{part} (GeV);Mean or median", nEBins - 1, JetEBorders);
        hSDE = new TH1D("hSDE", "Standard deviation of #DeltaE distribution;E^{part} (GeV);Standard deviation", nEBins - 1, JetEBorders);

        for (int nef = 0; nef < NEFbins; nef++)
        {
            hMedianENEF[nef] = new TH1D(Form("hMedianENEF%d", nef), Form("Mean and median of #DeltaE distribution. NEF bin%d", nef), nEBins - 1, JetEBorders);
            hMeanENEF[nef] = new TH1D(Form("hMeanENEF%d", nef), Form("Mean and median of #DeltaE distribution. NEF bin%d", nef), nEBins - 1, JetEBorders);
            hSDENEF[nef] = new TH1D(Form("hSDENEF%d", nef), Form("Standard deviation of #DeltaE distribution. NEF bin%d", nef), nEBins - 1, JetEBorders);
        }

        // energy loop
        for (int ipt = 0; ipt < nEBins - 1; ++ipt)
        {
            hjetRatioE[ipt] = new TH1D(Form("hjetRatioE_%d", ipt), Form("Jet-by-jet #Delta E distribution, E: %d - %d GeV;(E^{det}-E^{part})/E^{part};probability", int(JetEBorders[ipt]), int(JetEBorders[ipt + 1])), 50, -1.0, 1.0);

            // Fill deltaE histos, no eta cut
            jetTree->Draw(Form("(jetE-jetE_match)/jetE_match>>hjetRatioE_%d", ipt), Form("jetR==%d && jetE_match>=%d  && jetE_match<%d  && jetEta>3.4+%f && jetEta<5.5-%f && jetEta_match>3.4+%f && jetEta_match<5.5-%f && jet_distmatch<%f && jetParts_match > %d", int(Rvals[Rvalue] * 10), int(JetEBorders[ipt]), int(JetEBorders[ipt + 1]), Rvals[Rvalue], Rvals[Rvalue], Rvals[Rvalue], Rvals[Rvalue], Rvals[Rvalue] * 0.6, constMin), "goff");
            // if(hjetRatioE[ipt]->GetEntries()==0) hjetRatioE[ipt]->Fill(-2);
            // hjetRatioE[ipt]->Scale(normalizations[iNorm] / 500.0);

            mediansE[ipt] = GetMedian(hjetRatioE[ipt]);
            meansE[ipt] = hjetRatioE[ipt]->GetMean();
            SDE[ipt] = hjetRatioE[ipt]->GetStdDev();

            hMedianE->Fill(JetEBorders[ipt + 1] - ((JetEBorders[ipt + 1] - JetEBorders[ipt]) / 2), mediansE[ipt]);
            hMedianE->SetBinError(ipt + 1, hjetRatioE[ipt]->GetMeanError());
            hMeanE->Fill(JetEBorders[ipt + 1] - ((JetEBorders[ipt + 1] - JetEBorders[ipt]) / 2), meansE[ipt]);
            hMeanE->SetBinError(ipt + 1, hjetRatioE[ipt]->GetMeanError());
            hSDE->Fill(JetEBorders[ipt + 1] - ((JetEBorders[ipt + 1] - JetEBorders[ipt]) / 2), SDE[ipt]);
            hSDE->SetBinError(ipt + 1, hjetRatioE[ipt]->GetStdDevError());

            std::cout << "In jet E range " << JetEBorders[ipt] << " - " << JetEBorders[ipt + 1] << "GeV, number of jets accepted = " << hjetRatioE[ipt]->GetEntries() << std::endl;

            hjetRatioENEF[ipt][0] = new TH1D(Form("hjetRatioENEF_%d_0", ipt), Form("Jet-by-jet #Delta E distribution, NEF < 1/3, E: %d - %d GeV;(E^{det}-E^{part})/E^{part};probability", int(JetEBorders[ipt]), int(JetEBorders[ipt + 1])), 50, -1.0, 1.0);
            hjetRatioENEF[ipt][1] = new TH1D(Form("hjetRatioENEF_%d_1", ipt), Form("Jet-by-jet #Delta E distribution, NEF > 2/3, E: %d - %d GeV;(E^{det}-E^{part})/E^{part};probability", int(JetEBorders[ipt]), int(JetEBorders[ipt + 1])), 50, -1.0, 1.0);

            // Fill deltaE histos, no eta cut

            jetTree->Draw(Form("(jetE-jetE_match)/jetE_match>>hjetRatioENEF_%d_0", ipt), Form("matchedjetECALEnergyFrac <= 0.23 && jetR==%d && jetE_match>=%d  && jetE_match<%d  && jetEta>3.4+%f && jetEta<5.5-%f && jetEta_match>3.4+%f && jetEta_match<5.5-%f && jet_distmatch<%f && jetParts_match > %d", int(Rvals[Rvalue] * 10), int(JetEBorders[ipt]), int(JetEBorders[ipt + 1]), Rvals[Rvalue], Rvals[Rvalue], Rvals[Rvalue], Rvals[Rvalue], Rvals[Rvalue] * 0.6, constMin), "goff");
            jetTree->Draw(Form("(jetE-jetE_match)/jetE_match>>hjetRatioENEF_%d_1", ipt), Form("matchedjetECALEnergyFrac >= 0.63 && jetR==%d && jetE_match>=%d  && jetE_match<%d  && jetEta>3.4+%f && jetEta<5.5-%f && jetEta_match>3.4+%f && jetEta_match<5.5-%f && jet_distmatch<%f && jetParts_match > %d", int(Rvals[Rvalue] * 10), int(JetEBorders[ipt]), int(JetEBorders[ipt + 1]), Rvals[Rvalue], Rvals[Rvalue], Rvals[Rvalue], Rvals[Rvalue], Rvals[Rvalue] * 0.6, constMin), "goff");

            // jetTree->Draw(Form("(jetE-jetE_match)/jetE_match>>hjetRatioENEF_%d_0", ipt), Form("jetECALEnergyFrac <= 0.1 && jetR==%d && jetE_match>=%d  && jetE_match<%d  && jetEta>3.4+%f && jetEta<5.5-%f && jetEta_match>3.4+%f && jetEta_match<5.5-%f && jet_distmatch<%f && jetParts_match > %d", int(Rvals[Rvalue] * 10), int(JetEBorders[ipt]), int(JetEBorders[ipt + 1]), Rvals[Rvalue], Rvals[Rvalue], Rvals[Rvalue], Rvals[Rvalue], Rvals[Rvalue] * 0.6, constMin), "goff");
            // jetTree->Draw(Form("(jetE-jetE_match)/jetE_match>>hjetRatioENEF_%d_1", ipt), Form("jetECALEnergyFrac >= 0.573 && jetR==%d && jetE_match>=%d  && jetE_match<%d  && jetEta>3.4+%f && jetEta<5.5-%f && jetEta_match>3.4+%f && jetEta_match<5.5-%f && jet_distmatch<%f && jetParts_match > %d", int(Rvals[Rvalue] * 10), int(JetEBorders[ipt]), int(JetEBorders[ipt + 1]), Rvals[Rvalue], Rvals[Rvalue], Rvals[Rvalue], Rvals[Rvalue], Rvals[Rvalue] * 0.6, constMin), "goff");
            // jetTree->Draw(Form("(jetE-jetE_match)/jetE_match>>hjetRatioENEF_%d_0", ipt), Form("matchedjetECALEnergyFrac <= 0.333333 && jetECALEnergyFrac <= 0.333333 && jetR==%d && jetE_match>=%d  && jetE_match<%d  && jetEta>3.4+%f && jetEta<5.5-%f && jetEta_match>3.4+%f && jetEta_match<5.5-%f && jet_distmatch<%f && jetParts_match > %d", int(Rvals[Rvalue] * 10), int(JetEBorders[ipt]), int(JetEBorders[ipt + 1]), Rvals[Rvalue], Rvals[Rvalue], Rvals[Rvalue], Rvals[Rvalue], Rvals[Rvalue] * 0.6, constMin), "goff");
            // jetTree->Draw(Form("(jetE-jetE_match)/jetE_match>>hjetRatioENEF_%d_1", ipt), Form("matchedjetECALEnergyFrac >= 0.666666 && jetECALEnergyFrac >= 0.666666 && jetR==%d && jetE_match>=%d  && jetE_match<%d  && jetEta>3.4+%f && jetEta<5.5-%f && jetEta_match>3.4+%f && jetEta_match<5.5-%f && jet_distmatch<%f && jetParts_match > %d", int(Rvals[Rvalue] * 10), int(JetEBorders[ipt]), int(JetEBorders[ipt + 1]), Rvals[Rvalue], Rvals[Rvalue], Rvals[Rvalue], Rvals[Rvalue], Rvals[Rvalue] * 0.6, constMin), "goff");

            // jetTree->Draw(Form("(jetE-jetE_match)/jetE_match>>hjetRatioENEF_%d_0", ipt), Form("jetECALEnergyFrac <= 0.333333 && jetR==%d && jetE_match>=%d  && jetE_match<%d  && jetEta>3.4+%f && jetEta<5.5-%f && jetEta_match>3.4+%f && jetEta_match<5.5-%f && jet_distmatch<%f && jetParts_match > %d", int(Rvals[Rvalue] * 10), int(JetEBorders[ipt]), int(JetEBorders[ipt + 1]), Rvals[Rvalue], Rvals[Rvalue], Rvals[Rvalue], Rvals[Rvalue], Rvals[Rvalue] * 0.6, constMin), "goff");
            // jetTree->Draw(Form("(jetE-jetE_match)/jetE_match>>hjetRatioENEF_%d_1", ipt), Form("jetECALEnergyFrac >= 0.666666 && jetR==%d && jetE_match>=%d  && jetE_match<%d  && jetEta>3.4+%f && jetEta<5.5-%f && jetEta_match>3.4+%f && jetEta_match<5.5-%f && jet_distmatch<%f && jetParts_match > %d", int(Rvals[Rvalue] * 10), int(JetEBorders[ipt]), int(JetEBorders[ipt + 1]), Rvals[Rvalue], Rvals[Rvalue], Rvals[Rvalue], Rvals[Rvalue], Rvals[Rvalue] * 0.6, constMin), "goff");

            for (int nef = 0; nef < NEFbins; nef++)
            {
                // hjetRatioENEF[ipt][nef]->Scale(normalizations[iNorm] / 500.0);

                mediansENEF[ipt][nef] = GetMedian(hjetRatioENEF[ipt][nef]);
                meansENEF[ipt][nef] = hjetRatioENEF[ipt][nef]->GetMean();
                SDENEF[ipt][nef] = hjetRatioENEF[ipt][nef]->GetStdDev();

                hMedianENEF[nef]->Fill(JetEBorders[ipt + 1] - ((JetEBorders[ipt + 1] - JetEBorders[ipt]) / 2), mediansENEF[ipt][nef]);
                hMedianENEF[nef]->SetBinError(ipt + 1, hjetRatioENEF[ipt][nef]->GetMeanError());
                hMeanENEF[nef]->Fill(JetEBorders[ipt + 1] - ((JetEBorders[ipt + 1] - JetEBorders[ipt]) / 2), meansENEF[ipt][nef]);
                hMeanENEF[nef]->SetBinError(ipt + 1, hjetRatioENEF[ipt][nef]->GetMeanError());
                hSDENEF[nef]->Fill(JetEBorders[ipt + 1] - ((JetEBorders[ipt + 1] - JetEBorders[ipt]) / 2), SDENEF[ipt][nef]);
                hSDENEF[nef]->SetBinError(ipt + 1, hjetRatioENEF[ipt][nef]->GetStdDevError());

                std::cout << "In jet E range " << JetEBorders[ipt] << " - " << JetEBorders[ipt + 1] << "GeV, number of jets accepted = " << hjetRatioENEF[ipt][nef]->GetEntries() << std::endl;
            }
        }

        fout->Write();
        fout->Close();
    }

    gSystem->Exec(Form("hadd Data20230816/JES/EnMergedR%d_0Mass_20PercPart.root Data20230816/JES/En20230816_pythia8_JetJet_*_0Mass.root", int(Rvals[Rvalue] * 10)));
}

double GetMedian(TH1D *h)
{
    Double_t x, q;
    q = 0.5;              // 0.5 for "median"
    h->ComputeIntegral(); // just a precaution
    h->GetQuantiles(1, &x, &q);
    std::cout << "median = " << x << std::endl;
    return x;
}
