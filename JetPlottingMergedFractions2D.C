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

void JetPlottingMergedFractions2D()
{
    const int nNorm = 8;
    const float normalizations[nNorm] = {0.0611191, 0.00717001, 0.000558759, 0.000107936, 4.32163e-05, 9.57109e-06, 1.24606e-06, 6.01382e-08}; //{0.0610658,0.00716477,0.000557627,0.000107816,4.31694e-05,9.62255e-06,1.24904e-06,5.99517e-08};
    int pTHardRange[] = {5, 10, 20, 30, 40, 60, 100, 200, 1000};

    const float nFolders[nNorm] = {500.0, 500.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0};
    // const float nFolders[nNorm] = {500.0, 500.0, 500.0, 500.0, 500.0, 500.0, 500.0, 500.0};

    const int nR = 3;
    const float Rvals[nR] = {0.2, 0.4, 0.6}; // Cone radii
    // const double EtaMin[2] = {3.4,3.9};
    // const double EtaMax[2] = {5.0,5.5};

    const float etaMin = 3.4;
    const float etaMax = 5.5; // changed in June23, 3.4 to 5.5

    // add eta range binning here, which depends on the value of R.
    const int nEtaBins = 3;
    const double EtaBinBorders[nEtaBins][nR] = {{3.6, 4.5, 5.3}, {3.8, 4.5, 5.1}, {4.0, 4.5, 4.9}};
    TString etaRange[nEtaBins - 1] = {"3.4+R < #eta_{jet} < 4.5", "4.5 < #eta_{jet} < 5.5-R"};
    const int nEBins = 15;                                                                                                                              // 16;
    const int npTBins = 9;                                                                                                                              // 16;
    const double JetEBorders[nEBins] = {100.0, 200.0, 300.0, 400.0, 500.0, 600.0, 700.0, 800.0, 900.0, 1000.0, 1250.0, 1500.0, 1750.0, 2000.0, 3000.0}; //{100.0, 400.0, 800.0, 1200.0, 1600.0, 2000.0};
    const double JetpTBorders[npTBins] = {5.0, 10.0, 20.0, 40.0, 60.0, 80.0, 100.0, 150.0, 200.0};
    const int nConst = 5;
    const int ConstVals[nConst] = {1, 2, 3, 5, 10}; // constituent values to check

    int constMin = 0; // min number of constituents in matchd jet 2

    for (int iNorm = 0; iNorm < nNorm; iNorm++)
    {

        TFile *jetFile;
        jetFile = TFile::Open(Form("/Users/hadi/focalJetAnalysis/jetOutput/20230728_pythia8_JetJet_%d-%dGeV/Merged_0Mass.root", pTHardRange[iNorm], pTHardRange[iNorm + 1]));

        TFile *fout;
        fout = new TFile(Form("Data20230816/2D/2D_20230816_pythia8_JetJet_%d-%dGeV_Merged_Output_0Mass.root", pTHardRange[iNorm], pTHardRange[iNorm + 1]), "RECREATE");

        TTree *jetTree = (TTree *)jetFile->Get("jetTree");

        TH2F *hRespMatrix[nR];
        TH2F *hRespMatrixE[nR][nEBins - 1];
        TH2F *hRespMatrixpT[nR][npTBins - 1];

        for (int Rvalue = 0; Rvalue < nR; Rvalue++)
        {
            // histograms
            hRespMatrix[Rvalue] = new TH2F(Form("hRespMatrix_R%d", Rvalue), Form("Detector ECAL fraction vs Particle neutral energy fraction, R=%0.1f;#it{E}_{ECAL}/#it{E}_{jet}^{det} (GeV);#it{E}_{neutral}/#it{E}_{jet}^{part} (GeV)", Rvals[Rvalue]), 41, -0.0125, 1.0125, 41, -0.0125, 1.0125);

            jetTree->Draw(Form("matchedjetECALEnergyFrac:jetECALEnergyFrac>>hRespMatrix_R%d", Rvalue), Form("jetR==%d && jetEta>%f && jetEta<%f && jetEta_match>%f && jetEta_match<%f && jet_distmatch<%f && jetParts_match > %d&& jetParts > %d", int(Rvals[Rvalue] * 10), etaMin + Rvals[Rvalue], etaMax - Rvals[Rvalue], etaMin + Rvals[Rvalue], etaMax - Rvals[Rvalue], Rvals[Rvalue] * 0.6, constMin, constMin), "goff");

            // hRespMatrix[Rvalue]->Scale(normalizations[iNorm] / nFolders[iNorm]);
            // hRespMatrix[Rvalue]->Scale(1./hRespMatrix[Rvalue]->Integral(), "width"); //only perform this normalization at final stage of plotting

            for (int iE = 0; iE < nEBins - 1; ++iE)
            {
                hRespMatrixE[Rvalue][iE] = new TH2F(Form("hRespMatrixE%d_R%d", iE, Rvalue), Form("Detector ECAL fraction vs Particle neutral energy fraction, E_{jet}^{det}: %d - %d GeV, R=%0.1f;#it{E}_{ECAL}/#it{E}_{jet}^{det} (GeV);#it{E}_{neutral}/#it{E}_{jet}^{part} (GeV)", int(JetEBorders[iE]), int(JetEBorders[iE + 1]), Rvals[Rvalue]), 41, -0.0125, 1.0125, 41, -0.0125, 1.0125);
                jetTree->Draw(Form("matchedjetECALEnergyFrac:jetECALEnergyFrac>>hRespMatrixE%d_R%d", iE, Rvalue), Form("jetR==%d && jetE>=%d  && jetE<%d && jetEta>%f && jetEta<%f && jetEta_match>%f && jetEta_match<%f && jet_distmatch<%f && jetParts_match > %d&& jetParts > %d", int(Rvals[Rvalue] * 10), int(JetEBorders[iE]), int(JetEBorders[iE + 1]), etaMin + Rvals[Rvalue], etaMax - Rvals[Rvalue], etaMin + Rvals[Rvalue], etaMax - Rvals[Rvalue], Rvals[Rvalue] * 0.6, constMin, constMin), "goff");

                // hRespMatrixE[Rvalue][iE]->Scale(normalizations[iNorm] / nFolders[iNorm]);
                // hRespMatrixE[Rvalue][iE]->Scale(1./hRespMatrixE[Rvalue][iE]->Integral(), "width"); //only perform this normalization at final stage of plotting
            }

            for (int ipT = 0; ipT < npTBins - 1; ++ipT)
            {
                hRespMatrixpT[Rvalue][ipT] = new TH2F(Form("hRespMatrixpT%d_R%d", ipT, Rvalue), Form("Detector ECAL fraction vs Particle neutral energy fraction, pT_{jet}^{det}: %d - %d GeV, R=%0.1f;#it{E}_{ECAL}/#it{E}_{jet}^{det} (GeV);#it{E}_{neutral}/#it{E}_{jet}^{part} (GeV)", int(JetpTBorders[ipT]), int(JetpTBorders[ipT + 1]), Rvals[Rvalue]), 41, -0.0125, 1.0125, 41, -0.0125, 1.0125);
                jetTree->Draw(Form("matchedjetECALEnergyFrac:jetECALEnergyFrac>>hRespMatrixpT%d_R%d", ipT, Rvalue), Form("jetR==%d && jetpT>=%d  && jetpT<%d && jetEta>%f && jetEta<%f && jetEta_match>%f && jetEta_match<%f && jet_distmatch<%f && jetParts_match > %d&& jetParts > %d", int(Rvals[Rvalue] * 10), int(JetpTBorders[ipT]), int(JetpTBorders[ipT + 1]), etaMin + Rvals[Rvalue], etaMax - Rvals[Rvalue], etaMin + Rvals[Rvalue], etaMax - Rvals[Rvalue], Rvals[Rvalue] * 0.6, constMin, constMin), "goff");

                // hRespMatrixpT[Rvalue][ipT]->Scale(normalizations[iNorm] / nFolders[iNorm]);
                // hRespMatrixpT[Rvalue][ipT]->Scale(1./hRespMatrixpT[Rvalue][ipT]->Integral(), "width"); //only perform this normalization at final stage of plotting
            }
        }

        fout->Write();
        fout->Close();
    }

    gSystem->Exec("hadd Data20230816/2D/2D_Merged_0Mass.root Data20230816/2D/2D_20230816_pythia8_JetJet_*_0Mass.root");
}