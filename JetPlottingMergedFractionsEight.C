/* Macro for plotting some histograms from jet trees 
This macro loops through all three jet R values.
*/
//Fill fraction histos from trees
void JetPlottingMergedFractionsEight(int NormValue=0) //loops over R and saves to same file
{   

    const int nNorm = 8;
    const Float_t normalizations[nNorm] = {0.0610658,0.00716477,0.000557627,0.000107816,4.31694e-05,9.62255e-06,1.24904e-06,5.99517e-08};

    //int NormValue = 0;

    TFile *jetFile;
    
    //Read in the data
    //if(NormValue==0)jetFile = TFile::Open("/home/lmh/alice/fromPuhti/April2023/pTmin2_20230417_pythia8_JetJet_5-15GeV/Merged.root");
    if(NormValue==0)jetFile = TFile::Open("/home/lmh/alice/fromPuhti/July2023/ptmin10/fractions/fractions_ptmin10_20230722_pythia8_JetJet_5-10GeV/Merged.root");
    if(NormValue==1)jetFile = TFile::Open("/home/lmh/alice/fromPuhti/July2023/ptmin10/fractions/fractions_ptmin10_20230722_pythia8_JetJet_10-20GeV/Merged.root");
    if(NormValue==2)jetFile = TFile::Open("/home/lmh/alice/fromPuhti/July2023/ptmin10/fractions/fractions_ptmin10_20230722_pythia8_JetJet_20-30GeV/Merged.root");
    if(NormValue==3)jetFile = TFile::Open("/home/lmh/alice/fromPuhti/July2023/ptmin10/fractions/fractions_ptmin10_20230722_pythia8_JetJet_30-40GeV/Merged.root");
    if(NormValue==4)jetFile = TFile::Open("/home/lmh/alice/fromPuhti/July2023/ptmin10/fractions/fractions_ptmin10_20230722_pythia8_JetJet_40-60GeV/Merged.root");
    if(NormValue==5)jetFile = TFile::Open("/home/lmh/alice/fromPuhti/July2023/ptmin10/fractions/fractions_ptmin10_20230722_pythia8_JetJet_60-100GeV/Merged.root");
    if(NormValue==6)jetFile = TFile::Open("/home/lmh/alice/fromPuhti/July2023/ptmin10/fractions/fractions_ptmin10_20230722_pythia8_JetJet_100-200GeV/Merged.root");
    if(NormValue==7)jetFile = TFile::Open("/home/lmh/alice/fromPuhti/July2023/ptmin10/fractions/fractions_ptmin10_20230722_pythia8_JetJet_200-GeV/Merged.root");

    const Int_t nR = 3; 
    const Float_t Rvals[nR] = {0.2, 0.4, 0.6}; // Cone radii
    //const double EtaMin[2] = {3.4,3.9}; 
    //const double EtaMax[2] = {5.0,5.5};
    

    const Float_t etaMin = 3.4; const Float_t etaMax = 5.5; //changed in June23, 3.4 to 5.5  

    //add eta range binning here, which depends on the value of R.   
    const Int_t nEtaBins  = 3;
    const double EtaBinBorders[nR][nEtaBins] = {{3.6,4.5,5.3},{3.8,4.5,5.1},{4.0,4.5,4.9}};
    TString etaRange[nEtaBins-1] = {"3.4+R < #eta_{jet} < 4.5","4.5 < #eta_{jet} < 5.5-R"}; //different range for R=0.6, 0.4, 0.2. 
    const Int_t nEBins  = 6;//16;
    const Int_t npTBins  = 9;//16;
    const double JetEBorders[nEBins] = {100.0, 400.0, 800.0, 1200.0, 1600.0, 2000.0}; 
    const double JetpTBorders[npTBins] = {5.0, 10.0, 20.0, 40.0, 60.0, 80.0, 100.0, 150.0, 200.0}; 
    const Int_t nConst = 5;
    const Int_t ConstVals[nConst] = {1, 2, 3, 5, 10}; // constituent values to check

    const int nCol = 10;
    const int gcolors[nCol]={1,2,6,4,7,1,2,4,6,7};
    const int gmarkers[nCol]={4,4,4,8,25,21,8,21,25,4};

    // Create a ROOT file to save the histograms
    TFile *fout;

    int constMin = 0; //min number of constituents in matchd jet 1

    if(NormValue==0)fout = new TFile("JetJetOutput/July2023/ptmin10/fractions/20230722_pythia8_JetJet_5-10GeV_Merged_Output.root", "RECREATE");
    if(NormValue==1)fout = new TFile("JetJetOutput/July2023/ptmin10/fractions/20230722_pythia8_JetJet_10-20GeV_Merged_Output.root", "RECREATE");
    if(NormValue==2)fout = new TFile("JetJetOutput/July2023/ptmin10/fractions/20230722_pythia8_JetJet_20-30GeV_Merged_Output.root", "RECREATE");
    if(NormValue==3)fout = new TFile("JetJetOutput/July2023/ptmin10/fractions/20230722_pythia8_JetJet_30-40GeV_Merged_Output.root", "RECREATE");
    if(NormValue==4)fout = new TFile("JetJetOutput/July2023/ptmin10/fractions/20230722_pythia8_JetJet_40-60GeV_Merged_Output.root", "RECREATE");
    if(NormValue==5)fout = new TFile("JetJetOutput/July2023/ptmin10/fractions/20230722_pythia8_JetJet_60-100GeV_Merged_Output.root", "RECREATE");
    if(NormValue==6)fout = new TFile("JetJetOutput/July2023/ptmin10/fractions/20230722_pythia8_JetJet_100-200GeV_Merged_Output.root", "RECREATE");
    if(NormValue==7)fout = new TFile("JetJetOutput/July2023/ptmin10/fractions/20230722_pythia8_JetJet_200-GeV_Merged_Output.root", "RECREATE");

    TTree *jetTree = (TTree *)jetFile->Get("jetTree");
    TTree *TruthjetTree = (TTree *)jetFile->Get("truthjetTree");
    //TTree *PDGTree = (TTree *)jetFile->Get("inPDGTree");


    //histograms: 
    //ADD eta bins, ADD const bins
    //also make the same histos for pT
    //also add histograms for truth jet NEF and not just matched

    TH1F *hDetECAL[nR]; //ECALfracdet matched
    TH1F *hPartNeutralMatch[nR]; //neutralfrac matched
    TH1F *hDetECALAll[nR]; //ECALfracdet all not just matched
    TH1F *hPartNeutral[nR]; //neutralfrac all not just matched

    TH1F *hDetECALEta[nR][nEtaBins-1]; //ECALfracdet matched
    TH1F *hPartNeutralMatchEta[nR][nEtaBins-1]; //neutralfrac matched
    TH1F *hDetECALAllEta[nR][nEtaBins-1]; //ECALfracdet all not just matched
    TH1F *hPartNeutralEta[nR][nEtaBins-1]; //neutralfrac all not just matched

    TH1F *hDetECALEtaConst[nR][nEtaBins-1][nConst]; //ECALfracdet matched
    TH1F *hPartNeutralMatchEtaConst[nR][nEtaBins-1][nConst]; //neutralfrac matched
    TH1F *hDetECALAllEtaConst[nR][nEtaBins-1][nConst]; //ECALfracdet all not just matched
    TH1F *hPartNeutralEtaConst[nR][nEtaBins-1][nConst]; //neutralfrac all not just matched

    TH1F *hDetECALConst[nR][nConst]; //ECALfracdet matched
    TH1F *hPartNeutralMatchConst[nR][nConst]; //neutralfrac matched
    TH1F *hDetECALAllConst[nR][nConst]; //ECALfracdet all not just matched
    TH1F *hPartNeutralConst[nR][nConst]; //neutralfrac all not just matched



    TH1F *hDetECALE[nEBins-1][nR]; //ECALfracdet matched
    TH1F *hPartNeutralMatchE[nEBins-1][nR]; //neutralfrac matched
    TH1F *hDetECALAllE[nEBins-1][nR]; //ECALfracdet all not just matched
    TH1F *hPartNeutralE[nEBins-1][nR]; //neutralfrac all not just matched

    TH1F *hDetECALEEta[nEBins-1][nR][nEtaBins-1]; //ECALfracdet matched
    TH1F *hPartNeutralMatchEEta[nEBins-1][nR][nEtaBins-1]; //neutralfrac matched
    TH1F *hDetECALAllEEta[nEBins-1][nR][nEtaBins-1]; //ECALfracdet all not just matched
    TH1F *hPartNeutralEEta[nEBins-1][nR][nEtaBins-1]; //neutralfrac all not just matched

    TH1F *hDetECALEConst[nEBins-1][nR][nConst]; //ECALfracdet matched
    TH1F *hPartNeutralMatchEConst[nEBins-1][nR][nConst]; //neutralfrac matched
    TH1F *hDetECALAllEConst[nEBins-1][nR][nConst]; //ECALfracdet all not just matched
    TH1F *hPartNeutralEConst[nEBins-1][nR][nConst]; //neutralfrac all not just matched



    TH1F *hDetECALpT[npTBins-1][nR]; //ECALfracdet matched
    TH1F *hPartNeutralMatchpT[npTBins-1][nR]; //neutralfrac matched
    TH1F *hDetECALAllpT[npTBins-1][nR]; //ECALfracdet all not just matched
    TH1F *hPartNeutralpT[npTBins-1][nR]; //neutralfrac all not just matched

    TH1F *hDetECALpTEta[npTBins-1][nR][nEtaBins-1]; //ECALfracdet matched
    TH1F *hPartNeutralMatchpTEta[npTBins-1][nR][nEtaBins-1]; //neutralfrac matched
    TH1F *hDetECALAllpTEta[npTBins-1][nR][nEtaBins-1]; //ECALfracdet all not just matched
    TH1F *hPartNeutralpTEta[npTBins-1][nR][nEtaBins-1]; //neutralfrac all not just matched

    TH1F *hDetECALpTConst[npTBins-1][nR][nConst]; //ECALfracdet matched
    TH1F *hPartNeutralMatchpTConst[npTBins-1][nR][nConst]; //neutralfrac matched
    TH1F *hDetECALAllpTConst[npTBins-1][nR][nConst]; //ECALfracdet all not just matched
    TH1F *hPartNeutralpTConst[npTBins-1][nR][nConst]; //neutralfrac all not just matched

    for(int Rvalue=0; Rvalue<nR; Rvalue++){

        //create histograms
        hDetECAL[Rvalue] = new TH1F(Form("hDetECAL_R%d", Rvalue), Form("Matched detector level jet ECAL energy fraction distribution, R = %0.1f", Rvals[Rvalue]), 21, -0.025, 1.025); //ECALfracdet matched
        hPartNeutralMatch[Rvalue] = new TH1F(Form("hPartNeutralMatch_R%d", Rvalue), Form("Matched particle level jet neutral energy fraction distribution, R = %0.1f", Rvals[Rvalue]), 21, -0.025, 1.025); //neutralfrac matched
        hDetECALAll[Rvalue] = new TH1F(Form("hDetECALAll_R%d", Rvalue), Form("All detector level jet ECAL energy fraction distribution, R = %0.1f", Rvals[Rvalue]), 21, -0.025, 1.025); //ECALfracdet all not just matched
        hPartNeutral[Rvalue] = new TH1F(Form("hPartNeutral_R%d", Rvalue), Form("All particle level jet neutral energy fraction distribution, R = %0.1f", Rvals[Rvalue]), 21, -0.025, 1.025);; //neutralfrac all not just matched
        
        //begin filling the histograms
        jetTree->Draw(Form("jetECALEnergyFrac>>hDetECAL_R%d", Rvalue), Form("jetR==%d && jetEta>=%f && jetEta<%f && jetEta_match>=%f && jetEta_match<%f && jet_distmatch<%f && jetParts_match > %d && jetParts > %d", int(Rvals[Rvalue] * 10),etaMin+Rvals[Rvalue],etaMax-Rvals[Rvalue],etaMin+Rvals[Rvalue],etaMax-Rvals[Rvalue], Rvals[Rvalue]*0.6, constMin, constMin), "goff"); 

        jetTree->Draw(Form("matchedjetECALEnergyFrac>>hPartNeutralMatch_R%d", Rvalue), Form("jetR==%d && jetEta>=%f && jetEta<%f && jetEta_match>=%f && jetEta_match<%f && jet_distmatch<%f && jetParts_match > %d && jetParts > %d", int(Rvals[Rvalue] * 10), etaMin+Rvals[Rvalue],etaMax-Rvals[Rvalue],etaMin+Rvals[Rvalue],etaMax-Rvals[Rvalue], Rvals[Rvalue]*0.6, constMin, constMin), "goff");

        jetTree->Draw(Form("jetECALEnergyFrac>>hDetECALAll_R%d", Rvalue), Form("jetR==%d && jetEta>=%f && jetEta<%f && jetParts > %d", int(Rvals[Rvalue] * 10),etaMin+Rvals[Rvalue],etaMax-Rvals[Rvalue], constMin), "goff");

        TruthjetTree->Draw(Form("truthjetNeutralEnergyFrac>>hPartNeutral_R%d", Rvalue), Form("truthjetR==%d && truthjetEta>=%f && truthjetEta<%f && jetParts_truth > %d", int(Rvals[Rvalue] * 10),etaMin+Rvals[Rvalue],etaMax-Rvals[Rvalue], constMin), "goff"); 

        hDetECAL[Rvalue]->Scale(normalizations[NormValue]);
        hPartNeutralMatch[Rvalue]->Scale(normalizations[NormValue]);
        hDetECALAll[Rvalue]->Scale(normalizations[NormValue]);
        hPartNeutral[Rvalue]->Scale(normalizations[NormValue]);

        for (int ieta = 0; ieta < nEtaBins-1; ++ieta)
        {
            hDetECALEta[Rvalue][ieta] = new TH1F(Form("hDetECALEta_%d_R%d", ieta, Rvalue), Form("Matched detector level jet ECAL energy fraction distribution, %s, R = %0.1f", etaRange[ieta].Data(), Rvals[Rvalue]), 21, -0.025, 1.025); //ECALfracdet matched
            hPartNeutralMatchEta[Rvalue][ieta] = new TH1F(Form("hPartNeutralMatchEta_%d_R%d", ieta, Rvalue), Form("Matched particle level jet neutral energy fraction distribution, %s, R = %0.1f", etaRange[ieta].Data(), Rvals[Rvalue]), 21, -0.025, 1.025); //neutralfrac matched
            hDetECALAllEta[Rvalue][ieta] = new TH1F(Form("hDetECALAllEta_%d_R%d", ieta, Rvalue), Form("All detector level jet ECAL energy fraction distribution, %s, R = %0.1f", etaRange[ieta].Data(), Rvals[Rvalue]), 21, -0.025, 1.025); //ECALfracdet all not just matched
            hPartNeutralEta[Rvalue][ieta] = new TH1F(Form("hPartNeutralEta_%d_R%d", ieta, Rvalue), Form("All particle level jet neutral energy fraction distribution, %s, R = %0.1f", etaRange[ieta].Data(), Rvals[Rvalue]), 21, -0.025, 1.025);; //neutralfrac all not just matched

            jetTree->Draw(Form("jetECALEnergyFrac>>hDetECALEta_%d_R%d", ieta, Rvalue), Form("jetR==%d && jetEta>=%f && jetEta<%f && jetEta_match>=%f && jetEta_match<%f && jet_distmatch<%f && jetParts_match > %d && jetParts > %d", int(Rvals[Rvalue] * 10), EtaBinBorders[Rvalue][ieta], EtaBinBorders[Rvalue][ieta+1], EtaBinBorders[Rvalue][ieta], EtaBinBorders[Rvalue][ieta+1], Rvals[Rvalue]*0.6, constMin, constMin), "goff");

            jetTree->Draw(Form("matchedjetECALEnergyFrac>>hPartNeutralMatchEta_%d_R%d", ieta, Rvalue), Form("jetR==%d && jetEta>=%f && jetEta<%f && jetEta_match>=%f && jetEta_match<%f && jet_distmatch<%f && jetParts_match > %d && jetParts > %d", int(Rvals[Rvalue] * 10), EtaBinBorders[Rvalue][ieta], EtaBinBorders[Rvalue][ieta+1], EtaBinBorders[Rvalue][ieta], EtaBinBorders[Rvalue][ieta+1], Rvals[Rvalue]*0.6, constMin, constMin), "goff"); 

            jetTree->Draw(Form("jetECALEnergyFrac>>hDetECALAllEta_%d_R%d", ieta, Rvalue), Form("jetR==%d && jetEta>=%f && jetEta<%f && jetParts > %d", int(Rvals[Rvalue] * 10), EtaBinBorders[Rvalue][ieta], EtaBinBorders[Rvalue][ieta+1], constMin), "goff");

            TruthjetTree->Draw(Form("truthjetNeutralEnergyFrac>>hPartNeutralEta_%d_R%d", ieta, Rvalue), Form("truthjetR==%d && truthjetEta>=%f && truthjetEta<%f && jetParts_truth > %d", int(Rvals[Rvalue] * 10), EtaBinBorders[Rvalue][ieta], EtaBinBorders[Rvalue][ieta+1], constMin), "goff"); 


            if(hDetECALEta[Rvalue][ieta]->GetEntries()==0) hDetECALEta[Rvalue][ieta]->Fill(-2);
            if(hPartNeutralMatchEta[Rvalue][ieta]->GetEntries()==0) hPartNeutralMatchEta[Rvalue][ieta]->Fill(-2);
            if(hDetECALAllEta[Rvalue][ieta]->GetEntries()==0) hDetECALAllEta[Rvalue][ieta]->Fill(-2);
            if(hPartNeutralEta[Rvalue][ieta]->GetEntries()==0) hPartNeutralEta[Rvalue][ieta]->Fill(-2);

            hDetECALEta[Rvalue][ieta]->Scale(normalizations[NormValue]);
            hPartNeutralMatchEta[Rvalue][ieta]->Scale(normalizations[NormValue]);
            hDetECALAllEta[Rvalue][ieta]->Scale(normalizations[NormValue]);
            hPartNeutralEta[Rvalue][ieta]->Scale(normalizations[NormValue]);

            for (int iconst = 0; iconst < nConst; ++iconst)
            {

                hDetECALEtaConst[Rvalue][ieta][iconst] = new TH1F(Form("hDetECALEta_%d_Const_%d_R%d", ieta, iconst, Rvalue), Form("Matched detector level jet ECAL energy fraction distribution, %s, R = %0.1f", etaRange[ieta].Data(), Rvals[Rvalue]), 21, -0.025, 1.025); //ECALfracdet matched
                hPartNeutralMatchEtaConst[Rvalue][ieta][iconst] = new TH1F(Form("hPartNeutralMatchEta_%d_Const_%d_R%d", ieta, iconst, Rvalue), Form("Matched particle level jet neutral energy fraction distribution, %s, R = %0.1f", etaRange[ieta].Data(), Rvals[Rvalue]), 21, -0.025, 1.025); //neutralfrac matched
                hDetECALAllEtaConst[Rvalue][ieta][iconst] = new TH1F(Form("hDetECALAllEta_%d_Const_%d_R%d", ieta, iconst, Rvalue), Form("All detector level jet ECAL energy fraction distribution, %s, R = %0.1f", etaRange[ieta].Data(), Rvals[Rvalue]), 21, -0.025, 1.025); //ECALfracdet all not just matched
                hPartNeutralEtaConst[Rvalue][ieta][iconst] = new TH1F(Form("hPartNeutralEta_%d_Const_%d_R%d", ieta, iconst, Rvalue), Form("All particle level jet neutral energy fraction distribution, %s, R = %0.1f", etaRange[ieta].Data(), Rvals[Rvalue]), 21, -0.025, 1.025);; //neutralfrac all not just matched

                jetTree->Draw(Form("jetECALEnergyFrac>>hDetECALEta_%d_Const_%d_R%d", ieta, iconst, Rvalue), Form("jetR==%d && jetEta>=%f && jetEta<%f && jetEta_match>=%f && jetEta_match<%f && jet_distmatch<%f && jetParts == %d", int(Rvals[Rvalue] * 10), EtaBinBorders[Rvalue][ieta], EtaBinBorders[Rvalue][ieta+1], EtaBinBorders[Rvalue][ieta], EtaBinBorders[Rvalue][ieta+1], Rvals[Rvalue]*0.6, ConstVals[iconst]), "goff");

                jetTree->Draw(Form("matchedjetECALEnergyFrac>>hPartNeutralMatchEta_%d_Const_%d_R%d", ieta, iconst, Rvalue), Form("jetR==%d && jetEta>=%f && jetEta<%f && jetEta_match>=%f && jetEta_match<%f && jet_distmatch<%f && jetParts_match == %d", int(Rvals[Rvalue] * 10), EtaBinBorders[Rvalue][ieta], EtaBinBorders[Rvalue][ieta+1], EtaBinBorders[Rvalue][ieta], EtaBinBorders[Rvalue][ieta+1], Rvals[Rvalue]*0.6, ConstVals[iconst]), "goff"); 

                jetTree->Draw(Form("jetECALEnergyFrac>>hDetECALAllEta_%d_Const_%d_R%d", ieta, iconst, Rvalue), Form("jetR==%d && jetEta>=%f && jetEta<%f && jetParts == %d", int(Rvals[Rvalue] * 10), EtaBinBorders[Rvalue][ieta], EtaBinBorders[Rvalue][ieta+1], ConstVals[iconst]), "goff");

                TruthjetTree->Draw(Form("truthjetNeutralEnergyFrac>>hPartNeutralEta_%d_Const_%d_R%d", ieta, iconst, Rvalue), Form("truthjetR==%d && truthjetEta>=%f && truthjetEta<%f && jetParts_truth == %d", int(Rvals[Rvalue] * 10), EtaBinBorders[Rvalue][ieta], EtaBinBorders[Rvalue][ieta+1], ConstVals[iconst]), "goff"); 


                if(hDetECALEtaConst[Rvalue][ieta][iconst]->GetEntries()==0) hDetECALEtaConst[Rvalue][ieta][iconst]->Fill(-2);
                if(hPartNeutralMatchEtaConst[Rvalue][ieta][iconst]->GetEntries()==0) hPartNeutralMatchEtaConst[Rvalue][ieta][iconst]->Fill(-2);
                if(hDetECALAllEtaConst[Rvalue][ieta][iconst]->GetEntries()==0) hDetECALAllEtaConst[Rvalue][ieta][iconst]->Fill(-2);
                if(hPartNeutralEtaConst[Rvalue][ieta][iconst]->GetEntries()==0) hPartNeutralEtaConst[Rvalue][ieta][iconst]->Fill(-2);

                hDetECALEtaConst[Rvalue][ieta][iconst]->Scale(normalizations[NormValue]);
                hPartNeutralMatchEtaConst[Rvalue][ieta][iconst]->Scale(normalizations[NormValue]);
                hDetECALAllEtaConst[Rvalue][ieta][iconst]->Scale(normalizations[NormValue]);
                hPartNeutralEtaConst[Rvalue][ieta][iconst]->Scale(normalizations[NormValue]);
            }
        }

        for (int iconst = 0; iconst < nConst; ++iconst)
        {

            hDetECALConst[Rvalue][iconst] = new TH1F(Form("hDetECALConst_%d_R%d", iconst, Rvalue), Form("Matched detector level jet ECAL energy fraction distribution, constituent N = %d, R = %0.1f", ConstVals[iconst], Rvals[Rvalue]), 21, -0.025, 1.025);

            hPartNeutralMatchConst[Rvalue][iconst] = new TH1F(Form("hPartNeutralMatchConst_%d_R%d", iconst, Rvalue), Form("Matched particle level jet neutral constituent energy fraction distribution, constituent N = %d, R = %0.1f", ConstVals[iconst], Rvals[Rvalue]), 21, -0.025, 1.025);

            hDetECALAllConst[Rvalue][iconst] = new TH1F(Form("hDetECALAllConst_%d_R%d", iconst, Rvalue), Form("All detector level jet ECAL energy fraction distribution, constituent N = %d, R = %0.1f", ConstVals[iconst], Rvals[Rvalue]), 21, -0.025, 1.025);

            hPartNeutralConst[Rvalue][iconst] = new TH1F(Form("hPartNeutralConst_%d_R%d", iconst, Rvalue), Form("All particle level jet neutral constituent energy fraction distribution, constituent N = %d, R = %0.1f", ConstVals[iconst], Rvals[Rvalue]), 21, -0.025, 1.025);
            
            //Fill histos

            jetTree->Draw(Form("jetECALEnergyFrac>>hDetECALConst_%d_R%d", iconst, Rvalue), Form("jetR==%d && jetEta>=%f && jetEta<%f && jetEta_match>=%f && jetEta_match<%f && jet_distmatch<%f  && jetParts == %d", int(Rvals[Rvalue] * 10), etaMin+Rvals[Rvalue],etaMax-Rvals[Rvalue],etaMin+Rvals[Rvalue],etaMax-Rvals[Rvalue], Rvals[Rvalue]*0.6, ConstVals[iconst]), "goff");

            jetTree->Draw(Form("matchedjetECALEnergyFrac>>hPartNeutralMatchConst_%d_R%d", iconst, Rvalue), Form("jetR==%d && jetEta>=%f && jetEta<%f && jetEta_match>=%f && jetEta_match<%f && jet_distmatch<%f && jetParts_match == %d ", int(Rvals[Rvalue] * 10),etaMin+Rvals[Rvalue],etaMax-Rvals[Rvalue],etaMin+Rvals[Rvalue],etaMax-Rvals[Rvalue], Rvals[Rvalue]*0.6,  ConstVals[iconst]), "goff"); 

            jetTree->Draw(Form("jetECALEnergyFrac>>hDetECALAllConst_%d_R%d", iconst, Rvalue), Form("jetR==%d && jetEta>=%f && jetEta<%f && jetParts == %d", int(Rvals[Rvalue] * 10),etaMin+Rvals[Rvalue],etaMax-Rvals[Rvalue], ConstVals[iconst]), "goff");

            TruthjetTree->Draw(Form("truthjetNeutralEnergyFrac>>hPartNeutralConst_%d_R%d", iconst, Rvalue), Form("truthjetR==%d && truthjetEta>=%f && truthjetEta<%f && jetParts_truth == %d", int(Rvals[Rvalue] * 10), etaMin+Rvals[Rvalue],etaMax-Rvals[Rvalue], ConstVals[iconst]), "goff"); 

            if(hDetECALConst[Rvalue][iconst]->GetEntries()==0) hDetECALConst[Rvalue][iconst]->Fill(-2);
            if(hPartNeutralMatchConst[Rvalue][iconst]->GetEntries()==0) hPartNeutralMatchConst[Rvalue][iconst]->Fill(-2);
            if(hDetECALAllConst[Rvalue][iconst]->GetEntries()==0) hDetECALAllConst[Rvalue][iconst]->Fill(-2);
            if(hPartNeutralConst[Rvalue][iconst]->GetEntries()==0) hPartNeutralConst[Rvalue][iconst]->Fill(-2);

            hDetECALConst[Rvalue][iconst]->Scale(normalizations[NormValue]);
            hPartNeutralMatchConst[Rvalue][iconst]->Scale(normalizations[NormValue]);
            hDetECALAllConst[Rvalue][iconst]->Scale(normalizations[NormValue]);
            hPartNeutralConst[Rvalue][iconst]->Scale(normalizations[NormValue]);
        }

        //energy loop
        // ADD  a constituent loop 
        // an eta range loop as well
        for (int iE = 0; iE < nEBins-1; ++iE) 
        {
            // create histograms
            hDetECALE[iE][Rvalue] = new TH1F(Form("hDetECALE_%d_R%d", iE, Rvalue), Form("Matched detector level jet ECAL energy fraction distribution, E_{jet}^{det}: %d - %d GeV, R = %0.1f", int(JetEBorders[iE]), int(JetEBorders[iE+1]), Rvals[Rvalue]), 21, -0.025, 1.025);

            hPartNeutralMatchE[iE][Rvalue] = new TH1F(Form("hPartNeutralMatchE_%d_R%d", iE, Rvalue), Form("Matched particle level jet neutral constituent energy fraction distribution, E_{jet}^{part}: %d - %d GeV, R = %0.1f", int(JetEBorders[iE]), int(JetEBorders[iE+1]), Rvals[Rvalue]), 21, -0.025, 1.025);

            hDetECALAllE[iE][Rvalue] = new TH1F(Form("hDetECALAllE_%d_R%d", iE, Rvalue), Form("All detector level jet ECAL energy fraction distribution, E_{jet}^{det}: %d - %d GeV, R = %0.1f", int(JetEBorders[iE]), int(JetEBorders[iE+1]), Rvals[Rvalue]), 21, -0.025, 1.025);

            hPartNeutralE[iE][Rvalue] = new TH1F(Form("hPartNeutralE_%d_R%d", iE, Rvalue), Form("All particle level jet neutral constituent energy fraction distribution, E_{jet}^{part}: %d - %d GeV, R = %0.1f", int(JetEBorders[iE]), int(JetEBorders[iE+1]), Rvals[Rvalue]), 21, -0.025, 1.025);
            
            //Fill histos

            jetTree->Draw(Form("jetECALEnergyFrac>>hDetECALE_%d_R%d", iE, Rvalue), Form("jetR==%d && jetE>=%d  && jetE<%d  && jetEta>=%f && jetEta<%f && jetEta_match>=%f && jetEta_match<%f && jet_distmatch<%f && jetParts_match > %d && jetParts > %d", int(Rvals[Rvalue] * 10), int(JetEBorders[iE]), int(JetEBorders[iE+1]),etaMin+Rvals[Rvalue],etaMax-Rvals[Rvalue],etaMin+Rvals[Rvalue],etaMax-Rvals[Rvalue], Rvals[Rvalue]*0.6, constMin, constMin), "goff");

            jetTree->Draw(Form("matchedjetECALEnergyFrac>>hPartNeutralMatchE_%d_R%d", iE, Rvalue), Form("jetR==%d && jetE_match>=%d  && jetE_match<%d  && jetEta>=%f && jetEta<%f && jetEta_match>=%f && jetEta_match<%f && jet_distmatch<%f && jetParts_match > %d && jetParts > %d", int(Rvals[Rvalue] * 10), int(JetEBorders[iE]), int(JetEBorders[iE+1]),etaMin+Rvals[Rvalue],etaMax-Rvals[Rvalue],etaMin+Rvals[Rvalue],etaMax-Rvals[Rvalue], Rvals[Rvalue]*0.6, constMin, constMin), "goff"); 

            jetTree->Draw(Form("jetECALEnergyFrac>>hDetECALAllE_%d_R%d", iE, Rvalue), Form("jetR==%d && jetE>=%d  && jetE<%d  && jetEta>=%f && jetEta<%f && jetParts > %d", int(Rvals[Rvalue] * 10), int(JetEBorders[iE]), int(JetEBorders[iE+1]),etaMin+Rvals[Rvalue],etaMax-Rvals[Rvalue], constMin), "goff");

            TruthjetTree->Draw(Form("truthjetNeutralEnergyFrac>>hPartNeutralE_%d_R%d", iE, Rvalue), Form("truthjetR==%d && truthjetE>=%d  && truthjetE<%d  && truthjetEta>=%f && truthjetEta<%f && jetParts_truth > %d", int(Rvals[Rvalue] * 10), int(JetEBorders[iE]), int(JetEBorders[iE+1]),etaMin+Rvals[Rvalue],etaMax-Rvals[Rvalue], constMin), "goff"); 

            if(hDetECALE[iE][Rvalue]->GetEntries()==0) hDetECALE[iE][Rvalue]->Fill(-2);
            if(hPartNeutralMatchE[iE][Rvalue]->GetEntries()==0) hPartNeutralMatchE[iE][Rvalue]->Fill(-2);
            if(hDetECALAllE[iE][Rvalue]->GetEntries()==0) hDetECALAllE[iE][Rvalue]->Fill(-2);
            if(hPartNeutralE[iE][Rvalue]->GetEntries()==0) hPartNeutralE[iE][Rvalue]->Fill(-2);

            hDetECALE[iE][Rvalue]->Scale(normalizations[NormValue]);
            hPartNeutralMatchE[iE][Rvalue]->Scale(normalizations[NormValue]);
            hDetECALAllE[iE][Rvalue]->Scale(normalizations[NormValue]);
            hPartNeutralE[iE][Rvalue]->Scale(normalizations[NormValue]);

            //eta loop here
            for (int ieta = 0; ieta < nEtaBins-1; ++ieta)
            {
                /* code */
                // create histograms
                hDetECALEEta[iE][Rvalue][ieta] = new TH1F(Form("hDetECALE_%dEta_%d_R%d", iE, ieta, Rvalue), Form("Matched detector level jet ECAL energy fraction distribution, E_{jet}^{det}: %d - %d GeV, %s, R = %0.1f", int(JetEBorders[iE]), int(JetEBorders[iE+1]), etaRange[ieta].Data(), Rvals[Rvalue]), 21, -0.025, 1.025);

                hPartNeutralMatchEEta[iE][Rvalue][ieta] = new TH1F(Form("hPartNeutralMatchE_%dEta_%d_R%d", iE, ieta, Rvalue), Form("Matched particle level jet neutral constituent energy fraction distribution, E_{jet}^{part}: %d - %d GeV, %s, R = %0.1f", int(JetEBorders[iE]), int(JetEBorders[iE+1]), etaRange[ieta].Data(), Rvals[Rvalue]), 21, -0.025, 1.025);

                hDetECALAllEEta[iE][Rvalue][ieta] = new TH1F(Form("hDetECALAllE_%dEta_%d_R%d", iE, ieta, Rvalue), Form("All detector level jet ECAL energy fraction distribution, E_{jet}^{det}: %d - %d GeV, %s, R = %0.1f", int(JetEBorders[iE]), int(JetEBorders[iE+1]), etaRange[ieta].Data(), Rvals[Rvalue]), 21, -0.025, 1.025);

                hPartNeutralEEta[iE][Rvalue][ieta] = new TH1F(Form("hPartNeutralE_%dEta_%d_R%d", iE, ieta, Rvalue), Form("All particle level jet neutral constituent energy fraction distribution, E_{jet}^{part}: %d - %d GeV, %s, R = %0.1f", int(JetEBorders[iE]), int(JetEBorders[iE+1]), etaRange[ieta].Data(), Rvals[Rvalue]), 21, -0.025, 1.025);
                
                //Fill histos

                jetTree->Draw(Form("jetECALEnergyFrac>>hDetECALE_%dEta_%d_R%d", iE, ieta, Rvalue), Form("jetR==%d && jetE>=%d  && jetE<%d  && jetEta>=%f && jetEta<%f && jetEta_match>=%f && jetEta_match<%f && jet_distmatch<%f && jetParts_match > %d && jetParts > %d", int(Rvals[Rvalue] * 10), int(JetEBorders[iE]), int(JetEBorders[iE+1]), EtaBinBorders[Rvalue][ieta], EtaBinBorders[Rvalue][ieta+1], EtaBinBorders[Rvalue][ieta], EtaBinBorders[Rvalue][ieta+1], Rvals[Rvalue]*0.6, constMin, constMin), "goff");

                jetTree->Draw(Form("matchedjetECALEnergyFrac>>hPartNeutralMatchE_%dEta_%d_R%d", iE, ieta, Rvalue), Form("jetR==%d && jetE_match>=%d  && jetE_match<%d  && jetEta>=%f && jetEta<%f && jetEta_match>=%f && jetEta_match<%f && jet_distmatch<%f && jetParts_match > %d && jetParts > %d", int(Rvals[Rvalue] * 10), int(JetEBorders[iE]), int(JetEBorders[iE+1]), EtaBinBorders[Rvalue][ieta], EtaBinBorders[Rvalue][ieta+1], EtaBinBorders[Rvalue][ieta], EtaBinBorders[Rvalue][ieta+1], Rvals[Rvalue]*0.6, constMin, constMin), "goff"); 

                jetTree->Draw(Form("jetECALEnergyFrac>>hDetECALAllE_%dEta_%d_R%d", iE, ieta, Rvalue), Form("jetR==%d && jetE>=%d  && jetE<%d  && jetEta>=%f && jetEta<%f && jetParts > %d", int(Rvals[Rvalue] * 10), int(JetEBorders[iE]), int(JetEBorders[iE+1]), EtaBinBorders[Rvalue][ieta], EtaBinBorders[Rvalue][ieta+1], constMin), "goff");

                TruthjetTree->Draw(Form("truthjetNeutralEnergyFrac>>hPartNeutralE_%dEta_%d_R%d", iE, ieta, Rvalue), Form("truthjetR==%d && truthjetE>=%d  && truthjetE<%d  && truthjetEta>=%f && truthjetEta<%f && jetParts_truth > %d", int(Rvals[Rvalue] * 10), int(JetEBorders[iE]), int(JetEBorders[iE+1]), EtaBinBorders[Rvalue][ieta], EtaBinBorders[Rvalue][ieta+1], constMin), "goff"); 

                if(hDetECALEEta[iE][Rvalue][ieta]->GetEntries()==0) hDetECALEEta[iE][Rvalue][ieta]->Fill(-2);
                if(hPartNeutralMatchEEta[iE][Rvalue][ieta]->GetEntries()==0) hPartNeutralMatchEEta[iE][Rvalue][ieta]->Fill(-2);
                if(hDetECALAllEEta[iE][Rvalue][ieta]->GetEntries()==0) hDetECALAllEEta[iE][Rvalue][ieta]->Fill(-2);
                if(hPartNeutralEEta[iE][Rvalue][ieta]->GetEntries()==0) hPartNeutralEEta[iE][Rvalue][ieta]->Fill(-2);

                hDetECALEEta[iE][Rvalue][ieta]->Scale(normalizations[NormValue]);
                hPartNeutralMatchEEta[iE][Rvalue][ieta]->Scale(normalizations[NormValue]);
                hDetECALAllEEta[iE][Rvalue][ieta]->Scale(normalizations[NormValue]);
                hPartNeutralEEta[iE][Rvalue][ieta]->Scale(normalizations[NormValue]);


            }

            for (int iconst = 0; iconst < nConst; ++iconst)
            {   
                // create histograms
                hDetECALEConst[iE][Rvalue][iconst] = new TH1F(Form("hDetECALE_%dConst_%d_R%d", iE, iconst, Rvalue), Form("Matched detector level jet ECAL energy fraction distribution, E_{jet}^{det}: %d - %d GeV, constituent N = %d, R = %0.1f", int(JetEBorders[iE]), int(JetEBorders[iE+1]), ConstVals[iconst], Rvals[Rvalue]), 21, -0.025, 1.025);

                hPartNeutralMatchEConst[iE][Rvalue][iconst] = new TH1F(Form("hPartNeutralMatchE_%dConst_%d_R%d", iE, iconst, Rvalue), Form("Matched particle level jet neutral constituent energy fraction distribution, E_{jet}^{part}: %d - %d GeV, constituent N = %d, R = %0.1f", int(JetEBorders[iE]), int(JetEBorders[iE+1]), ConstVals[iconst], Rvals[Rvalue]), 21, -0.025, 1.025);

                hDetECALAllEConst[iE][Rvalue][iconst] = new TH1F(Form("hDetECALAllE_%dConst_%d_R%d", iE, iconst, Rvalue), Form("All detector level jet ECAL energy fraction distribution, E_{jet}^{det}: %d - %d GeV, constituent N = %d, R = %0.1f", int(JetEBorders[iE]), int(JetEBorders[iE+1]), ConstVals[iconst], Rvals[Rvalue]), 21, -0.025, 1.025);

                hPartNeutralEConst[iE][Rvalue][iconst] = new TH1F(Form("hPartNeutralE_%dConst_%d_R%d", iE, iconst, Rvalue), Form("All particle level jet neutral constituent energy fraction distribution, E_{jet}^{part}: %d - %d GeV, constituent N = %d, R = %0.1f", int(JetEBorders[iE]), int(JetEBorders[iE+1]), ConstVals[iconst], Rvals[Rvalue]), 21, -0.025, 1.025);
                
                //Fill histos

                jetTree->Draw(Form("jetECALEnergyFrac>>hDetECALE_%dConst_%d_R%d", iE, iconst, Rvalue), Form("jetR==%d && jetE>=%d  && jetE<%d  && jetEta>=%f && jetEta<%f && jetEta_match>=%f && jetEta_match<%f && jet_distmatch<%f  && jetParts == %d", int(Rvals[Rvalue] * 10), int(JetEBorders[iE]), int(JetEBorders[iE+1]),etaMin+Rvals[Rvalue],etaMax-Rvals[Rvalue],etaMin+Rvals[Rvalue],etaMax-Rvals[Rvalue], Rvals[Rvalue]*0.6, ConstVals[iconst]), "goff");

                jetTree->Draw(Form("matchedjetECALEnergyFrac>>hPartNeutralMatchE_%dConst_%d_R%d", iE, iconst, Rvalue), Form("jetR==%d && jetE_match>=%d  && jetE_match<%d  && jetEta>=%f && jetEta<%f && jetEta_match>=%f && jetEta_match<%f && jet_distmatch<%f && jetParts_match == %d ", int(Rvals[Rvalue] * 10), int(JetEBorders[iE]), int(JetEBorders[iE+1]),etaMin+Rvals[Rvalue],etaMax-Rvals[Rvalue],etaMin+Rvals[Rvalue],etaMax-Rvals[Rvalue], Rvals[Rvalue]*0.6,  ConstVals[iconst]), "goff"); 

                jetTree->Draw(Form("jetECALEnergyFrac>>hDetECALAllE_%dConst_%d_R%d", iE, iconst, Rvalue), Form("jetR==%d && jetE>=%d  && jetE<%d  && jetEta>=%f && jetEta<%f && jetParts == %d", int(Rvals[Rvalue] * 10), int(JetEBorders[iE]), int(JetEBorders[iE+1]),etaMin+Rvals[Rvalue],etaMax-Rvals[Rvalue], ConstVals[iconst]), "goff");

                TruthjetTree->Draw(Form("truthjetNeutralEnergyFrac>>hPartNeutralE_%dConst_%d_R%d", iE, iconst, Rvalue), Form("truthjetR==%d && truthjetE>=%d  && truthjetE<%d  && truthjetEta>=%f && truthjetEta<%f && jetParts_truth == %d", int(Rvals[Rvalue] * 10), int(JetEBorders[iE]), int(JetEBorders[iE+1]),etaMin+Rvals[Rvalue],etaMax-Rvals[Rvalue], ConstVals[iconst]), "goff"); 

                if(hDetECALEConst[iE][Rvalue][iconst]->GetEntries()==0) hDetECALEConst[iE][Rvalue][iconst]->Fill(-2);
                if(hPartNeutralMatchEConst[iE][Rvalue][iconst]->GetEntries()==0) hPartNeutralMatchEConst[iE][Rvalue][iconst]->Fill(-2);
                if(hDetECALAllEConst[iE][Rvalue][iconst]->GetEntries()==0) hDetECALAllEConst[iE][Rvalue][iconst]->Fill(-2);
                if(hPartNeutralEConst[iE][Rvalue][iconst]->GetEntries()==0) hPartNeutralEConst[iE][Rvalue][iconst]->Fill(-2);

                hDetECALEConst[iE][Rvalue][iconst]->Scale(normalizations[NormValue]);
                hPartNeutralMatchEConst[iE][Rvalue][iconst]->Scale(normalizations[NormValue]);
                hDetECALAllEConst[iE][Rvalue][iconst]->Scale(normalizations[NormValue]);
                hPartNeutralEConst[iE][Rvalue][iconst]->Scale(normalizations[NormValue]);
            }
        }

        //you can add a pT loop here
        for (int ipT = 0; ipT < npTBins-1; ++ipT) 
        {
            // create histograms
            hDetECALpT[ipT][Rvalue] = new TH1F(Form("hDetECALpT_%d_R%d", ipT, Rvalue), Form("Matched detector level jet ECAL energy fraction distribution, pT_{jet}^{det}: %d - %d GeV, R = %0.1f", int(JetpTBorders[ipT]), int(JetpTBorders[ipT+1]), Rvals[Rvalue]), 21, -0.025, 1.025);

            hPartNeutralMatchpT[ipT][Rvalue] = new TH1F(Form("hPartNeutralMatchpT_%d_R%d", ipT, Rvalue), Form("Matched particle level jet neutral constituent energy fraction distribution, pT_{jet}^{part}: %d - %d GeV, R = %0.1f", int(JetpTBorders[ipT]), int(JetpTBorders[ipT+1]), Rvals[Rvalue]), 21, -0.025, 1.025);

            hDetECALAllpT[ipT][Rvalue] = new TH1F(Form("hDetECALAllpT_%d_R%d", ipT, Rvalue), Form("All detector level jet ECAL energy fraction distribution, pT_{jet}^{det}: %d - %d GeV, R = %0.1f", int(JetpTBorders[ipT]), int(JetpTBorders[ipT+1]), Rvals[Rvalue]), 21, -0.025, 1.025);

            hPartNeutralpT[ipT][Rvalue] = new TH1F(Form("hPartNeutralpT_%d_R%d", ipT, Rvalue), Form("All particle level jet neutral constituent energy fraction distribution, pT_{jet}^{part}: %d - %d GeV, R = %0.1f", int(JetpTBorders[ipT]), int(JetpTBorders[ipT+1]), Rvals[Rvalue]), 21, -0.025, 1.025);
            
            //Fill histos

            jetTree->Draw(Form("jetECALEnergyFrac>>hDetECALpT_%d_R%d", ipT, Rvalue), Form("jetR==%d && jetpT>=%d  && jetpT<%d  && jetEta>=%f && jetEta<%f && jetEta_match>=%f && jetEta_match<%f && jet_distmatch<%f && jetParts_match > %d && jetParts > %d", int(Rvals[Rvalue] * 10), int(JetpTBorders[ipT]), int(JetpTBorders[ipT+1]),etaMin+Rvals[Rvalue],etaMax-Rvals[Rvalue],etaMin+Rvals[Rvalue],etaMax-Rvals[Rvalue], Rvals[Rvalue]*0.6, constMin, constMin), "goff");

            jetTree->Draw(Form("matchedjetECALEnergyFrac>>hPartNeutralMatchpT_%d_R%d", ipT, Rvalue), Form("jetR==%d && jetpT_match>=%d  && jetpT_match<%d  && jetEta>=%f && jetEta<%f && jetEta_match>=%f && jetEta_match<%f && jet_distmatch<%f && jetParts_match > %d && jetParts > %d", int(Rvals[Rvalue] * 10), int(JetpTBorders[ipT]), int(JetpTBorders[ipT+1]),etaMin+Rvals[Rvalue],etaMax-Rvals[Rvalue],etaMin+Rvals[Rvalue],etaMax-Rvals[Rvalue], Rvals[Rvalue]*0.6, constMin, constMin), "goff"); 

            jetTree->Draw(Form("jetECALEnergyFrac>>hDetECALAllpT_%d_R%d", ipT, Rvalue), Form("jetR==%d && jetpT>=%d  && jetpT<%d  && jetEta>=%f && jetEta<%f && jetParts > %d", int(Rvals[Rvalue] * 10), int(JetpTBorders[ipT]), int(JetpTBorders[ipT+1]),etaMin+Rvals[Rvalue],etaMax-Rvals[Rvalue], constMin), "goff");

            TruthjetTree->Draw(Form("truthjetNeutralEnergyFrac>>hPartNeutralpT_%d_R%d", ipT, Rvalue), Form("truthjetR==%d && truthjetpT>=%d  && truthjetpT<%d  && truthjetEta>=%f && truthjetEta<%f && jetParts_truth > %d", int(Rvals[Rvalue] * 10), int(JetpTBorders[ipT]), int(JetpTBorders[ipT+1]),etaMin+Rvals[Rvalue],etaMax-Rvals[Rvalue], constMin), "goff"); 

            if(hDetECALpT[ipT][Rvalue]->GetEntries()==0) hDetECALpT[ipT][Rvalue]->Fill(-2);
            if(hPartNeutralMatchpT[ipT][Rvalue]->GetEntries()==0) hPartNeutralMatchpT[ipT][Rvalue]->Fill(-2);
            if(hDetECALAllpT[ipT][Rvalue]->GetEntries()==0) hDetECALAllpT[ipT][Rvalue]->Fill(-2);
            if(hPartNeutralpT[ipT][Rvalue]->GetEntries()==0) hPartNeutralpT[ipT][Rvalue]->Fill(-2);

            hDetECALpT[ipT][Rvalue]->Scale(normalizations[NormValue]);
            hPartNeutralMatchpT[ipT][Rvalue]->Scale(normalizations[NormValue]);
            hDetECALAllpT[ipT][Rvalue]->Scale(normalizations[NormValue]);
            hPartNeutralpT[ipT][Rvalue]->Scale(normalizations[NormValue]);

            //eta loop here
            for (int ieta = 0; ieta < nEtaBins-1; ++ieta)
            {
                /* code */
                // create histograms
                hDetECALpTEta[ipT][Rvalue][ieta] = new TH1F(Form("hDetECALpT_%dEta_%d_R%d", ipT, ieta, Rvalue), Form("Matched detector level jet ECAL energy fraction distribution, pT_{jet}^{det}: %d - %d GeV, %s, R = %0.1f", int(JetpTBorders[ipT]), int(JetpTBorders[ipT+1]), etaRange[ieta].Data(), Rvals[Rvalue]), 21, -0.025, 1.025);

                hPartNeutralMatchpTEta[ipT][Rvalue][ieta] = new TH1F(Form("hPartNeutralMatchpT_%dEta_%d_R%d", ipT, ieta, Rvalue), Form("Matched particle level jet neutral constituent energy fraction distribution, pT_{jet}^{part}: %d - %d GeV, %s, R = %0.1f", int(JetpTBorders[ipT]), int(JetpTBorders[ipT+1]), etaRange[ieta].Data(), Rvals[Rvalue]), 21, -0.025, 1.025);

                hDetECALAllpTEta[ipT][Rvalue][ieta] = new TH1F(Form("hDetECALAllpT_%dEta_%d_R%d", ipT, ieta, Rvalue), Form("All detector level jet ECAL energy fraction distribution, pT_{jet}^{det}: %d - %d GeV, %s, R = %0.1f", int(JetpTBorders[ipT]), int(JetpTBorders[ipT+1]), etaRange[ieta].Data(), Rvals[Rvalue]), 21, -0.025, 1.025);

                hPartNeutralpTEta[ipT][Rvalue][ieta] = new TH1F(Form("hPartNeutralpT_%dEta_%d_R%d", ipT, ieta, Rvalue), Form("All particle level jet neutral constituent energy fraction distribution, pT_{jet}^{part}: %d - %d GeV, %s, R = %0.1f", int(JetpTBorders[ipT]), int(JetpTBorders[ipT+1]), etaRange[ieta].Data(), Rvals[Rvalue]), 21, -0.025, 1.025);
                
                //Fill histos

                jetTree->Draw(Form("jetECALEnergyFrac>>hDetECALpT_%dEta_%d_R%d", ipT, ieta, Rvalue), Form("jetR==%d && jetpT>=%d  && jetpT<%d  && jetEta>=%f && jetEta<%f && jetEta_match>=%f && jetEta_match<%f && jet_distmatch<%f && jetParts_match > %d && jetParts > %d", int(Rvals[Rvalue] * 10), int(JetpTBorders[ipT]), int(JetpTBorders[ipT+1]), EtaBinBorders[Rvalue][ieta], EtaBinBorders[Rvalue][ieta+1], EtaBinBorders[Rvalue][ieta], EtaBinBorders[Rvalue][ieta+1], Rvals[Rvalue]*0.6, constMin, constMin), "goff");

                jetTree->Draw(Form("matchedjetECALEnergyFrac>>hPartNeutralMatchpT_%dEta_%d_R%d", ipT, ieta, Rvalue), Form("jetR==%d && jetpT_match>=%d  && jetpT_match<%d  && jetEta>=%f && jetEta<%f && jetEta_match>=%f && jetEta_match<%f && jet_distmatch<%f && jetParts_match > %d && jetParts > %d", int(Rvals[Rvalue] * 10), int(JetpTBorders[ipT]), int(JetpTBorders[ipT+1]), EtaBinBorders[Rvalue][ieta], EtaBinBorders[Rvalue][ieta+1], EtaBinBorders[Rvalue][ieta], EtaBinBorders[Rvalue][ieta+1], Rvals[Rvalue]*0.6, constMin, constMin), "goff"); 

                jetTree->Draw(Form("jetECALEnergyFrac>>hDetECALAllpT_%dEta_%d_R%d", ipT, ieta, Rvalue), Form("jetR==%d && jetpT>=%d  && jetpT<%d  && jetEta>=%f && jetEta<%f && jetParts > %d", int(Rvals[Rvalue] * 10), int(JetpTBorders[ipT]), int(JetpTBorders[ipT+1]), EtaBinBorders[Rvalue][ieta], EtaBinBorders[Rvalue][ieta+1], constMin), "goff");

                TruthjetTree->Draw(Form("truthjetNeutralEnergyFrac>>hPartNeutralpT_%dEta_%d_R%d", ipT, ieta, Rvalue), Form("truthjetR==%d && truthjetE>=%d  && truthjetE<%d  && truthjetEta>=%f && truthjetEta<%f && jetParts_truth > %d", int(Rvals[Rvalue] * 10), int(JetpTBorders[ipT]), int(JetpTBorders[ipT+1]), EtaBinBorders[Rvalue][ieta], EtaBinBorders[Rvalue][ieta+1], constMin), "goff"); 

                if(hDetECALpTEta[ipT][Rvalue][ieta]->GetEntries()==0) hDetECALpTEta[ipT][Rvalue][ieta]->Fill(-2);
                if(hPartNeutralMatchpTEta[ipT][Rvalue][ieta]->GetEntries()==0) hPartNeutralMatchpTEta[ipT][Rvalue][ieta]->Fill(-2);
                if(hDetECALAllpTEta[ipT][Rvalue][ieta]->GetEntries()==0) hDetECALAllpTEta[ipT][Rvalue][ieta]->Fill(-2);
                if(hPartNeutralpTEta[ipT][Rvalue][ieta]->GetEntries()==0) hPartNeutralpTEta[ipT][Rvalue][ieta]->Fill(-2);

                hDetECALpTEta[ipT][Rvalue][ieta]->Scale(normalizations[NormValue]);
                hPartNeutralMatchpTEta[ipT][Rvalue][ieta]->Scale(normalizations[NormValue]);
                hDetECALAllpTEta[ipT][Rvalue][ieta]->Scale(normalizations[NormValue]);
                hPartNeutralpTEta[ipT][Rvalue][ieta]->Scale(normalizations[NormValue]);
            }

            for (int iconst = 0; iconst < nConst; ++iconst)
            {
                /* code */
                // create histograms
                hDetECALpTConst[ipT][Rvalue][iconst] = new TH1F(Form("hDetECALpT_%dConst_%d_R%d", ipT, iconst, Rvalue), Form("Matched detector level jet ECAL energy fraction distribution, pT_{jet}^{det}: %d - %d GeV, constituent N = %d, R = %0.1f", int(JetpTBorders[ipT]), int(JetpTBorders[ipT+1]), ConstVals[iconst], Rvals[Rvalue]), 21, -0.025, 1.025);

                hPartNeutralMatchpTConst[ipT][Rvalue][iconst] = new TH1F(Form("hPartNeutralMatchpT_%dConst_%d_R%d", ipT, iconst, Rvalue), Form("Matched particle level jet neutral constituent energy fraction distribution, pT_{jet}^{part}: %d - %d GeV, constituent N = %d, R = %0.1f", int(JetpTBorders[ipT]), int(JetpTBorders[ipT+1]), ConstVals[iconst], Rvals[Rvalue]), 21, -0.025, 1.025);

                hDetECALAllpTConst[ipT][Rvalue][iconst] = new TH1F(Form("hDetECALAllpT_%dConst_%d_R%d", ipT, iconst, Rvalue), Form("All detector level jet ECAL energy fraction distribution, pT_{jet}^{det}: %d - %d GeV, constituent N = %d, R = %0.1f", int(JetpTBorders[ipT]), int(JetpTBorders[ipT+1]), ConstVals[iconst], Rvals[Rvalue]), 21, -0.025, 1.025);

                hPartNeutralpTConst[ipT][Rvalue][iconst] = new TH1F(Form("hPartNeutralpT_%dConst_%d_R%d", ipT, iconst, Rvalue), Form("All particle level jet neutral constituent energy fraction distribution, pT_{jet}^{part}: %d - %d GeV, constituent N = %d, R = %0.1f", int(JetpTBorders[ipT]), int(JetpTBorders[ipT+1]), ConstVals[iconst], Rvals[Rvalue]), 21, -0.025, 1.025);
                
                //Fill histos

                jetTree->Draw(Form("jetECALEnergyFrac>>hDetECALpT_%dConst_%d_R%d", ipT, iconst, Rvalue), Form("jetR==%d && jetpT>=%d  && jetpT<%d  && jetEta>=%f && jetEta<%f && jetEta_match>=%f && jetEta_match<%f && jet_distmatch<%f  && jetParts == %d", int(Rvals[Rvalue] * 10), int(JetpTBorders[ipT]), int(JetpTBorders[ipT+1]),etaMin+Rvals[Rvalue],etaMax-Rvals[Rvalue],etaMin+Rvals[Rvalue],etaMax-Rvals[Rvalue], Rvals[Rvalue]*0.6, ConstVals[iconst]), "goff");

                jetTree->Draw(Form("matchedjetECALEnergyFrac>>hPartNeutralMatchpT_%dConst_%d_R%d", ipT, iconst, Rvalue), Form("jetR==%d && jetpT_match>=%d  && jetpT_match<%d  && jetEta>=%f && jetEta<%f && jetEta_match>=%f && jetEta_match<%f && jet_distmatch<%f && jetParts_match == %d ", int(Rvals[Rvalue] * 10), int(JetpTBorders[ipT]), int(JetpTBorders[ipT+1]),etaMin+Rvals[Rvalue],etaMax-Rvals[Rvalue],etaMin+Rvals[Rvalue],etaMax-Rvals[Rvalue], Rvals[Rvalue]*0.6,  ConstVals[iconst]), "goff"); 

                jetTree->Draw(Form("jetECALEnergyFrac>>hDetECALAllpT_%dConst_%d_R%d", ipT, iconst, Rvalue), Form("jetR==%d && jetpT>=%d  && jetpT<%d  && jetEta>=%f && jetEta<%f && jetParts == %d", int(Rvals[Rvalue] * 10), int(JetpTBorders[ipT]), int(JetpTBorders[ipT+1]),etaMin+Rvals[Rvalue],etaMax-Rvals[Rvalue], ConstVals[iconst]), "goff");

                TruthjetTree->Draw(Form("truthjetNeutralEnergyFrac>>hPartNeutralpT_%dConst_%d_R%d", ipT, iconst, Rvalue), Form("truthjetR==%d && truthjetE>=%d  && truthjetE<%d  && truthjetEta>=%f && truthjetEta<%f && jetParts_truth == %d", int(Rvals[Rvalue] * 10), int(JetpTBorders[ipT]), int(JetpTBorders[ipT+1]),etaMin+Rvals[Rvalue],etaMax-Rvals[Rvalue], ConstVals[iconst]), "goff"); 

                if(hDetECALpTConst[ipT][Rvalue][iconst]->GetEntries()==0) hDetECALpTConst[ipT][Rvalue][iconst]->Fill(-2);
                if(hPartNeutralMatchpTConst[ipT][Rvalue][iconst]->GetEntries()==0) hPartNeutralMatchpTConst[ipT][Rvalue][iconst]->Fill(-2);
                if(hDetECALAllpTConst[ipT][Rvalue][iconst]->GetEntries()==0) hDetECALAllpTConst[ipT][Rvalue][iconst]->Fill(-2);
                if(hPartNeutralpTConst[ipT][Rvalue][iconst]->GetEntries()==0) hPartNeutralpTConst[ipT][Rvalue][iconst]->Fill(-2);

                hDetECALpTConst[ipT][Rvalue][iconst]->Scale(normalizations[NormValue]);
                hPartNeutralMatchpTConst[ipT][Rvalue][iconst]->Scale(normalizations[NormValue]);
                hDetECALAllpTConst[ipT][Rvalue][iconst]->Scale(normalizations[NormValue]);
                hPartNeutralpTConst[ipT][Rvalue][iconst]->Scale(normalizations[NormValue]);
            }
        }
    }

    fout->Write();
    fout->Close();


}