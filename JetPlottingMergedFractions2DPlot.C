void JetPlottingMergedFractions2DPlot(){

//setting NormValue = 8 if merged file in use

    const double limitYconst = 0.2;
    const double limitYpT    = 0.2;
    const double limitYEta     = 0.2;
    const double limitYconstmin = 0.0;
    const double limitYpTmin    = 0.0;
    const double limitYEtamin     = 0.0;

    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    // Open the input file
    TFile *fin;
    if(NormValue==0)fin = new TFile("JetJetOutput/July2023/fractions/fractions/20230722_pythia8_JetJet_5-10GeV_Merged_Output.root", "READ");
    if(NormValue==1)fin = new TFile("JetJetOutput/July2023/fractions/fractions/20230722_pythia8_JetJet_10-20GeV_Merged_Output.root", "READ");
    if(NormValue==2)fin = new TFile("JetJetOutput/July2023/fractions/fractions/20230722_pythia8_JetJet_20-30GeV_Merged_Output.root", "READ");
    if(NormValue==3)fin = new TFile("JetJetOutput/July2023/fractions/fractions/20230722_pythia8_JetJet_30-40GeV_Merged_Output.root", "READ");
    if(NormValue==4)fin = new TFile("JetJetOutput/July2023/fractions/fractions/20230722_pythia8_JetJet_40-60GeV_Merged_Output.root", "READ");
    if(NormValue==5)fin = new TFile("JetJetOutput/July2023/fractions/fractions/20230722_pythia8_JetJet_60-100GeV_Merged_Output.root", "READ");
    if(NormValue==6)fin = new TFile("JetJetOutput/July2023/fractions/fractions/20230722_pythia8_JetJet_100-200GeV_Merged_Output.root", "READ");
    if(NormValue==7)fin = new TFile("JetJetOutput/July2023/fractions/fractions/20230722_pythia8_JetJet_200-GeV_Merged_Output.root", "READ");
    if(NormValue==8)fin = new TFile("JetJetOutput/July2023/fractions/fractions/Merged10-GeV.root", "READ");


    //Add normalization to probability density function: 1/Int, "width"
}