
void drawPDG()
{
	//TFile *jetFile = TFile::Open("/home/lmh/alice/fromPuhti/v19/20230127_PythiaMBTrig-2_34/Merged.root");
    //TFile *jetFile = TFile::Open("/home/lmh/alice/fromPuhti/v19/EnergyTests/20230127_PythiaMBTrig-2/006/analysisJets.root");
    TFile *jetFile = TFile::Open("/home/lmh/alice/fromPuhti/v19/20230127_PythiaMBTrig-2_232002/Merged.root");

    //TTree *jetTree = (TTree *)jetFile->Get("jetTree");
    //TTree *TruthjetTree = (TTree *)jetFile->Get("truthjetTree");
    TTree *PDGTree = (TTree *)jetFile->Get("inPDGTree");

    //TH1F *hPDG = (TH1F*)jetFile->Get("oneConstPDG_h"); //change to this if you want to study 1,2,3... const jets
    TH1F *hPDG = new TH1F("hPDG","hPDG", 12001, -6000.5, 6000.5); //8001, -4000.5, 4000.5
    TH1F *hEtaY = new TH1F("hEtaY","Particle level eta-y, in eta acceptance 3.4-5.8", 201, -1, 7.5);


    TCanvas *c1 = new TCanvas("c1", "PDG distribution of all incoming particles in eta acceptance 3.4-5.8");

    PDGTree->Draw("partPDG>>hPDG"); 


    TCanvas *c2 = new TCanvas("c2", "Eta-y distribution of all incoming particles (in eta acceptance 3.4-5.8)");
    c2->SetLogy();
    PDGTree->Draw("partEta-partRap>>hEtaY");

/*
    for (int i = 1; i < hPDG->GetNbinsX()+1; ++i)
    {
        if(hPDG->GetBinContent(i)!=0){
            cout << "PDG " << -1*((hPDG->GetNbinsX()-1)/2.0)+i-1 << " " << hPDG->GetBinContent(i) << endl;
        }
        
    }
    cout << hPDG->GetEntries() << endl;*/
    


}


