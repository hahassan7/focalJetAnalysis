Double_t EtaToTheta(Double_t eta)
{
  return 2*TMath::ATan(TMath::Exp(-eta))*180/TMath::Pi();
}

AliGenerator* GeneratorCustom(TString opt = "") {

    AliGenBox *gener = new AliGenBox(1);
    gener->SetPtRange(0., 10.);
    gener->SetPhiRange(0.0, 360.0);  // full polar angle around beam axis
    gener->SetEtaRange(3.0, 6.0);
    gener->SetOrigin(0,0,0);
    //vertex position
    gener->SetSigma(0,0,0);         //Sigma in (X,Y,Z) (cm) on IP position
    gener->SetPart(kGamma);
    //gener->SetPart(kPi0);
    //gener->SetPart(kPiPlus);
    //gener->SetPart(kElectron);
    return gener;
}

/**AliGenerator* GeneratorCustom(TString opt = "") {
    float energy = 14000;  // GeV
    AliGenerator * generator = 0x0;
    AliGenPythiaFOCAL *gener = new AliGenPythiaFOCAL(-1); 
    gener->SetMomentumRange(0,999999); 
    gener->SetThetaRange(0., 45.); 
    gener->SetYRange(-12,12); 
    gener->SetPtRange(0,1000); 
    gener->SetEnergyCMS(energy); // LHC energy 
    gener->SetSigma(0, 0, 5.3); // Sigma in (X,Y,Z) (cm) on IP position 
    gener->SetCutVertexZ(1.); // Truncate at 1 sigma 
    gener->SetVertexSmear(kPerEvent); // Smear per event 
    gener->SetTrackingFlag(1); // Particle transport 
    gener->SetProcess(kPyMb); // Min. bias events 
    gener->SetDecayPhotonInFOCAL(kTRUE);
    gener->SetCheckFOCAL(kTRUE);  
    gener->SetFOCALEta(3.5, 6.2);
    gener->SetTriggerParticleMinPt(2.0);
    generator = gener;
    return generator;
}**/
