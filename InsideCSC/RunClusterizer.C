void RunClusterizer(TString indir, TString clusteringOutputFileDir, Int_t jobstart, Int_t jobend)
{
    gSystem->Load("AliJBaseTrack_cxx.so");
    gSystem->Load("AliJHMRCluster_cxx.so");
    gSystem->Load("AliJHMREvent_cxx.so");
    gROOT->ProcessLine(Form(".x Clusterizer.C(\"%s\",\"%s\",%d,%d)", indir.Data(), clusteringOutputFileDir.Data(), jobstart, jobend));
}
