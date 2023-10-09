//void RunJetAnalysis(Int_t startFolder, Int_t endFolder, Int_t startEvent, Int_t endEvent, TString outputComment) {
void RunJetAnalysis(Int_t startFolder, Int_t endFolder, Int_t startEvent, Int_t endEvent) {
  /*const char *dataSampleTag="resutls_PythiaMBtrig_ceneknew";
  const char* simFolder="/eos/user/h/hahassan/FocalSim/PythiaMBtrig";

  const char *outputdir = "results_jets_PythiaMBtrig";
  const char *clustersFolder="../clustering/results_PythiaMBtrig";*/

  const char *dataSampleTag="MBPythia";
  const char* simFolder="/scratch/project_2003583/focal-full-sim/output/20230202_Spagetti_pythiaMBtrig-2";

 // const char *outputdir = Form("/scratch/project_2003583/focal-full-sim/output_jets/%s", outputComment.Data());
  const char *outputdir = "/scratch/project_2003583/focal-full-sim/output_jets/20230202_Spagetti_pythiaMBtrig-2";
  const char *clustersFolder="/scratch/project_2003583/focal-full-sim/output_clustering/20230202_Spagetti_pythiaMBtrig-2";

  gSystem->mkdir(outputdir, true);

  gROOT->ProcessLine(Form(".x AnalyzeJets.C(%d,%d,%d,%d,\"%s\",\"%s\",\"%s\",\"%s\")",startFolder, endFolder, startEvent, endEvent, simFolder, clustersFolder, dataSampleTag, outputdir));
}
