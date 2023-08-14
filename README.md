# focalJetAnalysis
Lots of focal jet analysis macros.

## In folder InsideCSC/ one finds the jet analysis macros run using the CSC Puhti supercomputer.
  - **AnalyzeJetsGrid.C** is the main macro, with the current version of the jet matching implementation. 
  - **AnalyzeJetsGridv1.C** uses the older version of jet matching, but can be useful for checking all the jet fractions, since not all the fractions are implemented in the main analysis macro (e.g. HCAL fraction, photon+electron fraction).
  - The output of the macros is merged manually with hadd to a file named Merged.root within each pThard bin folder

## The rest of the macros are run locally to fill and plot histograms from the jet output
  - In folder jetOutput/ one finds example data sets with which the plotting macros can be tested.
      - These are currently the jet analysis outputs from CSC merged into files where there is now 500k events in each file. No normalization has yet been done.

 ### Macros for plotting response matrices, JES, JER:
  - **JetPlottingMerged.C** runs the histogram filling, pT and E response matrix and , also two different eta bins.
    - This macro is run for all pThard bins separately, and the files are to be merged manually using hadd, separately for each value of R. 
    - The macro scales all the histograms with the corresponding factor normalizations[pThardbin], which are the cross section / ntrials given by pythia. **Now includes a factor 1/500 because of us running 500*1k events, change this if needed!**
    - **JetPlottingMergedFile.C** then calculates the mean and median and stdev (JES and JER) of the merged histograms, pT and E
    - Finally, **focalJetResolutions.C** and **focalJetResolutionsEnergy.C** plot the distributions, comparing between two R values and also two eta bins
      - NOTE: the first version of this plotting macro was provided by Friederike Bock!
  - Similarly,**JetPlottingMergedEta.C** runs the histo filling, but just in smaller eta bins
    - **JetPlottingMergedFileEta.C** calculates the mean, median, stdev
    - **focalJetResolutionsEta.C** then plots these histograms, only for three eta bins at a time, chosen with exampleBins[3]
    - The macros **focalJetResolutionsTest.C** and **focalJetResolutionsETest.C** draw the eta bins in one figure for just one value of R.
  - **PlotEtaFigs.C** is used to plot just a single deltaE distribution, can change to deltapT
  - **paperPlotsHeader.h** is used in these plotting macros (Not my own handwriting, provided by Friederike Bock)
  - For running the JES and JER plotting with NEF fractions, do the following:
      - Run **JetPlottingMergedEnergyJES.C** over all files, parameter 0-7, to fill the histograms (you can make changes here, like looking at the NEF instead of detector level ECAL fraction
      - Then merge all the files with **hadd**
      - Run over the merged file with **JetPlottingMergedFileEnergyJES.C** to create the JES and JER distributions
      - Then plot with plotting macro **focalJetResolutionsEnergyJES.C** 

  ### Macros for plotting fractions:
  - One can now use the same input file as for the other plotting, used to be a different input file
  - **JetPlottingMergedFractionsEight.C** goes through all values of R and fills fraction histograms for the given pthard bin data
  - Then there are three plotting macros that take in the same file, which the previous macro created
    -  **JetPlottingMergedFractionsEightPlot.C** plots the general fraction histograms, no selection on E or pT bin
    -  **JetPlottingMergedFractionsEightPlotEnergyBins.C** plots same fractions in jet energy bins
    -  **JetPlottingMergedFractionsEightPlotpTBins.C** plots same fractions in jet pT bins

    
### Macros for plotting jet spectra:
  - **SpectraFillHistos.C** fills the histos for detector level jet pT E and const spectra
  - **Spectra.C** plots them
  - **SpectraFillHistosPartLevel.C** fills histos for matched particle level jets
  - **SpectraPartLevel.C** plots them
  - **SpectraFillHistosPartLevelTruth.C** fills histos for all truth level particle jets
  - **SpectraPartLevelTruth.C** plots them
  - **drawPDG.C** plots the spectra of particle PDG codes and eta-y for all particles in the jets

### Old versions of macros, ignore these:
  - **JetPlotting.C** and **JetPlottingEnergy.C** are old versions of the code, not used anymore, saved just in case they might be needed later
  - **JetPlottingConst.C** was used to test differences in jets with different number of constituents
  - **JetPlottingMergedFractions.C** was used to plot the fractions from a separate file when there was only 5 pthard bins, and doesn't loop over R
  - This was followed with the plotting macro **JetPlottingMergedFileFractions.C**
  - **PlotDetResolutionFractions.C** used to be used to plot the fractions in a simple manner
  - **PlotDetResolution.C** was used to look at total energy differences
  - **NewPlotter.C** was some attempted version of the plotting...

