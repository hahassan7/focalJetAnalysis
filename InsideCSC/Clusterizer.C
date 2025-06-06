#include "Clusterizer.h"
#include "AliJHMREvent.h"

AliJHMREvent *event;

void Clusterizer(TString simFolder, TString clusteringOutputFileDir, Int_t jobstart, Int_t jobend)
{
    gSystem->Load("/alicesw/sw/slc7_x86-64/AliRoot/latest/lib/libpythia6.so");
    gSystem->Load("/alicesw/sw/slc7_x86-64/AliRoot/latest/lib/libAliPythia6.so");

    InitCombinedData(clusteringOutputFileDir);
    InitClusterizer();
    CreateClusterMaps();

    int nevtot = 0;
    int goodFiles = 0;
    for (Int_t ifolder = jobstart; ifolder <= jobend; ifolder++) {
        cout << "Processing folder " << Form("%03d", ifolder) << endl;

        TString inputFile = Form("%s/%03d/%s", simFolder.Data(), ifolder, "galice.root");
        int checkfileFocal = CheckFile(Form("%s/%03d/%s", simFolder.Data(), ifolder, "FOCAL.Hits.root"));
        int checkfileKine = CheckFile(Form("%s/%03d/%s", simFolder.Data(), ifolder, "Kinematics.root"));
        if (!checkfileFocal || !checkfileKine) continue;
        goodFiles++;

        TFile * outputFile = fFOCALCluster->CreateOutputFile(Form("%s/clusters_%s_%i.root",
                                                    clusteringOutputFileDir.Data(),
                                                    clusteringParametersTag.Data(),
                                                    ifolder));

        LoadClusterizerHit(inputFile);

        int nev = fRunLoader->GetNumberOfEvents();
        cout << "\tN EVENTS : " << nev << endl;
        for (int ievt = 0; ievt < nev; ievt++) {
            LoadEvent(simFolder, ifolder, ievt);
            GetKinematics(nevtot);
            Clusterize(ifolder, ievt);
            SaveCombinedOutput();
            nevtot++;
        } // end of event loop

        fFOCALCluster->CloseOutputFile();
    }
    outfile->cd();
    tree->Write("", TObject::kOverwrite);
    //    fFOCALCluster->AnaEnd();
    delete fFOCALCluster;
    /*
    cout << "Deleting maps" << endl;
    // causes crash; could try to narrow down which one?
    delete [] clusterMap;
    delete finalClusterMap;
    delete fineClusterMap;
    delete coarseClusterMap;
    */

    cout << "===========================DONE==========================" << endl;
    cout << goodFiles << "/" << jobend-jobstart+1 << " GOOD FILES AND CLUSTERIZED" << endl;
    cout << "=========================================================" << endl;

}

//******************************************************************************
//      FUNCTIONS
//******************************************************************************

int CheckFile(TString filename)
{
    if (gSystem->AccessPathName(filename.Data())) {
        cout << "!!!" << endl;
        cout << "FILE NOT FOUND : " << filename << endl;
        cout << "!!!" << endl;
        return 0;
    }

    TFile file(filename.Data());
    if (file.TestBit(TFile::kRecovered)) {
        cout << "!!!" << endl;
        cout << "FILE NOT WORKING : " << filename << endl;
        cout << "!!!" << endl;
        return 0;
    }
    return 1;
}

void InitCombinedData(TString clusteringOutputFileDir)
{
    outfile = TFile::Open(Form("%s/output.root", clusteringOutputFileDir.Data()), "RECREATE");
    tree = new TTree("Run", "MC and cluster data for FoCal");
    event = new AliJHMREvent();
    event->InitEvent();
    tree->Branch("Events", &event, 32000, 2);
}

void SaveCombinedOutput()
{
    tree->Fill();
    event->Clear();
}

void GetKinematics(int nevtot)
{
    TTree* treeH = fFOCALLoader->TreeH();
    TTree* treeK = fRunLoader->TreeK();
    AliStack *stack = fRunLoader->Stack();

    //int npart = stack->GetNtrack();
    int npart = stack->GetNprimary(); // Save only primaries
    for (int ipart = 0; ipart < npart; ipart++) {
        auto track = stack->Particle(ipart);
        int pdg = track->GetPdgCode();
        if (TDatabasePDG::Instance()->GetParticle(pdg)==NULL) {
            //std::cout << "Event " << nevtot << " : Particle with id " << pdg << " not found, skip" << std::endl;
            continue;
        }
        bool isPrimary = stack->IsPhysicalPrimary(ipart);
        Float_t charge = TDatabasePDG::Instance()->GetParticle(pdg)->Charge();
        Float_t px = track->Px();
        Float_t py = track->Py();
        Float_t pz = track->Pz();
        Float_t e = track->Energy();
        Int_t motherId = track->GetFirstMother();
        Int_t motherPdg = -1;
        if (motherId!=-1) {
            auto motherTrack = stack->Particle(motherId);
            motherPdg = motherTrack->GetPdgCode();
        }
        Float_t mass = TMath::Sqrt(e*e + track->P()*track->P());
        //if (pdg==22) cout << "PDG 22, isPrimary=" << isPrimary << ", mass=" << mass << ", mass(true)=" << track->GetMass() << ", status=" << track->GetStatusCode() << ", mother=" << motherPdg << ", isSecondaryFromWeakDecay=" << stack->IsSecondaryFromWeakDecay(ipart) << ", IsSecondaryFromMaterial=" << stack->IsSecondaryFromMaterial(ipart) << endl;
        event->AddTrack(pdg, px, py, pz, e, charge, ipart, motherId, motherPdg, isPrimary);
    }
}

void InitClusterizer()
{
    fFOCALCluster = new AliFOCALClusterizerv2();
    fFOCALCluster->InitGeometry(detectorGeometryFile);
    fFOCALCluster->InitParameters(clusterizerParametersFile);
    fFOCALCluster->SetDebugMode(clusteringDebugMode);

    //	  fFOCALCluster->SetOutput(outputFile);
    fFOCALCluster->SetLocalEnergyThreshold(eThreshold); //0.5MeV (1/2MIP)
    //	  fFOCALCluster->SetClusteringEnergyThreshold(0.004); //40MeV
    fFOCALCluster->SetDistance(dist);  // set the range you will look to make a cluster
    fFOCALCluster->SetPixelNoiseProb(pixelNoiseProb);
    fFOCALCluster->SetPadNoiseSigma(padNoiseSigma);

    // Set Logarithmic and power law weighting+parameters
    fFOCALCluster->SetMode(true, true);
    fFOCALCluster->SetW_0(logWeight, powWeight);
    fFOCALCluster->SetModeCorrectGW(false);

    geometry = fFOCALCluster->GetGeometry();
}

void CreateClusterMaps()
{
    // Create clusterMaps
    if (fFOCALCluster->GetDebugMode()) cout << "Creating cluster maps..." << endl;
    nSegments = geometry->GetVirtualNSegments();
    clusterMap = new ObjectMap* [nSegments];

    firstFineSeed   = -1;
    firstCoarseSeed = -1;

    for (Int_t segment = 0 ; segment < nSegments; segment++) {

        Bool_t isPixel = geometry->GetVirtualIsPixel(segment);
        Int_t nCols = (Int_t)(geometry->GetFOCALSizeX() / geometry->GetVirtualPadSize(segment));
        Int_t nRows = (Int_t)(geometry->GetFOCALSizeY() / geometry->GetVirtualPadSize(segment));

        clusterMap[segment] = new ObjectMap(nCols,nRows);

        if (!coarseClusterMap && !isPixel) {
            cout << "\t\tCOARSE" << endl;
            coarseClusterMap = new ObjectMap(nCols,nRows);
        }

        if (!fineClusterMap && isPixel) {
            cout << "\t\tFINE" << endl;
            fineClusterMap = new ObjectMap(nCols,nRows);
        }

        if (!finalClusterMap && isPixel) {
            finalClusterMap = new ObjectMap(nCols,nRows);
        }
    }
}

void LoadClusterizerHit(TString inputfile)
{
    if (fRunLoader) fRunLoader->Delete();
    fRunLoader = AliRunLoader::Open(inputfile.Data(), "READ");

    if (!fRunLoader) {
        cout << "ERROR : RunLoader not found for " << inputfile.Data() << endl;
    }


    if (!fRunLoader->GetAliRun()) fRunLoader->LoadgAlice();
    if (!fRunLoader->TreeE()) fRunLoader->LoadHeader();
    if (!fRunLoader->TreeK()) fRunLoader->LoadKinematics();

    gAlice = fRunLoader->GetAliRun();

    fFOCAL  = (AliFOCAL*)gAlice->GetDetector("FOCAL");
    fFOCALLoader = dynamic_cast<AliFOCALLoader*>(fRunLoader->GetLoader("FOCALLoader"));
    fFOCALLoader->LoadHits("READ");
}

void LoadEvent(TString simFolder, int ifolder, int ievt)
{
    actualDir = fFOCALCluster->SetOutputForEvent(Form("Event%i",ievt));

    Int_t ie = fRunLoader->GetEvent(ievt);
    if (ie) return;

    TTree* treeH = fFOCALLoader->TreeH();
    cout << "\tEvent: " << ievt << " with " << treeH->GetEntries() << " tracks" << endl;

    fFOCALCluster->SetInput(treeH, "HITS");
    fFOCALCluster->LoadHits2Cells(); // Creates fCells list

    actualDir->cd();

    TObjString * sPar;
    TParameter<Int_t> * iPar;

    sPar = new TObjString(Form("%s/%i",simFolder.Data(),ifolder));
    sPar->Write("SimFolder");

    sPar = new TObjString(simFolder.Data());
    sPar->Write("GeneralSimFolder");

    iPar = new TParameter<Int_t>("FolderNumber",ifolder);
    iPar->Write("FolderNumber");

    iPar = new TParameter<Int_t>("EventNumber",ievt);
    iPar->Write("EventNumber");
}

void Clusterize(int ifolder, int ievt)
{
    // Prepare for clustering
    fFOCALCluster->ClearClusterArrays();
    fFOCALCluster->SortCells();
    fFOCALCluster->CreateAndFillCellMap();
    fFOCALCluster->ResetSeeds();

    // Cluster individual segments - except the 4th - it needs external seeds
    Bool_t a[7] = {true,true,true,false,false,false,false};
    fFOCALCluster->MakeClustersPerSegment(a, 0);

    // Get preseeds from 2nd segment
    TObjArray * preSeeds = new TObjArray(10);
    TClonesArray * subClusters = fFOCALCluster->GetClusterItr();
    Int_t nClusters = subClusters->GetEntries();
    for (Int_t i = 0; i < nClusters; i++) {

        AliFOCALCluster * cluster = (AliFOCALCluster*) subClusters->UncheckedAt(i);

        if (cluster->Segment() != FINE_SEEDS[0]) continue;

        TVectorF * preSeed = new TVectorF(4);
        (*preSeed)[0] = cluster->X()/cluster->Z();
        (*preSeed)[1] = 0;
        (*preSeed)[2] = cluster->Y()/cluster->Z();
        (*preSeed)[3] = 0;

        preSeeds->AddLast(preSeed);
    }

    Bool_t a2[6] = {false,false,false,true,false,false};
    fFOCALCluster->MakeClustersPerSegment(a2, preSeeds);
    preSeeds->Delete();
    delete preSeeds;

    // Get preseeds from 3rd segment
    preSeeds = new TObjArray(10);
    nClusters = subClusters->GetEntries();
    for (Int_t i = 0; i < nClusters; i++) {

        AliFOCALCluster * cluster = (AliFOCALCluster*) subClusters->UncheckedAt(i);

        if (cluster->Segment() != COARSE_SEEDS[0]) continue;

        TVectorF * preSeed = new TVectorF(4);
        (*preSeed)[0] = cluster->X()/cluster->Z();
        (*preSeed)[1] = 0;
        (*preSeed)[2] = cluster->Y()/cluster->Z();
        (*preSeed)[3] = 0;

        preSeeds->AddLast(preSeed);
    }

    Bool_t a3[7] = {false,false,false,false,true,true,false};
    fFOCALCluster->MakeClustersPerSegment(a3,preSeeds);
    preSeeds->Delete();
    delete preSeeds;

    TClonesArray * clusters = fFOCALCluster->GetCluster();

    // #### Here, the clusterizer has done all it can and it is up to this script to sum-up
    // #### the individual sub-clusters to clusters
    // Sum-up sub-clusters from individual segments to clusters

    // START Combine algorithm

    // First create and initializer maps for fast-searching of clusters
    AliFOCALRinger ringer;
    ringer.Init(100);

    coarseClusterMap->ResetMap();
    fineClusterMap->ResetMap();
    //finalClusterMap->ResetMap(); // Not used here...

    Int_t segIndex [AliFOCALCluster::N_DEPTH_FIELDS];
    Float_t width1 [AliFOCALCluster::N_DEPTH_FIELDS];
    Float_t width2 [AliFOCALCluster::N_DEPTH_FIELDS];
    Float_t phi [AliFOCALCluster::N_DEPTH_FIELDS];
    Float_t segmentEnergy [AliFOCALCluster::N_DEPTH_FIELDS];
    Float_t seedE[AliFOCALCluster::N_DEPTH_FIELDS];
    Int_t nCells[AliFOCALCluster::N_DEPTH_FIELDS];

    for (Int_t segment = 0 ; segment < nSegments; segment++) {
        clusterMap[segment]->ResetMap();
    }

    Int_t nSubClusters = subClusters->GetEntries();

    if (fFOCALCluster->GetDebugMode()) cout << "Registering clusters into maps..." << endl;

    // Now register all sub-clusters in the map and reset their flag
    for (Int_t c = 0; c < nSubClusters; c++) {

        AliFOCALCluster * subCluster = (AliFOCALCluster*) subClusters->UncheckedAt(c);
        subCluster->SetFlag(true);

        Int_t segment = subCluster->Segment();

        if ((segment < 0) || (segment >=nSegments)) continue;

        Int_t col = (Int_t)((subCluster->X() + geometry->GetFOCALSizeX()/2)/geometry->GetVirtualPadSize(segment));
        Int_t row = (Int_t)((subCluster->Y() + geometry->GetFOCALSizeY()/2)/geometry->GetVirtualPadSize(segment));

        if (!clusterMap[segment]->Get(col,row)) {
            clusterMap[segment]->InsertAt(col,row,subCluster);
        } else {
            Bool_t found = false;
            ringer.SetRing(1);
            Int_t xCol,yCol;
            while (ringer.GetNext(xCol,yCol)) {
                if (!clusterMap[segment]->Get(col+xCol,row+yCol)) {
                    clusterMap[segment]->InsertAt(col+xCol,row+yCol,subCluster);
                    found = true;
                    break;
                }
            }
            if (!found && fFOCALCluster->GetDebugMode())
                cout << Form("WARNING sub-cluster at [%.2f,%.2f] not registered in the map, no space left in vicinity!",
            subCluster->X(),subCluster->Y());
        }
    }
    fFOCALCluster->MakeHCALClusters();

    // Now start combining clusters
    if (fFOCALCluster->GetDebugMode()) cout << "Combining clusters..." << endl;

    //Float_t finalZ = geometry->GetFOCALZ0();
    Float_t finalZ = geometry->GetVirtualSegmentZ(2);

    // For-loop over all seed segments
    for (Int_t seedSegment = 0; seedSegment < sizeof(ALL_SEEDS)/sizeof(ALL_SEEDS[0]); seedSegment++) {

        Int_t ALL_SEED = ALL_SEEDS[seedSegment];

        Bool_t isPixelSegment = geometry->GetVirtualIsPixel(ALL_SEED);
        Int_t finalSegment;
        Int_t calibSegment;
        ObjectMap * map = 0;
        if (!isPixelSegment) {
            finalSegment = -2;
            if (fFOCALCluster->GetDebugMode()) cout << "COARSE CLUSTERS:" << endl;
            map = coarseClusterMap;
            calibSegment = COARSE_CALIB;
            if (firstCoarseSeed == -1) firstCoarseSeed = ALL_SEED;
        } else {
            finalSegment = -3;
            if (fFOCALCluster->GetDebugMode()) cout << "FINE CLUSTERS:" << endl;
            map = fineClusterMap;
            calibSegment = FINE_CALIB;
            if (firstFineSeed == -1) firstFineSeed = ALL_SEED;
        }

        if (fFOCALCluster->GetDebugMode()) {
            Float_t * pars = fFOCALCluster->GetCalibrationPars(calibSegment);
            cout << Form("Calibration parameters for segment %i: p0=%.3e, p1=%.3e, p2=%.3e, p3=%.3e.",ALL_SEED,pars[0],pars[1],pars[2],pars[3]) << endl;
            delete pars;
        }

        for (Int_t c = 0; c < nSubClusters; c++) {

            AliFOCALCluster * subCluster = (AliFOCALCluster*) subClusters->UncheckedAt(c);

            // Skip already merged clusters
            if (!subCluster->Flag()) continue;

            Int_t segment = subCluster->Segment();

            if (segment != ALL_SEED)
            continue;

            // reset info fields
            for (Int_t i = 0; i < AliFOCALCluster::N_DEPTH_FIELDS; i++) {

                segIndex[i] = -1;
                width1[i] = -1;
                width2[i] = -1;
                phi[i] = - 2.0;
                segmentEnergy[i] = 0;
                seedE[i] = 0;
                nCells[i] = 0;
            }

            Float_t maxDistance = fFOCALCluster->GetMinRing(segment,0)*geometry->GetVirtualPadSize(segment);
            Float_t z1 = subCluster->Z();
            Float_t x1 = subCluster->X()*finalZ/z1;
            Float_t y1 = subCluster->Y()*finalZ/z1;

            // resulting cluster fields
            Float_t totalWeight = ENERGY_WEIGHTS[segment];
            Float_t energy = subCluster->E()*ENERGY_WEIGHTS[segment];
            Float_t x = x1*energy;
            Float_t y = y1*energy;

            if (fFOCALCluster->GetDebugMode())
            cout << Form("\t==========\n\tMerging: Seg=%i\tE = %.2f\tx = %.2f\ty = %.2f",segment,energy, x1, y1) << endl;

            if (segment < AliFOCALCluster::N_DEPTH_FIELDS) {
                segIndex[segment] = segment;
                width1[segment] = subCluster->GetWidth1(0);
                width2[segment] = subCluster->GetWidth2(0);
                phi[segment] = subCluster->GetPhi(0);
                segmentEnergy[segment] = subCluster->E();
                seedE[segment] = subCluster->GetSeedEnergy(0);
                nCells[segment] += subCluster->GetNcells(0);
                if (fFOCALCluster->GetDebugMode())
                cout << Form("\t\tAdding info: Segment=%i\tWidth1=%.2f\tWidth2=%.2f\tLongE=%.2f",
                segment, subCluster->GetWidth1(0), subCluster->GetWidth2(0), subCluster->E()) << endl;
            }

            // Set flag as the cluster to false, not to merge it again
            subCluster->SetFlag(false);

            // Now walk through neighbourhood in other coarse segments
            for (Int_t segment2 = 0; segment2 < nSegments; segment2++) {

                // Skip segments with of diffetent type (regarding pixels/pads)
                if (geometry->GetVirtualIsPixel(segment2) != isPixelSegment)
                continue;

                // skip seed segment
                if (segment2 == segment)
                continue;

                // Whether multiple clusters are merged from the same segment
                Bool_t merged = false;
                Bool_t multipleMerges = false;

                // Get seed col,row position for segment2
                Int_t col = (Int_t)((x1*geometry->GetVirtualSegmentZ(segment2)/finalZ+geometry->GetFOCALSizeX()/2)/geometry->GetVirtualPadSize(segment2));
                Int_t row = (Int_t)((y1*geometry->GetVirtualSegmentZ(segment2)/finalZ+geometry->GetFOCALSizeY()/2)/geometry->GetVirtualPadSize(segment2));

                Int_t xCol,yCol;
                for (Int_t ring = 0; ring <= 3; ring++){

                    ringer.SetRing(ring);
                    while(ringer.GetNext(xCol,yCol)) {
                        AliFOCALCluster * subCluster2 = (AliFOCALCluster*) clusterMap[segment2]->Get(col+xCol,row+yCol);

                        if (!subCluster2)
                        continue;

                        // Skip already merged clusters
                        if (!subCluster2->Flag())
                        continue;

                        // Get positions of sub-cluster2 recalculated for seed segment
                        Float_t z2 = subCluster2->Z();
                        Float_t x2 = subCluster2->X() * finalZ / z2;
                        Float_t y2 = subCluster2->Y() * finalZ / z2;
                        Float_t e2 = subCluster2->E();

                        // Merge the clusters, if they are withing dist
                        Float_t clusterDistance = TMath::Sqrt((TMath::Power((x1-x2),2) + TMath::Power((y1-y2),2)));
                        if (clusterDistance < maxDistance) {

                            if (merged && !multipleMerges)
                            multipleMerges = true;

                            if (!merged)
                            merged = true;

                            energy += e2*ENERGY_WEIGHTS[segment2];
                            totalWeight += ENERGY_WEIGHTS[segment2];
                            x += x2*e2*ENERGY_WEIGHTS[segment2];
                            y += y2*e2*ENERGY_WEIGHTS[segment2];
                            subCluster2->SetFlag(false); // So the clusters is not merged again

                            // Add info
                            if (segment2 < AliFOCALCluster::N_DEPTH_FIELDS) {
                                segIndex[segment2] = segment2;
                                width1[segment2] = subCluster2->GetWidth1(0);
                                width2[segment2] = subCluster2->GetWidth2(0);
                                phi[segment2] = subCluster2->GetPhi(0);
                                segmentEnergy[segment2] += subCluster2->E();
                                seedE[segment2] = TMath::Max(subCluster2->GetSeedEnergy(0),seedE[segment2]);
                                nCells[segment2] += subCluster2->GetNcells(0);
                            }

                            if (fFOCALCluster->GetDebugMode()) {
                                cout << Form("\t\t\tWith: Seg=%i\tE = %.2f\tx = %.2f\ty = %.2f",segment2, e2, x2, y2) << endl;
                                cout << Form("\t\t\t\tAdding info: Segment=%i\tWidth1=%.2f\tWidth2=%.2f\tLongE=%.2f",
                                segment2, subCluster2->GetWidth1(0), subCluster2->GetWidth2(0), subCluster2->E()) << endl;
                            }
                        } // end dist if
                    } // end ringer while
                } // end for rings

                // If there were multiple merges, the width paramteres are not correct
                if (multipleMerges) {
                    width1[segment2] = -2;
                    width2[segment2] = -2;
                    phi[segment2] = -3.0;
                }

            } // end segment for-loop

            // Calibrate the resulting
            Float_t * pars = fFOCALCluster->GetCalibrationPars(calibSegment);
            Float_t calibratedEnergy = Calibration(energy,calibSegment,pars,clusteringCalibrationRun,fFOCALCluster->GetDebugMode());
            delete pars;

            if (calibratedEnergy <= 0) {
                if (fFOCALCluster->GetDebugMode())
                cout << "\t\tCLUSTER_REJECTED (ZERO ENERGY): " << Form("Seg=%i\tE=%f\tx=%f\ty=%f", finalSegment, calibratedEnergy, x/energy, y/energy) << endl;
                continue;
            }

            // Create new cluster
            new((*subClusters)[subClusters->GetEntries()]) AliFOCALCluster(x/energy, y/energy, finalZ, calibratedEnergy, finalSegment);
            if (fFOCALCluster->GetDebugMode())
            cout << "\t\tCLUSTER_FOUND: " << Form("Seg=%i\tE=%f\tx=%f\ty=%f", finalSegment, calibratedEnergy, x/energy, y/energy) << endl;

            AliFOCALCluster * newCluster = (AliFOCALCluster*) subClusters->Last();

            // Set info for new cluster
            for (Int_t i = 0; i < AliFOCALCluster::N_DEPTH_FIELDS; i++) {

                newCluster->SetSegIndex(i,segIndex[i]);
                newCluster->SetWidth1(i,width1[i]);
                newCluster->SetWidth2(i,width2[i]);
                newCluster->SetPhi(i,phi[i]);
                newCluster->SetSegmentEnergy(i,segmentEnergy[i]);
                newCluster->SetSeedEnergy(i,seedE[i]);
                newCluster->SetNcells(i,nCells[i]);
            }

            Int_t col = (Int_t)((x/energy+geometry->GetFOCALSizeX()/2)/geometry->GetVirtualPadSize(ALL_SEED));
            Int_t row = (Int_t)((y/energy+geometry->GetFOCALSizeY()/2)/geometry->GetVirtualPadSize(ALL_SEED));
            // Register the cluster in the map
            Int_t xCol,yCol;
            if (!map->Get(col,row)) {
                if (!map->InsertAt(col,row,subClusters->Last()))
                cout << "Error inserting cluster in map at col " << col << " row " << row << endl;
                // keep info so that final cluster is registered at same coordinates
                xCol = col;
                yCol = row;
            } else {
                Bool_t found = false;
                ringer.SetRing(1);
                while(ringer.GetNext(xCol,yCol)) {
                    if (!map->Get(col+xCol,row+yCol)) {
                        map->InsertAt(col+xCol,row+yCol,subClusters->Last());
                        found = true;
                        // keep info so that final cluster is registered at same coordinates
                        xCol = col+xCol;
                        yCol = row+yCol;
                        break;
                    }
                }
                if (!found && fFOCALCluster->GetDebugMode())
                cout << Form("WARNING cluster at [%.2f,%.2f] of segment %i not registered in the map, no space left in vicinity!",
                x/energy,y/energy,finalSegment);
            } // end else register

            /*              NOT USED HERE
            // If it is a fine cluster, create also a final cluster and register it to the same spot
            if (isPixelSegment) {
            new((*clusters)[clusters->GetEntries()]) AliFOCALCluster(x/energy, y/energy, finalZ, 0, -1);
            finalClusterMap->InsertAt(xCol,yCol,clusters->Last());
        }
        */

        } // end loop sub-clusters
    } // end loop seed segments

    // Now, for each found cluster in fine layers, try to find closest coarse one
    if (fFOCALCluster->GetDebugMode()) cout << "FINAL CLUSTERS:" << endl;
    Int_t nSubClusters2 = subClusters->GetEntries();
    TObjArray * closestCoarseClusters = new TObjArray(nSubClusters2-nSubClusters);
    if (fFOCALCluster->GetDebugMode()) cout << "  ASSIGNING WEIGHTS:" << endl;

    // First assing weights to the nearest coarse cluster
    for (Int_t c = nSubClusters; c < nSubClusters2; c++) {

        AliFOCALCluster * subCluster = (AliFOCALCluster*) subClusters->UncheckedAt(c);

        // We want only clusters from fine layers
        if (subCluster->Segment() != -3) continue;

        Float_t x1 = subCluster->X();
        Float_t y1 = subCluster->Y();
        Float_t z1 = subCluster->Z();
        Float_t e1 = subCluster->E();

        if (fFOCALCluster->GetDebugMode())
            cout << Form("\t==========\n\tAssigning Weight: Seg=%i\tE = %.2f\tx = %.2f\ty = %.2f",subCluster->Segment(),e1, x1, y1) << endl;

        // corresponding col and row for fine segment
        Int_t col = (Int_t)((x1+geometry->GetFOCALSizeX()/2)/geometry->GetVirtualPadSize(firstCoarseSeed));
        Int_t row = (Int_t)((y1+geometry->GetFOCALSizeY()/2)/geometry->GetVirtualPadSize(firstCoarseSeed));

        // Now search neighbourhood for fine clusters
        Int_t xCol,yCol;

        // Set max dist : MaxRing in Coarse segments
        Int_t maxRing = fFOCALCluster->GetMaxRing(firstCoarseSeed,0);
        AliFOCALCluster * closest = 0;
        Float_t closestDistance = -1;

        for (Int_t ring = 0; ring <= maxRing; ring++){

            // This is to be sure that found cluster is really closest - we need to check also the next ring
            Bool_t wasFoundPreviously = false;
            if (closest) wasFoundPreviously = true;

            ringer.SetRing(ring);
            while (ringer.GetNext(xCol,yCol)) {
                AliFOCALCluster * subCluster2 = (AliFOCALCluster*) coarseClusterMap->Get(col+xCol,row+yCol);

                // skip empty places
                if (!subCluster2) continue;

                // Get positions of sub-cluster2 recalculated for seed segment
                Float_t z2 = subCluster2->Z();
                Float_t x2 = subCluster2->X();
                Float_t y2 = subCluster2->Y();
                Float_t e2 = subCluster2->E();

                if (!closest) {
                    closest = subCluster2;
                    closestDistance = TMath::Sqrt((TMath::Power((x1-x2),2) + TMath::Power((y1-y2),2)));
                } else {
                    Float_t actualDistance = TMath::Sqrt((TMath::Power((x1-x2),2) + TMath::Power((y1-y2),2)));
                    if (actualDistance < closestDistance) {
                        closest = subCluster2;
                        closestDistance = actualDistance;
                        if (wasFoundPreviously)
                        wasFoundPreviously = false;
                    }
                }

            } // end ringer while
            if (wasFoundPreviously) break;
        } // end for rings

        // Only if nearby coarse cluster was found, we discard any fine clusters that do not fulfill this condition
        if (closest) {
            if (fFOCALCluster->GetDebugMode())
                cout << Form("\t\t\tTo: Seg=%i\tE = %.2f\tx = %.2f\ty = %.2f",-2, closest->E(), closest->X(), closest->Y()) << endl;
            // Found closest cluster, register it in an object array
            closestCoarseClusters->AddAt(closest,c-nSubClusters);
            // and add the weight (= energy of the fine cluster) to it
            closest->AddWeight(e1);
        }
    } // end for assing weights

    if (fFOCALCluster->GetDebugMode()) {
        cout << "  CREATING NEW CLUSTERS FROM FINE ONES:" << endl;
    }
    // Second: Create final clusters from fine ones and give energy to them
    for (Int_t c = nSubClusters; c < nSubClusters2; c++) {

        AliFOCALCluster * subCluster2 = (AliFOCALCluster*) closestCoarseClusters->UncheckedAt(c-nSubClusters);

        // Skip if no coarse cluster was found for this one (should also skip all coarse clusters)
        if (!subCluster2) continue;

        AliFOCALCluster * subCluster = (AliFOCALCluster*) subClusters->UncheckedAt(c);

        // We want only clusters from fine layers (should not be necessary, but...)
        if (subCluster->Segment() != -3) {
            cout << "Wrong cluster segment ? : " << subCluster->Segment() << endl;
            continue;
        }

        Float_t x1 = subCluster->X();
        Float_t y1 = subCluster->Y();
        Float_t z1 = subCluster->Z();
        Float_t e1 = subCluster->E();

        // Calculate energy of the new cluster
        Float_t newEnergy = subCluster2->E() * (e1/subCluster2->Weight());
        // Create the new cluster
        new ((*clusters)[clusters->GetEntries()]) AliFOCALCluster(x1, y1, finalZ, newEnergy, -1);

        // report
        if (fFOCALCluster->GetDebugMode())
            cout << "\t\tSEMIFINAL_CLUSTER_FOUND_FROMFINE: " << Form("Seg=%i\tE=%f\tx=%f\ty=%f", -1, newEnergy, x1, y1) << endl;

        AliFOCALCluster * newCluster = (AliFOCALCluster*) clusters->Last();
        // Add info
        for (Int_t i = 0; i < AliFOCALCluster::N_DEPTH_FIELDS; i++) {

            if (subCluster->GetSegIndex(i) >= 0) {
                newCluster->SetSegIndex(i,subCluster->GetSegIndex(i));
                newCluster->SetWidth1(i,subCluster->GetWidth1(i));
                newCluster->SetWidth2(i,subCluster->GetWidth2(i));
                newCluster->SetPhi(i,subCluster->GetPhi(i));
                newCluster->SetSegmentEnergy(i,subCluster->GetSegmentEnergy(i));
                newCluster->SetSeedEnergy(i,subCluster->GetSeedEnergy(i));
                newCluster->SetNcells(i,subCluster->GetNcells(i));
                if (fFOCALCluster->GetDebugMode())
                    cout << Form("\t\t\tAdding info: Segment=%i\tWidth1=%.2f\tWidth2=%.2f\tLongE=%.2f",
                newCluster->GetSegIndex(i),newCluster->GetWidth1(i),
                newCluster->GetWidth2(i),newCluster->GetSegmentEnergy(i)) << endl;
            } else if (subCluster2->GetSegIndex(i) >= 0) {
                newCluster->SetSegIndex(i,subCluster2->GetSegIndex(i));
                newCluster->SetWidth1(i,subCluster2->GetWidth1(i));
                newCluster->SetWidth2(i,subCluster2->GetWidth2(i));
                newCluster->SetPhi(i,subCluster2->GetPhi(i));
                newCluster->SetSegmentEnergy(i,subCluster2->GetSegmentEnergy(i)*(e1/subCluster2->Weight()));
                newCluster->SetSeedEnergy(i,subCluster2->GetSeedEnergy(i));
                newCluster->SetNcells(i,subCluster2->GetNcells(i));
                if (fFOCALCluster->GetDebugMode())
                    cout << Form("\t\t\tAdding info (cl2): Segment=%i\tWidth1=%.2f\tWidth2=%.2f\tLongE=%.2f/%.2f",
                newCluster->GetSegIndex(i),newCluster->GetWidth1(i),
                newCluster->GetWidth2(i),newCluster->GetSegmentEnergy(i),
                subCluster->GetSegmentEnergy(i)) << endl;
            }
        }
        // Add HCAL info
        // Need to extrapolate to proper position?
        //cout << "Adding HCAL info; using x " << x1 << " y " << y1 << " z " << z1 << endl;
        Float_t eHCAL = fFOCALCluster->GetHCALEnergy(x1,y1);
        newCluster->SetHCALEnergy(eHCAL);
        eHCAL = fFOCALCluster->GetHCALIsoEnergy(x1,y1,0.2);
        newCluster->SetIsoEnergyR2(eHCAL);
        eHCAL = fFOCALCluster->GetHCALIsoEnergy(x1,y1,0.4);
        newCluster->SetIsoEnergyR4(eHCAL);
    }

    if (fFOCALCluster->GetDebugMode())
        cout << "  CREATING NEW CLUSTERS FROM LONELY COARSE ONES:" << endl;

    // Third: find coarse clusters with no fine one nearby and create new final cluster for each
    for (Int_t c = nSubClusters; c < nSubClusters2; c++) {

        AliFOCALCluster * subCluster = (AliFOCALCluster*) subClusters->UncheckedAt(c);

        // We want only clusters from coarse layers
        if (subCluster->Segment() != -2) continue;

        // We want only lonely coarse clusters now
        if (subCluster->Weight() != 0) continue;

        Float_t x1 = subCluster->X();
        Float_t y1 = subCluster->Y();
        Float_t z1 = subCluster->Z();
        Float_t e1 = subCluster->E();

        //create new cluster
        new ((*clusters)[clusters->GetEntries()]) AliFOCALCluster(x1, y1, finalZ, e1, -1);

        AliFOCALCluster * newCluster = (AliFOCALCluster*) clusters->Last();

        // Set info for new cluster
        for (Int_t i = 0; i < AliFOCALCluster::N_DEPTH_FIELDS; i++) {
            newCluster->SetSegIndex(i,subCluster->GetSegIndex(i));
            newCluster->SetWidth1(i,subCluster->GetWidth1(i));
            newCluster->SetWidth2(i,subCluster->GetWidth2(i));
            newCluster->SetPhi(i,subCluster->GetPhi(i));
            newCluster->SetSegmentEnergy(i,subCluster->GetSegmentEnergy(i));
            newCluster->SetSeedEnergy(i,subCluster->GetSeedEnergy(i));
            newCluster->SetNcells(i,subCluster->GetNcells(i));
        }

        // Add HCAL info
        // Need to extrapolate to proper position?
        //cout << "Adding HCAL info; using x " << x1 << " y " << y1 << " z " << z1 << endl;
        Float_t eHCAL = fFOCALCluster->GetHCALEnergy(x1,y1);
        newCluster->SetHCALEnergy(eHCAL);
        eHCAL = fFOCALCluster->GetHCALIsoEnergy(x1,y1,0.2);
        newCluster->SetIsoEnergyR2(eHCAL);
        eHCAL = fFOCALCluster->GetHCALIsoEnergy(x1,y1,0.4);
        newCluster->SetIsoEnergyR4(eHCAL);

        // report
        if (fFOCALCluster->GetDebugMode())
            cout << "\t\tSEMIFINAL_CLUSTER_FOUND_FROMCOARSE: " << Form("Seg=%i\tE=%f\tx=%f\ty=%f", -1, e1, x1, y1) << endl;

    } // end for final clusters

    if (fFOCALCluster->GetDebugMode())
        cout << "  REMOVE CLUSTERS WITH 0 ENERGY AND CREATE FINAL LIST" << endl;
    // Remove finalClusters 0 energy
    nClusters = clusters->GetEntries();
    for (Int_t c = 0; c < nClusters; c++) {

        AliFOCALCluster * cluster = (AliFOCALCluster*) clusters->UncheckedAt(c);
        if (cluster->E() == 0) {
            if (fFOCALCluster->GetDebugMode())
                cout << Form("\tRemoving final cluster [%.2f,%.2f] with 0 energy.",cluster->X(),cluster->Y()) << endl;
            clusters->Remove(cluster);
            continue;
        }

        if (cluster->E() < 0) {
            if (fFOCALCluster->GetDebugMode())
                cout << Form("\tRemoving final cluster [%.2f,%.2f] with negative (%.2e) energy.",
            cluster->X(),cluster->Y(),cluster->E()) << endl;
            clusters->Remove(cluster);
            continue;
        }
        if (fFOCALCluster->GetDebugMode())
            cout << Form("\tFINAL_CLUSTER_FOUND: Seg=%i\tE=%f\tx=%f\ty=%f", -1, cluster->E(), cluster->X(),cluster->Y()) << endl;

        double clustX = cluster->X();
        double clustY = cluster->Y();
        double clustZ = cluster->Z();
        double clustE = cluster->E();
        double clustEHCAL = cluster->GetHCALEnergy();
        float segE[AliFOCALCluster::N_DEPTH_FIELDS];
        float seedE[AliFOCALCluster::N_DEPTH_FIELDS];
        float width1[AliFOCALCluster::N_DEPTH_FIELDS];
        float width2[AliFOCALCluster::N_DEPTH_FIELDS];
        for (int i=0; i<AliFOCALCluster::N_DEPTH_FIELDS; i++) {
            segE[i] = cluster->GetSegmentEnergy(i);
            seedE[i] = cluster->GetSeedEnergy(i);
            width1[i] = cluster->GetWidth1(i);
            width2[i] = cluster->GetWidth2(i);
        }
        event->AddCluster(clustX, clustY, clustZ, clustE, clustEHCAL, segE, seedE, width1, width2);
    }
    clusters->Compress();

    // Clear up
    delete closestCoarseClusters;

    // END Combine algorithm

    // Clean up
    fFOCALCluster->DeleteCellMap();

    // Save clusters to output ttree
    fFOCALCluster->FillClusterTree();
    fFOCALCluster->SaveEventOutput(Form("Event%i",ievt));
    if (saveFigures)
        fFOCALCluster->SaveImage(Form("figs/event_%02i", ievt));
}

Float_t Calibration(Float_t signal, Int_t segment, Float_t * pars, bool calibrationRun, bool debug) {

    if (calibrationRun) return signal;

    //float result = TMath::Power((signal - pars[3])/pars[0],1./pars[2]) - pars[1];
    /*
    Float_t powerBase = signal+pars[1];

    Float_t result = 0;

    if (powerBase < 0) {
    if (debug)
    cout << Form("\t\tCalibration in segment %i: Negative power base (%.2e), returning 0",segment,powerBase) << endl;
    result = 0;
    }
    else
    result = pars[0]*TMath::Power(powerBase,pars[2]) + pars[3];
    */
    Float_t result = (signal - pars[1]*(signal/pars[0]*signal/pars[0]))/pars[0];

    if (result < 0) {
        if (debug)
            cout << Form("\t\tCalibration in segment %i: Negative value (%.2e), returning 0",segment,result) << endl;
        result = 0;
    }

    if (debug)
        cout << Form("\t\tCalibration in segment %i from %.2e to %.2e GeV",segment,signal,result) << endl;

    return result;
}
