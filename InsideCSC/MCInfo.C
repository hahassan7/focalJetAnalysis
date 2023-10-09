void MCInfo()
{

	// Parameteres for getting weighting function info
	Float_t infoParameterSigma = 1.0;
	Float_t infoParameterFitRadius = 5;
	Float_t infoParameterHistRadius = 5;

	const Float_t profileRadius = 5; // [cm] distance for collecting digits for shower profile

	Long_t totalEvents = 0;

	// Get geometry;
	AliFOCALGeometry *geometry = new AliFOCALGeometry();
	const Char_t *detectorGeometryFile = gSystem->ExpandPathName("$FOCAL/geometryFiles/geometry_Spaghetti.txt");
	geometry->Init(detectorGeometryFile);
	Int_t nSegments = geometry->GetVirtualNSegments();

	Float_t z_FOCAL_front = geometry->GetFOCALZ0() - geometry->GetFOCALSizeZ() / 2 - 0.1;
	// Create digitizer
	AliFOCALDigitizer *digitizer = new AliFOCALDigitizer();
	digitizer->SetGeometry(geometry);
	digitizer->Initialize("HD");

	// Event info
	Char_t folder[200];
	Int_t event;

	// Incident particle info
	Int_t pdgCode[4];
	Float_t e[4];
	Float_t pt[4];
	Float_t phi[4];
	Float_t theta[4];
	Float_t Vx[4], Vy[4], Vz[4];
	Float_t cVx[4], cVy[4], cVz[4];
	Int_t conv_flag[4];

	// Fit info
	Float_t *totalE = new Float_t[2 * nSegments];
	Float_t *maxE = new Float_t[2 * nSegments];
	Float_t *sigma = new Float_t[2 * nSegments];
	Float_t *ampl = new Float_t[2 * nSegments];
	Float_t *chi2 = new Float_t[2 * nSegments];
	Int_t *ndef = new Int_t[2 * nSegments];

	// profile histos
	TH2F **showProf2DEvent = new TH2F *[nSegments];
	TH1F **showProf1DEvent = new TH1F *[nSegments];
	for (Int_t i = 0; i < nSegments; i++)
	{
		Int_t nBins = 5;
		Float_t cellSize = 1;
		if (i == 1 || i == 3)
		{
			nBins = 50;
			cellSize = 0.1;
		}
		Float_t lim = 5 + 0.5 * cellSize;
		showProf2DEvent[i] = new TH2F(Form("prof2DEvent_%d", i), Form("Shower profile seg %d;x (cm);y (cm)", i), 2 * nBins + 1, -lim, lim, 2 * nBins + 1, -lim, lim);
		showProf1DEvent[i] = new TH1F(Form("profEvent_%d", i), Form("Shower profile seg %d;r (cm);dN/rdr", i), nBins + 1, -0.5 * cellSize, lim);
	}

	// Alice run loader
	AliRunLoader *fRunLoader = AliRunLoader::Open("galice.root");

	if (!fRunLoader->GetAliRun())
		fRunLoader->LoadgAlice();
	if (!fRunLoader->TreeE())
		fRunLoader->LoadHeader();
	if (!fRunLoader->TreeK())
		fRunLoader->LoadKinematics();

	gAlice = fRunLoader->GetAliRun();

	// Focal init
	AliFOCAL *fFOCAL = (AliFOCAL *)gAlice->GetDetector("FOCAL");
	AliFOCALLoader *fFOCALLoader = dynamic_cast<AliFOCALLoader *>(fRunLoader->GetLoader("FOCALLoader"));
	fFOCALLoader->LoadHits("READ");

	// Output ttree definition
	TFile *f = new TFile("MCInfoResults.root", "RECREATE");
	f->cd();
	TTree *tInfo = new TTree("MCInfo", "MonteCarlo and other information");

	tInfo->Branch("Event", &event, "Event/I");

	// Particles, mother + daughters
	tInfo->Branch("PdgCode", pdgCode, "PdgCode[4]/I");
	tInfo->Branch("Energy", e, "Energy[4]/F");
	tInfo->Branch("Pt", pt, "Pt[4]/F");
	tInfo->Branch("Phi", phi, "Phi[4]/F");
	tInfo->Branch("Theta", theta, "Theta[4]/F");

	tInfo->Branch("Vx", Vx, "Vx[4]/F");
	tInfo->Branch("Vy", Vy, "Vy[4]/F");
	tInfo->Branch("Vz", Vz, "Vz[4]/F");
	tInfo->Branch("conv_flag", conv_flag, "conv_flag[4]/I");
	tInfo->Branch("cVx", cVx, "cVx[4]/F");
	tInfo->Branch("cVy", cVy, "cVy[4]/F");
	tInfo->Branch("cVz", cVz, "cVz[4]/F");

	tInfo->Branch("TotalEnergy", totalE, "TotalEnergy[40]/F");
	tInfo->Branch("MaxEnergy", maxE, "MaxEnergy[40]/F");
	tInfo->Branch("Sigma", sigma, "Sigma[40]/F");
	tInfo->Branch("Amplitude", ampl, "Amplitude[40]/F");
	tInfo->Branch("Chi2", chi2, "Chi2[40]/F");
	tInfo->Branch("NDEF", ndef, "NDEF[40]/I");

	int nCols(0), nRows(0);
	geometry->GetVirtualNColRow(6, nCols, nRows);
	// TH2I* hDigitColRowMap = new TH2I("DigitColRowMap", "Digit Col vs Row;Col;Row", nCols, 0, nCols, nRows, 0, nRows);
	// TH2F* hDigitColRowMapEnergy = new TH2F("DigitColRowMapEnergy", "Digit Col vs Row Vs Energy;Col;Row", nCols, 0, nCols, nRows, 0, nRows);

	// TH1F *hECALClust = new TH1F("hECALClust", "ECAL clusters energy", 1000, 0, 1000);
	// hECALClust->GetXaxis()->SetTitle("Energy (GeV)");
	// TH1F *hHCALClust = new TH1F("hHCALClust", "HCAL clusters energy", 1000, 0, 1000);
	// hHCALClust->GetXaxis()->SetTitle("Energy (GeV)");
	// TH1F *hTotalClust = new TH1F("hTotalClust", "Total clusters energy", 1000, 0, 1000);
	// hTotalClust->GetXaxis()->SetTitle("Energy (GeV)");

	// TH2F *hHCALEnergy_IntegTot = new TH2F("hHCALEnergy_IntegTot", "Total clusters energy from integral Vs total;Integral;Total", 1000, 0, 1000, 1000, 0, 1000);

	// profile histos
	TProfile2D **showProf2D = new TProfile2D *[nSegments];
	TProfile **showProf1D = new TProfile *[nSegments];
	TProfile **showProf1D_dEdr = new TProfile *[nSegments];
	for (Int_t i = 0; i < nSegments; i++)
	{
		Int_t nBins = 5;
		Float_t cellSize = 1;
		if (i == 1 || i == 3)
		{
			nBins = 50;
			cellSize = 0.1;
		}
		Float_t lim = 5 + 0.5 * cellSize;
		showProf2D[i] = new TProfile2D(Form("prof2D_%d", i), Form("Shower profile seg %d;x (cm);y (cm)", i), 2 * nBins + 1, -lim, lim, 2 * nBins + 1, -lim, lim);
		showProf1D[i] = new TProfile(Form("prof_%d", i), Form("Shower profile seg %d;r (cm);dN/rdr", i), nBins + 1, -0.5 * cellSize, lim);
		showProf1D_dEdr[i] = new TProfile(Form("prof_dEdr_%d", i), Form("Shower profile seg %d;r (cm);dN/dr", i), nBins + 1, -0.5 * cellSize, lim);
	}

	// Loop over events in the folder
	for (Int_t ievt = 0; ievt < fRunLoader->GetNumberOfEvents(); ievt++)
	{

		totalEvents++;
		event = ievt;

		Int_t ie = fRunLoader->GetEvent(ievt);
		if (ie != 0)
		{
			cout << "Error reading event " << ievt << endl;
			continue;
		}
		TTree *treeH = fFOCALLoader->TreeH();
		if (!treeH)
		{
			std::cout << "TreeH is corrupt\n";
			break;
			continue;
		}
		cout << "Event: " << ievt << " with " << treeH->GetEntries() << " tracks" << endl;
		TTree *treeK = fRunLoader->TreeK();

		for (Int_t i = 0; i < 4; i++)
		{
			pdgCode[i] = 0;
			e[i] = 0;
			pt[i] = 0;
			phi[i] = 0;
			theta[i] = 0;
			;
			Vx[i] = 0;
			Vy[i] = 0;
			Vz[i] = 0;
			conv_flag[i] = 0;
			cVx[i] = 0;
			cVy[i] = 0;
			cVz[i] = 0;
		}

		AliStack *stack = fRunLoader->Stack();
		if (stack->GetNprimary() != 1)
			cout << "More than one primary found; this will lead to unexpected results... " << endl;

		auto primary = stack->Particle(0);
		pdgCode[0] = primary->GetPdgCode();
		e[0] = primary->Energy();
		pt[0] = primary->Pt();
		phi[0] = primary->Phi();
		theta[0] = primary->Theta();
		Vx[0] = primary->Vx();
		Vy[0] = primary->Vy();
		Vz[0] = primary->Vz();
		if (primary->GetFirstDaughter() <= 0 && primary->GetLastDaughter() <= 0)
			cout << "No daughter tracks... " << endl;
		else
		{
			Int_t iTrk = 1;
			for (Int_t i = primary->GetFirstDaughter(); i <= primary->GetLastDaughter(); i++)
			{
				TParticle *part = stack->Particle(i);
				if (iTrk < 4)
				{
					pdgCode[iTrk] = part->GetPdgCode();
					e[iTrk] = part->Energy();
					pt[iTrk] = part->Pt();
					phi[iTrk] = part->Phi();
					theta[iTrk] = part->Theta();
					Vx[iTrk] = part->Vx();
					Vy[iTrk] = part->Vy();
					Vz[iTrk] = part->Vz();
				}
				else
					cout << "WARNING: too many daughters" << endl;

				if (part->GetPdgCode() == 22 && part->GetFirstDaughter() >= 0)
				{
					// check conversion
					auto cdaughter = stack->Particle(part->GetFirstDaughter());
					if (cdaughter->Vz() < z_FOCAL_front)
					{
						conv_flag[iTrk] = 1;
					}
					cVx[iTrk] = cdaughter->Vx();
					cVy[iTrk] = cdaughter->Vy();
					cVz[iTrk] = cdaughter->Vz();
				}
				iTrk++;
			}
		}

		//          // Call the digitizer
		digitizer->Hits2Digits(treeH->GetBranch("FOCAL"));
		TClonesArray *digitsArray = digitizer->GetSegmentDigits();
		//
		//	        //
		TH2F **histograms = new TH2F *[nSegments];
		TH1D **histograms1D = new TH1D *[nSegments];
		// TH2F* histForHCAL;

		Int_t nCol, nRow;
		Float_t sizeX = geometry->GetFOCALSizeX();
		Float_t sizeY = geometry->GetFOCALSizeY();

		char name[20] = "Original";
		char name2[20] = "Projection";

		for (Int_t i = 0; i < nSegments; i++)
		{
			geometry->GetVirtualNColRow(i, nCol, nRow);
			Float_t x0, y0, z0, r;
			z0 = geometry->GetVirtualSegmentZ(i);
			r = TMath::Tan(theta[0]) * z0; // Phi = 0 && x = 1 && y = 0;
			x0 = TMath::Cos(phi[0]) * r;
			y0 = TMath::Sin(phi[0]) * r;
			histograms[i] = new TH2F(Form("%s_%i_%i", name, ievt, i), Form("%s_%i_%i", name, ievt, i), nCol, 0, nCol, nRow, 0, nRow);
			histograms[i]->SetStats(0);

			// if(i==6){
			// 	histForHCAL = new TH2F("DigitColRowMapEnergy", "Digit Col vs Row Vs Energy;Col;Row", nCol, 0, nCol, nRow, 0, nRow);
			// }
			//            histograms1D[i] = new TH1F(Form("%s_%i",name2,i),Form("%s_%i",name2,i),nCol,-1*sizeX/2,sizeX/2);
		}

		cout << "Filling histograms" << endl;

		Int_t *rowSeed = new Int_t[nSegments];
		Int_t *colSeed = new Int_t[nSegments];
		Float_t *ampSeed = new Float_t[nSegments];
		for (Int_t iSeg = 0; iSeg < nSegments; iSeg++)
		{
			rowSeed[iSeg] = -1;
			colSeed[iSeg] = -1;
			ampSeed[iSeg] = -1;
		}
		TObjArray digitsForProfile;

		double HCALEnergy(0.), ECALEnergy(0.);
		for (Int_t i = 0; i < digitsArray->GetEntries(); i++)
		{
			AliFOCALdigit *digit = (AliFOCALdigit *)digitsArray->UncheckedAt(i);
			Int_t col = digit->GetCol();
			Int_t row = digit->GetRow();
			Float_t x, y, z;
			Int_t segment = digit->GetSegment();

			Float_t x0, y0, z0, r;
			z0 = geometry->GetVirtualSegmentZ(segment);
			r = TMath::Tan(theta[0]) * z0; // Phi = 0 && x = 1 && y = 0;
			x0 = TMath::Cos(phi[0]) * r;
			y0 = TMath::Sin(phi[0]) * r;

			geometry->GetXYZFromColRowSeg(col, row, segment, x, y, z);
			Int_t energy = digit->GetAmp();
			histograms[segment]->Fill(col, row, (Float_t)energy);
			//		        histograms1D[segment]->Fill(x,(Float_t)energy);
			
			// prepare for profiles
			if (TMath::Abs(x - x0) < profileRadius &&
				TMath::Abs(y - y0) < profileRadius)
			{
				digitsForProfile.AddLast(digit);
				if (energy > ampSeed[segment])
				{
					rowSeed[segment] = row;
					colSeed[segment] = col;
					ampSeed[segment] = energy;
				}
			}
		}

		// Fill profiles
		//
		// First fill per-event profiles,then average. This matter for the exmpy bins
		//  (alternative would be to use TProfile, but fill empty cells as a zero in the prof)

		cout << "Reset per event hists" << endl;
		for (Int_t i = 0; i < nSegments; i++)
		{
			showProf2DEvent[i]->Reset();
			showProf1DEvent[i]->Reset();
		}
		cout << "Fill digits" << endl;
		for (Int_t iDigit = 0; iDigit < digitsForProfile.GetEntriesFast(); iDigit++)
		{
			AliFOCALdigit *digit = (AliFOCALdigit *)digitsForProfile.UncheckedAt(iDigit);
			Int_t col = digit->GetCol();
			Int_t row = digit->GetRow();
			Float_t x, y, z;
			Int_t segment = digit->GetSegment();
			geometry->GetXYZFromColRowSeg(col, row, segment, x, y, z);
			Float_t x0, y0, z0;
			geometry->GetXYZFromColRowSeg(colSeed[segment], rowSeed[segment], segment, x0, y0, z0);
			Int_t energy = digit->GetAmp();

			showProf2DEvent[segment]->Fill(x - x0, y - y0, energy / ampSeed[segment]);
			Float_t r = TMath::Sqrt(pow(x - x0, 2) + pow(y - y0, 2));
			showProf1DEvent[segment]->Fill(r, energy / ampSeed[segment]);
		}

		cout << "Add hists" << endl;
		for (Int_t iSeg = 0; iSeg < nSegments; iSeg++)
		{
			Int_t xbins = showProf2DEvent[iSeg]->GetXaxis()->GetNbins();
			Int_t ybins = showProf2DEvent[iSeg]->GetYaxis()->GetNbins();
			// could use this to count # bins per r bin
			for (Int_t ibin = 0; ibin < xbins; ibin++)
			{
				for (Int_t jbin = 0; jbin < ybins; jbin++)
				{
					Float_t x = showProf2DEvent[iSeg]->GetXaxis()->GetBinCenter(ibin + 1);
					Float_t y = showProf2DEvent[iSeg]->GetYaxis()->GetBinCenter(jbin + 1);
					showProf2D[iSeg]->Fill(x, y, showProf2DEvent[iSeg]->GetBinContent(ibin + 1, jbin + 1));
				}
			}
			xbins = showProf1DEvent[iSeg]->GetNbinsX();
			for (Int_t ibin = 0; ibin < xbins; ibin++)
			{
				Float_t r = showProf1DEvent[iSeg]->GetXaxis()->GetBinCenter(ibin + 1);
				// Would be better to use actual (discretised) area instead of 1/r
				if (r != 0)
					showProf1D[iSeg]->Fill(r, showProf1DEvent[iSeg]->GetBinContent(ibin + 1) / r);
				showProf1D_dEdr[iSeg]->Fill(r, showProf1DEvent[iSeg]->GetBinContent(ibin + 1));
			}
		}

		cout << "Fitting histograms" << endl;

		for (Int_t i = 0; i < nSegments; i++)
		{
			// Compute initial parameters
			Float_t x0, y0, z0, r;
			z0 = geometry->GetVirtualSegmentZ(i);
			r = TMath::Tan(theta[0]) * z0; // Phi = 0 && x = 1 && y = 0;
			x0 = TMath::Cos(phi[0]) * r;
			y0 = TMath::Sin(phi[0]) * r;
			Float_t amp = histograms[i]->GetMaximum();
			// Define function
			TF2 *fc = new TF2("cauchy", "[3]/(1+((x-[0])*(x-[0]) + (y-[1])*(y-[1]))/([2]*[2]))", x0 - infoParameterFitRadius, x0 + infoParameterFitRadius, y0 - infoParameterFitRadius, y0 + infoParameterFitRadius);
			fc->SetParameters(x0, y0, infoParameterSigma, 1);
			fc->SetParNames("x0", "y0", "#sigma", "I");

			geometry->GetVirtualNColRow(i, nCol, nRow);
			Int_t yBin = nRow * infoParameterHistRadius / sizeY / 2;
			histograms1D[i] = histograms[i]->ProjectionX(Form("%s_%i", name2, i), yBin, yBin);
			histograms1D[i]->SetStats(0);
			TF1 *fc1D = new TF1("cauchy1D", "[2]/(1+TMath::Power((x-[0]),2)/([1]*[1]))", x0 - infoParameterFitRadius, x0 + infoParameterFitRadius);
			fc1D->SetParameters(x0, infoParameterSigma, 1);
			fc1D->SetParNames("x0", "#sigma", "I");

			cout << Form("Initial parameters: %.2e, %.2e, %.2e, %i", x0, y0, infoParameterSigma, 1) << endl;

			// Fit the histogram
			if (pdgCode[0] == 22)
			{
				histograms[i]->Fit("cauchy", "RN", "goff");
				histograms1D[i]->Fit("cauchy1D", "R", "goff");
			}

			// Get Results
			Double_t pars[4];
			fc->GetParameters(pars);
			totalE[i * 2] = histograms[i]->Integral();
			maxE[i * 2] = amp;
			sigma[i * 2] = TMath::Abs(pars[2]);
			ampl[i * 2] = pars[3];
			chi2[i * 2] = fc->GetChisquare();
			ndef[i * 2] = fc->GetNDF();
			cout << Form("Final parameters: %.2e, %.2e, %.2e, %.2e; Chi/n = %.2e/%i",
						 pars[0], pars[1], pars[2], pars[3], fc->GetChisquare(), fc->GetNDF())
				 << endl;
			fc1D->GetParameters(pars);
			totalE[i * 2 + 1] = 0;
			maxE[i * 2 + 1] = histograms1D[i]->GetMaximum();
			sigma[i * 2 + 1] = TMath::Abs(pars[1]);
			ampl[i * 2 + 1] = pars[2];
			chi2[i * 2 + 1] = fc->GetChisquare();
			ndef[i * 2 + 1] = fc->GetNDF();
			cout << Form("Final parameters1D: %.2e, %.2e, %.2e; Chi/n = %.2e/%i",
						 pars[0], pars[1], pars[2], fc1D->GetChisquare(), fc1D->GetNDF())
				 << endl;

			cout << "Saving..." << endl;
			fc->Delete();
			fc1D->Delete();
			histograms[i]->Delete();
			histograms1D[i]->Delete();
		}
		// Fill the ttree
		tInfo->Fill();

		delete[] rowSeed;
		delete[] colSeed;
		delete[] ampSeed;

		delete[] histograms;
		delete[] histograms1D;

	} // end events

	f->cd();
	// hDigitColRowMap->Write();
	// hDigitColRowMapEnergy->Scale(1./fRunLoader->GetNumberOfEvents());
	// tInfo->Write();
	f->Write();
	f->Close();

	delete[] showProf2D;
	delete[] showProf1D;
	delete[] showProf1D_dEdr;

	fRunLoader->Delete();

	delete digitizer;
	delete[] maxE;
	delete[] sigma;
	delete[] ampl;
	delete[] chi2;
	delete[] ndef;
}
