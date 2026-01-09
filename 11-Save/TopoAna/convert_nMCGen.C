
// * This is a script to convert the int type nMCGen to double type nMCGenSelf
// in order to adapt the basf2 output file with int type nMCGen to the topoana.
// * The number of MCGen branches (NBRANCHES) is set to 200 in this script, if you use a 
// different number when running basf2, please modify NBRANCHES.
// * After running this script, please set the following variable in the *.card file:
// % TBranch name of the number of particles 
// to MCGenSelf.
#define NBRANCHES 200
double GetnMCGenFromMCGenMothIndex(double *MCGenMothIndices)
{
	for (int i = 0; i< NBRANCHES; ++i) {
		if (MCGenMothIndices[NBRANCHES-1-i]) return NBRANCHES - i;
	}
	std::cout << "Failed to get nMCGen from MCGenMothIndex! Exit!" << std::endl;
	exit(-1);
	return -1;
}
void convert_nMCGen(const char * inputFileName, const char * inputTreeName, const char * outputFileName)
{
	TFile * finput = TFile::Open(inputFileName);
	TTree * tinput = (TTree*)finput->Get(inputTreeName);
	TFile * foutput = TFile::Open(outputFileName, "recreate");
	tinput->SetBranchStatus("*", 1);
	int nMCGen;
	tinput->SetBranchAddress("nMCGen", &nMCGen);
	double MCGenPDGs[NBRANCHES];
	double MCGenMothIndices[NBRANCHES];
	for (int i = 0; i < NBRANCHES; ++i) {
			tinput->SetBranchAddress(Form("MCGenPDG_%d", i), &(MCGenPDGs[i]));
			tinput->SetBranchAddress(Form("MCGenMothIndex_%d", i), &(MCGenMothIndices[i]));
	}
	TTree * toutput = tinput->CloneTree(0);
	double nMCGenSelf;
	toutput->Branch("nMCGenSelf", &nMCGenSelf, "nMCGenSelf/D");
	for (int i = 0; i < tinput->GetEntries(); ++i) {
		tinput->GetEntry(i);
		if (nMCGen <= NBRANCHES) nMCGenSelf = nMCGen;
		else nMCGen = GetnMCGenFromMCGenMothIndex(MCGenMothIndices);
		toutput->Fill();
	}
	foutput->cd();
	toutput->Write(0,TObject::kOverwrite);
	foutput->Close();
	finput->Close();
}