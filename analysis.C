

const double speedOfLight = 299792458; // m/s
const double pionMass = 0.13957; // GeV /c^2

TFile* data;
TFile* fout;

TTree* recTree;

TH1D* hnM2;
TH1D* hA;
TH1D* hB;
TH1D* hC;
TH1D* hSqrt;
TH1D* hTofTime;
TH1D* hTofLength;
TH2F* hControlPlot;

Double_t nSigPPion,nSigPKaon,nSigPProton;
Double_t missingPt, invMassPion, invMassKaon, invMassProton;
Double_t deltaTOF, mSquared;
Double_t deltaDeltaTOFPion, deltaTOFPionExpected;
Double_t deltaDeltaTOFKaon, deltaTOFKaonExpected;
Double_t deltaDeltaTOFProton, deltaTOFProtonExpected;

Double_t dEdx1, dEdx2;
Double_t momentum1, momentum2;
Double_t tranMomenta1, tranMomenta2;
Double_t charge1, charge2;
Double_t TOFlength1, TOFlength2;
Double_t TOFtime1, TOFtime2;

Double_t vertexZ;

Double_t DcaXY1, DcaZ1, DcaXY2, DcaZ2;
Double_t NhitFit1, NhitsDEdx1, NhitFit2, NhitsDEdx2;
Double_t Eta1, Eta2;


Bool_t elastic, fourPiState;

void Init();
void ConnectInput();
void Make();

void analysis()
{
    TString input = "/home/truhlar/Desktop/STAR/CEP/Analysis/Data/MergeWO1.root";
    TString output = "/home/truhlar/Desktop/STAR/CEP/Analysis/Outputs/analOutput.root";

    data = TFile::Open(input, "read");
	if (!data){
      cout<<"Error: cannot open "<<input<<endl;
      return;
    }

    fout = new TFile(output,"RECREATE");
	Init(); // Preparing histograms 
	ConnectInput(); // Connecting input

////////////// Making histograms

	Long64_t nev = recTree->GetEntries();
	cout<<"Proccesing "<<nev<<" events"<<endl;
	//event loop
	//nev = 1000;
	for(Long64_t iev=0; iev<nev; ++iev) { //get the event
		recTree->GetEntry(iev); 
		Make();
	} 

    fout->Write();
    fout->Close();
}




void Make(){

// Calculate mSquared
	double speedOfLight2 = speedOfLight*speedOfLight;
	double speedOfLight4 = speedOfLight2*speedOfLight2;
	double length1Squared = TOFlength1*TOFlength1/(100*100); // convert from cm to m
	double length2Squared = TOFlength2*TOFlength2/(100*100); // convert from cm to m
	double deltaTime2 = ((TOFtime2 - TOFtime1)*(TOFtime2 - TOFtime1))/(pow(10.0,18.0)); // convert from ns to s
	double deltaTime4 = deltaTime2*deltaTime2;
	double oneOverMomentum1sq = 1/(momentum1*momentum1);
	double oneOverMomentum2sq = 1/(momentum2*momentum2);
	double cEq = -2*length1Squared*length2Squared + speedOfLight4*deltaTime4 + length2Squared*length2Squared + length1Squared*length1Squared -2*speedOfLight2*deltaTime2*(length2Squared + length1Squared);
	double bEq1 = -2*length1Squared*length2Squared*(oneOverMomentum1sq + oneOverMomentum2sq);
    double bEq2 = 2*length1Squared*length1Squared*oneOverMomentum1sq + 2*length2Squared*length2Squared*oneOverMomentum2sq;
    double bEq3 = -2*speedOfLight2*deltaTime2*(length1Squared*oneOverMomentum1sq + length2Squared*oneOverMomentum2sq);
    double bEq =  bEq1 + bEq2 + bEq3; 
	double aEq = -2*length1Squared*length2Squared*oneOverMomentum1sq*oneOverMomentum2sq + length1Squared*length1Squared*oneOverMomentum1sq*oneOverMomentum1sq + length2Squared*length2Squared*oneOverMomentum2sq*oneOverMomentum2sq;
	mSquared = (-bEq + sqrt(bEq*bEq-4*aEq*cEq)) / (2*aEq);

	if(TOFtime1 < 0 || TOFtime2 < 0 || TOFlength1 < 0 || TOFlength2 < 0)
		mSquared = -3.0;
    else{
        double controlEq = (TOFlength1/(100*speedOfLight))*sqrt(1 + pionMass*pionMass*oneOverMomentum1sq) + (TOFlength2/(100*speedOfLight))*sqrt(1 + pionMass*pionMass*oneOverMomentum2sq); 
    	hControlPlot->Fill(TOFtime2 - TOFtime1, controlEq*pow(10.0,9.0));
      //  cout<<bEq1<<" : "<<bEq2<<" : "<<bEq3<<endl;
        hnM2->Fill(mSquared);
        hA->Fill(aEq);
        hB->Fill(bEq);
        hC->Fill(cEq);
        hSqrt->Fill(sqrt(bEq*bEq-4*aEq*cEq));
        hTofTime->Fill(TOFtime1);
        hTofTime->Fill(TOFtime2);
        hTofLength->Fill(TOFlength1);
        hTofLength->Fill(TOFlength2);
    }


}


void Init(){

	hnM2 = new TH1D("mSquared", "mSquared", 200, -0.5, 1.5);
    hA = new TH1D("A", "A", 2000, -5.0, 45.0);
    hB = new TH1D("B", "B", 2000, -30.0, 30.0);
    hC = new TH1D("C", "C", 2000, -20.0, 20.0);
    hSqrt = new TH1D("SQRT", "SQRT", 2000, -5.0, 45.0);
    hTofTime = new TH1D("TOF time", "TOF time", 2000, -30.0, 30.0);
    hTofLength = new TH1D("TOF length", "TOF length", 200, 200.0, 500.0);
    hControlPlot = new TH2F("hControlPlot", "hControlPlot",100,-7,7,100,0,40);
}

void ConnectInput(){
	recTree = dynamic_cast<TTree*>( data->Get("recTree") );
	if (!recTree){
      cout<<"Error: cannot get recTree"<<endl;
      return;
   }

	recTree->SetBranchAddress("nSigPPion", &nSigPPion);
	recTree->SetBranchAddress("nSigPKaon", &nSigPKaon);
	recTree->SetBranchAddress("nSigPProton", &nSigPProton);
	recTree->SetBranchAddress("missingPt", &missingPt);
	recTree->SetBranchAddress("dEdx1", &dEdx1);
	recTree->SetBranchAddress("dEdx2", &dEdx2);
	recTree->SetBranchAddress("momentum1", &momentum1);
	recTree->SetBranchAddress("momentum2", &momentum2);
	recTree->SetBranchAddress("charge1", &charge1);
	recTree->SetBranchAddress("charge2", &charge2);    
	recTree->SetBranchAddress("mSquared", &mSquared); 
	recTree->SetBranchAddress("deltaTOF", &deltaTOF);
	recTree->SetBranchAddress("TOFtime1", &TOFtime1);
	recTree->SetBranchAddress("TOFtime2", &TOFtime2);
	recTree->SetBranchAddress("TOFlength1", &TOFlength1);
	recTree->SetBranchAddress("TOFlength2", &TOFlength2);
	recTree->SetBranchAddress("elastic", &elastic);
}


