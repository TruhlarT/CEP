#include "PID.h"


//_____________________________________________________________________________
PID::PID(string innam, string outnam)
{
	//constructor
	input = innam;
	output = outnam;

	speedOfLight = 299792458; // m/s
	pionMass = 0.13957; // GeV /c^2

	particleID[0] = TString("Pion");
	particleID[1] = TString("Kaon");
	particleID[2] = TString("Proton");
	doubleParticleID[0] = TString("PionVsKaon"); 
	doubleParticleID[1] = TString("PionsVsProton");
	doubleParticleID[2] = TString("KaonVsProton");

	cout << "PID::PID() called" << endl;

}//PID

//_____________________________________________________________________________
PID::~PID()
{
  //destructor

  cout << "PID::~PID() destructor called" << endl;

}//~PID


//_____________________________________________________________________________
void PID::PlotHistogram(){

	data = TFile::Open(input, "read");
	if (!data){
      cout<<"Error: cannot open "<<input<<endl;
      return;
   }

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
	Plot();

//////////////////// Output

	cout<<"Creating output file: "<<output<<endl;
	fout = new TFile(output,"UPDATE");

	nM2->Write("nM2");
	nM2->SaveAs("PID/nM2.png");
	nM2->Close();

	cdEdx->Write("cdEdx");
	cdEdx->SaveAs("PID/cdEdx.png");
	cdEdx->Close();

	cdEdxLog->Write("cdEdxLog");
	cdEdxLog->SaveAs("PID/cdEdxLog.png");
	cdEdxLog->Close();


	cDeltaTOF->Write("cDeltaTOF");
	cDeltaTOF->SaveAs("PID/cDeltaTOF.png");
	cDeltaTOF->Close();

	cDeltaDeltaTOF->Write("cDeltaDeltaTOF");
	cDeltaDeltaTOF->SaveAs("PID/cDeltaDeltaTOF.png");
	cDeltaDeltaTOF->Close();

	cdEdxBasicPID->Write("cdEdxBasicPID");
	cdEdxBasicPID->SaveAs("PID/cdEdxBasicPID.png");
	cdEdxBasicPID->Close();

	for(int i = 0; i < 3; ++i){

		nSigmaParticleVsM2[i]->Write("nSigma" + particleID[i] + "VsM2");
		nSigmaParticleVsM2[i]->SaveAs("PID/nSigma" + particleID[i] + "VsM2.png");
		nSigmaParticleVsM2[i]->Close();

		int i1 = i;
		int i2 = i+1;
		if(i==2){
			i1 = 0;
			i2 = 2;
		}
		nSigmaParticleVsOther[i]->Write("nSigma" + particleID[i1] + "Vs" + particleID[i2]);
		nSigmaParticleVsOther[i]->SaveAs("PID/nSigma" + particleID[i1] + "Vs" + particleID[i2] + ".png");
		nSigmaParticleVsOther[i]->Close();

		cChMomVsnSigma[i]->Write("cChMomVsnSigma"+particleID[i]);
		cChMomVsnSigma[i]->SaveAs("PID/cChMomVsnSigma"+particleID[i] + ".png");
		cChMomVsnSigma[i]->Close();

	}

   fout->Close();

}

void PID::Make(){
	hnSigma[0]->Fill(nSigPPion,nSigPKaon);
	hnSigma[1]->Fill(nSigPKaon,nSigPProton);
	hnSigma[2]->Fill(nSigPPion,nSigPProton);

	hdEdx->Fill(charge1*momentum1,dEdx1);
	hdEdx->Fill(charge2*momentum2,dEdx2);

	hdEdxLog->Fill(dEdx1);
	hdEdxLog->Fill(dEdx2);

	hChMomVsnSigma[0]->Fill(charge1*momentum1,nSigPPion);
	hChMomVsnSigma[1]->Fill(charge1*momentum1,nSigPKaon);
	hChMomVsnSigma[2]->Fill(charge1*momentum1,nSigPProton);

	hChMomVsnSigma[0]->Fill(charge2*momentum2,nSigPPion);
	hChMomVsnSigma[1]->Fill(charge2*momentum2,nSigPKaon);
	hChMomVsnSigma[2]->Fill(charge2*momentum2,nSigPProton);

// Calculate mSquared
	double speedOfLight2 = speedOfLight*speedOfLight;
	double speedOfLight4 = speedOfLight2*speedOfLight2;
	double length1Squared = TOFlength1*TOFlength1/(100*100*speedOfLight2); // convert from cm to m
	double length2Squared = TOFlength2*TOFlength2/(100*100*speedOfLight2); // convert from cm to m
	double deltaTime2 = (TOFtime2 - TOFtime1)*(TOFtime2 - TOFtime1)/(pow(10.0,18.0)); // convert from ns to s
	double deltaTime4 = deltaTime2*deltaTime2;
	double oneOverMomentum1sq = 1/(momentum1*momentum1);
	double oneOverMomentum2sq = 1/(momentum2*momentum2);
	double cEq = -2*length1Squared*length2Squared + speedOfLight4*deltaTime4 + length2Squared*length2Squared + length1Squared*length1Squared -2*speedOfLight2*deltaTime2*(length2Squared + length1Squared);
	double bEq = -2*length1Squared*length2Squared*speedOfLight2*(oneOverMomentum1sq + oneOverMomentum2sq) + 2*length1Squared*length1Squared*speedOfLight2*oneOverMomentum1sq + 2*length2Squared*length2Squared*speedOfLight2*oneOverMomentum2sq -2*speedOfLight4*deltaTime2*(length1Squared*oneOverMomentum1sq + length2Squared*oneOverMomentum2sq);
	double aEq = -2*length1Squared*length2Squared*speedOfLight4*oneOverMomentum1sq*oneOverMomentum2sq + speedOfLight4*(length1Squared*length1Squared*oneOverMomentum1sq*oneOverMomentum1sq + length2Squared*length2Squared*oneOverMomentum2sq*oneOverMomentum2sq);
	mSquared = (-bEq + sqrt(bEq*bEq-4*aEq*cEq)) / (2*aEq);

	if(TOFtime1 < 0 || TOFtime2 < 0 || TOFlength1 < 0 || TOFlength2 < 0)
		mSquared = -1.0;

//////////////////////

	hnM2->Fill(mSquared);
	hnSigmaVsM2[0]->Fill(nSigPPion,mSquared);
	hnSigmaVsM2[1]->Fill(nSigPKaon,mSquared);
	hnSigmaVsM2[2]->Fill(nSigPProton,mSquared);

/////////////// Calculate deltaDeltaTOF
	double expectedTime1 = (TOFlength1/ speedOfLight) * sqrt(1 + pow(pionMass/ momentum1, 2)) * pow(10,7); // in ns
	double expectedTime2 = (TOFlength2 / speedOfLight) * sqrt(1 + pow(pionMass/ momentum2,2)) * pow(10,7); // in ns
	double deltaTOF = TOFtime2 - TOFtime1;
	double deltaTOFexpected = expectedTime2 - expectedTime1;
	deltaDeltaTOF = deltaTOF - deltaTOFexpected;


   if(TOFtime1 < 0 || TOFtime2 < 0 || TOFlength1 < 0 || TOFlength2 < 0){
		deltaDeltaTOF = -999;
		deltaTOF = -999;   
		deltaTOFexpected = -999;  
   }
////////////////
   //cout<<deltaTOF <<" "<<deltaTOFexpected <<"  "<< deltaDeltaTOF<<endl;
	hdeltaDeltaTOF -> Fill(deltaDeltaTOF);
	hdeltaTOF -> Fill(deltaTOF);
	hdeltaTOFexpected -> Fill(deltaTOFexpected);

	//Basic PID
	if(nSigPProton < 3 && nSigPKaon > 3 && nSigPPion > 3){
		hdEdxBasicPID[2]->Fill(charge1*momentum1,dEdx1);
		hdEdxBasicPID[2]->Fill(charge2*momentum2,dEdx2);
	}else if( nSigPKaon < 3 && nSigPProton > 3 && nSigPPion > 3){
		hdEdxBasicPID[1]->Fill(charge1*momentum1,dEdx1);
		hdEdxBasicPID[1]->Fill(charge2*momentum2,dEdx2);	
	}else{
		hdEdxBasicPID[0]->Fill(charge1*momentum1,dEdx1);
		hdEdxBasicPID[0]->Fill(charge2*momentum2,dEdx2);
	}

}


void PID::Plot(){

	for(int i = 0; i < 3; ++i){

		int i1 = i;
		int i2 = i+1;
		if(i==2){
			i1 = 0;
			i2 = 2;
		}

		nSigmaParticleVsOther[i] = new TCanvas("nSigma" + particleID[i1] + "Vs" + particleID[i2],"nSigma" + particleID[i1] + "Vs" + particleID[i2],1200,800);
		gPad->SetMargin(0.1,0.09,0.14,0.02);
		gStyle->SetPalette(1);
		hnSigma[i]->SetStats(0);
		hnSigma[i]->SetTitle(" ; n#sigma^{pair}_{" + particleID[i1] + "}; n#sigma^{pair}_{" + particleID[i2] + "}");
		hnSigma[i]->GetXaxis()->SetTitleFont(42);
		hnSigma[i]->GetYaxis()->SetTitleFont(42);
		hnSigma[i]->GetXaxis()->SetLabelFont(42);
		hnSigma[i]->GetYaxis()->SetLabelFont(42);
		hnSigma[i]->GetXaxis()->SetTitleSize(.055);
		hnSigma[i]->GetYaxis()->SetTitleSize(.055);
		hnSigma[i]->GetXaxis()->SetTitleOffset(1);
		hnSigma[i]->GetYaxis()->SetTitleOffset(0.8);
		hnSigma[i]->Draw("colz");

		TPaveText *textPub = new TPaveText(0.68,0.75,0.75,0.97,"brNDC");
		textPub -> SetTextSize(0.055);
		textPub -> SetTextFont(42);
		textPub -> SetFillColor(0);	
		textPub -> AddText("p + p #rightarrow p + #pi^{+} + #pi^{-} + p");
		textPub -> AddText("#sqrt{s} = 510 GeV");
		textPub -> Draw("same");

		nSigmaParticleVsOther[i]->Update();

		nSigmaParticleVsM2[i] = new TCanvas("nSigma" + particleID[i] + "VsM2","nSigma" + particleID[i] + "VsM2",1200,800);
		gPad->SetMargin(0.1,0.09,0.14,0.02);
		gStyle->SetPalette(1);
		hnSigmaVsM2[i]->SetStats(0);
		hnSigmaVsM2[i]->SetTitle(" ; n#sigma^{pair}_{" + particleID[i] + "}; m^{2}_{TOF} [GeV^{2}/c^{4}]");
		hnSigmaVsM2[i]->GetXaxis()->SetTitleFont(42);
		hnSigmaVsM2[i]->GetYaxis()->SetTitleFont(42);
		hnSigmaVsM2[i]->GetXaxis()->SetLabelFont(42);
		hnSigmaVsM2[i]->GetYaxis()->SetLabelFont(42);
		hnSigmaVsM2[i]->GetXaxis()->SetTitleSize(.055);
		hnSigmaVsM2[i]->GetYaxis()->SetTitleSize(.055);
		hnSigmaVsM2[i]->GetXaxis()->SetTitleOffset(1);
		hnSigmaVsM2[i]->GetYaxis()->SetTitleOffset(0.8);
		hnSigmaVsM2[i]->Draw("colz");

		TPaveText *textPub1 = new TPaveText(0.68,0.75,0.75,0.97,"brNDC");
		textPub1 -> SetTextSize(0.055);
		textPub1 -> SetTextFont(42);
		textPub1 -> SetFillColor(0);	
		textPub1 -> AddText("p + p #rightarrow p + #pi^{+} + #pi^{-} + p");
		textPub1 -> AddText("#sqrt{s} = 510 GeV");
		textPub1 -> Draw("same");

		nSigmaParticleVsM2[i]->Update();

		cChMomVsnSigma[i] = new TCanvas("cChMomVsnSigma"+particleID[i],"cChMomVsnSigma"+particleID[i],1200,800);
		gPad->SetMargin(0.1,0.09,0.14,0.02);
		gStyle->SetPalette(1);
		hChMomVsnSigma[i]->SetStats(0);
		hChMomVsnSigma[i]->SetTitle(" ; #frac{q}{e} #times p [GeV/c] ; n#sigma^{pair}_{" + particleID[i] + "}");
		hChMomVsnSigma[i]->GetXaxis()->SetTitleFont(42);
		hChMomVsnSigma[i]->GetYaxis()->SetTitleFont(42);
		hChMomVsnSigma[i]->GetXaxis()->SetLabelFont(42);
		hChMomVsnSigma[i]->GetYaxis()->SetLabelFont(42);
		hChMomVsnSigma[i]->GetXaxis()->SetTitleSize(.055);
		hChMomVsnSigma[i]->GetYaxis()->SetTitleSize(.055);
		hChMomVsnSigma[i]->GetXaxis()->SetTitleOffset(1);
		hChMomVsnSigma[i]->GetYaxis()->SetTitleOffset(0.8);
		hChMomVsnSigma[i]->Draw("colz");

		TPaveText *textPub2 = new TPaveText(0.23,0.75,0.3,0.97,"brNDC");
		textPub2 -> SetTextSize(0.055);
		textPub2 -> SetTextFont(42);
		textPub2 -> SetFillColor(0);	
		textPub2 -> AddText("p + p #rightarrow p + #pi^{+} + #pi^{-} + p");
		textPub2 -> AddText("#sqrt{s} = 510 GeV");
		textPub2 -> Draw("same");

		cChMomVsnSigma[i]->Update();
	}


	cdEdx = new TCanvas("cdEdx","cdEdx",1200,800);
	gPad->SetMargin(0.1,0.09,0.14,0.02);
	gStyle->SetPalette(1);
	hdEdx->SetStats(0);
	hdEdx->SetTitle(" ; #frac{q}{e} #times p [GeV/c] ;dE/dx [keV/cm]");
	hdEdx->GetXaxis()->SetTitleFont(42);
	hdEdx->GetYaxis()->SetTitleFont(42);
	hdEdx->GetXaxis()->SetLabelFont(42);
	hdEdx->GetYaxis()->SetLabelFont(42);
	hdEdx->GetXaxis()->SetTitleSize(.055);
	hdEdx->GetYaxis()->SetTitleSize(.055);
	hdEdx->GetXaxis()->SetTitleOffset(1);
	hdEdx->GetYaxis()->SetTitleOffset(0.8);
	hdEdx->Draw("colz");

	TPaveText *textPub3 = new TPaveText(0.24,0.75,0.31,0.97,"brNDC");
	textPub3 -> SetTextSize(0.055);
	textPub3 -> SetTextFont(42);
	textPub3 -> SetFillColor(0);	
	textPub3 -> AddText("p + p #rightarrow p + #pi^{+} + #pi^{-} + p");
	textPub3 -> AddText("#sqrt{s} = 510 GeV");
	textPub3 -> Draw("same");

	cdEdx->Update();

	nM2 = new TCanvas("nM2","nM2",1200,800);
	gPad->SetMargin(0.9,0.02,0.14,0.02);;
	gStyle->SetPalette(1);
	hnM2->SetStats(0);
	hnM2->SetTitle(" ; m^{2}_{TOF} [GeV^{2}/c^{4}]; Number of events");
	hnM2->GetXaxis()->SetTitleFont(42);
	hnM2->GetYaxis()->SetTitleFont(42);
	hnM2->GetXaxis()->SetLabelFont(42);
	hnM2->GetYaxis()->SetLabelFont(42);
	hnM2->GetXaxis()->SetTitleSize(.055);
	hnM2->GetYaxis()->SetTitleSize(.055);
	hnM2->GetXaxis()->SetTitleOffset(1);
	hnM2->GetYaxis()->SetTitleOffset(0.8);
	hnM2->Draw("colz");

	TPaveText *textPub4 = new TPaveText(0.68,0.75,0.75,0.97,"brNDC");
	textPub4 -> SetTextSize(0.055);
	textPub4 -> SetTextFont(42);
	textPub4 -> SetFillColor(0);	
	textPub4 -> AddText("p + p #rightarrow p + #pi^{+} + #pi^{-} + p");
	textPub4 -> AddText("#sqrt{s} = 510 GeV");
	textPub4 -> Draw("same");

	nM2->Update();

	cDeltaDeltaTOF = new TCanvas("DeltaDeltaTOF","DeltaDeltaTOF",1200,800);
	gPad->SetMargin(0.9,0.02,0.14,0.02);
	gStyle->SetPalette(1);
	hdeltaDeltaTOF->SetStats(0);
	hdeltaDeltaTOF->SetTitle(" ; #Delta #Delta TOF [ns]; Number of events");
	hdeltaDeltaTOF->GetXaxis()->SetTitleFont(42);
	hdeltaDeltaTOF->GetYaxis()->SetTitleFont(42);
	hdeltaDeltaTOF->GetXaxis()->SetLabelFont(42);
	hdeltaDeltaTOF->GetYaxis()->SetLabelFont(42);
	hdeltaDeltaTOF->GetXaxis()->SetTitleSize(.055);
	hdeltaDeltaTOF->GetYaxis()->SetTitleSize(.055);
	hdeltaDeltaTOF->GetXaxis()->SetTitleOffset(1);
	hdeltaDeltaTOF->GetYaxis()->SetTitleOffset(0.8);
	hdeltaDeltaTOF->Draw("colz");

	TPaveText *textPub5 = new TPaveText(0.68,0.75,0.75,0.97,"brNDC");
	textPub5 -> SetTextSize(0.055);
	textPub5 -> SetTextFont(42);
	textPub5 -> SetFillColor(0);	
	textPub5 -> AddText("p + p #rightarrow p + #pi^{+} + #pi^{-} + p");
	textPub5 -> AddText("#sqrt{s} = 510 GeV");
	textPub5 -> Draw("same");

	cDeltaDeltaTOF->Update();

	cDeltaTOF = new TCanvas("DeltaTOF","DeltaTOF",1200,800);
	gPad->SetMargin(0.9,0.02,0.14,0.02);
	gStyle->SetPalette(1);
	hdeltaTOF->SetStats(0);
	hdeltaTOF->SetTitle(" ; #Delta TOF [ns]; Number of events");
	hdeltaTOF->GetXaxis()->SetTitleFont(42);
	hdeltaTOF->GetYaxis()->SetTitleFont(42);
	hdeltaTOF->GetXaxis()->SetLabelFont(42);
	hdeltaTOF->GetYaxis()->SetLabelFont(42);
	hdeltaTOF->GetXaxis()->SetTitleSize(.055);
	hdeltaTOF->GetYaxis()->SetTitleSize(.055);
	hdeltaTOF->GetXaxis()->SetTitleOffset(1);
	hdeltaTOF->GetYaxis()->SetTitleOffset(0.8);
	hdeltaTOF->Draw("same");
	hdeltaTOFexpected->SetLineColor(2);
	hdeltaTOFexpected->Draw("same");


	TPaveText *textPub6 = new TPaveText(0.68,0.75,0.75,0.97,"brNDC");
	textPub6 -> SetTextSize(0.055);
	textPub6 -> SetTextFont(42);
	textPub6 -> SetFillColor(0);	
	textPub6 -> AddText("p + p #rightarrow p + #pi^{+} + #pi^{-} + p");
	textPub6 -> AddText("#sqrt{s} = 510 GeV");
	textPub6 -> Draw("same");


	TLegend *legend = new TLegend(0.62,0.6,0.85,0.75,"","brNDC");
	legend -> AddEntry(hdeltaTOF, "Data", "l");
	legend -> AddEntry(hdeltaTOFexpected, "Pion assumption", "fl");
	legend -> SetTextFont(42);
	legend -> SetFillStyle(0);
	legend -> SetLineColor(0);
	legend -> SetTextSize(0.04);
	legend -> Draw("same");

	cDeltaTOF->Update();


	cdEdxBasicPID = new TCanvas("cdEdxBasicPID","cdEdxBasicPID",1200,800);
	gPad->SetMargin(0.9,0.02,0.14,0.02);
	//gStyle->SetPalette(1);
	hdEdxBasicPID[0]->SetStats(0);
	hdEdxBasicPID[0]->SetTitle(" ; #frac{q}{e} #times p [GeV/c] ;dE/dx [keV/cm]");
	hdEdxBasicPID[0]->GetXaxis()->SetTitleFont(42);
	hdEdxBasicPID[0]->GetYaxis()->SetTitleFont(42);
	hdEdxBasicPID[0]->GetXaxis()->SetLabelFont(42);
	hdEdxBasicPID[0]->GetYaxis()->SetLabelFont(42);
	hdEdxBasicPID[0]->GetXaxis()->SetTitleSize(.055);
	hdEdxBasicPID[0]->GetYaxis()->SetTitleSize(.055);
	hdEdxBasicPID[0]->GetXaxis()->SetTitleOffset(1);
	hdEdxBasicPID[0]->GetYaxis()->SetTitleOffset(0.8);
	hdEdxBasicPID[0]->SetMarkerColor(4);
	hdEdxBasicPID[0]->SetFillColor(4);
	hdEdxBasicPID[0]->SetMarkerStyle(29);
	hdEdxBasicPID[0]->SetMarkerSize(1);
	hdEdxBasicPID[0]->Draw("SCAT");
	gPad->Update();
	
	hdEdxBasicPID[1]->SetMarkerColor(2);
	hdEdxBasicPID[1]->SetFillColor(2);
	hdEdxBasicPID[1]->SetMarkerStyle(29);
	hdEdxBasicPID[1]->SetMarkerSize(1);
	hdEdxBasicPID[1]->Draw("same");
	gPad->Update();
	hdEdxBasicPID[2]->SetFillColor(8);
	hdEdxBasicPID[2]->SetMarkerColor(8);
	hdEdxBasicPID[2]->SetMarkerStyle(29);
	hdEdxBasicPID[2]->SetMarkerSize(1);
	hdEdxBasicPID[2]->Draw("SAME");
	cdEdxBasicPID->cd();
	gPad->Update();

	TPaveText *textPubPID = new TPaveText(0.24,0.75,0.31,0.97,"brNDC");
	textPubPID -> SetTextSize(0.055);
	textPubPID -> SetTextFont(42);
	textPubPID -> SetFillColor(0);	
	textPubPID -> AddText("p + p #rightarrow p + #pi^{+} + #pi^{-} + p");
	textPubPID -> AddText("#sqrt{s} = 510 GeV");
	textPubPID -> Draw("SAME");


	TLegend *legendPID = new TLegend(0.72,0.6,0.85,0.85,"","brNDC");
	legendPID -> AddEntry(hdEdxBasicPID[0], "#pi^{+} + #pi^{-}", "p");
	legendPID -> AddEntry(hdEdxBasicPID[1], "K^{+} + K^{-}", "p");
	legendPID -> AddEntry(hdEdxBasicPID[2], "p + #bar{p}", "p");
	legendPID -> SetTextFont(42);
	legendPID -> SetFillStyle(0);
	legendPID -> SetLineColor(0);
	legendPID -> SetTextSize(0.04);
	legendPID -> Draw("same");

	cdEdxBasicPID->Update();


	cdEdxLog = new TCanvas("dEdxLog","dEdxLog",1200,800);
	gPad->SetMargin(0.9,0.02,0.14,0.02);
	gStyle->SetPalette(1);
	gPad->SetLogy();
	hdEdxLog->SetStats(0);
	hdEdxLog->SetTitle(" ; dEdx; Number of events");
	hdEdxLog->GetXaxis()->SetTitleFont(42);
	hdEdxLog->GetYaxis()->SetTitleFont(42);
	hdEdxLog->GetXaxis()->SetLabelFont(42);
	hdEdxLog->GetYaxis()->SetLabelFont(42);
	hdEdxLog->GetXaxis()->SetTitleSize(.055);
	hdEdxLog->GetYaxis()->SetTitleSize(.055);
	hdEdxLog->GetXaxis()->SetTitleOffset(1);
	hdEdxLog->GetYaxis()->SetTitleOffset(0.8);
	hdEdxLog->Draw("same");

	TPaveText *textPub51 = new TPaveText(0.68,0.75,0.75,0.97,"brNDC");
	textPub51 -> SetTextSize(0.055);
	textPub51 -> SetTextFont(42);
	textPub51 -> SetFillColor(0);	
	textPub51 -> AddText("p + p #rightarrow p + #pi^{+} + #pi^{-} + p");
	textPub51 -> AddText("#sqrt{s} = 510 GeV");
	textPub51 -> Draw("same");

	cdEdxLog->Update();
}

void PID::Init(){

	for(int i = 0; i <3; ++i){
		hnSigma[i] = new TH2D("nSigma"+doubleParticleID[i],"nSigma "+doubleParticleID[i],100,0,50,100,0,50 );
		hnSigmaVsM2[i] = new TH2D("nSigma"+particleID[i]+"VsM2","nSigma "+particleID[i]+" Vs M2",100,0,50,100,0,1 );
		hChMomVsnSigma[i] = new TH2D("ChMomVsnSigma"+particleID[i],"charge momentum Vs nSigma"+particleID[i],80,-4,4,100,0,50 );
		hdEdxBasicPID[i] = new TH2D("dEdxBasicPID0"+particleID[i],"dEdx",90,-3,3,40,0,10 );
	}
			

	hdEdx = new TH2D("dEdx","dEdx",120,-3,3,80,0,10 );
	hnM2 = new TH1D("mSquared", "mSquared", 200, 0, 1.0);

	hdEdxLog = new TH1D("dEdxLog", "dEdxLog", 200, 0, 10);
	hdeltaDeltaTOF = new TH1D("deltaDeltaTOF", "deltaDeltaTOF", 1200, -10, 10);
	hdeltaTOF = new TH1D("deltaTOF", "deltaTOF", 200, -10, 10);
	hdeltaTOFexpected = new TH1D("deltaTOFexpected", "deltaTOFexpected", 200, -10, 10);
}

void PID::ConnectInput(){
	recTree = dynamic_cast<TTree*>( data->Get("recTree") );
	if (!recTree){
      cout<<"Error: cannot get recTree"<<endl;
      return;
   }

	recTree->SetBranchAddress("nSigPPion", &nSigPPion);
	recTree->SetBranchAddress("nSigPKaon", &nSigPKaon);
	recTree->SetBranchAddress("nSigPProton", &nSigPProton);
	recTree->SetBranchAddress("missingPt", &missingPt);
	recTree->SetBranchAddress("invMass", &invMass);
	recTree->SetBranchAddress("dEdx1", &dEdx1);
	recTree->SetBranchAddress("dEdx2", &dEdx2);
	recTree->SetBranchAddress("momentum1", &momentum1);
	recTree->SetBranchAddress("momentum2", &momentum2);
	recTree->SetBranchAddress("charge1", &charge1);
	recTree->SetBranchAddress("charge2", &charge2);
	recTree->SetBranchAddress("deltaDeltaTOF", &deltaDeltaTOF);
	recTree->SetBranchAddress("mSquared", &mSquared); 
	recTree->SetBranchAddress("deltaTOF", &deltaTOF);
	recTree->SetBranchAddress("TOFtime1", &TOFtime1);
	recTree->SetBranchAddress("TOFtime2", &TOFtime2);
	recTree->SetBranchAddress("TOFlength1", &TOFlength1);
	recTree->SetBranchAddress("TOFlength2", &TOFlength2);
	recTree->SetBranchAddress("elastic", &elastic);
}


