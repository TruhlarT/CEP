#include "PID.h"
#include "Plot.h"

//_____________________________________________________________________________
PID::PID(string innam, string outnam, bool text)
{
	//constructor
	input = innam;
	output = outnam;

	TEXT = text;
	cout << "PID::PID() called" << endl;

}//PID

//_____________________________________________________________________________
PID::~PID()
{
  //destructor

  cout << "PID::~PID() destructor called" << endl;

}//~PID


void PID::SetBichsel(TH1* hist, Int_t color, Float_t xMin, Float_t xMax, Int_t width){
	hist->Scale(1000000);
	hist->SetLineColor(color);	
	hist->SetLineWidth(width);
	hist->SetAxisRange(xMin,xMax);
	hist->Smooth();
}

void PID::DrawBichsel(){
	TFile* bichsel = new TFile("/home/truhlar/Desktop/STAR/CEP/Analysis/Data/bichselP.root" ,"r");	

	TH1D* hBichselPion = (TH1D*)bichsel -> Get("h_pi_70M");
	TH1D* hBichselKaon = (TH1D*)bichsel -> Get("h_k_70M");
	TH1D* hBichselProton = (TH1D*)bichsel -> Get("h_p_70M");
	TH1D* hBichselDeuteron = (TH1D*)bichsel -> Get("h_d_70M");

	TH1D* hBichselPionNegative = new TH1D("hBichselPionNegative", "h1", 1500, -14.5,0);
	TH1D* hBichselKaonNegative = new TH1D("hBichselKaonNegative", "h2", 1500, -14.5,0);
	TH1D* hBichselProtonNegative = new TH1D("hBichselProtonNegative", "h3", 1500, -14.5,0);
	TH1D* hBichselDeuteronNegative = new TH1D("hBichselDeuteronNegative", "h4", 1500, -14.5,0);
	for(int i = 1; i < 1500; i++){
		hBichselPionNegative->SetBinContent(i,hBichselPion->GetBinContent(1500-i));
		hBichselKaonNegative->SetBinContent(i,hBichselKaon->GetBinContent(1500-i));
		hBichselProtonNegative->SetBinContent(i,hBichselProton->GetBinContent(1500-i));
		hBichselDeuteronNegative->SetBinContent(i,hBichselDeuteron->GetBinContent(1500-i));
	}

	SetBichsel(hBichselPion,kBlack,0.09,3);
	hBichselPion->Draw("c SAME");

	SetBichsel(hBichselKaon,kRed,0.2,3);
	hBichselKaon->Draw("c SAME");

	SetBichsel(hBichselProton,kCyan,0.31,3);
	hBichselProton->Draw("c SAME");

	SetBichsel(hBichselDeuteron,kGreen,0.615,3);
	hBichselDeuteron->Draw("c SAME");

	SetBichsel(hBichselPionNegative,kBlack,-3,-0.095);
	hBichselPionNegative->Draw("c SAME");

	SetBichsel(hBichselKaonNegative,kRed,-3,-0.25);
	hBichselKaonNegative->Draw("c SAME");

	SetBichsel(hBichselProtonNegative,kCyan,-3,-0.315);
	hBichselProtonNegative->Draw("c SAME");

	SetBichsel(hBichselDeuteronNegative,kGreen,-3,-0.62);
   hBichselDeuteronNegative->Draw("c SAME");

   Plot tool;
   TLegend *leg = new TLegend(0.68,0.7,0.8,0.95);
   tool.SetLegendStyle(leg);
	leg->AddEntry(hBichselPion,"Pion","l");
	leg->AddEntry(hBichselKaon,"Kaon","l");
	leg->AddEntry(hBichselProton,"Proton","l");
	leg->AddEntry(hBichselDeuteron,"Deuteron","l");
	leg->Draw();

   fout->cd();
}

//_____________________________________________________________________________
void PID::PlotHistogram(){

	data = TFile::Open(input, "read");
	if (!data){
      cout<<"Error: cannot open "<<input<<endl;
      return;
   }

	TTree* tree = dynamic_cast<TTree*>( data->Get("recTree") );
	TTree* treeBack = dynamic_cast<TTree*>( data->Get("Background") );
	
	if (!tree || !treeBack){
      cout<<"Error: cannot open one of the TTree"<<endl;
      return;
   }

   cout<<"Creating output file: "<<output<<"StRP.root"<<endl;
	fout = new TFile(output +"StRP.root","UPDATE");

	TCanvas *cCanvas = new TCanvas("cCanvas","cCanvas",800,700);
	gPad->SetMargin(0.12,0.02,0.1,0.02); // (Float_t left, Float_t right, Float_t bottom, Float_t top)
	gPad->SetTickx();
	gPad->SetTicky();  
	gStyle->SetOptStat("");
	gStyle->SetPalette(1);
	gStyle->SetLineWidth(2);      //axis line
	gStyle->SetFrameLineWidth(2); //frame line

	tree->UseCurrentStyle();
	treeBack->UseCurrentStyle();

	Plot tool;

	TString variable;
// Plot mSquared 
	variable = "mSquared";
	tree->Draw(variable +">>" + variable +"Sig( 200, -0.5, 2.0)");
	tmpHist = (TH1F*)gPad->GetPrimitive(variable +"Sig");
	tmpHist->SetTitle(" ; m^{2}_{TOF} [GeV^{2}/c^{4}]; Number of events");
	//tmpHist->GetXaxis()->SetRangeUser(0,2.5);
	tool.SetGraphStyle(tmpHist);
	tool.SetMarkerStyle(tmpHist);
	tmpHist->Draw();
	tool.DrawText(tmpHist);
	tool.DrawTextStar(tmpHist);

	gPad->SetLogy();
	cCanvas->Update();
	cCanvas->SaveAs( output + "PID/" + variable + ".png");
	cCanvas->Write(variable);
	gPad->SetLogy(0);
//////////////////////////////////////////
	// Plot deltaTOF 
	variable = "deltaTOFexpected";
	tree->Draw(variable +">>" + variable +"Bcg(200, -10, 10)");
	tmpHist2 = (TH1F*)gPad->GetPrimitive(variable+"Bcg");
	tool.SetMarkerStyle(tmpHist2,2,20,1,2,1,1);
	variable = "deltaTOF";
	tree->Draw(variable +">>" + variable +"Sig(200, -10, 10)");
	tmpHist = (TH1F*)gPad->GetPrimitive(variable +"Sig");
	tmpHist->SetTitle(" ; #Delta TOF [ns]; Number of events");
	//tmpHist->GetXaxis()->SetRangeUser(0,2.5);
	tool.SetGraphStyle(tmpHist);
	tool.SetMarkerStyle(tmpHist);
	tmpHist->Draw("");
	tool.DrawText(tmpHist,0,true);
	tool.DrawTextStar(tmpHist);
	tmpHist2->Draw("SAME");

	TLegend* leg1 = new TLegend(0.72, 0.68, 0.9, 0.78);
	tool.SetLegendStyle(leg1);
	leg1 -> AddEntry(tmpHist, "Data", "l");
	leg1 -> AddEntry(tmpHist2, "#pi assumption", "fl");
	leg1->Draw("same");

	gPad->SetLogy();
	cCanvas->Update();
	cCanvas->SaveAs( output + "PID/" + variable + ".png");
	cCanvas->Write(variable);
	gPad->SetLogy(0);
//////////////////////////////////////////
	// Plot deltaDeltaTOF 
	variable = "deltaDeltaTOF";

	tree->Draw(variable +">>" + variable +"Sig(1200, -10, 10)");
	tmpHist = (TH1F*)gPad->GetPrimitive(variable +"Sig");
	tmpHist->SetTitle(" ; #Delta #Delta TOF [ns]; Number of events");
	//tmpHist->GetXaxis()->SetRangeUser(0,2.5);
	tool.SetGraphStyle(tmpHist);
	tool.SetMarkerStyle(tmpHist);
	tmpHist->Draw();
	tool.DrawText(tmpHist,0,true);
	tool.DrawTextStar(tmpHist);

	gPad->SetLogy();
	cCanvas->Update();
	cCanvas->SaveAs( output + "PID/" + variable + ".png");
	cCanvas->Write(variable);
	gPad->SetLogy(0);
//////////////////////////////////////////
// Plot nSigPPion 
	variable = "nSigPPion";

	tree->Draw(variable +">>" + variable +"Sig(50, 0, 7)");
	tmpHist = (TH1F*)gPad->GetPrimitive(variable +"Sig");
	tmpHist->SetTitle(" ; n#sigma^{pair}_{#pi} ; Number of events");
	//tmpHist->GetXaxis()->SetRangeUser(0,2.5);
	tool.SetGraphStyle(tmpHist);
	tool.SetMarkerStyle(tmpHist);
	tmpHist->GetYaxis()->SetTitleOffset(1.4);
	tmpHist->Draw();
	tool.DrawText(tmpHist,0,true,0.72,0.74,0.85,0.9);
	tool.DrawTextStar(tmpHist,2);

	TLine *left4 = new TLine(3,0,3,8000);
	tool.SetLineStyle(left4,10,1,4);
   left4->Draw("same");

	cCanvas->Update();
	cCanvas->SaveAs( output + "PID/" + variable + ".png");
	cCanvas->Write(variable);
//////////////////////////////////////////
//////////////////////////////////////////
	// Plot dEdx
	variable = "dEdx";
	tree->Draw(variable +"1>>" + variable +"Sig1(80,0,25)");
	tmpHist = (TH1F*)gPad->GetPrimitive(variable +"Sig1");
	tree->Draw(variable +"2>>" + variable +"Sig2(80,0,25)");
	tmpHist3 = (TH1F*)gPad->GetPrimitive(variable +"Sig2");
	tmpHist->Add(tmpHist3);
	tmpHist->SetTitle(" ; dEdx [keV/cm]; Number of tracks");;
	//tmpHist->GetXaxis()->SetRangeUser(0,2.5);
	tool.SetGraphStyle(tmpHist);
	tool.SetMarkerStyle(tmpHist);
	tmpHist->Draw();
	tool.DrawText(tmpHist);
	tool.DrawTextStar(tmpHist);

	gPad->SetLogy();
	cCanvas->Update();
	cCanvas->SaveAs( output + "PID/" + variable + ".png");
	cCanvas->Write(variable);
	gPad->SetLogy(0);

//////////////////////////////////////////
	TCanvas *cCanvas2D = new TCanvas("cCanvas2D","cCanvas2D",800,700);
	gPad->SetMargin(0.09,0.13,0.1,0.02); // (Float_t left, Float_t right, Float_t bottom, Float_t top)
	gStyle->SetPalette(1);
	gPad->SetTickx();
	gPad->SetTicky(); 
	TString variable2;
////////// Plot dEdx ////////////	
	variable = "dEdx";
	variable2 = "momentum";
	tree->Draw(variable + "1:charge1*tranMomenta1" + ">>" + variable + "Vs" + variable2 + "1Sig(120,-3,3,80,0,10)","","colz");
	tmp2DHist2 = (TH2F*)gPad->GetPrimitive(variable + "Vs" + variable2 + "1Sig");
	tree->Draw(variable + "2:charge2*tranMomenta2" + ">>" + variable + "Vs" + variable2 + "2Sig(120,-3,3,80,0,10)","","colz");
	tmp2DHist = (TH2F*)gPad->GetPrimitive(variable + "Vs" + variable2 + "2Sig");
	tmp2DHist->Add(tmp2DHist2);
	tmp2DHist->SetTitle(" ; #frac{q}{e} #times p_{T} [GeV/c] ;dE/dx [keV/cm]");
	tool.SetGraphStyle(tmp2DHist,4,20,1,4,1,1,0.9,0.8);
	tmp2DHist->Draw("colz");
	tool.DrawText(tmp2DHist,0,true,0.08,0.78,0.25,0.9,12);
	tool.DrawTextStar(tmp2DHist);

	cCanvas2D->Update();
	cCanvas2D->SaveAs( output + "PID/" + variable + "Vs" + variable2 + ".png");
	cCanvas2D->Write(variable + "Vs" + variable2);
 	////// Plot Bichsel////////////
	DrawBichsel();
	tool.DrawText(tmp2DHist,0,true,0.08,0.78,0.25,0.9,12);
	tool.DrawTextStar(tmp2DHist);
	cCanvas2D->Update();
	cCanvas2D->SaveAs( output + "PID/" + variable + "wBichsel" + ".png");
	cCanvas2D->Write(variable + "wBichsel");
/////////////////////////////////////////////////////
	TString particleID[] = {"Pion","Kaon","Proton"};

////////// Plot 2D histograms: ////////////
for(int i = 0; i < 3; ++i){
	int j1 = i;
	int j2 = i+1;
	if(i==2){
		j1 = 0;
		j2 = 2;
	}

	variable = "nSigP" + particleID[j1];
	variable2 = "nSigP" + particleID[j2];
	//tree->Draw(variable2 + ":" + variable + ">>" + variable + "Vs" + variable2 + "Sig(100,0,50,100,0,50)","deltaDeltaTOF < 1 && deltaDeltaTOF > -1 ","colz");
	tree->Draw(variable2 + ":" + variable + ">>" + variable + "Vs" + variable2 + "Sig(100,0,50,100,0,50)","","colz");
	tmp2DHist = (TH2F*)gPad->GetPrimitive(variable + "Vs" + variable2 + "Sig");
	tmp2DHist->SetTitle(" ; n#sigma^{pair}_{" + particleID[j1] + "}; n#sigma^{pair}_{" + particleID[j2] + "}");
	tool.SetGraphStyle(tmp2DHist,4,20,1,4,1,1,0.9,0.9);
	tmp2DHist->Draw("colz");
	if(i!=0){
		tool.DrawText(tmp2DHist,0,true,0.6,0.75,0.75,0.9);
		tool.DrawTextStar(tmp2DHist,1);
	}else{
		tool.DrawText(tmp2DHist,0,true,0.08,0.78,0.25,0.9,12);
		tool.DrawTextStar(tmp2DHist);	
	}

	cCanvas2D->Update();
	cCanvas2D->SaveAs( output + "PID/nSigP" + particleID[j1] + "Vs" + particleID[j2] + ".png");
	cCanvas2D->Write("nSigP" + particleID[j1] + "Vs" + particleID[j2]);

	variable = "nSigP" + particleID[i];
	variable2 = "mSquared";
	tree->Draw(variable2 + ":" + variable + ">>" + variable + "Vs" + variable2 + "Sig(100,0,50,100,0,1)","","colz");
	tmp2DHist = (TH2F*)gPad->GetPrimitive(variable + "Vs" + variable2 + "Sig");
	tmp2DHist->SetTitle(" ; n#sigma^{pair}_{" + particleID[i] + "}; m^{2}_{TOF} [GeV^{2}/c^{4}]");
	tool.SetGraphStyle(tmp2DHist,4,20,1,4,1,1,0.9,1.0);

	tmp2DHist->Draw("colz");
	tool.DrawText(tmp2DHist,0,true,0.6,0.75,0.75,0.9);
	tool.DrawTextStar(tmp2DHist,1);

	cCanvas2D->Update();
	cCanvas2D->SaveAs( output + "PID/" + variable + "Vs" + variable2 + ".png");
	cCanvas2D->Write(variable + "Vs" + variable2);

	variable = "momentum";
	variable2 = "nSigP" + particleID[i];
	tree->Draw(variable2 + ":" + "charge1*tranMomenta1" + ">>" + variable + "Vs" + variable2 + "Sig(80,-4,4,100,0,50)","","colz");
	tmp2DHist2 = (TH2F*)gPad->GetPrimitive(variable + "Vs" + variable2 + "Sig");
	tree->Draw(variable2 + ":" + "charge2*tranMomenta2" + ">>" + variable + "Vs" + variable2 + "2Sig(80,-4,4,100,0,50)","","colz");
	tmp2DHist = (TH2F*)gPad->GetPrimitive(variable + "Vs" + variable2 + "2Sig");
	tmp2DHist->Add(tmp2DHist2);
	tmp2DHist->SetTitle(" ; #frac{q}{e} #times p_{T} [GeV/c] ; n#sigma^{pair}_{" + particleID[i] + "}");
	tool.SetGraphStyle(tmp2DHist,4,20,1,4,1,1,0.9,0.9);
	tmp2DHist->Draw("colz");
	tool.DrawText(tmp2DHist,0,true,0.08,0.78,0.25,0.9,12);
	tool.DrawTextStar(tmp2DHist);

	cCanvas2D->Update();
	cCanvas2D->SaveAs( output + "PID/" + variable + "Vs" + variable2 + ".png");
	cCanvas2D->Write(variable + "Vs" + variable2);

}	
/////////////////////////////////////////////////////	

////////// Plot dEdx with basic PID ////////////	
	cCanvas->cd();
	variable = "dEdx";
	variable2 = "momentum";
	tree->Draw(variable + "1:charge1*tranMomenta1" + ">>" + variable + "Vs" + variable2 + "1Sig(120,-3,3,80,0,10)","nSigPProton < 3 && nSigPKaon > 3 && nSigPPion > 3","colz");
	tmp2DHist2 = (TH2F*)gPad->GetPrimitive(variable + "Vs" + variable2 + "1Sig");
	tree->Draw(variable + "2:charge2*tranMomenta2" + ">>" + variable + "Vs" + variable2 + "2Sig(120,-3,3,80,0,10)","nSigPProton < 3 && nSigPKaon > 3 && nSigPPion > 3","colz");
	tmp2DHist3 = (TH2F*)gPad->GetPrimitive(variable + "Vs" + variable2 + "2Sig");
	tmp2DHist3->Add(tmp2DHist2);
	tool.SetMarkerStyle(tmp2DHist3,8,29,1,8,1,1);
	tree->Draw(variable + "1:charge1*tranMomenta1" + ">>" + variable + "Vs" + variable2 + "1bSig(120,-3,3,80,0,10)","nSigPKaon < 3 && nSigPProton > 3 && nSigPPion > 3","colz");
	tmp2DHist4 = (TH2F*)gPad->GetPrimitive(variable + "Vs" + variable2 + "1bSig");
	tree->Draw(variable + "2:charge2*tranMomenta2" + ">>" + variable + "Vs" + variable2 + "2bSig(120,-3,3,80,0,10)","nSigPKaon < 3 && nSigPProton > 3 && nSigPPion > 3","colz");
	tmp2DHist5 = (TH2F*)gPad->GetPrimitive(variable + "Vs" + variable2 + "2bSig");
	tmp2DHist5->Add(tmp2DHist4);
	tool.SetMarkerStyle(tmp2DHist5,2,29,1,8,1,1);
	tree->Draw(variable + "1:charge1*tranMomenta1" + ">>" + variable + "Vs" + variable2 + "1cSig(120,-3,3,80,0,10)","(nSigPProton >= 3 || nSigPKaon <= 3 || nSigPPion <= 3) && (nSigPKaon >= 3 || nSigPProton >= 3 || nSigPPion <= 3) ","colz");
	tmp2DHist6 = (TH2F*)gPad->GetPrimitive(variable + "Vs" + variable2 + "1cSig");
	tree->Draw(variable + "2:charge2*tranMomenta2" + ">>" + variable + "Vs" + variable2 + "2cSig(120,-3,3,80,0,10)","(nSigPProton >= 3 || nSigPKaon <= 3 || nSigPPion <= 3) && (nSigPKaon >= 3 || nSigPProton >= 3 || nSigPPion <= 3) ","colz");
	tmp2DHist = (TH2F*)gPad->GetPrimitive(variable + "Vs" + variable2 + "2cSig");
	tmp2DHist->Add(tmp2DHist6);
	tmp2DHist->SetTitle(" ; #frac{q}{e} #times p_{T} [GeV/c] ;dE/dx [keV/cm]");
	tool.SetGraphStyle(tmp2DHist,4,20,1,4,1,1,0.9,1.0);
	tool.SetMarkerStyle(tmp2DHist,4,29,1,4,1,1);;
	tmp2DHist->Draw("SCAT");
	tool.DrawText(tmp2DHist);
	tool.DrawTextStar(tmp2DHist);
	tmp2DHist3->Draw("SAME");
	tmp2DHist5->Draw("SAME");

	TLegend *legendPID = new TLegend(0.72,0.55,0.85,0.80,"","brNDC");
	tool.SetLegendStyle(legendPID);
	legendPID -> AddEntry(tmp2DHist, "#pi^{+} + #pi^{-}", "p");
	legendPID -> AddEntry(tmp2DHist5, "K^{+} + K^{-}", "p");
	legendPID -> AddEntry(tmp2DHist3, "p + #bar{p}", "p");
	legendPID -> Draw("same");

	cCanvas->Update();
	cCanvas->SaveAs( output + "PID/" + variable + "Vs" + variable2 + "WithPID.png");
	cCanvas->Write(variable + "Vs" + variable2 + "WithPID");
/////////////////////////////////////////////////////
	cCanvas->Close();
	cCanvas2D->Close();

   fout->Close();

}


