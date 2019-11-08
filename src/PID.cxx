#include "PID.h"
#include "Plot.h"

//_____________________________________________________________________________
PID::PID(TFile* dataInput, TFile* fileOut, TString outnam, bool text, TString inputCuts)
{
	//constructor
	output = outnam;

	data = dataInput;
	fout = fileOut;

	cuts = inputCuts;
	if(inputCuts != "")
		cutsWithPrefix = " && " + inputCuts;
	else
		cutsWithPrefix="";

	TEXT = text;
	cout << "PID::PID() called" << endl;

}//PID

//_____________________________________________________________________________
PID::~PID()
{
  //destructor

  cout << "PID::~PID() destructor called" << endl;

}//~PID


void PID::SetBichsel(TH1* hist, Int_t color, Float_t xMin, Float_t xMax, Float_t width){
	hist->Scale(1000000,"nosw2");
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


	SetBichsel(hBichselPion,kBlack,0.085,3);
	hBichselPion->Draw("c SAME");

	SetBichsel(hBichselKaon,kRed,0.27,3);
	hBichselKaon->Draw("c SAME");

	SetBichsel(hBichselProton,kCyan,0.5,3);
	hBichselProton->Draw("c SAME");

	SetBichsel(hBichselDeuteron,kGreen,.97,3);
	hBichselDeuteron->Draw("c SAME");

	SetBichsel(hBichselPionNegative,kBlack,-3,-0.095);
	hBichselPionNegative->Draw("c SAME");

	SetBichsel(hBichselKaonNegative,kRed,-3,-0.3);
	hBichselKaonNegative->Draw("c SAME");

	SetBichsel(hBichselProtonNegative,kCyan,-3,-0.52);
	hBichselProtonNegative->Draw("c SAME");

	SetBichsel(hBichselDeuteronNegative,kGreen,-3,-1.1);
	hBichselDeuteronNegative->Draw("c SAME");

	Plot tool;
	TLegend *leg = new TLegend(0.68,0.7,0.8,0.95);
	tool.SetLegendStyle(leg);
	leg->AddEntry(hBichselPion,"Pion","l");
	leg->AddEntry(hBichselKaon,"Kaon","l");
	leg->AddEntry(hBichselProton,"Proton","l");
	leg->AddEntry(hBichselDeuteron,"Deuteron","l");
	leg->Draw();

	//bichsel->Close();
}

//_____________________________________________________________________________
void PID::PlotHistogram(){

	if (!data){
		cout<<"Error: cannot open input file"<<endl;
		return;
    }

	TTree* tree = dynamic_cast<TTree*>( data->Get("recTree") );
	TTree* treeBack = dynamic_cast<TTree*>( data->Get("Background") );
	
	if (!tree || !treeBack){
		cout<<"Error: cannot open one of the TTree"<<endl;
		return;
    }

	if (!fout){
		cout<<"Error: cannot open output file"<<endl;
		return;
    }


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

	tree->Draw(variable +">>" + variable + "Pion( 200, -0.5, 1.5)", "nSigPairPion < 3 " + cutsWithPrefix);
	tmpHist2 = (TH1F*)gPad->GetPrimitive(variable + "Pion");
	tool.SetMarkerStyle(tmpHist2,2,20,1,3,1,1);
	tree->Draw(variable +">>" + variable + "Kaon( 200, -0.5, 1.5)", "nSigPairKaon < 3 && nSigPairPion > 3 && nSigPairProton > 3 " + cutsWithPrefix);
	tmpHist3 = (TH1F*)gPad->GetPrimitive(variable + "Kaon");
	tool.SetMarkerStyle(tmpHist3,2,20,1,2,1,1);
	tree->Draw(variable +">>" + variable + "Proton( 200, -0.5, 1.5)", "nSigPairProton < 3 && nSigPairKaon > 3 && nSigPairPion > 3" + cutsWithPrefix);
	tmpHist4 = (TH1F*)gPad->GetPrimitive(variable + "Proton");
	tool.SetMarkerStyle(tmpHist4,2,20,1,1,1,1);

	tree->Draw(variable +">>" + variable +"Sig( 200, -0.5, 1.5)",cuts);
	tmpHist = (TH1F*)gPad->GetPrimitive(variable +"Sig");
	tmpHist->SetTitle(" ; m^{2}_{TOF} [GeV^{2}/c^{4}]; Number of events");
	//tmpHist->GetXaxis()->SetRangeUser(0,2.5);
	tool.SetGraphStyle(tmpHist);
	tool.SetMarkerStyle(tmpHist);
	tmpHist->Draw();
	//tool.DrawText(tmpHist);
	//tool.DrawTextStar(tmpHist);
    tmpHist2->Draw("SAME");
	tmpHist3->Draw("SAME");
    tmpHist4->Draw("SAME");

	TLine *mLine = new TLine(0.2,0,0.2,8000);
	tool.SetLineStyle(mLine,10,1,4);
    mLine->Draw("same");

    mLine = new TLine(0.32,0,0.32,8000);
	tool.SetLineStyle(mLine,10,1,4);
    mLine->Draw("same");

    mLine = new TLine(0.7,0,0.7,8000);
	tool.SetLineStyle(mLine,10,1,4);
    mLine->Draw("same");

    mLine = new TLine(1.1,0,1.1,8000);
	tool.SetLineStyle(mLine,10,1,4);
    mLine->Draw("same");

	TPaveText *textPub = new TPaveText(0.72,0.8,0.84,0.95,"brNDC");
	textPub -> SetTextSize(0.04);
	textPub -> SetFillColor(0);
	textPub -> SetTextFont(42);
	textPub -> AddText("Pair PID based on dE/dx");
	textPub -> Draw("same");

	leg1 = new TLegend(0.72, 0.65, 0.9, 0.8);
	tool.SetLegendStyle(leg1);
	leg1 -> AddEntry(tmpHist, "All pairs", "l");
    leg1 -> AddEntry(tmpHist2, "#pi^{+} #pi^{-}", "fl");
	leg1 -> AddEntry(tmpHist3, "K^{+}K^{-}", "fl");
    leg1 -> AddEntry(tmpHist4, "p#bar{p}", "fl");
	leg1->Draw("same");

	gPad->SetLogy();
	cCanvas->Update();
	//cCanvas->SaveAs( output + "PID/" + variable + ".png");
	cCanvas->Write(variable);
	gPad->SetLogy(0);
//////////////////////////////////////////
	// Plot deltaTOF kaon
	variable = "deltaTOFExpectedKaon";
	tree->Draw(variable +">>" + variable +"Assump(200, -10, 10)",cuts);
	tmpHist2 = (TH1F*)gPad->GetPrimitive(variable+"Assump");
	tool.SetMarkerStyle(tmpHist2,2,20,1,2,1,1);

    variable = "deltaTOFExpectedPion";
    tree->Draw(variable +">>" + variable +"Assump(200, -10, 10)",cuts);
    tmpHist3 = (TH1F*)gPad->GetPrimitive(variable+"Assump");
    tool.SetMarkerStyle(tmpHist3,2,20,1,3,1,1);

    variable = "deltaTOFExpectedProton";
    tree->Draw(variable +">>" + variable +"Assump(200, -10, 10)",cuts);
    tmpHist4 = (TH1F*)gPad->GetPrimitive(variable+"Assump");
    tool.SetMarkerStyle(tmpHist4,2,20,1,1,1,1);

	variable = "deltaTOF";
	tree->Draw(variable +">>" + variable +"Sig(200, -10, 10)",cuts);
	tmpHist = (TH1F*)gPad->GetPrimitive(variable +"Sig");
	tmpHist->SetTitle(" ; #Delta TOF [ns]; Number of events");
	//tmpHist->GetXaxis()->SetRangeUser(0,2.5);
	tool.SetGraphStyle(tmpHist);
	tool.SetMarkerStyle(tmpHist);
	tmpHist->Draw("");
	tool.DrawText(tmpHist,0,true);
	tool.DrawTextStar(tmpHist);
	tmpHist3->Draw("SAME");
    tmpHist2->Draw("SAME");
    tmpHist4->Draw("SAME");

	leg1 = new TLegend(0.72, 0.58, 0.9, 0.78);
	tool.SetLegendStyle(leg1);
	leg1 -> AddEntry(tmpHist, "Data", "l");
    leg1 -> AddEntry(tmpHist2, "#pi assumption", "fl");
	leg1 -> AddEntry(tmpHist3, "K assumption", "fl");
    leg1 -> AddEntry(tmpHist4, "p assumption", "fl");
	leg1->Draw("same");

	gPad->SetLogy();
	cCanvas->Update();
	//cCanvas->SaveAs( output + "PID/" + variable + ".png");
	cCanvas->Write(variable);
	gPad->SetLogy(0);	
//////////////////////////////////////////
	// Plot deltaDeltaTOFPion 
	variable = "deltaDeltaTOFPion";

	tree->Draw(variable +">>" + variable +"Sig(1200, -10, 10)",cuts);
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
	//cCanvas->SaveAs( output + "PID/" + variable + ".png");
	cCanvas->Write(variable);
	gPad->SetLogy(0);
//////////////////////////////////////////
	// Plot deltaDeltaTOFPion 
	variable = "deltaDeltaTOFKaon";

	tree->Draw(variable +">>" + variable +"Sig(1200, -10, 10)",cuts);
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
	//cCanvas->SaveAs( output + "PID/" + variable + ".png");
	cCanvas->Write(variable);
	gPad->SetLogy(0);
//////////////////////////////////////////
	// Plot deltaDeltaTOFPion 
	variable = "deltaDeltaTOFProton";

	tree->Draw(variable +">>" + variable +"Sig(1200, -10, 10)",cuts);
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
	//cCanvas->SaveAs( output + "PID/" + variable + ".png");
	cCanvas->Write(variable);
	gPad->SetLogy(0);
//////////////////////////////////////////
// Plot nSigPairPion 
	variable = "nSigPairPion";

	tree->Draw(variable +">>" + variable +"Sig(50, 0, 7)",cuts);
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
	//cCanvas->SaveAs( output + "PID/" + variable + ".png");
	cCanvas->Write(variable);
//////////////////////////////////////////
//////////////////////////////////////////
	// Plot dEdx
	variable = "dEdx";
	tree->Draw(variable +"1>>" + variable +"Sig1(80,0,25)",cuts);
	tmpHist = (TH1F*)gPad->GetPrimitive(variable +"Sig1");
	tree->Draw(variable +"0>>" + variable +"Sig2(80,0,25)",cuts);
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
	//cCanvas->SaveAs( output + "PID/" + variable + ".png");
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
	tree->Draw(variable + "1:charge1*transMomentum1" + ">>" + variable + "Vs" + variable2 + "1Sig(120,-3,3,80,0,10)",cuts,"colz");
	tmp2DHist2 = (TH2F*)gPad->GetPrimitive(variable + "Vs" + variable2 + "1Sig");
	tree->Draw(variable + "0:charge0*transMomentum0" + ">>" + variable + "Vs" + variable2 + "2Sig(120,-3,3,80,0,10)",cuts,"colz");
	tmp2DHist = (TH2F*)gPad->GetPrimitive(variable + "Vs" + variable2 + "2Sig");
	tmp2DHist->Add(tmp2DHist2);
	tmp2DHist->SetTitle(" ; #frac{q}{e} #times p_{T} [GeV/c] ;dE/dx [keV/cm]");
	tool.SetGraphStyle(tmp2DHist,4,20,1,4,1,1,0.9,0.8);
	tmp2DHist->Draw("colz");
	tool.DrawText(tmp2DHist,0,true,0.08,0.78,0.25,0.9,12);
	tool.DrawTextStar(tmp2DHist);

	cCanvas2D->Update();
	//cCanvas2D->SaveAs( output + "PID/" + variable + "Vs" + variable2 + ".png");
	cCanvas2D->Write(variable + "Vs" + variable2);
 	

 	////// Plot Bichsel////////////
	TDirectory* currentDir = TDirectory::CurrentDirectory();
	DrawBichsel();
	tool.DrawText(tmp2DHist,0,true,0.08,0.78,0.25,0.9,12);
	tool.DrawTextStar(tmp2DHist);
	currentDir->cd();

	cCanvas2D->Update();
	//cCanvas2D->SaveAs( output + "PID/" + variable + "wBichsel" + ".png");
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

	variable = "nSigPair" + particleID[j2];
	variable2 = "nSigPair" + particleID[j1];
	//tree->Draw(variable2 + ":" + variable + ">>" + variable + "Vs" + variable2 + "Sig(100,0,50,100,0,50)","deltaDeltaTOF < 1 && deltaDeltaTOF > -1 ","colz");
	tree->Draw(variable2 + ":" + variable + ">>" + variable + "Vs" + variable2 + "Sig(100,0,35,100,0,35)",cuts,"colz");
	tmp2DHist = (TH2F*)gPad->GetPrimitive(variable + "Vs" + variable2 + "Sig");
	tmp2DHist->SetTitle(" ; n#sigma^{pair}_{" + particleID[j2] + "}; n#sigma^{pair}_{" + particleID[j1] + "}");
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
	//cCanvas2D->SaveAs( output + "PID/nSigPair" + particleID[j1] + "Vs" + particleID[j2] + ".png");
	cCanvas2D->Write("nSigPair" + particleID[j2] + "Vs" + particleID[j1]);

	variable = "nSigPair" + particleID[i];
	variable2 = "mSquared";
	tree->Draw(variable2 + ":" + variable + ">>" + variable + "Vs" + variable2 + "Sig(100,0,16,100,-0.5,1.5)",cuts,"colz");
	tmp2DHist = (TH2F*)gPad->GetPrimitive(variable + "Vs" + variable2 + "Sig");
	tmp2DHist->SetTitle(" ; n#sigma^{pair}_{" + particleID[i] + "}; m^{2}_{TOF} [GeV^{2}/c^{4}]");
	tool.SetGraphStyle(tmp2DHist,4,20,1,4,1,1,0.9,1.0);

	tmp2DHist->Draw("colz");
	tool.DrawText(tmp2DHist,0,true,0.6,0.75,0.75,0.9);
	tool.DrawTextStar(tmp2DHist,1);

	cCanvas2D->Update();
	//cCanvas2D->SaveAs( output + "PID/" + variable + "Vs" + variable2 + ".png");
	cCanvas2D->Write(variable2 + "Vs" + variable);

	variable = "momentum";
	variable2 = "nSigPair" + particleID[i];
	tree->Draw(variable2 + ":" + "charge1*transMomentum1" + ">>" + variable + "Vs" + variable2 + "Sig(80,-4,4,100,0,50)",cuts,"colz");
	tmp2DHist2 = (TH2F*)gPad->GetPrimitive(variable + "Vs" + variable2 + "Sig");
	tree->Draw(variable2 + ":" + "charge0*transMomentum0" + ">>" + variable + "Vs" + variable2 + "2Sig(80,-4,4,100,0,50)",cuts,"colz");
	tmp2DHist = (TH2F*)gPad->GetPrimitive(variable + "Vs" + variable2 + "2Sig");
	tmp2DHist->Add(tmp2DHist2);
	tmp2DHist->SetTitle(" ; #frac{q}{e} #times p_{T} [GeV/c] ; n#sigma^{pair}_{" + particleID[i] + "}");
	tool.SetGraphStyle(tmp2DHist,4,20,1,4,1,1,0.9,0.9);
	tmp2DHist->Draw("colz");
	tool.DrawText(tmp2DHist,0,true,0.08,0.78,0.25,0.9,12);
	tool.DrawTextStar(tmp2DHist);

	cCanvas2D->Update();
	//cCanvas2D->SaveAs( output + "PID/" + variable + "Vs" + variable2 + ".png");
	cCanvas2D->Write(variable2 + "Vs" + variable);

}	
/////////////////////////////////////////////////////	

////////// Plot dEdx with basic PID ////////////	
	cCanvas->cd();
	variable = "dEdx";
	variable2 = "momentum";
	tree->Draw(variable + "1:charge1*transMomentum1" + ">>" + variable + "Vs" + variable2 + "1Sig(120,-3,3,80,0,10)","nSigPairProton < 3" + cutsWithPrefix,"colz");
	tmp2DHist2 = (TH2F*)gPad->GetPrimitive(variable + "Vs" + variable2 + "1Sig");
	tree->Draw(variable + "0:charge0*transMomentum0" + ">>" + variable + "Vs" + variable2 + "2Sig(120,-3,3,80,0,10)","nSigPairProton < 3" + cutsWithPrefix,"colz");
	tmp2DHist3 = (TH2F*)gPad->GetPrimitive(variable + "Vs" + variable2 + "2Sig");
	tmp2DHist3->Add(tmp2DHist2);
	tool.SetMarkerStyle(tmp2DHist3,8,29,1,8,1,1);
	tree->Draw(variable + "1:charge1*transMomentum1" + ">>" + variable + "Vs" + variable2 + "1bSig(120,-3,3,80,0,10)","nSigPairKaon < 3" + cutsWithPrefix,"colz");
	tmp2DHist4 = (TH2F*)gPad->GetPrimitive(variable + "Vs" + variable2 + "1bSig");
	tree->Draw(variable + "0:charge0*transMomentum0" + ">>" + variable + "Vs" + variable2 + "2bSig(120,-3,3,80,0,10)","nSigPairKaon < 3 " + cutsWithPrefix,"colz");
	tmp2DHist5 = (TH2F*)gPad->GetPrimitive(variable + "Vs" + variable2 + "2bSig");
	tmp2DHist5->Add(tmp2DHist4);
	tool.SetMarkerStyle(tmp2DHist5,2,29,1,8,1,1);
	tree->Draw(variable + "1:charge1*transMomentum1" + ">>" + variable + "Vs" + variable2 + "1cSig(120,-3,3,80,0,10)","nSigPairPion < 3" + cutsWithPrefix,"colz");
	tmp2DHist6 = (TH2F*)gPad->GetPrimitive(variable + "Vs" + variable2 + "1cSig");
	tree->Draw(variable + "0:charge0*transMomentum0" + ">>" + variable + "Vs" + variable2 + "2cSig(120,-3,3,80,0,10)"," nSigPairPion < 3" + cutsWithPrefix,"colz");
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
	//cCanvas->SaveAs( output + "PID/" + variable + "Vs" + variable2 + "WithPID.png");
	cCanvas->Write(variable + "Vs" + variable2 + "WithPID");
/////////////////////////////////////////////////////
	cCanvas->Close();
	cCanvas2D->Close();

}


