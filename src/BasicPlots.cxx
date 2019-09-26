#include "BasicPlots.h"
#include "Plot.h"


//_____________________________________________________________________________
BasicPlots::BasicPlots(string innam, string outnam, bool text)
{
	//constructor
	TEXT = text;
	bcground = true;
	input = innam;
	output = outnam;

	siz = 0.045;
	font = 42;

	cout << "BasicPlots::BasicPlots() called" << endl;

}//BasicPlots

//_____________________________________________________________________________
BasicPlots::~BasicPlots()
{
  //destructor

  cout << "BasicPlots::~BasicPlots() destructor called" << endl;

}//~BasicPlots


//_____________________________________________________________________________
void BasicPlots::PlotHistogram() {

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
	fout = new TFile(output +"StRP.root","RECREATE");

	TCanvas *cCanvas = new TCanvas("cCanvas","cCanvas",800,700);
	gPad->SetMargin(0.9,0.02,0.105,0.02); // (Float_t left, Float_t right, Float_t bottom, Float_t top)
	gPad->SetTickx();
	gPad->SetTicky();  
	gStyle->SetOptStat("");
	gStyle->SetPalette(1);
	gStyle->SetLineWidth(2);      //axis line
	gStyle->SetFrameLineWidth(2); //frame line
	//gStyle->SetErrorX(0); // suppres X errors

	tree->UseCurrentStyle();
	treeBack->UseCurrentStyle();

	Plot tool;

// Plot Inv Mass Inelastic + Elastic
	treeBack->Draw("invMass>>invMassBackground(50,0,2.5)","nSigPPion<3");
	TH1F *tmpHist2 = (TH1F*)gPad->GetPrimitive("invMassBackground");
	tool.SetMarkerStyle(tmpHist2,2,20,1,2,1,1);

	tree->Draw("invMass>>invMassSignal(50,0,2.5)","nSigPPion<3");
	TH1F *tmpHist = (TH1F*)gPad->GetPrimitive("invMassSignal");
	tmpHist->SetTitle(" ; m(#pi^{+}#pi^{-}) [GeV/c^{2}]; Number of events");
	//tmpHist->GetXaxis()->SetRangeUser(0,2.5);
	tool.SetGraphStyle(tmpHist);
	tool.SetMarkerStyle(tmpHist);
	tmpHist->Draw("E");
	tool.DrawText(tmpHist, true);
	tmpHist2->Draw("ESAME");

	TLegend* leg1 = new TLegend(0.58, 0.7, 0.78, 0.82);
	tool.SetLegendStyle(leg1);
	leg1->AddEntry(tmpHist,"In+El (unlike-sign pairs)","p");
	leg1->AddEntry(tmpHist2,"In+El (like-sign pairs)","p");
	leg1->Draw("same");


	if(TEXT){
		TPaveText *textPub3 = new TPaveText(0.6,0.6,0.7,0.7,"brNDC");
		tool.SetTextStyle(textPub3);
		textPub3 -> AddText("Drop at about 1 GeV");
		textPub3 -> Draw("same");

		TPaveText *textPub4 = new TPaveText(0.75,0.5,0.8,0.6,"brNDC");
		tool.SetTextStyle(textPub4);
		textPub4 -> AddText("Peak at about 1270 MeV");
		textPub4 -> Draw("same");	
	}


	cCanvas->Update();
	cCanvas->SaveAs(output + "BasicPlots/invMass.png");
	cCanvas->Write("invMass");
//////////////////////////////////////////
// Plot Inv Mass Inelastic
	treeBack->Draw("invMass>>invMassBackground(50,0,2.5)","nSigPPion<3 && !elastic");
	tmpHist2 = (TH1F*)gPad->GetPrimitive("invMassBackground");
	tool.SetMarkerStyle(tmpHist2,2,20,1,2,1,1);

	tree->Draw("invMass>>invMassSignal(50,0,2.5)","nSigPPion<3 && !elastic");
	tmpHist = (TH1F*)gPad->GetPrimitive("invMassSignal");
	tmpHist->SetTitle(" ; m(#pi^{+}#pi^{-}) [GeV/c^{2}];Number of track pairs");
	tmpHist->GetXaxis()->SetRangeUser(0,2.5);
	tool.SetGraphStyle(tmpHist);
	tool.SetMarkerStyle(tmpHist);
	tmpHist->Draw("E");
	tool.DrawText(tmpHist, true);
	tmpHist2->Draw("ESAME");

	leg1 = new TLegend(0.58, 0.7, 0.78, 0.82);
	tool.SetLegendStyle(leg1);
	leg1->AddEntry(tmpHist,"In (unlike-sign pairs)","p");
	leg1->AddEntry(tmpHist2,"In (like-sign pairs)","p");
	leg1->Draw("same");

	cCanvas->Update();
	cCanvas->SaveAs(output + "BasicPlots/invMassIn.png");
	cCanvas->Write("invMassIn");
//////////////////////////////////////////
// Plot Inv Mass Elastic
	treeBack->Draw("invMass>>invMassBackground(50,0,2.5)","nSigPPion<3 && elastic");
	tmpHist2 = (TH1F*)gPad->GetPrimitive("invMassBackground");
	tool.SetMarkerStyle(tmpHist2,2,20,1,2,1,1);

	tree->Draw("invMass>>invMassSignal(50,0,2.5)","nSigPPion<3 && elastic");
	tmpHist = (TH1F*)gPad->GetPrimitive("invMassSignal");
	tmpHist->SetTitle(" ; m(#pi^{+}#pi^{-}) [GeV/c^{2}];Number of track pairs");
	tmpHist->GetXaxis()->SetRangeUser(0,2.5);
	tool.SetGraphStyle(tmpHist);
	tool.SetMarkerStyle(tmpHist);
	tmpHist->Draw("E");
	tool.DrawText(tmpHist, true);
	tmpHist2->Draw("ESAME");

	leg1 = new TLegend(0.58, 0.7, 0.78, 0.82);
	tool.SetLegendStyle(leg1);
	leg1->AddEntry(tmpHist,"El (unlike-sign pairs)","p");
	leg1->AddEntry(tmpHist2,"El (like-sign pairs)","p");
	leg1->Draw("same");

	cCanvas->Update();
	cCanvas->SaveAs(output + "BasicPlots/invMassEl.png");
	cCanvas->Write("invMassEl");
//////////////////////////////////////////

	hCuts = (TH1F*)data -> Get("AnalysisFlow");
	hMissingPtTPC = (TH1D*)data -> Get("All/MissingPt_TPC2t_Combi");
	hMissingPtTOF = (TH1D*)data -> Get("All/MissingPt_TOF2trk_Combi");
	hMissingPtQ0 = (TH1D*)data -> Get("All/MissingPt_Q0_Combi");
	hMissingPtExc = (TH1D*)data -> Get("All/MissingPt_Excl_Combi");

// //////////////////////////////////////////////////////////
// Plot Cuts Flow
	TString Labels[] = { TString("All"), TString("CPT trigger"), TString("El+InEl"), TString("2 TPC-TOF tracks"), 
	                  TString("Same vertex"), TString("TotCharge 0"), TString("p_{T}^{miss} < 0.1 GeV/c"), TString(""),};
	//gPad->SetMargin(0.9,0.02,0.1,0.02); // (Float_t left, Float_t right, Float_t bottom, Float_t top)
	gPad->SetLogy();
	gPad->SetTickx(0);
	gPad->SetTicky(0);

	TH1F* hCutsFlow = new TH1F("cuts", "cuts", 8, 1, 9);
	int i = 1;
	for(int tb=1; tb<8; ++tb){ 
		if(tb==3){
			hCutsFlow->GetXaxis()->SetBinLabel(tb, Labels[tb-1]);
			hCutsFlow->SetBinContent(tb,hCuts->GetBinContent(i)+hCuts->GetBinContent((++i)++));
			cout<<hCutsFlow->GetBinContent(tb)<<endl;
			continue;
		}
		hCutsFlow->GetXaxis()->SetBinLabel(tb, Labels[tb-1]);
		hCutsFlow->SetBinContent(tb,hCuts->GetBinContent(i));
		cout<<hCutsFlow->GetBinContent(tb)<<endl;
		i++;
	}
	hCutsFlow->GetXaxis()->SetBinLabel(8, "PID");
	tree->Draw("nSigPPion>>nSigPPion","nSigPPion<3");
	hCutsFlow->SetBinContent(8,((TH1F*)gPad->GetPrimitive("nSigPPion"))->GetEntries());
	hCutsFlow->SetTitle("; ; Number of events");
	tool.SetGraphStyle(hCutsFlow);
	tool.SetMarkerStyle(hCutsFlow);
	hCutsFlow->GetYaxis()->SetRangeUser(800,300000000);
	hCutsFlow->GetXaxis()->SetLabelSize(0.045);
	hCutsFlow->Draw();

	TPaveText *textPub = new TPaveText(0.75,0.8,0.85,0.96,"brNDC");
	tool.SetTextStyle(textPub);
	textPub -> AddText("p + p #rightarrow p + #pi^{+} + #pi^{-} + p");
	textPub -> AddText("#sqrt{s} = 510 GeV");
	textPub -> AddText("Cuts flow");
	textPub -> Draw("same");

	TPaveText *textSTAR = new TPaveText(0.15,0.91,0.2,0.97,"brNDC"); //for text "star"
	//TPaveText *textSTAR = new TPaveText(0.20,0.92,0.27,0.97,"brNDC");
	textSTAR -> SetTextSize(0.045);
	textSTAR -> SetFillColor(0);
	textSTAR -> SetTextFont(62);
	textSTAR -> AddText("STAR");
	//textSTAR -> AddText("THIS THESIS");
	textSTAR -> Draw("same");
	
	cCanvas->Update();
	cCanvas->SaveAs(output + "BasicPlots/Cuts.png");
	cCanvas->Write("CutsFlow");
//////////////////////////////////////////
////////////// 
	//TCanvas *cMissingPt = new TCanvas("cMissingPt","cMissingPt",1200,800);
	TCanvas *cMissingPt = new TCanvas("cMissingPt","cMissingPt",800,700);
	gPad->SetMargin(0.9,0.02,0.1,0.02); // (Float_t left, Float_t right, Float_t bottom, Float_t top)
	//gPad->SetMargin(0.9,0.02,0.14,0.02);
	gPad->SetLogy();
	hMissingPtTPC->SetStats(0);
	hMissingPtTPC->SetTitle(" ; p_{T}^{miss} [GeV/c];Number of events");
	tool.SetGraphStyle(hMissingPtTPC);
	//tool.SetMarkerStyle(hMissingPtTPC);
	hMissingPtTPC->SetMinimum(10);

	hMissingPtTPC->Draw("same");	
	hMissingPtQ0->SetFillColor(4);
	hMissingPtQ0->SetLineColor(4);
	hMissingPtQ0->SetFillStyle(1001);
	hMissingPtQ0->Draw("same");
	hMissingPtExc->SetFillColorAlpha(2, 0.5);
	hMissingPtExc->SetLineColor(2);
	//hMissingPtExc->SetFillColor(2);
	hMissingPtExc->SetFillStyle(1001);
	hMissingPtExc->Draw("same");
	cMissingPt->cd();


	leg1 = new TLegend(0.69,0.62,0.81,0.82);
	tool.SetLegendStyle(leg1);
	leg1 -> AddEntry(hMissingPtTPC, "2 TPC-TOF tracks", "l");
	leg1 -> AddEntry(hMissingPtQ0, "Total charge 0", "fl");
	leg1 -> AddEntry(hMissingPtExc, "Exclusive", "fl");
	leg1->Draw("same");

	TPaveText *textPub1 = new TPaveText(0.25,0.75,0.35,0.75,"brNDC");
	textPub1 -> SetTextSize(0.04);
	textPub1 -> SetFillColor(0);
	textPub1 -> SetTextFont(42);
	if(TEXT)
		textPub1 -> AddText("Exclusive peak");
	textPub1 -> Draw("same");

	TPaveText *textPub2 = new TPaveText(0.75,0.85,0.85,0.96,"brNDC");
	tool.SetTextStyle(textPub2);
	textPub2 -> AddText("p + p #rightarrow p + #pi^{+} + #pi^{-} + p");
	textPub2 -> AddText("#sqrt{s} = 510 GeV");
	textPub2 -> Draw("same");

	//textSTAR = new TPaveText(0.15,0.91,0.2,0.97,"brNDC"); //for text "star"
	textSTAR = new TPaveText(0.21,0.91,0.28,0.97,"brNDC");
	textSTAR -> SetTextSize(0.045);
	textSTAR -> SetFillColor(0);
	textSTAR -> SetTextFont(62);
	//textSTAR -> AddText("STAR");
	textSTAR -> AddText("THIS THESIS");
	textSTAR -> Draw("same");

	TLine *left = new TLine(0.1,0,0.1,600000);
	tool.SetLineStyle(left,10,1,4);
   left->Draw("same");

	cMissingPt->Update();
	cMissingPt-> Write("hMissingPt");
	cMissingPt->SaveAs(output + "BasicPlots/cMissingPt.png");
	cMissingPt->Close();



   cCanvas->Close();
   fout->Close();

}//BasicPlots::PlotHistogram

