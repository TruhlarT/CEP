#include "trackQuality.h"
#include "Plot.h"

//_____________________________________________________________________________
trackQuality::trackQuality(string innam, string outnam, bool text)
{
	//constructor
	input = innam;
	output = outnam;

	TEXT = text;

	cout << "trackQuality::trackQuality() called" << endl;

}//trackQuality

//_____________________________________________________________________________
trackQuality::~trackQuality()
{
  //destructor

  cout << "trackQuality::~trackQuality() destructor called" << endl;

}//~trackQuality


//_____________________________________________________________________________
void trackQuality::PlotHistogram(){

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
	gPad->SetMargin(0.9,0.02,0.1,0.02); // (Float_t left, Float_t right, Float_t bottom, Float_t top)
	gPad->SetTickx();
	gPad->SetTicky();  
	gStyle->SetOptStat("");
	gStyle->SetPalette(1);
	gStyle->SetLineWidth(2);      //axis line
	gStyle->SetFrameLineWidth(2); //frame line

	tree->UseCurrentStyle();
	treeBack->UseCurrentStyle();

	Plot tool;


// Plot z vertex
	treeBack->Draw("vertexZ>>vertexZBcg","nSigPPion<3");
	tmpHist2 = (TH1F*)gPad->GetPrimitive("vertexZBcg");
	tool.SetMarkerStyle(tmpHist2,2,20,1,2,1,1);

	tree->Draw("vertexZ>>vertexZSig","nSigPPion<3");
	tmpHist = (TH1F*)gPad->GetPrimitive("vertexZSig");
	tmpHist->SetTitle(" ; Z_{vrtx} [cm]; Number of events");
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
		TLine *left4 = new TLine(-80,0,-80,45);
		tool.SetLineStyle(left4,10,1,4);
	   left4->Draw("same");

		TLine *left5 = new TLine(80,0,80,45);
		tool.SetLineStyle(left5,10,1,4);
	   left5->Draw("same");		
	}
	cCanvas->Update();
	cCanvas->SaveAs( output + "trackQuality/zVertex.png");
	cCanvas->Write("zVertex");
//////////////////////////////////////////
// Plot DCAXY vertex
	treeBack->Draw("DcaXY1>>DcaXY1Bcg(100,0,3.5)","nSigPPion<3");
	tmpHist2 = (TH1F*)gPad->GetPrimitive("DcaXY1Bcg");
	treeBack->Draw("DcaXY2>>DcaXY2Bcg(100,0,3.5)","nSigPPion<3");
	tmpHist = (TH1F*)gPad->GetPrimitive("DcaXY2Bcg");
	tmpHist2->Add(tmpHist);
	tool.SetMarkerStyle(tmpHist2,2,20,1,2,1,1);

	tree->Draw("DcaXY1>>DcaXY1Sig(100,0,3.5)","nSigPPion<3");
	tmpHist = (TH1F*)gPad->GetPrimitive("DcaXY1Sig");
	tree->Draw("DcaXY2>>DcaXY2Sig(100,0,3.5)","nSigPPion<3");
	tmpHist3 = (TH1F*)gPad->GetPrimitive("DcaXY2Sig");
	tmpHist->Add(tmpHist3);
	tmpHist->SetTitle(" ; DCA_{xy} [cm]; Number of events");
	tool.SetGraphStyle(tmpHist);
	tool.SetMarkerStyle(tmpHist);
	tmpHist->Draw("E");
	tool.DrawText(tmpHist);
	tmpHist2->Draw("ESAME");

	leg1 = new TLegend(0.62, 0.7, 0.82, 0.82);
	tool.SetLegendStyle(leg1);
	leg1->AddEntry(tmpHist,"In+El (unlike-sign pairs)","p");
	leg1->AddEntry(tmpHist2,"In+El (like-sign pairs)","p");
	leg1->Draw("same");
	if(TEXT){
		TLine *left6 = new TLine(1.5,0,1.5,70);
		tool.SetLineStyle(left6,10,1,4);
	   left6->Draw("same");
	}
	cCanvas->Update();
	cCanvas->SaveAs( output + "trackQuality/DcaXY.png");
	cCanvas->Write("DcaXY");
//////////////////////////////////////////
// Plot DCAZ vertex
	treeBack->Draw("DcaZ1>>DcaZ1Bcg(100,-1.5,1.5)","nSigPPion<3");
	tmpHist2 = (TH1F*)gPad->GetPrimitive("DcaZ1Bcg");
	treeBack->Draw("DcaZ2>>DcaZ2Bcg(100,-1.5,1.5)","nSigPPion<3");
	tmpHist = (TH1F*)gPad->GetPrimitive("DcaZ2Bcg");
	tmpHist2->Add(tmpHist);
	tool.SetMarkerStyle(tmpHist2,2,20,1,2,1,1);

	tree->Draw("DcaZ1>>DcaZ1Sig(100,-1.5,1.5)","nSigPPion<3");
	tmpHist = (TH1F*)gPad->GetPrimitive("DcaZ1Sig");
	tree->Draw("DcaZ2>>DcaZ2Sig(100,-1.5,1.5)","nSigPPion<3");
	tmpHist3 = (TH1F*)gPad->GetPrimitive("DcaZ2Sig");
	tmpHist->Add(tmpHist3);
	tmpHist->SetTitle(" ; DCA_{z} [cm]; Number of events");
	tool.SetGraphStyle(tmpHist);
	tool.SetMarkerStyle(tmpHist);
	tmpHist->Draw("E");
	tool.DrawText(tmpHist);
	tmpHist2->Draw("ESAME");

	leg1 = new TLegend(0.62, 0.7, 0.82, 0.82);
	tool.SetLegendStyle(leg1);
	leg1->AddEntry(tmpHist,"In+El (unlike-sign pairs)","p");
	leg1->AddEntry(tmpHist2,"In+El (like-sign pairs)","p");
	leg1->Draw("same");
	if(TEXT){
		TLine *left7 = new TLine(1,0,1,500);
		tool.SetLineStyle(left7,10,1,4);
	   left7->Draw("same");

		TLine *left8 = new TLine(-1,0,-1,500);
		tool.SetLineStyle(left8,10,1,4);
	   left8->Draw("same");
	}
	cCanvas->Update();
	cCanvas->SaveAs( output + "trackQuality/DcaZ.png");
	cCanvas->Write("DcaZ");
//////////////////////////////////////////
// Plot NhitsDEdx vertex
	treeBack->Draw("NhitsDEdx1>>NhitsDEdx1Bcg(51,0,50)","nSigPPion<3");
	tmpHist2 = (TH1F*)gPad->GetPrimitive("NhitsDEdx1Bcg");
	treeBack->Draw("NhitsDEdx2>>NhitsDEdx2Bcg(51,0,50)","nSigPPion<3");
	tmpHist = (TH1F*)gPad->GetPrimitive("NhitsDEdx2Bcg");
	tmpHist2->Add(tmpHist);
	tool.SetMarkerStyle(tmpHist2,2,20,1,2,1,1);

	tree->Draw("NhitsDEdx1>>NhitsDEdx1Sig(61,0,60)","nSigPPion<3");
	tmpHist = (TH1F*)gPad->GetPrimitive("NhitsDEdx1Sig");
	tree->Draw("NhitsDEdx2>>NhitsDEdx2Sig(61,0,60)","nSigPPion<3");
	tmpHist3 = (TH1F*)gPad->GetPrimitive("NhitsDEdx2Sig");
	tmpHist->Add(tmpHist3);
	tmpHist->SetTitle(" ; N^{dE/dx}_{hits} ; Number of events");
	tool.SetGraphStyle(tmpHist);
	tool.SetMarkerStyle(tmpHist);
	tmpHist->Draw("E");
	tool.DrawText(tmpHist,.25,0.77,.33,.91);
	tmpHist2->Draw("ESAME");

	leg1 = new TLegend(0.13,0.65,0.2,0.77);
	tool.SetLegendStyle(leg1);
	leg1->AddEntry(tmpHist,"In+El (unlike-sign pairs)","p");
	leg1->AddEntry(tmpHist2,"In+El (like-sign pairs)","p");
	leg1->Draw("same");

	if(TEXT){
		TLine *left9 = new TLine(15,0,15,120);
		tool.SetLineStyle(left9,10,1,4);
	   left9->Draw("same");
	}

	cCanvas->Update();
	cCanvas->SaveAs( output + "trackQuality/NhitsDEdx.png");
	cCanvas->Write("NhitsDEdx");
//////////////////////////////////////////
// Plot Eta vertex
	treeBack->Draw("Eta1>>Eta1Bcg(100,-2,3.5)","nSigPPion<3");
	tmpHist2 = (TH1F*)gPad->GetPrimitive("Eta1Bcg");
	treeBack->Draw("Eta2>>Eta2Bcg(100,-2,3.5)","nSigPPion<3");
	tmpHist = (TH1F*)gPad->GetPrimitive("Eta2Bcg");
	tmpHist2->Add(tmpHist);
	tool.SetMarkerStyle(tmpHist2,2,20,1,2,1,1);

	tree->Draw("Eta1>>Eta1Sig(100,-2,3.5)","nSigPPion<3");
	tmpHist = (TH1F*)gPad->GetPrimitive("Eta1Sig");
	tree->Draw("Eta2>>Eta2Sig(100,-2,3.5)","nSigPPion<3");
	tmpHist3 = (TH1F*)gPad->GetPrimitive("Eta2Sig");
	tmpHist->Add(tmpHist3);
	tmpHist->SetTitle(" ; #eta ; Number of tracks");
	tool.SetGraphStyle(tmpHist);
	tool.SetMarkerStyle(tmpHist);
	tmpHist->Draw("E");
	tool.DrawText(tmpHist);
	tmpHist2->Draw("ESAME");

	leg1 = new TLegend(0.58, 0.7, 0.78, 0.82);
	tool.SetLegendStyle(leg1);
	leg1->AddEntry(tmpHist,"In+El (unlike-sign pairs)","p");
	leg1->AddEntry(tmpHist2,"In+El (like-sign pairs)","p");
	leg1->Draw("same");
	if(TEXT){
		TLine *left2 = new TLine(-0.7,0,-0.7,100);
		tool.SetLineStyle(left2,10,1,4);
	   left2->Draw("same");

		TLine *left3 = new TLine(0.7,0,0.7,80);
		tool.SetLineStyle(left3,10,1,4);
	   left3->Draw("same");
	}
	cCanvas->Update();
	cCanvas->SaveAs( output + "trackQuality/Eta.png");
	cCanvas->Write("Eta");
//////////////////////////////////////////

	TCanvas *cCanvas2D = new TCanvas("cCanvas2D","cCanvas2D",800,700);
	gPad->SetMargin(0.1,0.09,0.1,0.02); // (Float_t left, Float_t right, Float_t bottom, Float_t top)
	gStyle->SetPalette(1);


// Plot XYEastCor vertex
	tree->Draw("yCorrelationsEast:xCorrelationsEast>>XYEastCorSig(100,-1,1.2,100,-1,1.2)","nSigPPion<3","colz");
	tmp2DHist = (TH2F*)gPad->GetPrimitive("XYEastCorSig");
	tmp2DHist->SetTitle(" ; p_{x} [GeV/c]; p_{y} [GeV/c]");
	tool.SetGraphStyle(tmp2DHist);
	tmp2DHist->Draw("colz");
	tool.DrawText(tmp2DHist,true,0.7,0.82,0.77,0.96);
	if(TEXT){
		const Int_t n = 200;
		Double_t x[n], y[n];
		Double_t tmp;
		for(int i = 0; i < n; ++i){
			x[i] = -0.036 + (0.152*i)/n;
			tmp = (x[i] + 0.3)*(x[i]+0.3) - 0.46*0.46;
			y[i] = sqrt(abs(tmp));
			if(tmp<0)
				y[i] = - y[i];
		}
		TGraph* gr = new TGraph(n,x,y);
		gr->SetLineWidth(4);
		gr->Draw("same");

		for(int i = 0; i < n; ++i){
			x[i] = -0.22 + (0.185*i)/n;
			tmp = (x[i] + -0.15)*(x[i]-0.15) - 0.42*0.42;
			y[i] = sqrt(abs(tmp));
			if(tmp<0)
				y[i] = - y[i];
		}
		TGraph* gr1 = new TGraph(n,x,y);
		gr1->SetLineWidth(4);
		gr1->Draw("same");

		TLine *left02 = new TLine(-0.219,-0.2,0.116,-0.2);
		tool.SetLineStyle(left02);
	   left02->Draw("same");


	   TLine *left01 = new TLine(-0.219,0.2,0.116,0.2);
		tool.SetLineStyle(left01);
	   left01->Draw("same");

	   for(int i = 0; i < n; ++i){
			x[i] = -0.036 + (0.152*i)/n;
			tmp = (x[i] + 0.3)*(x[i]+0.3) - 0.46*0.46;
			y[i] = sqrt(abs(tmp));
		}
		TGraph* gr2 = new TGraph(n,x,y);
		gr2->SetLineWidth(4);
		gr2->Draw("same");

		for(int i = 0; i < n; ++i){
			x[i] = -0.22 + (0.185*i)/n;
			tmp = (x[i] + -0.15)*(x[i]-0.15) - 0.42*0.42;
			y[i] = sqrt(abs(tmp));
		}
		TGraph* gr3 = new TGraph(n,x,y);
		gr3->SetLineWidth(4);
		gr3->Draw("same");
	}
	cCanvas2D->Update();
	cCanvas2D->SaveAs( output + "trackQuality/XYEastCor.png");
	cCanvas2D->Write("XYEastCor");
//////////////////////////////////////////

// Plot XYWestCor vertex
	tree->Draw("yCorrelationsWest:xCorrelationsWest>>XYWestCorSig(100,-1,1.2,100,-1,1.2)","nSigPPion<3","colz");
	tmp2DHist = (TH2F*)gPad->GetPrimitive("XYWestCorSig");
	tmp2DHist->SetTitle(" ; p_{x} [GeV/c]; p_{y} [GeV/c]");
	tool.SetGraphStyle(tmp2DHist);
	tmp2DHist->Draw("colz");
	tool.DrawText(tmp2DHist,true,0.7,0.82,0.77,0.96);

	if(TEXT){
		const Int_t n = 200;
		Double_t x[n], y[n];
		Double_t tmp;
		for(int i = 0; i < n; ++i){
			x[i] = -0.036 + (0.152*i)/n;
			tmp = (x[i] + 0.3)*(x[i]+0.3) - 0.46*0.46;
			y[i] = sqrt(abs(tmp));
			if(tmp<0)
				y[i] = - y[i];
		}
		TGraph* gr = new TGraph(n,x,y);
		gr->SetLineWidth(4);
		gr->Draw("same");

		for(int i = 0; i < n; ++i){
			x[i] = -0.22 + (0.185*i)/n;
			tmp = (x[i] + -0.15)*(x[i]-0.15) - 0.42*0.42;
			y[i] = sqrt(abs(tmp));
			if(tmp<0)
				y[i] = - y[i];
		}
		TGraph* gr1 = new TGraph(n,x,y);
		gr1->SetLineWidth(4);
		gr1->Draw("same");

		TLine *left02 = new TLine(-0.219,-0.2,0.116,-0.2);
		tool.SetLineStyle(left02);
	   left02->Draw("same");


	   TLine *left01 = new TLine(-0.219,0.2,0.116,0.2);
		tool.SetLineStyle(left01);
	   left01->Draw("same");

	   for(int i = 0; i < n; ++i){
			x[i] = -0.036 + (0.152*i)/n;
			tmp = (x[i] + 0.3)*(x[i]+0.3) - 0.46*0.46;
			y[i] = sqrt(abs(tmp));
		}
		TGraph* gr2 = new TGraph(n,x,y);
		gr2->SetLineWidth(4);
		gr2->Draw("same");

		for(int i = 0; i < n; ++i){
			x[i] = -0.22 + (0.185*i)/n;
			tmp = (x[i] + -0.15)*(x[i]-0.15) - 0.42*0.42;
			y[i] = sqrt(abs(tmp));
		}
		TGraph* gr3 = new TGraph(n,x,y);
		gr3->SetLineWidth(4);
		gr3->Draw("same");
	}
	cCanvas2D->Update();
	cCanvas2D->SaveAs( output + "trackQuality/hXYWestCor.png");
	cCanvas2D->Write("hXYWestCor");
//////////////////////////////////////////
// Plot ZvrtxVsEta vertex
	tree->Draw("Eta2:vertexZ>>ZvrtxVsEta2Sig(100,-200,200,100,-1,1.1)","nSigPPion<3","colz");
	tmp2DHist2 = (TH2F*)gPad->GetPrimitive("ZvrtxVsEta2Sig");
	tree->Draw("Eta1:vertexZ>>ZvrtxVsEtaSig(100,-200,200,100,-1,1.1)","nSigPPion<3","colz");
	tmp2DHist = (TH2F*)gPad->GetPrimitive("ZvrtxVsEtaSig");
	tmp2DHist->Add(tmp2DHist2);
	tmp2DHist->SetTitle(" ; Z_{vrtx} [cm]; #eta");
	tool.SetGraphStyle(tmp2DHist);
	tmp2DHist->Draw("colz");
	tool.DrawText(tmp2DHist,true,0.7,0.82,0.77,0.96);
	if(TEXT){
		TLine *left = new TLine(-80,-0.7,-80,0.7);
		tool.SetLineStyle(left);
	   left->Draw("same");

		TLine *left1 = new TLine(80,-0.7,80,0.7);
		tool.SetLineStyle(left1);
	   left1->Draw("same");

		TLine *line = new TLine(-80,0.7,80,0.7);
		tool.SetLineStyle(line);
	   line->Draw("same");

		TLine *line2 = new TLine(-80,-0.7,80,-0.7);
		tool.SetLineStyle(line2);
	   line2->Draw("same");	
	}
	cCanvas2D->Update();
	cCanvas2D->SaveAs( output + "trackQuality/ZvrtxVsEta.png");
	cCanvas2D->Write("ZvrtxVsEta");
//////////////////////////////////////////

	cCanvas->Close();
	cCanvas2D->Close();
   fout->Close();
}

