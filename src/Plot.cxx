#include "Plot.h"


//_____________________________________________________________________________
Plot::Plot()
{

	siz = 0.045;
	font = 42;

	cout << "Plot::Plot() called" << endl;

}//Plot

//_____________________________________________________________________________
Plot::~Plot()
{
  //destructor

  cout << "Plot::~Plot() destructor called" << endl;

}//~Plot


void Plot::SetGraphStyle(TH1* hist, Int_t markColor, Int_t markStyle, Int_t markSize, 
						Int_t lineColor, Int_t lineStyle, Float_t lineWidth, Float_t xOffSet, Float_t yOffSet)
{
	hist->GetXaxis()->SetTitleFont(font);
	hist->GetYaxis()->SetTitleFont(font);
	hist->GetXaxis()->SetLabelFont(font);
	hist->GetYaxis()->SetLabelFont(font);
	hist->GetXaxis()->SetLabelSize(siz+0.005);
	hist->GetYaxis()->SetLabelSize(siz+0.005);
	hist->GetXaxis()->SetTitleSize(siz+0.005);
	hist->GetYaxis()->SetTitleSize(siz+0.005);
	hist->GetXaxis()->SetTitleOffset(xOffSet);
	hist->GetYaxis()->SetTitleOffset(yOffSet);

}//SetGraphStyle


void Plot::SetMarkerStyle(TH1* hist, Int_t markColor, Int_t markStyle, Int_t markSize, 
						Int_t lineColor, Int_t lineStyle, Float_t lineWidth)
{
	hist->SetMarkerColor(markColor);
	hist->SetMarkerSize(markSize);
	hist->SetMarkerStyle(markStyle);
	hist->SetLineColor(lineColor);
	hist->SetLineStyle(lineStyle);
	hist->SetLineWidth(lineWidth);

}//SetMarkerStyle


void Plot::DrawText(TH1* hist, Int_t state, Int_t pair, Float_t xMin, Float_t yMin, Float_t xMax, Float_t yMax, Int_t align, Bool_t data)
{
	TString centralState;
	switch(state){
		case 0:
			centralState = "X";
			xMax+=0.02;
			xMin+=0.02;
			break;
		case 1:
			centralState = "#pi^{+} #pi^{-}";
			break;
		case 2:
			centralState = "K^{+} K^{-}";
			break;
		case 3:
			centralState = "p #bar{p}";
			break;
		case 4:
			centralState = "#pi^{+} #pi^{+} #pi^{-} #pi^{-}";
			break;
	}
	if(pair == 2)
		xMin-=0.06;

	TPaveText *textPub = new TPaveText(xMin,yMin,xMax,yMax,"brNDC");
	textPub -> SetTextSize(siz-0.005);
	textPub -> SetTextAlign(align);
	textPub -> SetFillColor(0);
	textPub -> SetTextFont(font);
	textPub -> AddText("p + p #rightarrow p + " + centralState +" + p");
	if(data) textPub -> AddText("#sqrt{s} = 510 GeV");
	int NentriesEl = hist->Integral();
	TString tileIdStrEl; tileIdStrEl.Form("%i h_{cand}^{Exc}",NentriesEl);
	if(pair == 1)
		tileIdStrEl.Form("%i Excl. events",NentriesEl);
	if(pair == 2)
		tileIdStrEl.Form("%i Corrected excl. events",NentriesEl);
	textPub -> AddText(tileIdStrEl);
	textPub -> Draw("same");

}

void Plot::DrawTextStar(TH1* hist, Int_t position, Bool_t star)
{
	TPaveText *textSTAR;

	switch(position){
		case 0:
			textSTAR = new TPaveText(0.18,0.91,0.28,0.97,"brNDC");
			break;
		case 1:
			textSTAR = new TPaveText(0.65,0.89,0.85,0.95,"brNDC");
			break;
		case 2:
			textSTAR = new TPaveText(0.75,0.89,0.9,0.95,"brNDC");
			break;
		case 3:
			textSTAR = new TPaveText(0.22,0.89,0.33,0.95,"brNDC");
			break;
	}
	
	textSTAR -> SetTextSize(siz-0.005);
	textSTAR -> SetFillColor(0);
	textSTAR -> SetTextFont(62);
	if(star)
		textSTAR->AddText("STAR Internal");
		//textSTAR->AddText("Request for Preliminary");
	else
		textSTAR -> AddText("THIS THESIS");
	textSTAR -> Draw("same");
}

void Plot::SetLegendStyle(TLegend* leg1)
{
	leg1->SetFillStyle(0);
	leg1->SetBorderSize(0);
	leg1->SetTextSize(siz-0.005);
	leg1->SetTextFont(font);
}

void Plot::SetTextStyle(TPaveText* text)
{
	text->SetFillStyle(0);
	text->SetFillColor(0);
	text->SetBorderSize(0);
	text->SetTextSize(siz-0.005);
	text->SetTextFont(font);
}


void Plot::SetLineStyle(TLine* line, Int_t style, Int_t color, Int_t width)
{
	line->SetLineStyle(style);
	line->SetLineColor(color);
	line->SetLineWidth(width);
}
