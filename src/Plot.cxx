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
						Int_t lineColor, Int_t lineStyle, Float_t lineWidth, Float_t xOffSet, Float_t yOffSet){
	hist->GetXaxis()->SetTitleFont(font);
	hist->GetYaxis()->SetTitleFont(font);
	hist->GetXaxis()->SetLabelFont(font);
	hist->GetYaxis()->SetLabelFont(font);
	hist->GetXaxis()->SetTitleSize(siz);
	hist->GetYaxis()->SetTitleSize(siz);
	hist->GetXaxis()->SetTitleOffset(xOffSet);
	hist->GetYaxis()->SetTitleOffset(yOffSet);

}//SetGraphStyle


void Plot::SetMarkerStyle(TH1* hist, Int_t markColor, Int_t markStyle, Int_t markSize, 
						Int_t lineColor, Int_t lineStyle, Float_t lineWidth){
	hist->SetMarkerColor(markColor);
	hist->SetMarkerSize(markSize);
	hist->SetMarkerStyle(markStyle);
	hist->SetLineColor(lineColor);
	hist->SetLineStyle(lineStyle);
	hist->SetLineWidth(lineWidth);

}//SetGraphStyle

void Plot::DrawText(TH1* hist, Bool_t pair, Float_t xMin, Float_t yMin, Float_t xMax, Float_t yMax){
	TPaveText *textPub = new TPaveText(xMin,yMin,xMax,yMax,"brNDC");
	textPub -> SetTextSize(siz-0.005);
	textPub -> SetFillColor(0);
	textPub -> SetTextFont(font);
	textPub -> AddText("p + p #rightarrow p' + #pi^{+} + #pi^{-} + p'");
	textPub -> AddText("#sqrt{s} = 510 GeV");
	int NentriesEl = hist->GetEntries();
	TString tileIdStrEl; tileIdStrEl.Form("%i Exclusive #pi candidates",NentriesEl);
	if(pair)
		tileIdStrEl.Form("%i Exclusive events",NentriesEl);
	textPub -> AddText(tileIdStrEl);
	textPub -> Draw("same");


	TPaveText *textSTAR = new TPaveText(0.15,0.91,0.2,0.97,"brNDC"); //for text "star"
	//TPaveText *textSTAR = new TPaveText(0.2,0.91,0.27,0.97,"brNDC");
	textSTAR -> SetTextSize(siz);
	textSTAR -> SetFillColor(0);
	textSTAR -> SetTextFont(62);
	//textSTAR -> AddText("THIS THESIS");
	textSTAR->AddText("STAR");
	textSTAR -> Draw("same");
}

void Plot::SetLegendStyle(TLegend* leg1){
	leg1->SetFillStyle(0);
	leg1->SetBorderSize(0);
	leg1->SetTextSize(siz-0.005);
	leg1->SetTextFont(font);
}

void Plot::SetTextStyle(TPaveText* text){
	text->SetFillStyle(0);
	text->SetFillColor(0);
	text->SetBorderSize(0);
	text->SetTextSize(siz-0.005);
	text->SetTextFont(font);
}


void Plot::SetLineStyle(TLine* line, Int_t style, Int_t color, Int_t width){
	line->SetLineStyle(style);
	line->SetLineColor(color);
	line->SetLineWidth(width);
}
