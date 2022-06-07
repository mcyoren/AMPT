//UU
#include <iostream>
#include <TGraph.h>
#include <TGraphErrors.h>
#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TF1.h"
#include "TProfile.h"
#include "TLatex.h"
#include "TMarker.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TNtuple.h"
#include "TMultiGraph.h"
#include "TImage.h"
#include "TFile.h"
#include "TLine.h"
using namespace std;

void Format_Graph(TGraphErrors *g, int Mstyle, double Msize, int Mcolor, int Lwidth, int Lcolor, double_t Malpha, double_t Lalpha)
{
    g->SetMarkerStyle(Mstyle);
    g->SetMarkerSize(Msize);
    g->SetMarkerColorAlpha(Mcolor,Malpha);
    if(Mcolor==3) g->SetMarkerColorAlpha(kGreen+2,Malpha);
    if(Mcolor==5) g->SetMarkerColorAlpha(kOrange+7,Malpha);
    g->SetLineWidth(Lwidth);
    g->SetLineColorAlpha(Lcolor,Lalpha);
    if(Lcolor==3) g->SetLineColorAlpha(kGreen+2,Lalpha);
    if(Lcolor==5) g->SetLineColorAlpha(kOrange+7,Lalpha);

}
void Format_Canvas(TCanvas *c2, int divide_x, int divide_y, bool if_logy)
{
    c2->Divide(divide_x,divide_y);
    TPad *c2_1; char name1[200];
    for (int centr = 0; centr < divide_x*divide_y; centr++) {
        c2->cd(centr + 1);
        c2->SetFillColor(0);
        c2->SetBorderMode(0);
        c2->SetBorderSize(2);
        if (if_logy) c2->SetLogy();
        c2->SetFrameBorderMode(0);
        sprintf(name1,"%s_%d",c2->GetName(), centr+1);

        c2_1= (TPad*) c2->GetListOfPrimitives()->FindObject(name1);
        c2_1->SetTickx();
        c2_1->SetTicky();
    }

}

void Format_Pad(double_t left, double_t right, double_t min, double_t max, char title_x[100], char title_y[100], double_t offset_x, double_t offset_y, double_t Tsize, double_t Lsize, char title[100])
{


    TH1F *second = new TH1F("", "", 100, left, right);

    second->SetMinimum(min);
    second->SetMaximum(max);
    second->SetStats(0);
    second->GetXaxis()->SetTitle(title_x);
    second->GetXaxis()->SetLabelFont(42);
    second->GetXaxis()->SetTitleFont(42);
    second->GetXaxis()->SetLabelSize(Lsize);
    second->GetXaxis()->SetTickSize(0.025);
    second->GetYaxis()->SetTickSize(0.015);
    second->GetYaxis()->SetLabelSize(Lsize);
    second->GetXaxis()->SetTitleSize(Tsize);
    second->GetXaxis()->SetTitleOffset(offset_x);
    second->GetYaxis()->SetTitle(title_y);
    second->GetYaxis()->SetLabelFont(42);
    second->GetYaxis()->SetTitleOffset(1.5);
    second->GetYaxis()->SetTitleFont(42);
    second->GetYaxis()->SetTitleOffset(offset_y);
    second->GetYaxis()->SetTitleSize(Tsize*1.);
    second->SetTitle(title);
    second->SetTitleSize(Tsize);
    second->SetTitleFont(42);

    second->Draw();
}

void Format_Function(TF1 *f, int Mstyle, int Msize, int Mcolor, int Lwidth, int Lcolor, double_t Malpha, double_t Lalpha)
{
    f->SetMarkerStyle(Mstyle);
    f->SetMarkerSize(Msize);
    f->SetMarkerColorAlpha(Mcolor,Malpha);
    f->SetLineWidth(Lwidth);
    f->SetLineColorAlpha(Lcolor,Lalpha);

}

void Format_Hist1(TH1 *h1, int Mstyle, int Msize, int Mcolor, int Lwidth, int Lcolor, double_t Malpha, double_t Lalpha)
{
    h1->SetMarkerStyle(Mstyle);
    h1->SetMarkerSize(Msize);
    h1->SetMarkerColorAlpha(Mcolor,Malpha);
    h1->SetLineWidth(Lwidth);
    h1->SetLineColorAlpha(Lcolor,Lalpha);

}

void Format_Hist2(TH2 *h2, int Mstyle, int Msize, int Mcolor, int Lwidth, int Lcolor, double_t Malpha, double_t Lalpha)
{
    h2->SetMarkerStyle(Mstyle);
    h2->SetMarkerSize(Msize);
    h2->SetMarkerColorAlpha(Mcolor,Malpha);
    h2->SetLineWidth(Lwidth);
    h2->SetLineColorAlpha(Lcolor,Lalpha);

}

void Format_Hist3(TH3 *h3, int Mstyle, int Msize, int Mcolor, int Lwidth, int Lcolor, double_t Malpha, double_t Lalpha)
{
    h3->SetMarkerStyle(Mstyle);
    h3->SetMarkerSize(Msize);
    h3->SetMarkerColorAlpha(Mcolor,Malpha);
    h3->SetLineWidth(Lwidth);
    h3->SetLineColorAlpha(Lcolor,Lalpha);

}

void Format_Profile(TProfile *p, int Mstyle, int Msize, int Mcolor, int Lwidth, int Lcolor, double_t Malpha, double_t Lalpha)
{
    p->SetMarkerStyle(Mstyle);
    p->SetMarkerSize(Msize);
    p->SetMarkerColorAlpha(Mcolor,Malpha);
    p->SetLineWidth(Lwidth);
    p->SetLineColorAlpha(Lcolor,Lalpha);

}

void Format_Text(TLatex *tex00, double size, double font, double color)
{
    tex00->SetTextSize(size);
    tex00->SetTextColor(color);
    tex00->SetTextFont(font);
}

void DrawLogo(float x0 = 0., float y0 = 0., float x1 = 0., float y1 = 0.) {

    TImage* img = TImage::Open("input/logo_prelim.png");

    if (!img) {
        printf("Could not create an image... exit\n");
        return;
    }

    TPad* l = new TPad("l", "l", x0, y0, x1, y1);
    //gPad->cd(0);
    l->SetFillColorAlpha(0,0);
    l->Draw("same");
    l->cd();
    img->Draw();
}

void DrawTypeC(double syst, double abc, double width, double Lwidth, int Lcolor, double Lalpha)
{
    double xxx[2]={abc,1e30}, yyy[2] = {1,1e30}, del_yyy[2]={syst,0}, del_xxx[2]={width,0};
    TGraphErrors *gr_typeC = new TGraphErrors(1,xxx,yyy,del_xxx,del_yyy);
    gr_typeC->SetMarkerSize(0);
    if(Lalpha==0)
    {
        gr_typeC->SetLineWidth(Lwidth);
        gr_typeC->SetLineColor(Lcolor);
        gr_typeC->SetFillStyle(0);
    }else
    {
        gr_typeC->SetLineWidth(0);
        gr_typeC->SetFillColorAlpha(Lcolor,Lalpha);
    }
    gr_typeC->Draw("same 2P");
}
