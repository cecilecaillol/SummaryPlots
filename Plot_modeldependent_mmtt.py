from RecoLuminosity.LumiDB import argparse
import math
import os
from HttStyles import GetStyleHtt
from HttStyles import MakeCanvas
import ROOT
import numpy as np
from array import array
from BR import get_total_width
from BR import gamma_quarks
from BR import gamma_mu
from BR import gamma_tau
from BR import gamma_photon
from BR import gamma_gg

def add_lumi():
    lowX=0.56
    lowY=0.835
    lumi  = ROOT.TPaveText(lowX, lowY+0.06, lowX+0.30, lowY+0.16, "NDC")
    lumi.SetBorderSize(   0 )
    lumi.SetFillStyle(    0 )
    lumi.SetTextAlign(   12 )
    lumi.SetTextColor(    1 )
    lumi.SetTextSize(0.045)
    lumi.AddText("19.7 fb^{-1} (8 TeV)")
    return lumi

def add_CMS():
    lowX=0.18
    lowY=0.835
    lumi  = ROOT.TPaveText(lowX, lowY+0.06, lowX+0.15, lowY+0.16, "NDC")
    lumi.SetTextFont(61)
    lumi.SetTextSize(0.045)
    lumi.SetBorderSize(   0 )
    lumi.SetFillStyle(    0 )
    lumi.SetTextAlign(   12 )
    lumi.SetTextColor(    1 )
    lumi.AddText("CMS Preliminary")
    return lumi 

def add_Preliminary():
    lowX=0.21
    lowY=0.70
    lumi  = ROOT.TPaveText(lowX, lowY+0.06, lowX+0.15, lowY+0.16, "NDC")
    lumi.SetTextSize(0.03)
    lumi.SetBorderSize(   0 )
    lumi.SetFillStyle(    0 )
    lumi.SetTextAlign(   12 )
    lumi.SetTextColor(    1 )
    lumi.SetTextFont(52)
    lumi.AddText("Preliminary")
    return lumi 

def make_legend():
    output = ROOT.TLegend(0.42, 0.68, 0.89, 0.90, "", "brNDC")
    output.SetLineWidth(0)
    output.SetLineStyle(0)
    output.SetFillStyle(0)
    output.SetBorderSize(0)
    output.SetTextFont(62)
    return output

def add_model(model):
    lowX=0.21
    lowY=0.20
    if (args.model==4):
       lowY=0.50
    lumi  = ROOT.TPaveText(lowX, lowY+0.06, lowX+0.15, lowY+0.16, "NDC")
    lumi.SetTextSize(0.05)
    lumi.SetBorderSize(   0 )
    lumi.SetFillStyle(    0 )
    lumi.SetTextAlign(   12 )
    lumi.SetTextColor(    1 )
    lumi.SetTextFont(62)
    lumi.AddText("2HDM+S type-"+str(model))
    return lumi

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('--model', type=int, default='1', help="Which type of 2HDM?")

    args = parser.parse_args()

    style1=GetStyleHtt()
    #style1.SetPalette(1)
    style1.cd()

    binbeta=1000
    minbeta=1.45
    maxbeta=5
    if (args.model==4):
       minbeta=0.28
       maxbeta=0.9

    x_mmtt1, y_mmtt1 = np.loadtxt('mmtt.txt', unpack=True)
    x_mmtt=array("d",x_mmtt1)
    y_mmtt=array("d",y_mmtt1)
    hist=ROOT.TH2F("hist","hist",len(x_mmtt)-1,x_mmtt[0],x_mmtt[len(x_mmtt)-1],binbeta,minbeta,maxbeta)
    for b in range(0,binbeta+1):
	tanbeta=0.001+minbeta+1.0*b*(maxbeta-minbeta)/(binbeta)
        for i in range(0,len(x_mmtt)):
	   #for b in range(0,binbeta):
	   #tanbeta=minbeta+1.0*b*(maxbeta-minbeta)/binbeta
           width=get_total_width(args.model,float(x_mmtt[i]),tanbeta)
           BRtt=gamma_tau(tanbeta,float(x_mmtt[i]),args.model)/width
	   y=y_mmtt[i]/(BRtt*BRtt)
	   print x_mmtt[i],tanbeta,y
           hist.Fill(x_mmtt[i],1.0*tanbeta,y)
    #hist.SetContour(500)

    canvas = MakeCanvas("asdf","asdf",800,800)
    canvas.SetRightMargin(0.2)
    canvas.SetLeftMargin(0.15)
    canvas.cd()
    canvas.SetLogz()
    hist.GetXaxis().SetTitle("m_{a} (GeV)")
    hist.GetYaxis().SetTitle("tan#beta")
    hist.GetZaxis().SetTitle("#frac{#sigma(h)}{#sigma_{SM}} #times B(h#rightarrowaa)")
    if (args.model==3):
       hist.GetZaxis().SetRangeUser(0.04,1.6)
    if (args.model==2):
       hist.GetZaxis().SetRangeUser(5,100)
    if (args.model==4):
       hist.GetZaxis().SetRangeUser(0.30,20)
    hist.Draw("colz")
    lumiBlurb1=add_CMS()
    lumiBlurb1.Draw("same")
    #lumiBlurb2=add_Preliminary()
    #lumiBlurb2.Draw("same")
    lumiBlurb=add_lumi()
    lumiBlurb.Draw("same")
    lumi=add_model(args.model)
    lumi.Draw("same")
    canvas.SaveAs('plots/plot_BRaa_mmtt_Type'+str(args.model)+'.png')
    canvas.SaveAs('plots/plot_BRaa_mmtt_Type'+str(args.model)+'.pdf')


