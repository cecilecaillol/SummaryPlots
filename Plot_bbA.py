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
    lowX=0.64
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
    lowX=0.21
    lowY=0.75
    lumi  = ROOT.TPaveText(lowX, lowY+0.06, lowX+0.15, lowY+0.16, "NDC")
    lumi.SetTextFont(61)
    lumi.SetTextSize(0.07)
    lumi.SetBorderSize(   0 )
    lumi.SetFillStyle(    0 )
    lumi.SetTextAlign(   12 )
    lumi.SetTextColor(    1 )
    lumi.AddText("CMS")
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
    output = ROOT.TLegend(0.42, 0.78, 0.93, 0.90, "", "brNDC")
    output.SetLineWidth(0)
    output.SetLineStyle(0)
    output.SetFillStyle(0)
    output.SetBorderSize(0)
    output.SetTextFont(62)
    return output

def add_model(model,tanbeta):
    lowX=0.21
    lowY=0.30
    lumi  = ROOT.TPaveText(lowX, lowY+0.06, lowX+0.15, lowY+0.16, "NDC")
    lumi.SetTextSize(0.04)
    lumi.SetBorderSize(   0 )
    lumi.SetFillStyle(    0 )
    lumi.SetTextAlign(   12 )
    lumi.SetTextColor(    1 )
    lumi.SetTextFont(62)
    if model<3:
       lumi.AddText("2HDM+S type-1/2")
    else:
        lumi.AddText("2HDM+S type-"+str(model))
        lumi.AddText("tan #beta = "+str(tanbeta))
    return lumi

if __name__ == "__main__":

    style1=GetStyleHtt()
    style1.cd()

    x_att1, y_att1 = np.loadtxt('att.txt', unpack=True)
    x_att=array("d",x_att1)
    y_att=array("d",y_att1)
    gatt = ROOT.TGraph(len(x_att), x_att,y_att)

    x_amm, y_amm = np.loadtxt('amm.txt', unpack=True)
    for i in range(0,len(x_amm)):
        BRmm=gamma_mu(1,float(x_amm[i]),1)
        BRtt=gamma_tau(1,float(x_amm[i]),1)
        #y_amm[i]=y_amm[i]*BRtt/BRmm
        y_amm[i]=100*y_amm[i]
    gamm = ROOT.TGraph(len(x_amm), x_amm.flatten('C'),y_amm.flatten('C'))

    canvas = MakeCanvas("asdf","asdf",800,800)
    canvas.cd()
    canvas.SetLogy()
    gamm.Draw("AC")
    gamm.SetLineColor(ROOT.EColor.kAzure-1)
    gamm.SetLineStyle(7)
    gamm.GetXaxis().SetRangeUser(25,80);
    gamm.GetXaxis().SetLimits(25,80);
    gamm.SetMinimum(5.0);
    gamm.SetMaximum(200);
    gamm.GetXaxis().SetTitle("m_{a} (GeV)");
    gamm.GetYaxis().SetTitle("95% CL on #sigma(bbA) #times B(a#rightarrow#tau#tau)")
    gamm.Draw("AC")
    canvas.Update();
    gatt.SetLineColor(ROOT.EColor.kGreen-1)
    gatt.Draw("same")
    canvas.Update();

    legend=make_legend()
    legend.AddEntry(gamm,"bbA#rightarrowbb#mu#mu, HIG-15-009","l")
    legend.AddEntry(gatt,"bbA#rightarrowbb#tau#tau, HIG-14-033","l")
    legend.Draw("same")

    lumiBlurb1=add_CMS()
    lumiBlurb1.Draw("same")
    lumiBlurb2=add_Preliminary()
    lumiBlurb2.Draw("same")
    lumiBlurb=add_lumi()
    lumiBlurb.Draw("same")
    #lumi=add_model(args.model,args.tanbeta)
    #lumi.Draw("same")
    canvas.SaveAs('plots/plot_bbA.png')
    canvas.SaveAs('plots/plot_bbA.pdf')


