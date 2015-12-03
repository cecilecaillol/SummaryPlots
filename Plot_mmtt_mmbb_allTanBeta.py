from RecoLuminosity.LumiDB import argparse
import math
import os
from HttStyles import GetStyleHtt
from HttStyles import MakeCanvas
from BR import get_total_width
from BR import gamma_quarks
from BR import gamma_mu
from BR import gamma_tau
from BR import gamma_photon
from BR import gamma_gg
import ROOT
import numpy as np
from array import array

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
    output = ROOT.TLegend(0.50, 0.67, 0.94, 0.80, "", "brNDC")
    output.SetLineWidth(0)
    output.SetLineStyle(0)
    output.SetFillStyle(0)
    output.SetBorderSize(0)
    output.SetTextFont(62)
    return output

def add_model(model,ma):
    lowX=0.64
    lowY=0.29
    lumi  = ROOT.TPaveText(lowX, lowY+0.06, lowX+0.15, lowY+0.16, "NDC")
    lumi.SetTextSize(0.04)
    lumi.SetBorderSize(   0 )
    lumi.SetFillStyle(    0 )
    lumi.SetTextAlign(   12 )
    lumi.SetTextColor(    1 )
    lumi.SetTextFont(62)
    #lumi.AddText("2HDM+S type "+str(model)+", tan#beta = "+str(tanbeta))
    lumi.AddText("2HDM+S type-"+str(model))
    lumi.AddText("m_{a} = "+str(int(ma))+" GeV")
    return lumi

if __name__ == "__main__":

    style1=GetStyleHtt()
    style1.cd()

    parser = argparse.ArgumentParser()
    parser.add_argument('--model', type=int, default='1', help="Which type of 2HDM?")
    parser.add_argument('--tanbeta', type=float, default='1', help="Which tan beta?")
    parser.add_argument('--ma', type=float, default='1', help="Which pseudoscalar mass?")

    args = parser.parse_args()

    limit_mmtt=1
    limit_mmbb=1

    x1, ymmtt1 = np.loadtxt('mmtt.txt', unpack=True)
    for i in range(0,len(x1)):
       if float(x1[i])==args.ma:
	  limit_mmtt1=ymmtt1[i]

    x1, ymmbb1 = np.loadtxt('mmbb.txt', unpack=True)
    for i in range(0,len(x1)):
       if float(x1[i])==args.ma:
          limit_mmbb1=ymmbb1[i]


    mintanbeta=0.1
    maxtanbeta=5
    n=100
    step=(maxtanbeta-mintanbeta)/n
    BRmm=1
    BRtt=1
    BRbb=1
    a_ymmbb=[]
    a_x=[]
    a_ymmtt=[]
    for i in range(0,100):
       tanbeta=mintanbeta+step*i
       width=get_total_width(args.model,float(args.ma),tanbeta)
       BRmm=gamma_mu(tanbeta,float(args.ma),args.model)/width
       BRtt=gamma_tau(tanbeta,args.ma,args.model)/width
       BRbb=gamma_quarks(tanbeta,args.ma,args.model,6)/width
       a_ymmbb.append(limit_mmbb1*0.00017/(2*BRbb*BRmm))
       a_ymmtt.append(limit_mmtt1/(BRtt*BRtt))
       #print limit_mmbb1,limit_mmtt1
       #print BRbb,BRmm,BRtt
       #print limit_mmbb1*0.00017/(2*BRbb*BRmm),limit_mmtt1/(BRtt*BRtt)
       a_x.append(tanbeta)

    x = array("d", a_x)
    ymmbb = array("d", a_ymmbb)
    ymmtt = array("d", a_ymmtt)

    gmmtt = ROOT.TGraph(len(x),x,ymmtt)
    gmmbb = ROOT.TGraph(len(x),x,ymmbb)

    ymax=100000
    gmmtt_shade = ROOT.TGraph(len(x))
    for i in range(0,len(x)):
      gmmtt_shade.SetPoint(i,x[i],ymax)
      gmmtt_shade.SetPoint(len(x)+i,x[len(x)-i-1],ymmtt[len(x)-i-1])
    gmmtt_shade.SetFillStyle(3001)
    adapt1=ROOT.gROOT.GetColor(ROOT.EColor.kPink+7)
    new_idx1=ROOT.gROOT.GetListOfColors().GetSize() + 1
    trans1=ROOT.TColor(new_idx1, adapt1.GetRed(), adapt1.GetGreen(),adapt1.GetBlue(), "",0.3)
    gmmtt_shade.SetFillColor(new_idx1)

    gmmbb_shade = ROOT.TGraph(len(x))
    for i in range(0,len(x)):
      gmmbb_shade.SetPoint(i,x[i],ymax)
      gmmbb_shade.SetPoint(len(x)+i,x[len(x)-i-1],ymmbb[len(x)-i-1])
    gmmbb_shade.SetFillStyle(3001)
    adapt2=ROOT.gROOT.GetColor(ROOT.EColor.kGreen-3)
    new_idx2=ROOT.gROOT.GetListOfColors().GetSize() + 1
    trans2=ROOT.TColor(new_idx2, adapt2.GetRed(), adapt2.GetGreen(),adapt2.GetBlue(), "",0.3)
    gmmbb_shade.SetFillColor(new_idx2)

    canvas = MakeCanvas("asdf","asdf",800,800)
    canvas.cd()
    canvas.SetLogy()
    gmmtt.Draw("AC")
    gmmtt.GetXaxis().SetRangeUser(0.1,5);
    gmmtt.GetXaxis().SetLimits(0.1,5);
    gmmtt.SetLineColor(ROOT.EColor.kPink+7)
    gmmtt.SetMinimum(0.005);
    gmmtt.SetMaximum(1000);
    if args.model==4:
       gmmtt.SetMinimum(0.1);
       gmmtt.SetMaximum(100000);
    gmmtt.GetXaxis().SetTitle("tan #beta");
    gmmtt.GetYaxis().SetTitle("95% CL on #frac{#sigma(h)}{#sigma_{SM}} #times BR(h#rightarrow aa)")
    gmmtt.Draw("AC")
    gmmtt_shade.Draw("fsame")
    canvas.Update();
    gmmbb.SetLineColor(ROOT.EColor.kGreen-3)
    gmmbb.Draw("lsame")
    gmmbb_shade.Draw("fsame")
    canvas.Update();

    legend=make_legend()
    legend.AddEntry(gmmtt,"h#rightarrow#mu#mu#tau#tau, HIG-15-011","l")
    legend.AddEntry(gmmbb,"h#rightarrow#mu#mubb, HIG-14-041","l")
    legend.Draw("same")

    lumiBlurb1=add_CMS()
    lumiBlurb1.Draw("same")
    lumiBlurb2=add_Preliminary()
    lumiBlurb2.Draw("same")
    lumiBlurb=add_lumi()
    lumiBlurb.Draw("same")
    lumi=add_model(args.model,args.ma)
    lumi.Draw("same")
    canvas.SaveAs('plots/plotType'+str(args.model)+'_ma'+str(int(args.ma))+'.png')
    canvas.SaveAs('plots/plotType'+str(args.model)+'_ma'+str(int(args.ma))+'.pdf')

