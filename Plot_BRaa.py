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
    output = ROOT.TLegend(0.52, 0.68, 0.95, 0.90, "", "brNDC")
    output.SetLineWidth(0)
    output.SetLineStyle(0)
    output.SetFillStyle(0)
    output.SetBorderSize(0)
    output.SetTextFont(62)
    return output

def add_model(model,tanbeta):
    lowX=0.21
    lowY=0.10
    lumi  = ROOT.TPaveText(lowX, lowY+0.06, lowX+0.15, lowY+0.16, "NDC")
    lumi.SetTextSize(0.04)
    lumi.SetBorderSize(   0 )
    lumi.SetFillStyle(    0 )
    lumi.SetTextAlign(   12 )
    lumi.SetTextColor(    1 )
    lumi.SetTextFont(62)
    if model<2:
       lumi.AddText("2HDM+S type-1")
    else:
        lumi.AddText("2HDM+S type-"+str(model))
        lumi.AddText("tan #beta = "+str(tanbeta))
    return lumi

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('--model', type=int, default='1', help="Which type of 2HDM?")
    parser.add_argument('--tanbeta', type=float, default='1', help="Which tan beta?")

    args = parser.parse_args()

    style1=GetStyleHtt()
    style1.cd()

    #### h->aa->mmtt ####
    x_mmtt1, y_mmtt1 = np.loadtxt('mmtt.txt', unpack=True)
    x_mmtt=array("d",x_mmtt1)
    y_mmtt=array("d",y_mmtt1)
    for i in range(0,len(x_mmtt)):
        width=get_total_width(args.model,float(x_mmtt[i]),args.tanbeta)
        BRtt=gamma_tau(args.tanbeta,float(x_mmtt[i]),args.model)/width
        y_mmtt[i]=y_mmtt[i]/(BRtt*BRtt)
    gmmtt = ROOT.TGraph(len(x_mmtt), x_mmtt,y_mmtt)

    #### h->aa->mmbb ####
    x_mmbb, y_mmbb = np.loadtxt('mmbb.txt', unpack=True)
    for i in range(0,len(x_mmbb)):
        width=get_total_width(args.model,float(x_mmbb[i]),args.tanbeta)
        BRmm=gamma_mu(args.tanbeta,float(x_mmbb[i]),args.model)/width
        BRbb=gamma_quarks(args.tanbeta,float(x_mmbb[i]),args.model,6)/width
	y_mmbb[i]=y_mmbb[i]*0.00017
        y_mmbb[i]=y_mmbb[i]/(2*BRmm*BRbb)
    gmmbb = ROOT.TGraph(len(x_mmbb), x_mmbb.flatten('C'),y_mmbb.flatten('C'))

    #### h->aa->tttt (HIG-14-019) ####
    x_tttt1, y_tttt1 = np.loadtxt('tttt1.txt', unpack=True)
    for i in range(0,len(x_tttt1)):
        width=get_total_width(args.model,float(x_tttt1[i]),args.tanbeta)
        BRtt=gamma_tau(args.tanbeta,float(x_tttt1[i]),args.model)/width
        y_tttt1[i]=y_tttt1[i]/(19.3*BRtt*BRtt)
    gtttt1 = ROOT.TGraph(len(x_tttt1), x_tttt1.flatten('C'),y_tttt1.flatten('C'))

    #### h->aa->tttt (HIG-14-022) ####
    x_tttt2, y_tttt2 = np.loadtxt('tttt2.txt', unpack=True)
    for i in range(0,len(x_tttt2)):
        width=get_total_width(args.model,float(x_tttt2[i]),args.tanbeta)
        BRtt=gamma_tau(args.tanbeta,float(x_tttt2[i]),args.model)/width
        y_tttt2[i]=y_tttt2[i]/(BRtt*BRtt)
    gtttt2 = ROOT.TGraph(len(x_tttt2), x_tttt2.flatten('C'),y_tttt2.flatten('C'))

    #### h->aa->mmmm ####
    x_mmmm, y_mmmm = np.loadtxt('mmmm.txt', unpack=True)
    for i in range(0,len(x_mmmm)):
        width=get_total_width(args.model,float(x_mmmm[i]),args.tanbeta)
        BRmm=gamma_mu(args.tanbeta,float(x_mmmm[i]),args.model)/width
        y_mmmm[i]=y_mmmm[i]/(19300*BRmm*BRmm)
    gmmmm = ROOT.TGraph(len(x_mmmm), x_mmmm.flatten('C'),y_mmmm.flatten('C'))

    #### Dummy histogram for plotting purposes ####
    x, y = np.loadtxt('dummy.txt', unpack=True)
    gx = ROOT.TGraph(len(x), x.flatten('C'),y.flatten('C'))

    #### Shaded areas above curves ####
    ymax=100000000000000000
    gmmtt_shade = ROOT.TGraph(len(x_mmtt))
    for i in range(0,len(x_mmtt)):
      gmmtt_shade.SetPoint(i,x_mmtt[i],ymax)
      gmmtt_shade.SetPoint(len(x_mmtt)+i,x_mmtt[len(x_mmtt)-i-1],y_mmtt[len(x_mmtt)-i-1])
    gmmtt_shade.SetFillStyle(3001)
    adapt1=ROOT.gROOT.GetColor(ROOT.EColor.kPink+7)
    new_idx1=ROOT.gROOT.GetListOfColors().GetSize() + 1
    trans1=ROOT.TColor(new_idx1, adapt1.GetRed(), adapt1.GetGreen(),adapt1.GetBlue(), "",0.3)
    gmmtt_shade.SetFillColor(new_idx1)

    gmmbb_shade = ROOT.TGraph(len(x_mmbb))
    for i in range(0,len(x_mmbb)):
      gmmbb_shade.SetPoint(i,x_mmbb[i],ymax)
      gmmbb_shade.SetPoint(len(x_mmbb)+i,x_mmbb[len(x_mmbb)-i-1],y_mmbb[len(x_mmbb)-i-1])
    gmmbb_shade.SetFillStyle(3001)
    adapt2=ROOT.gROOT.GetColor(ROOT.EColor.kGreen-3)
    new_idx2=ROOT.gROOT.GetListOfColors().GetSize() + 1
    trans2=ROOT.TColor(new_idx2, adapt2.GetRed(), adapt2.GetGreen(),adapt2.GetBlue(), "",0.3)
    gmmbb_shade.SetFillColor(new_idx2)

    gtttt1_shade = ROOT.TGraph(len(x_tttt1))
    for i in range(0,len(x_tttt1)):
      gtttt1_shade.SetPoint(i,x_tttt1[i],ymax)
      gtttt1_shade.SetPoint(len(x_tttt1)+i,x_tttt1[len(x_tttt1)-i-1],y_tttt1[len(x_tttt1)-i-1])
    gtttt1_shade.SetFillStyle(3001)
    adapt3=ROOT.gROOT.GetColor(ROOT.EColor.kBlue-3)
    new_idx3=ROOT.gROOT.GetListOfColors().GetSize() + 1
    trans3=ROOT.TColor(new_idx3, adapt3.GetRed(), adapt3.GetGreen(),adapt3.GetBlue(), "",0.3)
    gtttt1_shade.SetFillColor(new_idx3)

    gtttt2_shade = ROOT.TGraph(len(x_tttt2))
    for i in range(0,len(x_tttt2)):
      gtttt2_shade.SetPoint(i,x_tttt2[i],ymax)
      gtttt2_shade.SetPoint(len(x_tttt2)+i,x_tttt2[len(x_tttt2)-i-1],y_tttt2[len(x_tttt2)-i-1])
    gtttt2_shade.SetFillStyle(3001)
    adapt4=ROOT.gROOT.GetColor(ROOT.EColor.kOrange-3)
    new_idx4=ROOT.gROOT.GetListOfColors().GetSize() + 1
    trans4=ROOT.TColor(new_idx4, adapt4.GetRed(), adapt4.GetGreen(),adapt4.GetBlue(), "",0.3)
    gtttt2_shade.SetFillColor(new_idx4)

    gmmmm_shade = ROOT.TGraph(len(x_mmmm))
    for i in range(0,len(x_mmmm)):
      gmmmm_shade.SetPoint(i,x_mmmm[i],ymax)
      gmmmm_shade.SetPoint(len(x_mmmm)+i,x_mmmm[len(x_mmmm)-i-1],y_mmmm[len(x_mmmm)-i-1])
    gmmmm_shade.SetFillStyle(3001)
    adapt5=ROOT.gROOT.GetColor(ROOT.EColor.kViolet-3)
    new_idx5=ROOT.gROOT.GetListOfColors().GetSize() + 1
    trans5=ROOT.TColor(new_idx5, adapt5.GetRed(), adapt5.GetGreen(),adapt5.GetBlue(), "",0.3)
    gmmmm_shade.SetFillColor(new_idx5)
  
    #### Vertical shaded areas for not valid computations ####
    shade1=ROOT.TGraph(2)
    shade1.SetPoint(0,3,ymax)
    shade1.SetPoint(2,5.2,0.00000001)
    shade1.SetPoint(1,5.2,ymax)
    shade1.SetPoint(3,3,0.000000001)
    shade1.SetFillStyle(3012)
    shade1.SetFillColor(ROOT.EColor.kGray+1)

    shade2=ROOT.TGraph(2)
    shade2.SetPoint(0,9.5,ymax)
    shade2.SetPoint(2,11,0.0000000001)
    shade2.SetPoint(1,11,ymax)
    shade2.SetPoint(3,9.5,0.000000001)
    shade2.SetFillStyle(3012)
    shade2.SetFillColor(ROOT.EColor.kGray+1)

    adapt1=ROOT.gROOT.GetColor(ROOT.EColor.kPink+7)

    #### Plotting ####
    canvas = MakeCanvas("asdf","asdf",800,800)
    canvas.cd()
    canvas.SetLogy()
    canvas.SetLogx()
    gx.Draw("AC")
    gx.GetXaxis().SetRangeUser(1,65)
    gx.GetXaxis().SetLimits(1,65)
    maximum=max(max(y_mmmm),max(y_tttt1),max(y_tttt2),max(y_mmtt),max(y_mmbb))
    minimum=min(min(y_mmmm),min(y_tttt1),min(y_tttt2),min(y_mmtt),min(y_mmbb))
    gx.SetMinimum(minimum/100)
    if minimum>100:
      gx.SetMinimum(1)
    gx.SetMaximum(maximum*100)
    if maximum>1000000:
      gx.SetMaximum(1000000)
    if max(y_tttt2)>10000:
      gx.SetMaximum(10000*max(y_tttt2))
    #gx.SetMinimum(0.03);
    #gx.SetMaximum(500);
    #if args.model<3:
    #   gx.SetMinimum(0.005)
    #   gx.SetMaximum(10000)
    gx.GetXaxis().SetTitle("m_{a} (GeV)");
    gx.GetYaxis().SetTitle("95% CL on #frac{#sigma(h)}{#sigma_{SM}} #times B(h#rightarrow aa)")
    gx.Draw("AC")
    canvas.Update();
    gmmtt.SetLineColor(ROOT.EColor.kPink+7)
    gmmbb.SetLineColor(ROOT.EColor.kGreen-3)
    shade1.Draw("fsame")
    shade2.Draw("fsame")
    gmmbb_shade.Draw("fsame")
    gmmbb.Draw("lsame")
    gmmtt_shade.Draw("fsame")
    gmmtt.Draw("lsame")

    gtttt1.SetLineColor(ROOT.EColor.kBlue-3)
    gtttt1_shade.Draw("fsame")
    gtttt1.Draw("lsame")
    gtttt2.SetLineColor(ROOT.EColor.kOrange-3)
    gtttt2_shade.Draw("fsame")
    gtttt2.Draw("lsame")
    gmmmm.SetLineColor(ROOT.EColor.kViolet-3)
    gmmmm_shade.Draw("fsame")
    gmmmm.Draw("lsame")
    line = ROOT.TLine(1,1,62.5,1)
    line.SetLineStyle(2)
    line.SetLineColor(ROOT.EColor.kRed)
    line.Draw("lsame");
    canvas.Update();

    legend=make_legend()
    legend.AddEntry(gmmtt,"h#rightarrowaa#rightarrow#mu#mu#tau#tau, HIG-15-011","l")
    legend.AddEntry(gmmbb,"h#rightarrowaa#rightarrow#mu#mubb, HIG-14-041","l")
    legend.AddEntry(gtttt1,"h#rightarrowaa#rightarrow#tau#tau#tau#tau, HIG-14-019","l")
    legend.AddEntry(gtttt2,"h#rightarrowaa#rightarrow#tau#tau#tau#tau, HIG-14-022","l")
    legend.AddEntry(gmmmm,"h#rightarrowaa#rightarrow#mu#mu#mu#mu, HIG-13-010","l")
    legend.AddEntry(shade1,"B(a#rightarrow XX) values not valid","f") 
    legend.Draw("same")

    lumiBlurb1=add_CMS()
    lumiBlurb1.Draw("same")
    lumiBlurb2=add_Preliminary()
    lumiBlurb2.Draw("same")
    lumiBlurb=add_lumi()
    lumiBlurb.Draw("same")
    lumi=add_model(args.model,args.tanbeta)
    lumi.Draw("same")
    postfix=""
    if (args.model>1):
        postfix="_tanbeta"+str(int(args.tanbeta))
    canvas.SaveAs('plots/plot_BRaa_Type'+str(args.model)+postfix+'.png')
    canvas.SaveAs('plots/plot_BRaa_Type'+str(args.model)+postfix+'.pdf')


