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
    output = ROOT.TLegend(0.42, 0.68, 0.93, 0.90, "", "brNDC")
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
        #lumi.AddText("tan #beta = "+str(tanbeta))
    return lumi

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('--model', type=int, default='1', help="Which type of 2HDM?")
    parser.add_argument('--tanbeta', type=float, default='1', help="Which tan beta?")
    parser.add_argument('--br', default='tau', help="As a function of which BR?")

    args = parser.parse_args()

    style1=GetStyleHtt()
    style1.cd()

    x_mmtt1, y_mmtt1 = np.loadtxt('mmtt.txt', unpack=True)
    x_mmtt=array("d",x_mmtt1)
    y_mmtt=array("d",y_mmtt1)
    for i in range(0,len(x_mmtt)):
        width=get_total_width(args.model,float(x_mmtt[i]),args.tanbeta)
        BRmm=gamma_mu(args.tanbeta,float(x_mmtt[i]),args.model)/width
        BRtt=gamma_tau(args.tanbeta,float(x_mmtt[i]),args.model)/width
        BRbb=gamma_quarks(args.tanbeta,float(x_mmtt[i]),args.model,6)/width
        y_mmtt[i]=y_mmtt[i]
	if args.br=="mu":
           y_mmtt[i]=y_mmtt[i]*BRmm*BRmm/(BRtt*BRtt)
    gmmtt = ROOT.TGraph(len(x_mmtt), x_mmtt,y_mmtt)

    x_mmbb, y_mmbb = np.loadtxt('mmbb.txt', unpack=True)
    for i in range(0,len(x_mmbb)):
        width=get_total_width(args.model,float(x_mmbb[i]),args.tanbeta)
        BRmm=gamma_mu(args.tanbeta,float(x_mmbb[i]),args.model)/width
        BRtt=gamma_tau(args.tanbeta,float(x_mmbb[i]),args.model)/width
        BRbb=gamma_quarks(args.tanbeta,float(x_mmbb[i]),args.model,6)/width
	y_mmbb[i]=y_mmbb[i]*0.00017
        if args.br=="mu":
           y_mmbb[i]=y_mmbb[i]*BRmm*BRmm/(2*BRmm*BRbb)
        if args.br=="tau":
           y_mmbb[i]=y_mmbb[i]*BRtt*BRtt/(2*BRmm*BRbb)
    gmmbb = ROOT.TGraph(len(x_mmbb), x_mmbb.flatten('C'),y_mmbb.flatten('C'))

    x_mmbb4, y_mmbb4 = np.loadtxt('mmbb.txt', unpack=True)
    for i in range(0,len(x_mmbb4)):
        width=get_total_width(args.model,float(x_mmbb4[i]),4)
        BRmm=gamma_mu(4,float(x_mmbb4[i]),args.model)/width
        BRtt=gamma_tau(4,float(x_mmbb4[i]),args.model)/width
        BRbb=gamma_quarks(4,float(x_mmbb4[i]),args.model,6)/width
        y_mmbb4[i]=y_mmbb4[i]*0.00017
        if args.br=="mu":
           y_mmbb4[i]=y_mmbb4[i]*BRmm*BRmm/(2*BRmm*BRbb)
        if args.br=="tau":
           y_mmbb4[i]=y_mmbb4[i]*BRtt*BRtt/(2*BRmm*BRbb)
    gmmbb4 = ROOT.TGraph(len(x_mmbb4), x_mmbb4.flatten('C'),y_mmbb4.flatten('C'))

    x_mmbb2, y_mmbb2 = np.loadtxt('mmbb.txt', unpack=True)
    for i in range(0,len(x_mmbb)):
        width=get_total_width(args.model,float(x_mmbb2[i]),2)
        BRmm=gamma_mu(2,float(x_mmbb2[i]),args.model)/width
        BRtt=gamma_tau(2,float(x_mmbb2[i]),args.model)/width
        BRbb=gamma_quarks(2,float(x_mmbb2[i]),args.model,6)/width
        y_mmbb2[i]=y_mmbb2[i]*0.00017
        if args.br=="mu":
           y_mmbb2[i]=y_mmbb[i]*BRmm*BRmm/(2*BRmm*BRbb)
        if args.br=="tau":
           y_mmbb2[i]=y_mmbb[i]*BRtt*BRtt/(2*BRmm*BRbb)
    gmmbb2 = ROOT.TGraph(len(x_mmbb2), x_mmbb2.flatten('C'),y_mmbb2.flatten('C'))

    x_tttt1, y_tttt1 = np.loadtxt('tttt1.txt', unpack=True)
    for i in range(0,len(x_tttt1)):
        width=get_total_width(args.model,float(x_tttt1[i]),args.tanbeta)
        BRmm=gamma_mu(args.tanbeta,float(x_tttt1[i]),args.model)/width
        BRtt=gamma_tau(args.tanbeta,float(x_tttt1[i]),args.model)/width
        BRbb=gamma_quarks(args.tanbeta,float(x_tttt1[i]),args.model,6)/width
        y_tttt1[i]=y_tttt1[i]/19.3
        if args.br=="mu":
           y_tttt1[i]=y_tttt1[i]*BRmm*BRmm/(BRtt*BRtt)
    gtttt1 = ROOT.TGraph(len(x_tttt1), x_tttt1.flatten('C'),y_tttt1.flatten('C'))

    x_tttt2, y_tttt2 = np.loadtxt('tttt2.txt', unpack=True)
    for i in range(0,len(x_tttt2)):
        y_tttt2[i]=y_tttt2[i]
        width=get_total_width(args.model,float(x_tttt2[i]),args.tanbeta)
        BRmm=gamma_mu(args.tanbeta,float(x_tttt2[i]),args.model)/width
        BRtt=gamma_tau(args.tanbeta,float(x_tttt2[i]),args.model)/width
        BRbb=gamma_quarks(args.tanbeta,float(x_tttt2[i]),args.model,6)/width
        if args.br=="mu":
           y_tttt2[i]=y_tttt2[i]*BRmm*BRmm/(BRtt*BRtt)
    gtttt2 = ROOT.TGraph(len(x_tttt2), x_tttt2.flatten('C'),y_tttt2.flatten('C'))

    x_mmmm, y_mmmm = np.loadtxt('mmmm.txt', unpack=True)
    for i in range(0,len(x_mmmm)):
        y_mmmm[i]=y_mmmm[i]/19300
    gmmmm = ROOT.TGraph(len(x_mmmm), x_mmmm.flatten('C'),y_mmmm.flatten('C'))

    x, y = np.loadtxt('dummy.txt', unpack=True)
    gx = ROOT.TGraph(len(x), x.flatten('C'),y.flatten('C'))

    ymax=100
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

    adapt1=ROOT.gROOT.GetColor(ROOT.EColor.kPink+7)

    canvas = MakeCanvas("asdf","asdf",800,800)
    canvas.cd()
    canvas.SetLogy()
    canvas.SetLogx()
    gx.Draw("AC")
    gx.GetXaxis().SetRangeUser(2,65);
    gx.GetXaxis().SetLimits(2,65);
    gx.SetMinimum(0.001);
    gx.SetMaximum(50);
    if (args.model==4):
       gx.SetMinimum(0.00001)
       gx.SetMaximum(100)
    if args.br=="mu":
       gx.GetXaxis().SetRangeUser(1,63);
       gx.GetXaxis().SetLimits(1,63);
       gx.SetMinimum(0.00000001);
       gx.SetMaximum(0.1);
    gx.GetXaxis().SetTitle("m_{a} (GeV)")
    gx.GetYaxis().SetTitleSize(0.048)
    gx.GetYaxis().SetTitleOffset(1.6)
    gx.GetYaxis().SetTitle("95% CL on #frac{#sigma(h)}{#sigma_{SM}} #times B(h#rightarrow aa) #times B(a#rightarrow #tau#tau)^{2}")
    if args.br=="mu":
	gx.GetYaxis().SetTitle("95% CL on #frac{#sigma(h)}{#sigma_{SM}} #times B(h#rightarrow aa) #times B(a#rightarrow #mu#mu)^{2}")
    gx.Draw("AC")
    canvas.Update();
    gmmtt.SetLineColor(ROOT.EColor.kPink+7)
    #gmmtt_shade.Draw("fsame")
    #gmmtt.Draw("lsame")
    #gmmtt.Draw("AC")
    #gmmtt.GetXaxis().SetRangeUser(0,63);
    #gmmtt.GetXaxis().SetLimits(0,63);
    #gmmtt.SetMinimum(0.001);
    #gmmtt.SetMaximum(1);
    #gmmtt.GetXaxis().SetTitle("m_{a} (GeV)");
    #gmmtt.GetYaxis().SetTitle("95% CL on BR(h#rightarrow aa) #times BR(a#rightarrow #tau#tau)^{2}")
    #gmmtt.Draw("AC")
    #canvas.Update();
    gmmbb.SetLineColor(ROOT.EColor.kGreen-3)
    gmmbb_shade.Draw("fsame")
    if args.model>2:
       gmmbb4.SetLineColor(ROOT.EColor.kGreen+3)
       gmmbb4.SetLineStyle(2)
       gmmbb4.Draw("lsame")
       gmmbb2.SetLineStyle(2)
       gmmbb2.SetLineColor(ROOT.EColor.kGreen)
       gmmbb2.Draw("lsame")
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
    if args.br=="mu":
       gmmmm_shade.Draw("fsame")
       gmmmm.Draw("lsame")
    canvas.Update();

    legend=make_legend()
    legend.AddEntry(gmmtt,"h#rightarrowaa#rightarrow#mu#mu#tau#tau, HIG-15-011","l")
    if args.model<3:
      legend.AddEntry(gmmbb,"h#rightarrowaa#rightarrow#mu#mubb, HIG-14-041","l")
    else:
      legend.AddEntry(gmmbb,"h#rightarrowaa#rightarrow#mu#mubb, HIG-14-041, tan#beta=1","l")
      legend.AddEntry(gmmbb2,"h#rightarrowaa#rightarrow#mu#mubb, HIG-14-041, tan#beta=2","l")
      legend.AddEntry(gmmbb4,"h#rightarrowaa#rightarrow#mu#mubb, HIG-14-041, tan#beta=4","l")
    legend.AddEntry(gtttt1,"h#rightarrowaa#rightarrow#tau#tau#tau#tau, HIG-14-019","l")
    legend.AddEntry(gtttt2,"h#rightarrowaa#rightarrow#tau#tau#tau#tau, HIG-14-022","l")
    if args.br=="mu":
       legend.AddEntry(gmmmm,"h#rightarrowaa#rightarrow#mu#mu#mu#mu, HIG-13-010","l")
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
    if (args.model>2):
        postfix="_tanbeta"+str(int(args.tanbeta))
    canvas.SaveAs('plots/plotType'+str(args.model)+'_BR'+args.br+postfix+'.png')
    canvas.SaveAs('plots/plotType'+str(args.model)+'_BR'+args.br+postfix+'.pdf')


