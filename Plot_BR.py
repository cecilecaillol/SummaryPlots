from RecoLuminosity.LumiDB import argparse
import math
from HttStyles import GetStyleHtt
from HttStyles import MakeCanvas
import ROOT
import numpy as np
from array import array


def get_factor(model,particles,tanbeta):
   if model==1:
	if particles==1:#up-type quarks
	    return (1.0/tanbeta)
	elif particles==2: #down-type quarks
	    return (-1.0/tanbeta)
	elif particles==3: #leptons
	    return (-1.0/tanbeta)
   elif model==2:
        if particles==1:#up-type quarks
            return (1.0/tanbeta)
        elif particles==2: #down-type quarks
            return (tanbeta)
        elif particles==3: #leptons
            return (tanbeta)
   elif model==3:
        if particles==1:#up-type quarks
            return (1.0/tanbeta)
        elif particles==2: #down-type quarks
            return (-1.0/tanbeta)
        elif particles==3: #leptons
            return (tanbeta)
   elif model==4:
        if particles==1:#up-type quarks
            return (1.0/tanbeta)
        elif particles==2: #down-type quarks
            return (tanbeta)
        elif particles==3: #leptons
            return (-1.0/tanbeta)

def gamma_tau(tanbeta,ma,model):
   G_mu=1.16637*0.00001
   mtau=1.77682
   xb=get_factor(model,3,tanbeta) 
   gamma=0
   if ma>2*mtau:
      gamma=(G_mu/(4*(2**0.5)*math.pi))*xb*xb*ma*mtau*mtau*(1-4*((mtau*mtau)/(ma*ma)))**0.5
   return gamma

def gamma_mu(tanbeta,ma,model):
   G_mu=1.16637*0.00001
   mmu=0.1056584
   xb=get_factor(model,3,tanbeta) 
   gamma=(G_mu/(4*(2**0.5)*math.pi))*xb*xb*ma*mmu*mmu*(1-4*((mmu*mmu)/(ma*ma)))**0.5
   return gamma

def running_alpha(scale):
   Nf=6
   if scale<172.5:
      Nf=5
   if scale<4.75:
      Nf=4
   if scale<1.29:
      Nf=3
   if scale<0.130:
      Nf=2
   if scale<0.0048:
      Nf=1
   if scale<0.0023:
      Nf=0
   beta0=11-(2.0/3)*Nf
   beta1=51-(19.0/3)*Nf
   beta2=2857-(5033.0/9)*Nf+(325.0/7)*Nf*Nf
   Lambda=0.216
   if scale<0.216:
      scale=0.46
   #scale=5
   lmu=math.log((1.0*scale*scale)/(Lambda*Lambda))
   run=((4*math.pi)/(beta0*lmu))*(1-(2*(beta1/(beta0*beta0))*(math.log(lmu)/lmu))+((4*beta1*beta1)/(beta0*beta0*beta0*beta0*lmu*lmu))*((math.log(lmu)-0.5)*(math.log(lmu)-0.5)+((beta2*beta0)/(8*beta1*beta1))-(5.0/4)))
   return run

def running_mass(scale,mq):
   Nf=6
   if scale<172.5:
      Nf=5
   if scale<4.75:
      Nf=4
   if scale<1.29:
      Nf=3
   if scale<0.130:
      Nf=2
   if scale<0.0048:
      Nf=1
   if scale<0.0023:
      Nf=0
   #if (mq<1):
   #   running=mq
   #else:
   running=1
   if mq>0.216:
      running=mq*(1-(4.0/3)*(running_alpha(mq)/math.pi)+(1.0414*Nf-14.3323)*((running_alpha(mq)/math.pi)*(running_alpha(mq)/math.pi))+(-0.65269*Nf*Nf+26.9239*Nf-198.7068)*(running_alpha(mq)/math.pi)*(running_alpha(mq)/math.pi)*(running_alpha(mq)/math.pi))
   if mq==0.130:
      running=0.2
   if mq<0.12:
      running=mq
   return running

def gamma_quarks(tanbeta,ma,model,quark):
   G_mu=1.16637*0.00001
   xb=1
   if (quark==1 or quark==3 or quark==5):
      xb=get_factor(model,1,tanbeta)
   else:
      xb=get_factor(model,2,tanbeta)
   mbpole=4.75#4.88
   mtpole=172.5#178
   mcpole=1.29#1.64
   mupole=0.0025
   mdpole=0.005
   mspole=0.130
   #############
   #mbpole=4.66#4.88
   #mtpole=173.21#178
   #mcpole=1.275#1.64
   #mupole=0.0023
   #mdpole=0.0048
   #mspole=0.130
   #############
   mq=1
   if quark==1:
     mq=mupole
   elif quark==2:
     mq=mdpole
   elif quark==3:
     mq=mcpole
   elif quark==4:
     mq=mspole
   elif quark==5:
     mq=mtpole
   elif quark==6:
     mq=mbpole
   scale=ma #renormalization scale strong coupling constant
   Nf=5 #number of active light quarks
   if ma<4.75:
      Nf=4
   elif ma<1.29:
      Nf=3
   qcdcorrections=1+5.67*(running_alpha(scale)/math.pi)+(35.94-1.35*Nf)*(running_alpha(scale)/math.pi)*(running_alpha(scale)/math.pi)+(running_alpha(scale)/math.pi)*(running_alpha(scale)/math.pi)*(3.83-math.log(1.0*ma*ma/(mtpole*mtpole))+(1.0/6)*math.log(running_mass(ma,mq)*running_mass(ma,mq)/(ma*ma))*math.log(running_mass(ma,mq)*running_mass(ma,mq)/(ma*ma)))
   gamma=0
   #print qcdcorrections
   if (2*mq)<ma:
      gamma=(3*G_mu/(4*(2**0.5)*math.pi))*xb*xb*ma*running_mass(ma,mq)*running_mass(ma,mq)*((1-(4.0*mq*mq)/(ma*ma))**0.5)*qcdcorrections
      #gamma=(3*G_mu/(4*(2**0.5)*math.pi))*xb*xb*ma*running_mass(ma,mq)*running_mass(ma,mq)*qcdcorrections
   return gamma

def get_f(x):
   if x<=1:
      return math.asin((1.0*x)**0.5)*math.asin((1.0*x)**0.5)
   else:
      return (-1.0/4)*(math.log((1+(1-1.0/x)**0.5)/(1-(1-1.0/x)**0.5))-1j*math.pi)*(math.log((1+(1-1.0/x)**0.5)/(1-(1-1.0/x)**0.5))-1j*math.pi)

def get_A12(x):
   return (2.0*get_f(x)/x)

def gamma_photon(tanbeta,ma,model):
   G_mu=1.16637*0.00001
   alpha=1.0/137
   Nc=3
   qu=2.0/3
   qd=-1.0/3
   mt=172.5
   mb=4.75
   xbu=get_factor(model,1,tanbeta)
   xbd=get_factor(model,2,tanbeta)
   termt=Nc*qu*qu*xbu*get_A12((1.0*ma*ma)/(4*mt*mt))
   termb=Nc*qd*qd*xbd*get_A12((1.0*ma*ma)/(4*mb*mb))
   gamma=(G_mu*alpha*alpha*ma*ma*ma/(128*(2**0.5)*math.pi*math.pi*math.pi))*abs(termt+termb)*abs(termt+termb)
   return gamma

def gamma_gg(tanbeta,ma,model):
   G_mu=1.16637*0.00001
   Nf=6
   if ma<172.5:
      Nf=5
   if ma<4.75:
      Nf=4
   if ma<1.29:
      Nf=3
   if ma<0.130:
      Nf=2
   if ma<0.0048:#0.005:
      Nf=1
   if ma<0.0023:#0.0025:
      Nf=0
   mt=172.5
   mb=4.75
   mc=1.29
   xbu=get_factor(model,1,tanbeta)
   xbd=get_factor(model,2,tanbeta)
   termt=1.0*xbu*get_A12((1.0*ma*ma)/(4*mt*mt))
   termb=1.0*xbd*get_A12((1.0*ma*ma)/(4*mb*mb))
   termc=1.0*xbu*get_A12((1.0*ma*ma)/(4*mc*mc))
   gamma=(G_mu*running_alpha(ma)*running_alpha(ma)*ma*ma*ma/(36*(2**0.5)*math.pi*math.pi*math.pi))*abs(3.0*(termt+termb+termc)/4)*abs(3.0*(termt+termb+termc)/4)*(1+((97.0/4)-(7.0/6)*Nf)*(running_alpha(ma)/math.pi))
   return gamma

def get_total_width(model,ma,tanbeta):
   return (gamma_quarks(tanbeta,ma,model,1)+gamma_quarks(tanbeta,ma,model,2)+gamma_quarks(tanbeta,ma,model,3)+gamma_quarks(tanbeta,ma,model,4)+gamma_quarks(tanbeta,ma,model,5)+gamma_quarks(tanbeta,ma,model,6)+gamma_tau(tanbeta,ma,model)+gamma_mu(tanbeta,ma,model)+gamma_gg(tanbeta,ma,model)+gamma_photon(tanbeta,ma,model))

def make_legend():
   output = ROOT.TLegend(0.50, 0.73, 0.89, 0.90, "", "brNDC")
   output.SetNColumns(3)
   output.SetLineWidth(0)
   output.SetLineStyle(0)
   output.SetFillStyle(0)
   output.SetBorderSize(0)
   output.SetTextFont(62)
   return output

def add_model(model,tanbeta):
    lowX=0.21
    lowY=0.20
    lumi  = ROOT.TPaveText(lowX, lowY+0.06, lowX+0.15, lowY+0.16, "NDC")
    lumi.SetTextSize(0.04)
    lumi.SetBorderSize(   0 )
    lumi.SetFillStyle(    0 )
    lumi.SetTextAlign(   12 )
    lumi.SetTextColor(    1 )
    lumi.SetTextFont(62)
    #lumi.AddText("2HDM+S type "+str(model)+", tan#beta = "+str(tanbeta))
    if model!=1:
       lumi.AddText("2HDM+S type-"+str(model)+",")
    if model==1:
       lumi.AddText("2HDM+S type-"+str(model))
    if model!=1:
       lumi.AddText("tan#beta = "+str(tanbeta))
    return lumi

parser = argparse.ArgumentParser()

parser.add_argument('--model', type=int, default='1', help="Which type of 2HDM?")
parser.add_argument('--tanbeta', type=float, default='1', help="Which tan beta?")
parser.add_argument('--ma', type=int, default='1', help="Which pseudoscalar mass?")
args = parser.parse_args()

model=args.model
ma=args.ma
tanbeta=args.tanbeta

style=GetStyleHtt()
style.cd()

#model=1
#ma=40
#tanbeta=5
width=get_total_width(model,ma,tanbeta)
print "running alpha ",running_alpha(ma)
print "up ",gamma_quarks(tanbeta,ma,model,1)/width
print "down ",gamma_quarks(tanbeta,ma,model,2)/width
print "charm ",gamma_quarks(tanbeta,ma,model,3)/width
print "strange ",gamma_quarks(tanbeta,ma,model,4)/width
print "top ",gamma_quarks(tanbeta,ma,model,5)/width
print "bottom ",gamma_quarks(tanbeta,ma,model,6)/width
print "tau ",gamma_tau(tanbeta,ma,model)/width
print "muon ",gamma_mu(tanbeta,ma,model)/width
print "gluon ",gamma_gg(tanbeta,ma,model)/width
print "photon ",gamma_photon(tanbeta,ma,model)/width

ptau=[]
pmu=[]
pu=[]
pd=[]
pb=[]
ps=[]
pc=[]
pg=[]
pp=[]
xtau=[]
xmu=[]
xu=[]
xd=[]
xb=[]
xs=[]
xc=[]
xg=[]
xp=[]

xa=[]

#for m in range(1,64):
for m in [float(j) / 10 for j in range(10, 625, 1)]:
   width=get_total_width(model,m,tanbeta)
   if gamma_tau(tanbeta,m,model)/width > 0.00001:
      ptau.append(gamma_tau(tanbeta,m,model)/width)
      xtau.append(m)
   if gamma_mu(tanbeta,m,model)/width > 0.00001:
      pmu.append(gamma_mu(tanbeta,m,model)/width)
      xmu.append(m)
   if (gamma_quarks(tanbeta,m,model,6)/width)>0.00001:
      pb.append(gamma_quarks(tanbeta,m,model,6)/width)
      xb.append(m)
   if gamma_quarks(tanbeta,m,model,1)/width > 0.00001:
      pu.append(gamma_quarks(tanbeta,m,model,1)/width)
      xu.append(m)
   if gamma_quarks(tanbeta,m,model,2)/width:
      pd.append(gamma_quarks(tanbeta,m,model,2)/width)
      xd.append(m)
   if gamma_quarks(tanbeta,m,model,3)/width > 0.00001:
      pc.append(gamma_quarks(tanbeta,m,model,3)/width)
      xc.append(m)
   if gamma_quarks(tanbeta,m,model,4)/width > 0.00001:
      ps.append(gamma_quarks(tanbeta,m,model,4)/width)
      xs.append(m)
   if (gamma_gg(tanbeta,m,model)/width) > 0.00001:
      pg.append(gamma_gg(tanbeta,m,model)/width)
      xg.append(m)
   if gamma_photon(tanbeta,m,model)/width > 0.00001:
      pp.append(gamma_photon(tanbeta,m,model)/width)
      xp.append(m)
   xa.append(m)


if len(xb) == 0:
  xb.append(10)
  pb.append(0)
if len(xc) == 0:
  xc.append(10)
  pc.append(0)
if len(xs) == 0:
  xs.append(10)
  ps.append(0)
if len(xg) == 0:
  xg.append(10)
  pg.append(0)
if len(xp) == 0:
  xp.append(10)
  pp.append(0)
if len(xtau) == 0:
  xtau.append(10)
  ptau.append(0)
if len(xmu) == 0:
  xmu.append(10)
  pmu.append(0)
if len(xd) == 0:
  xd.append(10)
  pd.append(0)
if len(xu) == 0:
  xu.append(10)
  pu.append(0)

x = array("d", xa)
xxtau = array("d",xtau)
xxmu = array("d",xmu)
xxb = array("d",xb)
xxs = array("d",xs)
xxc = array("d",xc)
xxd = array("d",xd)
xxu = array("d",xu)
xxg = array("d",xg)
xxp = array("d",xp)
ytau = array("d",ptau)  
gtau = ROOT.TGraph(len(xxtau),xxtau,ytau) 
ymu = array("d",pmu)        
gmu = ROOT.TGraph(len(xxmu),xxmu,ymu)
yu = array("d",pu)        
gu = ROOT.TGraph(len(xxu),xxu,yu)
yd = array("d",pd)        
gd = ROOT.TGraph(len(xxd),xxd,yd)
ys = array("d",ps)        
gs = ROOT.TGraph(len(xxs),xxs,ys)
yc = array("d",pc)        
gc = ROOT.TGraph(len(xxc),xxc,yc)
yb = array("d",pb)        
gb = ROOT.TGraph(len(xxb),xxb,yb)
yg = array("d",pg)        
gg = ROOT.TGraph(len(xxg),xxg,yg)
yp = array("d",pp)        
gp = ROOT.TGraph(len(xxp),xxp,yp)

canvas = MakeCanvas("asdf","asdf",800,800)
canvas.cd()
canvas.SetLogy()
canvas.SetLogx()
canvas.SetTicks(1,1)
gtau.Draw("AC")
gtau.GetXaxis().SetRangeUser(1,63);
gtau.GetXaxis().SetLimits(1,63);
gtau.SetLineColor(ROOT.EColor.kGreen+2)
gtau.SetLineWidth(4)
gtau.SetMinimum(0.0001);
gtau.SetMaximum(1.1);
gtau.GetXaxis().SetTitle("m_{a} (GeV)");
gtau.GetYaxis().SetTitle("B(a#rightarrow XX)")
gtau.GetYaxis().SetTitleOffset(1.4)
gtau.SetTitle("")
gtau.Draw("AC")
canvas.Update();
gmu.SetLineColor(ROOT.EColor.kPink-3)
gmu.SetLineWidth(4)
gmu.Draw("same")
gu.SetLineColor(ROOT.EColor.kYellow)
gu.SetLineWidth(4)
gu.Draw("same")
gd.SetLineColor(ROOT.EColor.kGreen-7)
gd.SetLineWidth(4)
gd.Draw("same")
gc.SetLineColor(ROOT.EColor.kAzure+5)
gc.SetLineWidth(4)
gc.Draw("same")
gs.SetLineColor(ROOT.EColor.kViolet-3)
gs.SetLineWidth(4)
gs.Draw("same")
gb.SetLineColor(ROOT.EColor.kBlue+2)
gb.SetLineWidth(4)
gb.Draw("same")
gg.SetLineColor(ROOT.EColor.kOrange-3)
gg.SetLineWidth(4)
gg.Draw("same")
gp.SetLineColor(ROOT.EColor.kRed-3)
gp.SetLineWidth(4)
gp.Draw("same")
canvas.Update();
legend=make_legend()
legend.AddEntry(gtau,"#tau#tau","l")
legend.AddEntry(gmu,"#mu#mu","l")
legend.AddEntry(gu,"uu","l")
legend.AddEntry(gd,"dd","l")
legend.AddEntry(gc,"cc","l")
legend.AddEntry(gs,"ss","l")
legend.AddEntry(gb,"bb","l")
legend.AddEntry(gg,"gg","l")
legend.AddEntry(gp,"#gamma#gamma","l")
legend.Draw("same")
lumi=add_model(model,tanbeta)
lumi.Draw("same")
canvas.SaveAs('plots/BR_2HDMS_type'+str(model)+'_tanbeta'+str(int(tanbeta))+'.png')
canvas.SaveAs('plots/BR_2HDMS_type'+str(model)+'_tanbeta'+str(int(tanbeta))+'.pdf')

