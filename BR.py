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
   G_mu=1.16637*0.00001 #Fermi constant
   mtau=1.77682
   xb=get_factor(model,3,tanbeta) 
   gamma=0
   if ma>2*mtau:
      gamma=(G_mu/(4*(2**0.5)*math.pi))*xb*xb*ma*mtau*mtau*(1-4*((mtau*mtau)/(ma*ma)))**0.5
   return gamma

def gamma_mu(tanbeta,ma,model):
   G_mu=1.16637*0.00001 #Fermi constant
   mmu=0.10566
   xb=get_factor(model,3,tanbeta) 
   gamma=(G_mu/(4*(2**0.5)*math.pi))*xb*xb*ma*mmu*mmu*(1-4*((mmu*mmu)/(ma*ma)))**0.5
   return gamma

def running_alpha(scale): #running of strong coupling constant at a given scale
   Nf=6 #number of active flavors
   if scale<172.5:
      Nf=5
   if scale<4.75:
      Nf=4
   if scale<1.29:
      Nf=3
   if scale<0.130:
      Nf=2
   if scale<0.005:
      Nf=1
   if scale<0.0025:
      Nf=0
   beta0=11-(2.0/3)*Nf
   beta1=51-(19.0/3)*Nf
   beta2=2857-(5033.0/9)*Nf+(325.0/7)*Nf*Nf
   Lambda=0.216 #QCD scale
   if scale<0.216:
      scale=0.46
   #scale=5
   lmu=math.log((1.0*scale*scale)/(Lambda*Lambda))
   run=((4*math.pi)/(beta0*lmu))*(1-(2*(beta1/(beta0*beta0))*(math.log(lmu)/lmu))+((4*beta1*beta1)/(beta0*beta0*beta0*beta0*lmu*lmu))*((math.log(lmu)-0.5)*(math.log(lmu)-0.5)+((beta2*beta0)/(8*beta1*beta1))-(5.0/4)))
   return run

def running_mass(scale,mq): #quark running mass at a given scale
   Nf=6 #number of active flavors
   if scale<172.5:
      Nf=5
   if scale<4.75:
      Nf=4
   if scale<1.29:
      Nf=3
   if scale<0.130:
      Nf=2
   if scale<0.005:
      Nf=1
   if scale<0.0025:
      Nf=0
   #if (mq<1):
   #   running=mq
   #else:
   running=1
   if mq>0.216:
      running=mq*(1-(4.0/3)*(running_alpha(mq)/math.pi)+(1.0414*Nf-14.3323)*((running_alpha(mq)/math.pi)*(running_alpha(mq)/math.pi))+(-0.65269*Nf*Nf+26.9239*Nf-198.7068)*(running_alpha(mq)/math.pi)*(running_alpha(mq)/math.pi)*(running_alpha(mq)/math.pi))
   if mq==0.130: #special case for c quark
      running=0.2
   if mq<0.12: #no correction for up and down quarks
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
      return ((math.asin(1.0*x**0.5))*(math.asin(1.0*x**0.5)))
   else:
      return (-1.0/4)*(math.log((1+(1-1.0/x)**0.5)/(1-(1-1.0/x)**0.5))-1j*math.pi)*(math.log((1+(1-1.0/x)**0.5)/(1-(1-1.0/x)**0.5))-1j*math.pi)

def get_A12(x):
   return (2.0*get_f(x)/x)

def gamma_photon(tanbeta,ma,model):#decay to photons through top and bottom quark loops
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

def gamma_gg(tanbeta,ma,model): #decay to gluons through t-, b- and c-loops
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
   if ma<0.005:
      Nf=1
   if ma<0.0025:
      Nf=0
   mt=172.5
   mb=4.75
   mc=1.29
   xbu=get_factor(model,1,tanbeta)
   xbd=get_factor(model,2,tanbeta)
   termt=xbu*get_A12((1.0*ma*ma)/(4*mt*mt))
   termb=xbd*get_A12((1.0*ma*ma)/(4*mb*mb))
   termc=xbu*get_A12((1.0*ma*ma)/(4*mc*mc))
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


