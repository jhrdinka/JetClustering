import ROOT
import sys
import os
from ROOT import TFile, TH1F, TGraphErrors, TMultiGraph, TLegend, TF1, gROOT

c = ROOT.TCanvas("", "", 600, 600) 
ROOT.gROOT.SetBatch(True)

font = 132

colors = []
colors.append(ROOT.kRed);
colors.append(ROOT.kBlue-3);
colors.append(ROOT.kGreen+2);
colors.append(ROOT.kOrange-3);
colors.append(ROOT.kYellow+2);
colors.append(ROOT.kMagenta+1);
colors.append(ROOT.kGray+1);

styles = []
styles.append(1)
styles.append(2)
styles.append(3)
styles.append(7)
styles.append(9)
styles.append(10)
styles.append(20)

#____________________________________________________
def myStyle():
    ROOT.gStyle.SetFrameBorderMode(0)
    ROOT.gStyle.SetCanvasBorderMode(0)
    ROOT.gStyle.SetPadBorderMode(0)

    ROOT.gStyle.SetFrameFillColor(0)
    ROOT.gStyle.SetPadColor(0)
    ROOT.gStyle.SetCanvasColor(0)
    ROOT.gStyle.SetTitleColor(1)
    ROOT.gStyle.SetStatColor(0)

    # set the paper & margin sizes
    ROOT.gStyle.SetPaperSize(20,26)
    ROOT.gStyle.SetPadTopMargin(0.10)
    ROOT.gStyle.SetPadRightMargin(0.03)
    ROOT.gStyle.SetPadBottomMargin(0.13)
    ROOT.gStyle.SetPadLeftMargin(0.125)
    ROOT.gStyle.SetPadTickX(1)
    ROOT.gStyle.SetPadTickY(1)
    
    ROOT.gStyle.SetTextFont(42) #132
    ROOT.gStyle.SetTextSize(0.09)
    ROOT.gStyle.SetLabelFont(42,"xyz")
    ROOT.gStyle.SetTitleFont(42,"xyz")
    ROOT.gStyle.SetLabelSize(0.045,"xyz") #0.035
    ROOT.gStyle.SetTitleSize(0.045,"xyz")
    ROOT.gStyle.SetTitleOffset(1.15,"y")
    
    # use bold lines and markers
    ROOT.gStyle.SetMarkerStyle(8)
    ROOT.gStyle.SetHistLineWidth(2)
    ROOT.gStyle.SetLineWidth(1)
    #ROOT.gStyle.SetLineStyleString(2,"[12 12]") // postscript dashes

    # do not display any of the standard histogram decorations
    ROOT.gStyle.SetOptTitle(1)
    ROOT.gStyle.SetOptStat(0) #("m")
    ROOT.gStyle.SetOptFit(0)
    
    #ROOT.gStyle.SetPalette(1,0)
    ROOT.gStyle.cd()
    ROOT.gROOT.ForceStyle()
#_____________________________________________________________________________________________________

def printResoHisto(histo,f, pt, eta, r, pu):
    print 'draw histogram with name: ',histo.GetName() 
    ROOT.gStyle.SetOptFit(1)
    c = ROOT.TCanvas(histo.GetName(), "", 600, 600)
    ROOT.gPad.SetLeftMargin(0.20) ; 
    ROOT.gPad.SetRightMargin(0.05) ; 
    ROOT.gPad.SetBottomMargin(0.20) ; 
    ROOT.gStyle.SetOptStat(0000000);
    ROOT.gStyle.SetTextFont(132);
    ROOT.gStyle.SetOptStat(0000000)
    c.cd()
    histo.SetLineWidth(3)
    histo.Draw()
    f.SetLineWidth(3)
    f.Draw('same')

    histo.SetMaximum(1.6*histo.GetMaximum())
    histo.SetMinimum(0.)
    histo.GetXaxis().SetRangeUser(0., 1.5)
   
    histo.GetYaxis().SetLabelFont(132)
    histo.GetYaxis().SetTitleFont(132)
    histo.GetXaxis().SetLabelFont(132)
    histo.GetXaxis().SetTitleFont(132)
   
    histo.GetXaxis().SetTitleOffset(1.2)
    histo.GetYaxis().SetTitleOffset(1.2)
    histo.GetXaxis().SetLabelOffset(0.02)
    histo.GetYaxis().SetLabelOffset(0.02)
    histo.GetXaxis().SetTitleSize(0.06)
    histo.GetYaxis().SetTitleSize(0.06)
    histo.GetXaxis().SetLabelSize(0.06)
    histo.GetYaxis().SetLabelSize(0.06)
    histo.GetXaxis().SetNdivisions(505)
    histo.GetYaxis().SetNdivisions(505)
   
   
    ptext = ROOT.TPaveText(0.21,0.62,0.50,.87)
    ptext.AddText("{0} pile-up".format(pu))
    ptext.AddText("anti-k_{{T}}, R = {0}".format(r))
    ptext.AddText('{0} < |#eta| < {1} '.format(eta[0], eta[1]))
    ptext.AddText("p_{{T}} = {0} GeV".format(pt))
    
    ptext.SetTextFont(132)
    ptext.SetTextSize(0.035) 
    ptext.SetFillColor(0)
    ptext.SetFillStyle(0)
    ptext.SetLineColor(0)
    ptext.SetBorderSize(1)
    ptext.Paint('NDC')
    ptext.Draw()
    
    '''Tleft = ROOT.TLatex(0.23, 0.92, lt) 
    Tleft.SetNDC(ROOT.kTRUE) 
    Tleft.SetTextSize(0.044) 
    Tleft.SetTextFont(132) '''
    
    ps = histo.FindObject('stats')
    #ps.SetX1NDC(0.60)
    #ps.SetX2NDC(0.95)
    #ps.SetY1NDC(0.70)
    #ps.SetY2NDC(0.90)
    c.Modified()
    c.Update()
    
    histo.GetXaxis().SetTitle('p_{T}^{reco} / p_{T}^{gen}')
    
    basename = os.path.basename(histo.GetName())
    if not os.path.exists('plots'):
        os.makedirs('plots')
    plotname = 'plots/'+basename+'.png'
    print plotname

    c.Print(plotname)
        
#_____________________________________________________________________________________________________
def drawMultiGraph(mg, name, lt, rt, pdir, xmin, xmax, ymin, ymax, logx, logy, bl, f):
    #myStyle()
    ROOT.gStyle.SetOptStat(0000000)
    
    canvas = ROOT.TCanvas(name, '', 600,600) 
    
    ROOT.gPad.SetLeftMargin(0.20) ; 
    ROOT.gPad.SetRightMargin(0.05) ; 
    ROOT.gPad.SetBottomMargin(0.20) ; 
    ROOT.gStyle.SetOptStat(0000000);
    ROOT.gStyle.SetOptFit(000000)
    ROOT.gStyle.SetTextFont(132);
   
    Tright = ROOT.TLatex(0.23, 0.92, rt) 
    Tright.SetNDC(ROOT.kTRUE) 
    Tright.SetTextSize(0.044) 
    Tright.SetTextFont(132) 

    Tleft = ROOT.TPaveText(0.26,0.65, 0.50,0.89, 'NDC')
    #Tleft.SetNDC(ROOT.kTRUE) 
    for txt in lt:
       Tleft.AddText(txt)
    Tleft.SetTextFont(132)
    #Tleft.SetTextSize(0.044) 
    Tleft.SetFillColor(0)
    Tleft.SetFillStyle(0)
    Tleft.SetLineColor(0)
    Tleft.SetBorderSize(1)
    Tleft.Paint('NDC')
    Tleft.Draw()

    canvas.cd(0)
    
    if logx: ROOT.gPad.SetLogx()
    if logy: ROOT.gPad.SetLogy()

    # fit stuff
    if f:
        mg.Draw("AP")
        '''
        strA = '{:.2f}'.format(f.GetParameter(0)/100.)
        strB = '{:.0f}'.format(f.GetParameter(1))
        '''
        
        strA = '{0:.0f}'.format(f.GetParameter(0))
        strB = '{0:.1f}'.format(f.GetParameter(1))
        
       #fittext = '#frac{#sigma(p_{T})}{p_{T}} = #frac{A [GeV]}{p_{T}} #oplus #frac{B%}{#sqrt{p_{T}}} #oplus C%'
        fittext = '#frac{#sigma(p_{T})}{p_{T}} = #frac{A%}{#sqrt{p_{T}}} #oplus B%'
        
        fittext = fittext.replace('A', strA)
        fittext = fittext.replace('B', strB)
       #fittext = fittext.replace('C', strC)
        
        Tfit = ROOT.TLatex(0.4, 0.55, fittext)
        Tfit.SetNDC(ROOT.kTRUE) 
        Tfit.SetTextSize(0.04) 
        Tfit.SetTextFont(132) 
        Tfit.SetTextColor(ROOT.kRed+1) 
        
        f.SetLineColor(ROOT.kRed+1)
        f.SetLineWidth(3)
        f.Draw('same')
        Tfit.Draw('same') 
        
    else:
        mg.Draw("ALP")
        
#    mg.GetHistogram().GetYaxis().SetLabelFont(132)
#    mg.GetHistogram().GetYaxis().SetTitleFont(132)
#    mg.GetHistogram().GetYaxis().SetLabelOffset(0.015)
#    mg.GetHistogram().GetYaxis().CenterTitle()
#    mg.GetHistogram().GetYaxis().SetNdivisions(505)
#    mg.GetHistogram().GetXaxis().SetNdivisions(505)
#    mg.GetHistogram().GetYaxis().SetTitleOffset(1.4)
#    
#    mg.GetHistogram().GetXaxis().SetTitleFont(132)
#    mg.GetHistogram().GetXaxis().SetLabelFont(132)
#    mg.GetHistogram().GetXaxis().SetLabelOffset(0.02)
#    mg.GetHistogram().GetXaxis().SetTitleOffset(1.5)
#    mg.GetHistogram().GetXaxis().SetTitleSize(0.06)
#    mg.GetHistogram().GetYaxis().SetTitleSize(0.06)
#    mg.GetHistogram().GetXaxis().SetLabelSize(0.06)
#    mg.GetHistogram().GetYaxis().SetLabelSize(0.06)
#    mg.GetHistogram().GetXaxis().SetRangeUser(xmin, xmax)
#    mg.SetMinimum(ymin)
#    mg.SetMaximum(ymax)
        
    
    if bl:
        leg = canvas.BuildLegend(0.62,0.68,0.92,0.88)
        leg.SetTextFont(132) 
        leg.SetFillColor(0)
        leg.SetFillStyle(0)
        leg.SetLineColor(0)
        leg.Draw() 
    
    Tleft.Draw() 
    Tright.Draw() 
#    mg.GetXaxis().SetRangeUser(xmin, xmax)

    canvas.Print('{0}/{1}.png'.format(pdir, name), 'png')

#_______________________________________________________________________

fileName = sys.argv[1]

#ptvals = {10., 20., 30.,50., 75., 100., 150., 200., 300., 500., 750., 1000., 1500., 2000., 5000.};
ptbins = [10., 20., 30., 50., 75., 100., 150., 200., 300., 500., 750., 1000., 1500., 2000., 5000., 7500.]
pts = [20,50,100,200,500,1000,2000]

mus = dict()
sigs = dict()

algo = 'ak4'

r = 0.4
eta = (0.0,1.3)
pu = 0

hfile = ROOT.TFile(fileName)

# retrieve histograms from ROOT files
num = 0
for i, pt in enumerate(ptbins):
    print i
    if i+1 == len(ptbins):
        break
    num = i
    histname = 'reco/reso/reco_{0}_{1}'.format(int(ptbins[i]),int(ptbins[i+1]))
    
    pt = 0.5*(ptbins[i]+ptbins[i+1])
    print pt
    
    print histname    
    hist = hfile.Get(histname)
    
    if hist.GetEntries() == 0:
        print 'Stops at ',pt 
        break
    histo = hist.Clone()

    histo.SetName(histname)
    
    fname = '{0}'.format(histname)
    
    if histo.Integral() > 0:
        histo.Scale(1/histo.Integral())
        
    x0 = histo.GetXaxis().GetBinCenter(histo.GetMaximumBin())      
    d = histo.GetRMS()
        
    # now perform gaussian fit in [x_max_sigm, x_max_sigp]
    fit = ROOT.TF1(fname, 'gaus',0.0, 2.0)

    s = 1.0
    histo.Fit(fname, 'Q', '', x0 - s*d, x0 + s*d)
    printResoHisto(histo, fit, pt, eta, r, pu)

    key = pt

    mus[key]  = (fit.GetParameter(1), fit.GetParError(1))
    sigs[key] = (fit.GetParameter(2)*100, fit.GetParError(2)*100)

    print '---------------------------------------------'
    print 'fitting results for {0} < pt < {1}'.format(ptbins[i],ptbins[i+1])
    print '---------------------------------------------'
    print ''

    mus[key]  = (fit.GetParameter(1), fit.GetParError(1))
    sigs[key] = (fit.GetParameter(2)*100 / fit.GetParameter(1), fit.GetParError(2)*100 / fit.GetParameter(1))

## Check of gen particle pt
#genhistname = 'gen/reso/gen_{0}_{1}'.format(int(ptbins[i]),int(ptbins[i+1]))
#genhist = hfile.Get(genhistname)
#genhisto = genhist.Clone()
#genfname = '{0}'.format(genhistname)
#genx0 = genhisto.GetXaxis().GetBinCenter(genhisto.GetMaximumBin())
#gend = genhisto.GetRMS()
## now perform gaussian fit in [x_max_sigm, x_max_sigp]                                                                                                                                                                                
#genfit = ROOT.TF1(genfname, 'gaus',0.0, 2.0)
#genhisto.SetName(genhistname)
#
#genhisto.Fit(genfname, 'Q', '', genx0 - s*gend, genx0 + s*gend)
#printResoHisto(genhisto, genfit, pt, eta, r, pu)
#

    print 'mu  = ', '{0:.2f} +/- {1:.2f}'.format(mus[key][0], mus[key][1])
    print 'res = ', '{0:.2f} +/- {1:.2f}'.format(sigs[key][0], sigs[key][1])
#    print 'mean of genParticles = ',genhisto.GetMean()
#    print 'mu   of genParticles = ',genfit.GetParameter(1)

# sigs[25]= (20., 0.1)

# now do multigraphs

#lt = 'QCD (uds) jets, #sqrt{{s}} = 100 TeV, {} PU'.format(pu)
rt = 'FCC-hh simulation'
lt = []
lt.append('#sqrt{s} = 100 TeV')
lt.append('QCD (uds) jets')
lt.append("anti-k_{{T}}, R = {0}".format(r))
lt.append("0 < |#eta| < 1.3".format(r))

#_____________________  VS PT__________________________

mg_reso_pt = ROOT.TMultiGraph()
mg_reso_pt.SetTitle(";p_{T} [GeV]; resolution (%)")

mg_resp_pt = ROOT.TMultiGraph()
mg_resp_pt.SetTitle(";p_{T} [GeV]; response")


reso_pt = ROOT.TGraphErrors()
reso_pt.SetLineColor(colors[0])
reso_pt.SetLineStyle(styles[0])
reso_pt.SetFillColor(0)
reso_pt.SetLineWidth(4)
reso_pt.SetMarkerColor(colors[0])

resp_pt = ROOT.TGraphErrors()
resp_pt.SetLineColor(colors[1])
resp_pt.SetLineStyle(styles[1])
resp_pt.SetFillColor(0)
resp_pt.SetLineWidth(4)
resp_pt.SetMarkerColor(colors[1])

for i, pt in enumerate(ptbins):
    
    if i+1 == len(ptbins):
        break
    if i >= num:
        break
    
    pt =  0.5*(ptbins[i]+ptbins[i+1])
    key = pt
    
    title = '{0} < |#eta| < {1} '.format(eta[0], eta[1])
    reso_pt.SetTitle(title)
    reso_pt.SetPoint(i,pt,sigs[key][0])
    reso_pt.SetPointError(i,0.,sigs[key][1])
#    pt = 0.5*(ptbins[i]+ptbins[i+1])
    
    resp_pt.SetTitle(title)
    resp_pt.SetPoint(i,pt,mus[key][0])
    resp_pt.SetPointError(i,0.,mus[key][1])
    
# fit pt resolution
freso = ROOT.TF1('reso', 'sqrt(pow([0]/sqrt(x),2) + pow([1],2))',0.0, 5000.0)
#freso.SetParLimits(0,0,100)
    
'''freso.SetParLimits(0,0,50)
freso.SetParLimits(1,60,200)
'''

reso_pt.Fit(freso, '', '', 10.0, 5000.0)

mg_reso_pt.Add(reso_pt)
mg_resp_pt.Add(resp_pt)

basename = os.path.basename(fileName)
basename = os.path.splitext(basename)[0]
print basename
name_reso = '{0}_reso_pt_{1}_{2}pu'.format(basename,algo,pu)
name_resp = '{0}_resp_pt_{1}_{2}pu'.format(basename,algo,pu)

#raw_input ("Press Enter to continue.. ")
drawMultiGraph(mg_reso_pt, name_reso, lt, rt, 'plots', 10, 5000., 0., 50., True, False, False, freso)   
freso = None
drawMultiGraph(mg_resp_pt, name_resp, lt, rt, 'plots', 10., 5000., 0.5, 1.2, True, False, True, freso)   

