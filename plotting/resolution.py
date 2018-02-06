import ROOT, sys, os
from ROOT import TFile, TH1F, TGraphErrors, TMultiGraph, TLegend, TF1,gROOT

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
   ptext.AddText("{} pile-up".format(pu))
   ptext.AddText("anti-k_{{T}}, R = {}".format(r))
   ptext.AddText('{} < |#eta| < {} '.format(eta[0], eta[1]))
   ptext.AddText("p_{{T}} = {} GeV".format(pt))
   
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

   ROOT.gStyle.SetOptFit()
   ps = histo.GetListOfFunctions().FindObject("stats")
   ps.SetX1NDC(0.60)
   ps.SetX2NDC(0.95)
   ps.SetY1NDC(0.70)
   ps.SetY2NDC(0.90)
   c.Modified()
   c.Update()

   histo.GetXaxis().SetTitle('p_{T}^{reco} / p_{T}^{gen}')

   basename = os.path.basename(histo.GetName())
   if not os.path.exists('plots'):
       os.makedirs('plots')
       plotname = 'plots/'+basename+'.png'

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

       strA = '{:.0f}'.format(f.GetParameter(0))
       strB = '{:.0f}'.format(f.GetParameter(1))
       fittext = '#frac{#sigma(p_{T})}{p_{T}} = #frac{A%}{#sqrt{p_{T}}} #oplus B%'

       fittext = fittext.replace('A', strA)
       fittext = fittext.replace('B', strB)

       Tfit = ROOT.TLatex(0.6, 0.40, fittext)
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

    mg.GetYaxis().SetLabelFont(132)
    mg.GetYaxis().SetTitleFont(132)
    mg.GetYaxis().SetLabelOffset(0.015)
    mg.GetYaxis().CenterTitle()
    mg.GetYaxis().SetNdivisions(505)
    mg.GetXaxis().SetNdivisions(505)
    mg.GetYaxis().SetTitleOffset(1.4)

    mg.GetXaxis().SetTitleFont(132)
    mg.GetXaxis().SetLabelFont(132)
    mg.GetXaxis().SetLabelOffset(0.02)
    mg.GetXaxis().SetTitleOffset(1.5)
    mg.GetXaxis().SetTitleSize(0.06)
    mg.GetYaxis().SetTitleSize(0.06)
    mg.GetXaxis().SetLabelSize(0.06)
    mg.GetYaxis().SetLabelSize(0.06)
    mg.GetXaxis().SetRangeUser(xmin, xmax)
    mg.SetMinimum(ymin)
    mg.SetMaximum(ymax)

   
    if bl:
        leg = canvas.BuildLegend(0.62,0.68,0.92,0.88)
        leg.SetTextFont(132) 
        leg.SetFillColor(0)
        leg.SetFillStyle(0)
        leg.SetLineColor(0)
        leg.Draw() 
    
    Tleft.Draw() 
    Tright.Draw() 

    canvas.Print('{}/{}.png'.format(pdir, name), 'png')

#_______________________________________________________________________

fileName = sys.argv[1]

#ptbins = [10., 20., 30.,50., 75., 100., 150., 200., 300., 500., 750., 1000., 1500., 2000., 5000.]
ptbins = [10., 20., 30.,50., 75., 100., 150., 200., 300., 500., 750., 1000., 1500.]

mus = dict()
sigs = dict()

algo = 'ak4'

r = 0.4
eta = (0.0,1.3)
pu = 0

hfile = ROOT.TFile(fileName)

# retrieve histograms from ROOT files
i = 0
for pt in ptbins:
    
    if i == len(ptbins)-1:
       break
    
    histname = 'reco/reso/reco_{:.0f}_{:.0f}'.format(ptbins[i],ptbins[i+1])

    hist = hfile.Get(histname)
    histo = hist.Clone()
    histo.SetName(histname)

    fname = 'f{}'.format(histname)

    if histo.Integral() > 0:
        histo.Scale(1/histo.Integral())

    x0 = histo.GetXaxis().GetBinCenter(histo.GetMaximumBin())      
    d = histo.GetRMS()

    # now perform gaussian fit in [x_max_sigm, x_max_sigp]
    f = ROOT.TF1(fname, 'gaus',0.0, 2.0)

    s = 1.0
    histo.Fit(fname, 'Q', '', x0 - s*d, x0 + s*d)
    printResoHisto(histo, f, pt, eta, r, pu)

    key = pt

    mus[key]  = (f.GetParameter(1), f.GetParError(1))
    sigs[key] = (f.GetParameter(2)*100, f.GetParError(2)*100)

    print '---------------------------------------------'
    print 'fitting results for {} < pt < {}'.format(ptbins[i],ptbins[i+1])
    print '---------------------------------------------'
    print ''

    mus[key]  = (f.GetParameter(1), f.GetParError(1))
    sigs[key] = (f.GetParameter(2)*100 / f.GetParameter(1), f.GetParError(2)*100 / f.GetParameter(1))

    print 'mu  = ', '{:.2f} +/- {:.2f}'.format(mus[key][0], mus[key][1])
    print 'res = ', '{:.2f} +/- {:.2f}'.format(sigs[key][0], sigs[key][1])
    print ''
    print ''
    i += 1

# now do multigraphs

#lt = 'QCD (uds) jets, #sqrt{{s}} = 100 TeV, {} PU'.format(pu)
rt = 'FCC-hh simulation'
lt = []
lt.append('#sqrt{s} = 100 TeV')
lt.append('QCD (uds) jets')
lt.append("anti-k_{{T}}, R = {}".format(r))
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

i = 0
for pt in ptbins:
   key = pt
   if i == len(ptbins)-1:
      break

   title = '{} < |#eta| < {} '.format(eta[0], eta[1])
   reso_pt.SetTitle(title)
   reso_pt.SetPoint(i,pt,sigs[key][0])
   reso_pt.SetPointError(i,0.,sigs[key][1])

   resp_pt.SetTitle(title)
   resp_pt.SetPoint(i,pt,mus[key][0])
   resp_pt.SetPointError(i,0.,mus[key][1])

   i += 1

# fit pt resolution
#freso = ROOT.TF1('reso', 'sqrt([0]^2/x^2 + [1]^2/x + [2]^2)',0.0, 10000.0)
freso = ROOT.TF1('reso', 'sqrt([1]^2/x + [2]^2)',0.0, 10000.0)
#freso.SetParLimits(0,0,100)
freso.SetParLimits(0,0,200)
freso.SetParLimits(1,0,5)
reso_pt.Fit(freso, '', '', 20.0, 5000.0)
#freso.SetLineWidth(3)
#freso.SetLineColor(ROOT.kBlack)
#freso.Draw('same')
#f.Draw('same')

mg_reso_pt.Add(reso_pt)
mg_resp_pt.Add(resp_pt)

if algo == 'ak4': r = 0.4
if algo == 'ak2': r = 0.2
if algo == 'ak1': r = 0.1

basename = os.path.basename(fileName)
basename = os.path.splitext(basename)[0]

name_reso = '{}_reso_pt_{}_{}pu'.format(basename,algo,pu)
name_resp = '{}_resp_pt_{}_{}pu'.format(basename,algo,pu)

drawMultiGraph(mg_reso_pt, name_reso, lt, rt, 'plots', 10, 5000., 0., 40., True, False, False, freso)   
freso = None
drawMultiGraph(mg_resp_pt, name_resp, lt, rt, 'plots', 10., 5000., 0.5, 1.2, True, False, True, freso)   

