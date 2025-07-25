#!/usr/bin/python3

import os, sys, math, json, itertools
import numpy as np
from array import array
from collections import OrderedDict
from optparse import OptionParser
import ROOT

def main(argv = None):
    
    if argv == None:
        argv = sys.argv[1:]
        
    usage = "usage: %prog [options]\n Convert pyhf json to combine datacard"
    
    parser = OptionParser(usage)
    parser.add_option("--input", default='datacard.json', help="Input pyhf json file [default: %default]")
    parser.add_option("--output", default='datacard', help="Output combine datacard file name [default: %default]")
    parser.add_option("--normshape", action='store_true', help="Enable histosys split into normalization and shape components")
    
    (options, args) = parser.parse_args(sys.argv[1:])
    
    return options

if __name__ == '__main__':
    
    options = main()
    
    d = json.load(open(options.input))
    
    fr = ROOT.TFile(options.output+'.root', "RECREATE")
    
    # BB-lite check
    bbl = False
    for ch in d['channels']:
        for s in ch['samples']:
            for m in s['modifiers']:
                if 'prop' in m['name'] or m['type'] == 'staterror' or m['type'] == 'shapesys':
                    bbl = True
                    break
            if bbl: break
        if bbl: break
        
    # Get nuisances
    mods, modsstat, modsnormshape, nuis, samples = {}, [], {}, [], []
    for ch in d['channels']:
        for s in ch['samples']:
            found = False
            for samp in samples:
                if s['name'] == samp:
                    found = True
                    break
            if not found: samples.append(s['name'])
            for m in s['modifiers']:
                if m['type'] not in ['normfactor', 'staterror', 'lumi', 'shapesys'] and 'prop' not in m['name']:
                    if m['name'] not in modsnormshape.keys(): modsnormshape[m['name']] = [m['type']]
                    elif m['type'] not in modsnormshape[m['name']]: modsnormshape[m['name']].append(m['type'])
                    if m['name'] not in mods.keys():
                        mods[m['name']] = [m]
                        nuis.append(m['name'])
                    else:
                        foundsys = False
                        for me in mods[m['name']]:
                            if me['type'] == m['type']:
                                foundsys = True
                                break
                        if not foundsys:
                            mods[m['name']].append(m)
                elif m['name'] not in modsstat and m['type'] in ['staterror', 'shapesys']:
                    modsstat.append(m['name'])
    nuis = [s.replace('_mergedns', '') for i, s in enumerate(nuis)]
    nuis = list(set(nuis))
                    
    # Create shape file
    poi = []
    for ch in d['channels']:
        h = {}
        fr.cd()
        sd = fr.mkdir(ch['name']);
        sd.cd()
        for s in ch['samples']:
            hname = s['name']
            hnamep = ch['name']+'_'+hname
            data = s['data']
            nb = len(data)
            bins = array('f', list(np.arange(nb+1)))
            h[hnamep] = ROOT.TH1D(hname, hname, nb, bins)
            for m in s['modifiers']:
                if m['type'] in ['histosys'] or (m['type'] in ['normsys'] and len(modsnormshape[m['name']]) == 2):
                    hsys = hnamep+'_'+m['name']
                    hsysname = hname+'_'+m['name']
                    hsysname = hsysname.replace('_mergedns', '')

                    hsname = hsys.replace('_mergedns', '')
                    if hsname+'Up' not in h:
                        h[hsname+'Up'] = ROOT.TH1D(hsysname+'Up', hsysname+'Up', nb, bins)
                        h[hsname+'Down'] = ROOT.TH1D(hsysname+'Down', hsysname+'Down', nb, bins)                    
            for i in range(len(data)):
                h[hnamep].SetBinContent(i+1, data[i])
                for m in s['modifiers']:
                    if m['type'] == 'normfactor':
                        poidesc = [s['name'], m['name'], ch['name']]
                        if poidesc not in poi: poi.append(poidesc)
                    elif 'prop' in m['name'] or 'staterror' in m['type']:
                        h[hnamep].SetBinError(i+1, m['data'][i])
                    elif m['type'] in ['histosys']:
                        hsys = (hnamep+'_'+m['name']).replace('_mergedns', '')
                        vup = m['data']['hi_data'][i]
                        vdown = m['data']['lo_data'][i]
                        if 'mergedns' in m['name']:
                            for mn in s['modifiers']:
                                if mn['type'] == 'normsys' and mn['name'] == m['name']:
                                    vup *= mn['data']['hi']
                                    vdown *= mn['data']['lo']
                                    break
                        h[hsys+'Up'].SetBinContent(i+1, vup)
                        h[hsys+'Down'].SetBinContent(i+1, vdown)
                    elif m['type'] in ['normsys']:
                        hsys = (hnamep+'_'+m['name']).replace('_mergedns', '')
                        found = False
                        for mn in s['modifiers']:
                            if mn['type'] == 'histosys' and mn['name'] == m['name']:
                                found = True
                                break
                        if (not found) and (len(modsnormshape[m['name']]) == 2):
                            vup = data[i]*m['data']['hi']
                            vdown = data[i]*m['data']['lo']
                            if vdown < 1E-7: vdown = 1E-7
                            h[hsys+'Up'].SetBinContent(i+1, vup)
                            h[hsys+'Down'].SetBinContent(i+1, vdown)
                        
        chsamp = [s['name'] for s in ch['samples']]
        for s in samples:
            if s not in chsamp:
                hname = s
                hnamep = ch['name']+'_'+s
                h[hnamep] = ROOT.TH1D(hname, hname, nb, bins)

        for hk in h.keys():
            h[hk].SetDirectory(sd)
            h[hk].Write()

        for obs in d['observations']:
            if obs['name'] == ch['name']:
                hname = 'data_obs'
                hnamep = ch['name']+'_'+hname
                data = obs['data']
                nb = len(data)
                h[hnamep] = ROOT.TH1D(hnamep, hnamep, nb, array('f', list(np.arange(nb+1))))
                for i in range(len(data)):
                    h[hnamep].SetBinContent(i+1, data[i])
                    h[hnamep].SetBinError(i+1, math.sqrt(data[i]))
                h[hnamep].SetDirectory(sd)
                h[hnamep].Write(hname)
             
    fr.Close()
    
    # Create datacard
    
    chans = []
    for ich, ch in enumerate(d['channels']): chans.append(ch['name'])
    nchan = len(chans)
    poisig = ''
    normf = []
    for p in poi:
        if (('r_' in p[1]) and not 'atlas-' in options.input) \
        or (('XS' in p[1] or 'mu_tttt' in p[1] or 'mu_SIG' in p[1]) and 'atlas-' in options.input):
            poisig = p[0]
        else: normf.append(p)
    samples.remove(poisig)
    samples = [poisig]+samples
    nsamp = len(set(samples))
    dc = 'imax '+str(nchan)+' number of bins\n'
    dc += 'jmax '+str(nsamp-1)+' number of processes minus 1\n'
    dc += 'kmax '+str(len(nuis))+' number of nuisance parameters\n'
    dc += '------------------------------------\n'
    for ch in chans:
        dc += 'shapes * '+ch+' '+options.output.split('/')[-1]+'.root '+ch+'/$PROCESS '+ch+'/$PROCESS_$SYSTEMATIC\n'
    dc += '------------------------------------\n'
    dc += 'bin          '+' '.join(chans)+'\n'
    dc += 'observation  '+' '.join(np.repeat('-1 ', nchan))+'\n'
    dc += '------------------------------------\n'
    procbin, procname, procsamp, rate, samplenames, smap = [], [], [], [], [], {}
    it = 0
    for ich, ch in enumerate(d['channels']):
        smap[ch['name']] = {}
        for i, s in enumerate(samples):
            if ich == 0: samplenames.append(s)
            for sn in ch['samples']:
                if sn['name'] == s and sum(sn['data']) > 0:
                    smap[ch['name']][s] = it
                    it += 1
                    procbin.append(ch['name'])
                    procname.append(s)
                    procsamp.append(str(i))
                    rate.append(str(-1))
                    break
    dc += 'bin          '+' '.join(procbin)+'\n'
    dc += 'process      '+' '.join(procname)+'\n'
    dc += 'process      '+' '.join(procsamp)+'\n'
    dc += 'rate         '+' '.join(rate)+'\n'
    dc += '------------------------------------\n'

    sysd = []
    modsc = mods.copy()
    for m1 in mods.keys():
        if '_mergedns' in m1:
            sname = m1.replace('_mergedns', '')
            if sname in modsc.keys():
                modsc.pop(sname)
    
    for m in mods.keys():
        s = mods[m]
        for se in s:
#            if m in sysd: continue
            sysl = ''
            if se['type'] not in ['histosys', 'normsys']: continue
            sname = se['name'].replace('_mergedns', '')
            if sname in sysd and 'atlas-' not in options.input: continue
            for si in modsnormshape[se['name']]:
                if si == 'histosys':
                    sysl += sname+' shape '
                    for ch in d['channels']:
                        for samp in samples:
                            foundsamp = False
                            sampf, ndata = None, 0
                            for sampp in ch['samples']:
                                if samp == sampp['name']:
                                    foundsamp = True
                                    sampf = sampp
                                    ndata = sum(sampp['data'])
                                    break
                            if (not foundsamp) and (ndata > 0):
                                sysl += ' - '
                                continue
                            found = False
                            if sampf:
                                for sysn in sampp['modifiers']:
                                    if sysn['name'].replace('_mergedns', '') == sname:
                                        sysl += ' 1.0 '
                                        found = True
                                        break
                                if (not found) and (ndata > 0):
                                    sysl += ' - '
                elif si == 'normsys' and 'mergedns' not in se['name']:
                    sysl += se['name']+' lnN '
                    for ch in d['channels']:
                        for samp in samplenames:
                            found, ndata = False, 0
                            for sa in ch['samples']:
                                if sa['name'] == samp:
                                    ndata = sum(sa['data'])
                                    for sysn in sa['modifiers']:
                                        if sysn['name'] == se['name'] and sysn['type'] == si:
                                            if abs(sysn['data']['lo']-1./sysn['data']['hi']) < 1E-7:
                                                sysl += ' '+str(sysn['data']['hi'])
                                            else:
                                                if sysn['data']['lo'] < 1E-7: sysn['data']['lo'] = 1E-7
                                                sysl += ' '+str(sysn['data']['lo'])+'/'+str(sysn['data']['hi'])
                                            found = True
                                            break
                            if (not found) and (ndata > 0): sysl += ' - '
                if sysl != '':
                    sysl += '\n'
                    dc += sysl
                    sysd.append(sname)
                
    # add normsys to histosys
    dcmod = dc
    lines = dcmod.split('\n')
    histsplit = []
    histnorm = {}
    shapes = []
    for m in modsnormshape.keys():
        if len(modsnormshape[m]) < 2: continue
        for iline, line in enumerate(lines):
            dl = line.split()
            if len(dl) > 0 and m == dl[0] and 'shape' in dl[1]:
                if m in shapes:
                    lines[iline] = ''
                    continue
                for iline2 in range(len(lines)):
                    if iline == iline2: continue
                    d2 = lines[iline2].split()
                    if len(d2) > 0 and m == d2[0] and 'ln' in d2[1]:
                        shapes.append(m)
                        histnorm[m] = d2[2:]
                        if '_splitns' not in dl[0]:
                            histsplit.append(dl[0])
                            dl[0] += '_splitns'
                        lines[iline] = ' '.join(dl)
                        lines[iline2] = ''
                
    linescl = []
    lines = OrderedDict.fromkeys((line for line in lines if line))
    for l in lines:
        if l != '':
            linescl.append(l)
    dcmod = '\n'.join(linescl)
    
    # update shape file
    fr = ROOT.TFile(options.output+'.root', "UPDATE")
    chanlist = fr.GetListOfKeys()
    for ch in chanlist:
        chname = ch.ReadObj().GetName()
        for h in histsplit:
            dirr = fr.GetDirectory(chname)
            if dirr:
                for p in samples:
                    if dirr.GetListOfKeys().Contains(p+'_'+h+'Up'):
                        dirr.cd()
                        sf = histnorm[h][smap[chname][p]]                        
                        hup = dirr.Get(p+'_'+h+'Up')
                        hdown = dirr.Get(p+'_'+h+'Down')
                        hup.SetName(p+'_'+h+'_splitnsUp')
                        hdown.SetName(p+'_'+h+'_splitnsDown')
                        hup.SetDirectory(dirr)
                        hdown.SetDirectory(dirr)
                        if sf != '-' and len(modsnormshape[h]) == 2:
                            nln, nsh = 0, 0
                            for chh in d['channels']:
                                if chh['name'] == chname:
                                    for pp in chh['samples']:
                                        if pp['name'] == p:
                                            for mm in pp['modifiers']:
                                                if mm['name'] == h:
                                                    if mm['type'] == 'histosys': nsh += 1
                                                    elif mm['type'] == 'normsys': nln += 1
                            if nln+nsh == 2:
                                sfv = sf.split('/')
                                if len(sfv) == 2:
                                    hup.Scale(float(sfv[1]))
                                    hdown.Scale(float(sfv[0]))
                                else:
                                    hup.Scale(float(sfv[0]))
                                    hdown.Scale(1./float(sfv[0]))
                        hup.Write()
                        hdown.Write()
                        dirr.Delete(p+'_'+h+'Up;*')
                        dirr.Delete(p+'_'+h+'Down;*')
    fr.Close()
            
    normf.sort()
    normf = list(k for k, _ in itertools.groupby(normf))
    for p in normf:
        dcmod += '\n'+p[1]+' rateParam '+p[2]+' '+p[0]+' 1 [0.,10.]'
        
    hasStat = bool(len(modsstat) > 0)
    dcmod += '\n------------------------------------\n'
    for ch in chans:
        if hasStat:
            if bbl: dcmod += ch+' autoMCStats 0 1 2\n'
            else: dcmod += ch+' autoMCStats 0 100000 2\n'
        elif not options.normshape: dcmod += ch+' autoMCStats -1 1 2\n'
            
    with open(options.output+'.txt', 'w') as f:
        f.write(dcmod)
        f.close()