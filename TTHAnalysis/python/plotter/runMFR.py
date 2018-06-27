import optparse, subprocess, ROOT, datetime, math, array, copy, os
import numpy as np


## =====================================================
## =====================================================
## ==== FIRST SOME UTILITY FUNCTIONS ===================
## =====================================================
## =====================================================

doPUreweighting = True
doPUandSF = False

def submitFRrecursive(ODIR, name, cmd, dryRun=False):
    outdir=ODIR+"/jobs/"
    if not os.path.isdir(outdir): 
        os.system('mkdir -p '+outdir)
    os.system('cp ${{HOME}}/index.php {od}/../'.format(od=outdir))
    os.system('cp ${{HOME}}/index.php {od}/../../'.format(od=outdir))
    os.system('cp ${{HOME}}/resubFRs.py {od}/../../'.format(od=outdir))
    srcfile = outdir+name+".sh"
    logfile = outdir+name+".log"
    srcfile_op = open(srcfile,"w")
    srcfile_op.write("#! /bin/sh\n")
    srcfile_op.write("ulimit -c 0\n")
    srcfile_op.write("cd {cmssw};\neval $(scramv1 runtime -sh);\ncd {d};\n".format( 
            d = os.getcwd(), cmssw = os.environ['CMSSW_BASE']))
    srcfile_op.write(cmd+'\n')
    os.system("chmod a+x "+srcfile)
    bsubcmd = "bsub -q 1nd -o {logfile} {srcfile}\n".format(d=os.getcwd(), logfile=logfile, srcfile=srcfile)
    if dryRun: 
        print "[DRY-RUN]: ", bsubcmd
    else: os.system(bsubcmd)

def printAggressive(s):
    print '='.join('' for i in range(len(s)+1))
    print s
    print '='.join('' for i in range(len(s)+1))

def readScaleFactor(path, process, reterr = False):
    infile = open(path,'r')
    lines = infile.readlines()
    
    for line in lines:
        if 'Process {proc} scaled by'.format(proc=process) in line:
            scale = float(line.split()[4])
            scaleerr = float(line.split()[-1])
    if  not reterr:
        return scale
    else:
        return scale, scaleerr

def readFakerate(path, process):
    infile = open(path,'r')
    lines = infile.readlines()
    index = 999
    for ind,line in enumerate(lines):
        if process in line and '===' in line:
            index = ind
    frs = []; errs = []
    for il, line in enumerate(lines):
        if il < index+3: continue
        if len(line)==1: break
        frs.append(float(line.split()[2]))
        down = float(line.split()[3].replace('--','-'))
        up   = float(line.split()[4].replace('++','+'))
        errs.append( (abs(down)+abs(up))/2.)
    print('this is frs:', frs)
    return frs, errs

def runCards(trees, friends, targetdir, fmca, fcut, fsyst, plotbin, enabledcuts, disabledcuts, processes, scaleprocesses, extraopts = ''):

    if not type(trees)==list: trees = [trees]
    treestring = ' '.join(' -P '+ t for t in list(trees))

    cmd  = ' makeShapeCardsSusy.py --s2v -f -j 6 -l {lumi} --od {td} {trees} {fmca} {fcut}'.format(lumi=lumi, td=targetdir, trees=treestring, fmca=fmca, fcut=fcut)
    cmd += ' {plotvar} {binning}'.format(plotvar=plotbin.split()[0], binning=plotbin.split()[1])
    if friends:
        if not type(friends)==list: friends = [friends]
        for f in friends:
            cmd += ' -F Friends {friends}/tree_Friend_{{cname}}.root'.format(friends=f)
    cmd += ' -W new_puwts2016(nTrueInt) ' 
    cmd += ' -p '+','.join(processes)
    cmd += ''.join(' -E ^'+cut for cut in enabledcuts )
    cmd += ''.join(' -X ^'+cut for cut in disabledcuts)
    if scaleprocesses:
        for proc,scale in scaleprocesses.items():
            cmd += ' --scale-process {proc} {scale} '.format(proc=proc, scale=scale)
    cmd += ' {fsyst} '.format(fsyst=fsyst)
    if extraopts:
        cmd += ' '+extraopts

#example command
#python makeShapeCardsSusy.py --s2v -P /afs/cern.ch/work/e/efascion/DPStrees/TREES_110816_2muss/ --Fs /afs/cern.ch/work/e/efascion/public/friendsForDPS_110816/ -l 12.9 dps-ww/final_mca.txt dps-ww/cutfinal.txt finalMVA_DPS 10,0.,1.0  --od dps-ww/cards -p DPSWW,WZ,ZZ,WWW,WpWpJJ,Wjets  -W 0.8874 --asimov dps-ww/syst.txt
    print '============================================================================================='
    print 'running: python', cmd
    print '============================================================================================='
    subprocess.call(['python']+cmd.split())#+['/dev/null'],stderr=subprocess.PIPE)
    

def runefficiencies(trees, friends, targetdir, fmca, fcut, ftight, fxvar, enabledcuts, disabledcuts, scaleprocesses, compareprocesses, showratio, extraopts = ''):
    
    if not type(trees)==list: trees = [trees]
    treestring = ' '.join(' -P '+ t for t in list(trees))
    cmd  = ' mcEfficiencies.py --s2v -f -j 6 -l {lumi} -o {td} {trees} {fmca} {fcut} {ftight} {fxvar}'.format(lumi=lumi, td=targetdir, trees=treestring, fmca=fmca, fcut=fcut, ftight=ftight, fxvar=fxvar)
    if friends:
        cmd += ' -F Friends {friends}/tree_Friend_{{cname}}.root'.format(friends=friends)
    cmd += ' --groupBy cut '
    if doPUreweighting: cmd += ' -W new_puwts2016(nTrueInt)'

    if doPUandSF and not '-W ' in extraopts: cmd += ' -W puWeight*LepGood_effSF[0] '

    cmd += ''.join(' -E ^'+cut for cut in enabledcuts )
    cmd += ''.join(' -X ^'+cut for cut in disabledcuts)
    cmd += ' --compare {procs}'.format(procs=(','.join(compareprocesses)  ))
    if scaleprocesses:
        for proc,scale in scaleprocesses.items():
            cmd += ' --scale-process {proc} {scale} '.format(proc=proc, scale=scale)
    showrat   = ''
    if showratio:
        showrat = ' --showRatio '
    cmd += showrat
    if extraopts:
        cmd += ' '+extraopts

    print 'running: python', cmd
    subprocess.call(['python']+cmd.split())#+['/dev/null'],stderr=subprocess.PIPE)


def runplots(trees, friends, targetdir, fmca, fcut, fplots, enabledcuts, disabledcuts, processes, scaleprocesses, fitdataprocess, plotlist, showratio, extraopts = '', invertedcuts = []):
    
    if not type(trees)==list: trees = [trees]
    treestring = ' '.join(' -P '+ t for t in list(trees))
    cmd  = ' mcPlots.py --s2v -f -j 6 -l {lumi} --pdir {td} {trees} {fmca} {fcut} {fplots}'.format(lumi=lumi, td=targetdir, trees=treestring, fmca=fmca, fcut=fcut, fplots=fplots)
    if friends:
        if not type(friends)==list: friends = [friends]
        for f in friends:
            cmd += ' -F Friends {friends}/tree_Friend_{{cname}}.root'.format(friends=f)
    cmd += ''.join(' -E ^'+cut for cut in enabledcuts )
    cmd += ''.join(' -X ^'+cut for cut in disabledcuts)
    cmd += ' --sP '+','.join(plot for plot in plotlist)
    cmd += ' -p '+','.join(processes)
    if invertedcuts:
        cmd += ''.join(' -I ^'+cut for cut in invertedcuts )
    if doPUandSF and not '-W ' in extraopts: cmd += ' -W puWeight*LepGood_effSF[0] '
    if doPUreweighting: cmd += ' -W new_puwts2016(nTrueInt)'
    cmd += ' -o '+targetdir+'/'+'_AND_'.join(plot for plot in plotlist)+'.root'
    if fitdataprocess:
        cmd+= ' --fitData '
        cmd+= ''.join(' --flp '+proc for proc in fitdataprocess)
    if scaleprocesses:
        for proc,scale in scaleprocesses.items():
            cmd += ' --scale-process {proc} {scale} '.format(proc=proc, scale=scale)
    showrat   = ''
    if showratio:
        showrat = ' --showRatio '
    cmd += showrat
    if extraopts:
        cmd += ' '+extraopts

    print 'running: python', cmd
    subprocess.call(['python']+cmd.split())#+['/dev/null'],stderr=subprocess.PIPE)


## =====================================================
## =====================================================
## ==== USEFUL FUNCTIONS FROM HERE ON OUT ==============
## =====================================================
## =====================================================


def simplePlot():
    print '=========================================='
    print 'running simple plots'
    print '=========================================='
    trees     = ['/eos/user/m/mdunser/dps-13TeV-combination/TREES_latest/']
    friends   = '/eos/user/m/mdunser/dps-13TeV-combination/TREES_latest/friends_jet_pu_lepSF/'
    targetdir = '/eos/user/e/eabram/...'

    fmca      = 'dpsww13TeV/dps2016/simple/mca_simple.txt'
    fcut      = 'dpsww13TeV/dps2016/simple/cuts_simple.txt'
    fplots    = 'dpsww13TeV/dps2016/simple/plots.txt'

    enable    = []
    disable   = ['trigger2mu']
    processes = ['WWherw', 'WZ', 'WWCUETPM8']
    fittodata = []
    scalethem = {}
    extraopts = '--plotmode=norm'
    makeplots = ['pt1']
    showratio = True
    runplots(trees, friends, targetdir, fmca, fcut, fplots, enable, disable, processes, scalethem, fittodata, makeplots, showratio, extraopts)
    


if __name__ == '__main__':
    parser = optparse.OptionParser(usage='usage: %prog [opts] ', version='%prog 1.0')
    parser.add_option('--pf'        , '--postfix'    , dest='postfix'      , type='string'       , default=''    , help='postfix for running each module')
    parser.add_option('-d'          , '--date'       , dest='date'         , type='string'       , default=''    , help='run with specified date instead of today')
    parser.add_option('-l'          , '--lumi'       , dest='lumi'         , type='float'        , default=0.    , help='change lumi by hand')
    parser.add_option('--simple'    ,                  dest='simple'       , action='store_true' , default=False , help='make simple plot')
    (opts, args) = parser.parse_args()

    global date, postfix, lumi, date
    postfix = opts.postfix
    lumi = 36.0 if not opts.lumi else opts.lumi
    date = datetime.date.today().isoformat()
    if opts.date:
        date = opts.date

    if opts.simple:
        print 'making simple plots'
        simplePlot()
