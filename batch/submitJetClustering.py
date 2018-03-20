#!/usr/bin/env python
import os, sys, subprocess
from optparse import OptionParser

def main():
    parser = OptionParser()

    parser.add_option ('-i','--input', help='input directory containing hgcal ntuple trees (most likely on eos)',
                       dest='input',
                       default='')

    parser.add_option ('-o','--output', help='output directory (default will have same name as input dir)',
                       dest='output',
                       default='jetsPt30_0PU')

    parser.add_option ('-q','--queue', help='lsf queue',
                       dest='queue',
                       default='1nd')

    parser.add_option ('-n', '--nev',  help='max number of event per file. default is 1000',
                       dest='nev',
                       default='1000')

    parser.add_option ('--njobs',  help='max number of jobs. Max is number of ntuple files',
                       dest='njobs',
                       default='10')

    parser.add_option("--algorithm", help="choose algorithm : antiKt, antiKt_cluster, or simpleCone",
                      dest="algorithm",
                      default='antiKt')

    parser.add_option("-c","--collect", help="collects jobs for given algorithm",
		      dest="collect", 
		      default='')

    parser.add_option("-e","--collectPt", help="collects jobs for given pt",
		      dest="collectPt",
		      default='')

    parser.add_option("-a","--allPts", help="collects jobs for all pts",
		      dest="allPts", action="store_true", 
		      default=False)

    parser.add_option("-t","--topoCluster", help="topoCluster for jet building",
                      dest="topoCluster", action="store_true",
                      default=False)

    parser.add_option("--coneCheck", help="check simple cone sum around jet axis",
                      dest="coneCheck", action="store_true",
                      default=False)

    (options, args) = parser.parse_args()
    input_dir     = options.input
    algo = options.algorithm
    check = "0"
    if options.coneCheck:
        algo = algo + '_coneCheck'
        check = "1"
        
    output_dir    = '/eos/experiment/fcc/users/c/cneubuse/JetClustering/'+algo+'/'+options.output
    max_events    = int(options.nev)
    queue         = options.queue
    collect       = options.collect

    if collect and options.collectPt:
        en = int(options.collectPt)
        print 'Collecting jobs for process: '+str(collect)
        input_dir = '/eos/experiment/fcc/users/c/cneubuse/JetClustering/'
        # /eos/experiment/fcc/users/c/cneubuse/JetClustering/antiKt/200GeVljets/out/
        print 'Find files in ', input_dir+collect+'/'+str(en)+'GeVljets'
        hadd_dir = input_dir+collect+'/'+str(en)+'GeVljets/out/'
        basename = input_dir+collect+'/'+str(en)+'GeVljets'
        basename = os.path.basename(basename)
        print basename
        outfile = hadd_dir + basename + '.root'
        hadd_files = hadd_dir + basename + '_*.root'
        cmd ='hadd -f '+ outfile + ' ' + hadd_files
        os.system(cmd)
        sys.exit('Collection of jobs done.')
       
    if collect and options.allPts:
        input_dir = '/eos/experiment/fcc/users/c/cneubuse/JetClustering/'
        energies = [20,50,100,200,500,1000,2000,5000,10000]
        outfile = input_dir+collect+'/all_'+collect+'.root'
        print outfile
        for energy in energies:
            dir_out = str(input_dir)+collect+'/'+str(energy)+'GeVljets/out/'
            print 'directory in which to find the merged files of energy {0}GeV : {1}'.format(energy,dir_out)
            os.system('cp '+dir_out+str(energy)+'GeVljets.root '+input_dir+collect+'/')
        hadd_files = input_dir+collect+'/*ljets.root'
        os.system('hadd -f '+ outfile + ' ' + hadd_files)
        sys.exit('All energies have been collected, saved as: '+outfile)

    # first create output dir
    if not os.path.exists(output_dir):
       os.makedirs(output_dir)
       os.makedirs(output_dir+'/std/')
       os.makedirs(output_dir+'/out/')
    else:
       sys.exit('Output dir: "'+output_dir+'" exists.')

    # find list of input files
    if options.topoCluster:
        cmd = "find {0} -type f -name '*_ntuple.root'".format(options.input)
    else :
        cmd = "find {0} -type f -name '*.root'".format(options.input)
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
    lst = process.communicate()[0]
    list_of_files = lst.splitlines()

    # just send one job per ntup file
    nfiles = len(list_of_files)
    njobs = min(nfiles, int(options.njobs))

    for job in xrange(njobs):
       
       print "output dir  : ", options.output
       print "output base : ", os.path.basename(options.output)
       basename = options.output + '_'+str(job)
       basename = os.path.basename(basename)
       print "output name : ",basename
       currentDir = os.getcwd()
       inputFile = list_of_files[job]
       
       if 'cms' in inputFile:
          runtype = 'CMS'
       else:
          runtype = 'FCC'
       
       outputFile = output_dir+'/out/'+basename+'.root'
       cmd = 'bsub -o '+output_dir+'/std/'+basename +'.out -e '+output_dir+'/std/'+basename +'.err -q '+queue
       cmd +=' -J '+basename+' "'+currentDir+'/submitJets.sh '+inputFile+' '+outputFile+' '+str(max_events)+' '+runtype+' '+str(check)+'"'
       
       print cmd
       # submitting jobs
       os.system(cmd)


#_______________________________________________________________________________________
if __name__ == "__main__":
    main()
