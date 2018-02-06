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

    parser.add_option("-c","--collect", help="collects jobs for given output name",
		      dest="collect", action="store_true", 
		      default=False)

    (options, args) = parser.parse_args()
    input_dir     = options.input
    output_dir    = options.output
    max_events    = int(options.nev)
    queue         = options.queue
    collect       = options.collect

    if collect:
       print 'Collecting jobs for process: '+output_dir
       hadd_dir = output_dir +'/out/'
       outfile = hadd_dir + output_dir + '.root'
       hadd_files = hadd_dir + output_dir + '_*.root'
       cmd ='hadd -f '+ outfile + ' ' + hadd_files
       os.system(cmd)
       cmd = 'rm -rf '+hadd_files+' '+output_dir+'/std'
       os.system(cmd)
       sys.exit('Collection of jobs done.')

    # first create output dir
    if not os.path.exists(output_dir):
       os.makedirs(output_dir)
       os.makedirs(output_dir+'/std/')
       os.makedirs(output_dir+'/out/')
    else:
       sys.exit('Output dir: "'+output_dir+'" exists.')

    # find list of input files
    cmd = "find {} -type f -name '*.root'".format(options.input)
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
    lst = process.communicate()[0]
    list_of_files = lst.splitlines()

    # just send one job per ntup file
    nfiles = len(list_of_files)
    njobs = min(nfiles, int(options.njobs))

    for job in xrange(njobs):
       
       basename = output_dir + '_'+str(job)
       basename = os.path.basename(basename)
       currentDir = os.getcwd()
       inputFile = list_of_files[job]
       
       if 'cms' in inputFile:
          runtype = 'CMS'
       else:
          runtype = 'FCC'
       
       outputFile = currentDir+'/'+output_dir+'/out/'+basename+'.root'
       cmd = 'bsub -o '+output_dir+'/std/'+basename +'.out -e '+output_dir+'/std/'+basename +'.err -q '+queue
       cmd +=' -J '+basename+' "submitJets.sh '+currentDir+' '+inputFile+' '+outputFile+' '+str(max_events)+' '+runtype+'"'
       
       #print cmd
       # submitting jobs
       os.system(cmd)


#_______________________________________________________________________________________
if __name__ == "__main__":
    main()
