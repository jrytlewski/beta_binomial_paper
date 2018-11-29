# -*- coding: utf-8 -*-
"""
@author: eyusko, updated by thether
@edit: rridder 09/25/17
@edit: sjohnson
"""

from differentialAbundanceBatch_2017_09_22 import main, readConfig, error_log, trainCoefs, find
import sys, os, optparse, traceback
import multiprocessing as mp


 #####      #####         ####         ####  ####      ####              ##########     #############     ##############         ####         ##################  ##############  ##############  ####  ####            ##############     ##########
 #######  #######        ######        ####  #####     ####            ##############   ###############   ##############        ######        ##################  ##############  ##############  ####  ####            ##############   ####      ####
 ################       ########       ####  ######    ####            ####       ####  ####        ####  ####                 ########              ####         ####            ####            ####  ####            ####            ####        ####
 #### ###### ####      ####  ####      ####  #######   ####            ####             ####       ####   ####                ####  ####             ####         ####            ####            ####  ####            ####            ######
 ####  ####  ####     ####    ####     ####  ########  ####            ####             #############     #########          ####    ####            ####         #########       #########       ####  ####            #########          #######
 ####        ####    ##############    ####  #### #### ####            ####             ###########       #########         ##############           ####         #########       #########       ####  ####            #########            #########
 ####        ####   ################   ####  ####  ########            ####             ####     ####     ####             ################          ####         ####            ####            ####  ####            ####                      ######
 ####        ####  ####          ####  ####  ####   #######            ####       ####  ####      ####    ####            ####          ####         ####         ####            ####            ####  ####            ####            ####        ####
 ####        ####  ####          ####  ####  ####    ######            ##############   ####       ####   ##############  ####          ####         ####         ##############  ####            ####  ##############  ##############   ####      ####
 ####        ####  ####          ####  ####  ####     #####              ##########     ####        ####  ##############  ####          ####         ####         ##############  ####            ####  ##############  ##############     ##########

def createfiles(outputDir, configFile, sampleloc, outBase, batchfile, replicateSamplesFile = None, parallel = False):
    if not os.path.exists(outputDir):
        os.makedirs(outputDir)
    outputDir +='/'

    #read config file
    items = readConfig(outputDir, configFile) #items = {} dictionary containing all config details. errors printed to error.log

    if items['method'] == 'betabinomial':
        coefs = trainCoefs(replicateSamplesFile, outBase, items)
    else:
        coefs = None

    stuff=[]
    batchData={}
    if os.path.isfile(batchfile):
        headers=['subject_id', 'sample_a', 'sample_b', "morisita_index", 'tcr_overlap',
                 'jaccard_index', 'bhattacharyya_distance', 'clones_shared', 'percent_shared', 'percent_shared_sample_a',
                 'percent_shared_sample_b', 'number_expanded', 'maximum_expansion(%)', 'number_contracted', 'maximum_contraction(%)',
                 'maximum_fold_change']

        items['batch']=True
        if items['reverse']=='yes': error_log(outputDir, 'Samples may not be reversed if a file is uploaded for batch analysis.')

        batchPairs = readAnalystBatchFile(outputDir, batchfile, sampleloc)
        if len(batchPairs)<=0:
            error_log(outputDir, 'No sample pairing were retrieved from the batch file.')

        if parallel:
            batchData = parallelBatch(outputDir, batchPairs, items, coefs)
        else:
            for pair in batchPairs:
                individualDic = one_pair(outputDir, pair, batchPairs[pair], items, coefs)
                if type(individualDic)==dict:
                    batchData[individualDic["pair_id"]]=individualDic
                else: #only returns ID
                    pair = individualDic
                    print ("Fail on ", pair, batchPairs[pair])
                    batchData[pair] = {"pair_id": pair, 'sample_a': batchPairs[pair][0][0], 'sample_b': batchPairs[pair][1][0]}
    else:
        print('batch file could not be opened')

    if len(batchData)>0:
        with open(outputDir+'/BatchResults.tsv', 'w') as batchwrite: #so that any data is written no matter whether an error pops up
            batchwrite.write('\t'.join(headers)+'\n')
            for ID in batchData:
                line_to_write = []
                for i in headers:
                    try:
                        line_to_write.append(str(batchData[ID][i])) # create a list of all data in the line to be written
                    except KeyError:
                        line_to_write.append('NA')
                batchwrite.write('\t'.join(line_to_write)) # add \t between all data and write
                batchwrite.write('\n')
    else:
        error_log(outputDir, 'No batch file was generated due to errors that prevented any comparison statistics from being calculated.')
        print('No pairings were successfully processed.')


 #############     #############         ########        ##########     ##############     ##########        ##########                  ##########            ####         #####      #####  #############     ####            ##############     ##########
 ###############   ###############     ############    ##############   ##############   ####      ####    ####      ####              ####      ####         ######        #######  #######  ###############   ####            ##############   ####      ####
 ####        ####  ####        ####   ####      ####   ####       ####  ####            ####        ####  ####        ####            ####        ####       ########       ################  ####        ####  ####            ####            ####        ####
 ####       ####   ####       ####   ####        ####  ####             ####            ######            ######                      ######                ####  ####      #### ###### ####  ####       ####   ####            ####            ######
 #############     #############     ####        ####  ####             #########          #######           #######                     #######           ####    ####     ####  ####  ####  #############     ####            #########          #######
 ###########       ###########       ####        ####  ####             #########            #########         #########                   #########      ##############    ####        ####  ###########       ####            #########            #########
 ####              ####     ####     ####        ####  ####             ####                      ######            ######                      ######   ################   ####        ####  ####              ####            ####                      ######
 ####              ####      ####     ####      ####   ####       ####  ####            ####        ####  ####        ####            ####        ####  ####          ####  ####        ####  ####              ####            ####            ####        ####
 ####              ####       ####     ############    ##############   ##############   ####      ####    ####      ####              ####      ####   ####          ####  ####        ####  ####              ##############  ##############   ####      ####
 ####              ####        ####      ########        ##########     ##############     ##########        ##########                  ##########     ####          ####  ####        ####  ####              ##############  ##############     ##########

def parallelBatch(outputDir, batchPairs, items, coefs):
    data = {}
    '''
    for pairID in batchPairs:
        data[pairID]=one_pair(outputDir, pairID, batchPairs[pairID], items)
    '''
    pool = mp.Pool(10) #test run uses about 1 GB per process... instance has 4 GB, so to be on safe side, we process only 2 at a time
    tasks = [pool.apply_async(one_pair, (outputDir, pairID, batchPairs[pairID], items, coefs,)) for pairID in batchPairs]
    for process in tasks:
        try:
            individualDic = process.get()
            if type(individualDic)==dict:
                data[individualDic["pair_id"]]=individualDic
            else: #only returns ID
                pair = individualDic
                print ("Fail on ", pair, batchPairs[pair])
                data[pair] = {"pair_id": pair, 'sample_a': batchPairs[pair][0][0], 'sample_b': batchPairs[pair][1][0]}

        except Exception as ex:
            print("Issue in parallel processing of process %s.\n\tException: %s\n\tError traceback should be shown. This should never be printed." %(str(process), ex)) #this should never happen because there is a try except around the entirety of one_pair function

    pool.close()
    pool.join()

    return data

# in this case, ID is the unique pair ID (subject ID_sampleA_VS_sampleB), while subject_id is the selected ID corresponding to the subject
def one_pair(outputDir, ID, samplesPair, items, coefs):
    output = outputDir+'/' + ID + '/'

    try: # all in a try for parallel processing
        if not os.path.exists(output): #make a subDirectory for all the differential abundance files generated for one match
            os.makedirs(output)
        subject_id = samplesPair[-1]
        del samplesPair[-1]
        individualDic=main(output, ID, samplesPair, items, coefs)
        if type(individualDic)==str: #minTotalError
            #os.rmdir(output)
            error_log(outputDir, individualDic)
            individualDic={} #gets ID and sample names
            # return ID
        '''
        rootdir = output
        with zipfile.ZipFile(output+'.zip', 'a') as zf:
            for subdir, dirs, files in os.walk(rootdir):
                for file in files:
                    zf.write(os.path.join(subdir,file), file) #the second argument is the name, so this places only the file in the zipped folder, not the entire path
                    os.remove(os.path.join(subdir,file))
            os.rmdir(rootdir)
        '''
        individualDic["pair_id"]=ID # unique ID corresponding to only this pair
        individualDic['sample_a']=samplesPair[0][0]
        individualDic['sample_b']=samplesPair[1][0]
        individualDic['subject_id']=subject_id # subject ID which can be repeated
        return individualDic
    except Exception as ex:
        error_log(outputDir, 'Issue in processing of %s pairing:\n\tException: %s' %(ID, ex))
        traceback.print_exception(*sys.exc_info()) #print to differential-abundance.log not error
        return ID


 #############            ####         ####  #############     ####  ####      ####    ##########               ##############  ####  ####            ##############
 ###############         ######        ####  ###############   ####  #####     ####  ##############             ##############  ####  ####            ##############
 ####        ####       ########       ####  ####        ####  ####  ######    ####  ####       ####            ####            ####  ####            ####
 ####       ####       ####  ####      ####  ####       ####   ####  #######   ####  ####                       ####            ####  ####            ####
 #############        ####    ####     ####  #############     ####  ########  ####  ####                       #########       ####  ####            #########
 ###########         ##############    ####  ###########       ####  #### #### ####  ####    #######            #########       ####  ####            #########
 ####               ################   ####  ####     ####     ####  ####  ########  ####    #######            ####            ####  ####            ####
 ####              ####          ####  ####  ####      ####    ####  ####   #######  ####       ####            ####            ####  ####            ####
 ####              ####          ####  ####  ####       ####   ####  ####    ######  ##############             ####            ####  ##############  ##############
 ####              ####          ####  ####  ####        ####  ####  ####     #####    ##########               ####            ####  ##############  ##############

def readAnalystBatchFile(outputDir, batchfile, sampleloc):
    batchPairs={}
    userinfo = [x.rstrip().split('\t') for x in open(batchfile, 'rU').readlines()]
    fields = userinfo[0]
    # Subject ID required in internal version, not Analyzer version
    sampi = fields.index('subject_id') # index of column containing subject_id
    s1i = fields.index('TSV 1') # index of column containing tsv 1
    s2i = fields.index('TSV 2') # index of column containing tsv 2
    for line in userinfo[1:]:
        if len(line)<=0: continue
        try:
            s1path=line[s1i].rstrip('\n\t\r') #sample 1 name
            s2path=line[s2i].rstrip('\n\t\r') #sample 2 name
            subject_id=line[sampi].rstrip('\n\t\r') #given subject ID
        except:
            if len(''.join(line))==0: continue
            error_log(outputDir, 'The line: '+'\t'.join(line) +' was not formatted correctly, and therefore omitted for analysis.')
            continue
        if ".tsv.gz" in s1path:
            s1name = s1path[:-7].split('/\\')[-1] #to get only the name without directories or extensions
            s2name = s2path[:-7].split('/\\')[-1]
        elif ".tsv" in s1path:
            s1name = s1path[:-4].split('/\\')[-1]
            s2name = s2path[:-4].split('/\\')[-1]
        else:
            sys.exit('include file extensions for pairings!')

        s1fullpath = sampleloc+s1path
        s2fullpath = sampleloc+s2path

        ID = subject_id + '_' + s1name + '_VS_' + s2name # make unique ID: combination of sample ID and pairings to allow for repeated subject IDs
        if os.path.isfile(s1fullpath) and os.path.isfile(s2fullpath):
            batchPairs[ID] = [[s1name, sampleloc+s1path], [s2name, sampleloc+s2path], subject_id]
        else:
            error_log(outputDir, 'paths of %s and %s were not valid. Their comparison is excluded.' %(sampleloc+s1path, sampleloc+s2path))
            print('paths of %s and %s were not valid. Their comparison is excluded.' %(sampleloc+s1path, sampleloc+s2path))
    return batchPairs


if __name__ == "__main__":
    # Read in the arguments
    parser = optparse.OptionParser()
    parser.add_option('--tsvDir', type='str', dest='tsvDir')
    parser.add_option('--outDir', type='str', dest='outDir')
    parser.add_option('--config', type='str', default="./configuration.ini", dest='configfile')
    parser.add_option('--batchfile', type='str', dest='batchfile')
    parser.add_option('--train', dest = 'replicateSamplesFile', default = '', help = 'if using betabinomial, filepath to the folder contianing training TSV')
    parser.add_option('--basename', dest = 'outBase', default = 'diff_ab_summary', help = 'output filename')
    parser.add_option('--parallel', dest='parallel', action = 'store_true', default=False, help = 'run in parallel') # stores False if not specified

    (options, args) = parser.parse_args()

    # We want to have a set of default behaviors such that the user only needs to
    # specify the sample (union) file and, if input and output directory tsv arguments
    # are blank, assume they are in the same directory.

    batchFile=options.batchfile
    if batchFile is None:
        print "Please provide a batch file. See --help for help."
        exit (2)

    if options.tsvDir is None:
        print "Using default sample file location."
        indx = find(batchFile, '/')
        sampleloc = batchFile[:indx[-1]+1]
        print "%r" % sampleloc
    # since you can only add '/' to a string, must ensure tsvDir is not None (is a string) first
    else:
        sampleloc=str(options.tsvDir)+'/'

    outputDir = options.outDir
    if outputDir is None:
        outputDir = sampleloc
        print "Results will be sent to %r" % outputDir

    configFile=options.configfile
    if configFile is None:
        print "Please provide a configuration file. See -help for help."
        exit (3)

    # To do: read in config file and make sure params are for bb are correct. More on this soon.


    createfiles(outputDir, configFile, sampleloc, options.outBase, batchFile, options.replicateSamplesFile, options.parallel)
