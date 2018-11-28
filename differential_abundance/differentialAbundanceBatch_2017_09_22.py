# -*- coding: utf-8 -*-
"""
@author: wdewitt 2015
@author: rridder 2016/07/20 - calculations, graphics, configs...
@edit: rridder 2016/08/05 - batch mode
@edit: eyusko 2016 - betabinomial
@edit: jrytlewski 2016/11/29 - source and minTotal calculations
@edit: rridder 2017/09/22 - analyzer standardization, outputs, batchsettings

Must be called from the runBatch script.

Dependencies:
Python 2.7
Matplotlib 1.5.1
Scipy 17.1
Numpy 1.11.1
"""

import sys, os, math
from scipy.stats import binom, binom_test, fisher_exact
import scipy, matplotlib
matplotlib.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
from copy import deepcopy
import gzip, glob
from scipy.special import betaln, gammaln, digamma, polygamma, logit, expit
from scipy.optimize import minimize


   ##########         ########      ####      ####  ##############  ####    ##########
 ##############     ############    #####     ####  ##############  ####  ##############
 ####       ####   ####      ####   ######    ####  ####            ####  ####       ####
 ####             ####        ####  #######   ####  ####            ####  ####
 ####             ####        ####  ########  ####  #########       ####  ####
 ####             ####        ####  #### #### ####  #########       ####  ####    #######
 ####             ####        ####  ####  ########  ####            ####  ####    #######
 ####       ####   ####      ####   ####   #######  ####            ####  ####       ####
 ##############     ############    ####    ######  ####            ####  ##############
   ##########         ########      ####     #####  ####            ####    ##########


def readConfig(outputDir, configFile): #one function to check all the variables coming from the config file
    '''
    ARGUMENTS: OutputDirectory for error reporting, Configuration File (with full location)
    RETURNS: dictionary of settings
    PURPOSE: One function to read the entire config file and make any changes to defaults.
             Ensures all settings are acceptable for the rest of the code.
    '''
    items={}
    if sys.version_info<(3, 0):
        import ConfigParser as ConfigParser #in python 2.7 it is capitalized ConfigParser
    elif sys.version_info>=(3, 0):
        import configparser as ConfigParser #in python 3 it is lowercase configparser
    config = ConfigParser.RawConfigParser()
    try: config.read(configFile) #make sure file is readable, and return the name as a list
    except:
        # all default options if config file is not readable
        items={'source':'repertoire', 'count': 'nucleotide', 'minTotal': 10, 'verbose': False, 'correction': 'bh',
               'scatter': ['frequency', 'symlog', 'extra', 'frequency', 'on', 'spwo', 4.5, True],
               'hist': ['linear', 'on'], 'alpha': 0.01, 'productiveOnly': True, 'alternative': 'two-sided',
               'method': 'binomial', 'reverse': 'no', 'format': 'pdf', 'batch_count_minTot_seqs':False}
        error_log(outputDir, '[union] was missing from the configuration file, so all settings were returned to default')
        return items
    try:
        count = ''.join(config.get('union', 'count').split()).lower() # get value of count from section union, remove space from 'amino acid' if it's the input
        if (count not in ['nucleotide', 'aminoacid']): # it will yield an error before exiting the function... unless there is some way we can do some sort of raw_input for users to edit their input
            error_log(outputDir, 'count must be "nucleotide" or "amino acid", not %s, used default, nucleotide.' %(count))
            count='nucleotide'
        if count == 'aminoacid': count ='aminoAcid' # correct capitalization for column matching
    except:
        count='nucleotide'
        error_log(outputDir, 'count could not be retrieved from the configuration settings, used default, nucleotide.')
    try:
        minTotal = config.get('union', 'minTotal')
        try:
            assert('-' != minTotal[0]) # fail if a negative value is given... this may not be written as minTotal>=0 because it is retrieved as a string
            minTotal=float(minTotal)
            if minTotal != int(minTotal): # ex: if 1.7 is input, minTotal must be rounded to 1
                error_log(outputDir, 'minTotal was %s (not a valid integer), used rounded value of %s.' %(minTotal, int(minTotal)))
            minTotal=int(minTotal)
        except: # if minTotal is negative
                error_log(outputDir, 'minTotal was %s (not a valid integer), used default, 10.' %(str(minTotal)))
                minTotal=10
    except:
        minTotal=10
        error_log(outputDir, 'minTotal could not be retrieved from the configuration settings, used default, 10.')
    try:
        productiveOnly=config.getboolean('union', 'productiveOnly')
        if count=='aminoAcid' and productiveOnly==False:
            error_log(outputDir, 'to count Amino Acids, productiveOnly must be True. productiveOnly was corrected to True.')
            productiveOnly=True
    except:
        try: error_log(outputDir, 'productiveOnly was %s (not a boolean value), used default, True.' %(config.get('union', 'productiveOnly')))
        except: error_log(outputDir, 'productiveOnly could not be retrieved from the configuration settings, used default, True.')
        productiveOnly=True
    try:
        source = config.get('abundance', 'source').lower() # get source from the abundance section
        if (source not in ['repertoire', 'nucleated']): # it will yeild an error before exiting the function... unless there is some way we can do some sort of raw_input for users to edit their input
            error_log(outputDir, 'source must be "repertoire" or "nucleated", not %s, used default, repertoire.' %(source))
            source='repertoire'
    except:
        source='repertoire'
        error_log(outputDir, 'source could not be retrieved from the configuration settings, used default, repertoire.')
    try:
        method = config.get('p value', 'method').lower() #get the comparison method
        if method not in ['binomial', 'fisher', 'betabinomial']:
            error_log(outputDir, 'method must be either "binomial", "fisher", or "betabinomial", not %s, used default, "binomial"' %(method))
            method='binomial'
    except:
        method='binomial'
        error_log(outputDir, 'method could not be retrieved from the configuration settings, used default, binomial.')
    if method == 'betabinomial' and source == 'nucleated':
        error_log(outputDir, "For method='betabinomial', source must be equal to 'repertoire'. Aborting!")
        sys.exit("For method='betabinomial', source must be equal to 'repertoire'. Aborting!")
    try:
        alternative = config.get('p value', 'alternative').lower() #get other parameters
        if alternative not in ['greater', 'less', 'two-sided']:
            error_log(outputDir, 'alternative must be "less", "greater", or "two-sided", not %s, used default, "two-sided"' %(alternative))
            alternative='two-sided'
    except:
        alternative='two-sided'
        error_log(outputDir, 'alternative could not be retrieved from the configuration settings, used default, two-sided.')
    try:
        correction = config.get('significance', 'correction').lower()
        if correction not in ['none', 'bh', 'by', 'bonferroni']:
            error_log(outputDir, 'correction must be "none", "bh" (Benjamini-Hochberg), "by" (Benjamini-Yukutieli), or "bonferroni" not %s, used default, "bh"' %(correction))
            correction = 'bh'
    except:
        correction='bh'
        error_log(outputDir, 'correction could not be retrieved from the configuration settings, used default, bh.')
    try:
        alpha = config.get('significance', 'alpha')
        try:
            alpha=float(alpha)
            assert(alpha<=1 and alpha>=0) #Note: alpha may be equal to 0 or 1
        except:
            error_log(outputDir, 'invalid entry of %s for alpha, used default, .01' %(alpha))
            alpha=.01
    except:
        alpha=.01
        error_log(outputDir, 'alpha could not be retrieved from the configuration settings, used default, .01.')
    try: verbose = config.getboolean('output', 'verbose')
    except:
        try:
            verbose=config.get('output', 'verbose').lower()
            if verbose not in ['extra', 'scatter']:
                error_log(outputDir, 'verbose was changed to False because %s is not a valid input.' %verbose)
                verbose=False
        except:
            error_log(outputDir, 'verbose could not be retrieved from the configuration settings, used default, False.')
            verbose=False
    try:
        reverse = config.get('output', 'reverse samples').lower()
        if reverse not in ['yes', 'no']:
            error_log(outputDir, 'reverse samples was changed to no because %s is not a valid input.' %reverse)
            reverse='no'
    except:
        error_log(outputDir, 'reverse samples could not be retrieved from the configuration settings, used default, no.')
        reverse='no'
    try:
        form = config.get('output', 'graphics format').lower()
        #print (plt.gcf().canvas.get_supported_filetypes())
        if form not in ['eps', 'jpg', 'jpeg', 'pdf', 'png', 'ps', 'raw', 'svg', 'svgz']: #use plt.gcf().canvas.get_supported_filetypes() for formats supported by backend.
            error_log(outputDir, 'graphics format was changed to pdf because %s is not a supported file format.' %form)
            form='pdf'
    except:
        error_log(outputDir, 'graphics format could not be retrieved from the configuration settings, used default, pdf.')
        form='pdf'
    try:
        batch_count_minTot_seqs=config.getboolean('batch', 'count minTotal excluded sequences')
    except:
        try: error_log(outputDir, 'count minTotal excluded sequences was %s (not a boolean value), used default, False.' %config.get('batch', 'count minTotal excluded sequences'))
        except: error_log(outputDir, 'count minTotal excluded sequences could not be retrieved from the configuration settings, used default, False.')
        batch_count_minTot_seqs=False
    try:
        hbase=config.get('histogram', 'base').lower()
        if hbase not in ['log', 'linear']:
            error_log(outputDir, 'histogram base was changed to linear because %s is not a valid input.' %hbase)
            hbase = 'linear'
    except:
        error_log(outputDir, 'histogram base could not be retrieved from the configuration settings, used default, linear.')
        hbase='linear'
    try:
        labels=config.get('histogram', 'labels').lower()
        if labels not in ['on', 'off']:
            error_log(outputDir, 'histogram labels was changed to on because %s is not a valid input.' %labels)
            labels = 'on'
    except:
        error_log(outputDir, 'histogram labels could not be retrieved from the configuration settings, used default, on.')
        labels='on'
    try:
        scale = config.get('scatter graph', 'scale').lower()
        if scale not in ['abundance', 'frequency']:
            error_log(outputDir, 'scatter scale was changed to frequency because %s is not a valid input.' %scale)
            scale= 'frequency'
    except:
        error_log(outputDir, 'scatter scale could not be retrieved from the configuration settings, used default, frequency.')
        scale='frequency'
    try:
        sbase = config.get('scatter graph', 'scatter base').lower()
        if sbase=='log': sbase='symlog'
        elif sbase != 'linear':
            error_log(outputDir, 'scatter base was changed to log because %s is not a valid input.' %sbase)
            sbase='symlog'
    except:
        error_log(outputDir, 'scatter base could not be retrieved from the configuration settings, used default, log.')
        sbase='symlog'
    try:
        sverbose = config.get('scatter graph', 'scatter verbose').lower()
        if sverbose not in ['true', 'false', 'extra']:
            error_log(outputDir, 'scatter verbose was changed to extra because %s is not a valid input.' %sverbose)
            sverbose='extra'
        if sverbose=='true': sverbose=True
        elif sverbose=='false': sverbose=False
    except:
            sverbose='extra'
            error_log(outputDir, 'scatter verbose could not be retrieved from the configuration settings, used default, extra.')
    try:
        line = config.get('scatter graph', 'line of equality comparison').lower()
        if line not in ['none', 'abundance', 'frequency']:
            error_log(outputDir, 'line of equality comparison was changed to frequency because %s is not a valid input.' %line)
            line='frequency'
    except:
        error_log(outputDir, 'line of equality comparison could not be retrieved from the configuration settings, used default, frequency.')
        line='frequency'
    pointopts=['transparentpoints', 'transparentpointswithoutline', 'solidpoints', 'solidpointswithoutline', 'whitepointswithtransparentoutline', 'whitepointswithwidthoutline']
    pointoptsshort=['tp', 'tpwo', 'sp', 'spwo', 'wpwto', 'wpwwo']
    try:
        points = ''.join(config.get('scatter graph', 'points').split()).lower()
        if points in pointopts:
            i=pointopts.index(points)
            points=pointoptsshort[i]
        if points not in pointoptsshort:
            error_log(outputDir, 'scatter points was changed to spwo because %s is not a valid input.' %points)
            points='spwo'
    except:
        error_log(outputDir, 'scatter points could not be retrieved from the configuration settings, used default, spwo.')
        points='spwo'
    try:
        minTotLine = config.get('scatter graph', 'threshold for statistical comparison line').lower()
        if minTotLine not in ['on', 'off']:
            error_log(outputDir, 'threshold for statistical comparison line was changed to on because %s is not a valid input.' %minTotLine)
            minTotLine='on'
    except:
        error_log(outputDir, 'threshold for statistical comparison line could not be retrieved from the configuration settings, used default, on.')
        minTotLine='on'
    try:
        size = config.get('scatter graph', 'size')
        size=float(size)
        if size>10 or size<1:
            error_log(outputDir, 'scatter size was changed to 4.5 because %s is not a valid input.' %size)
            size=4.5
    except:
        size=4.5
        error_log(outputDir, 'scatter size could not be retrieved as a valid number, used default, 4.5')
    try: lverbose = config.getboolean('scatter graph', 'legend verbose')
    except:
        error_log(outputDir, 'legend verbose could not be retrieved as a boolean value, used default, True')
        lverbose=True

    #asserts guarantee config variables for the rest of the code:
    assert(source in ['repertoire', 'nucleated'])
    assert(count in ['nucleotide', 'aminoAcid'])
    assert(method in ['binomial', 'fisher', 'betabinomial'])
    assert(alternative in ['greater', 'less', 'two-sided'])
    assert(correction in ['none', 'bh', 'by', 'bonferroni'])
    assert(type(productiveOnly)==bool)
    assert(type(alpha)==float)
    assert(alpha<=1 and alpha>=0) # alpha may equal 1
    assert(verbose in [True, False, 'extra', 'scatter'])
    assert(hbase in ['log', 'linear'])
    assert(scale in ['frequency', 'abundance'])
    assert(sbase in ['symlog', 'linear'])
    assert(sverbose in [True, False, 'extra'])
    assert(1<=size<=10)
    assert(points in pointoptsshort)
    assert(line in ['none', 'frequency', 'abundance'])
    assert(minTotLine in ['on', 'off'])
    assert(reverse in ['yes', 'no'])
    assert(type(lverbose)==bool)
    items['source']=source
    items['count']=count
    items['minTotal']= minTotal
    items['productiveOnly']=productiveOnly
    items['method']=method
    items['alternative']=alternative
    items['correction']=correction
    items['alpha']=alpha
    items['verbose']=verbose
    items['hist']=[hbase, labels]
    items['scatter']=[scale, sbase, sverbose, line, minTotLine, points, size, lverbose]
    items['reverse']=reverse
    items['format']=form
    items['batch_count_minTot_seqs']=batch_count_minTot_seqs
    return items


def returnConfig(outbase, items):
    '''
    ARGUMENTS: output file location, dictionary of settings
    RETURNS: NONE
    PURPOSE: Writes a configuration file to the output directory.
             It is corrected for any changes that were made to the settings throughout the code.
             Occurs at the very end of main in order to ensure any adjusted settings (especially minTotal) are reflected.
    '''
    f=open(outbase+'ConfigurationFile.txt', 'w')
    if items['scatter'][1]=='symlog': items['scatter'][1]='log'
    scale, base, verbose, line, minTotLine, points, size, lverbose=items['scatter']
    f.write('# %s\n' %outbase)
    f.write('# Reconstructed Configuration File with corrected inputs for further testing\n\n')
    f.write('[union]\ncount = %s\nminTotal = %s\nproductiveOnly = %s\n\n' %(items['count'], items['minTotal'], items['productiveOnly']))
    #f.write('[abundance]\nsource = %s\n\n' %items(['source']))
    f.write('[p value]\nmethod = %s\nalternative = %s\n\n' %(items['method'], items['alternative']))
    f.write('[significance]\ncorrection = %s\nalpha = %s\n\n' %(items['correction'], items['alpha']))
    f.write('[output]\nverbose = %s\nreverse samples = %s\ngraphics format = %s\n\n' %(items['verbose'], items['reverse'], items['format']))
    f.write('[histogram]\nbase = %s\nlabels = %s\n\n' %(items['hist'][0], items['hist'][1]))
    f.write('[scatter graph]\nscale = %s\nscatter base = %s\nscatter verbose = %s\npoints = %s\nsize = %s\nline of equality comparison = %s\nthreshold for statistical comparison line = %s\nlegend verbose = %s\n\n'
            %(scale, base, verbose, points, size, line, minTotLine, lverbose))
    f.write('# see error.log for a list of any corrections.')
    f.close()


 ##################     ##########     ####          ####            #############     ##############         ####         ###########
 ##################   ####      ####   ####          ####            ###############   ##############        ######        ##############
        ####         ####        ####  ####          ####            ####        ####  ####                 ########       ####      #####
        ####         ######             ####        ####             ####       ####   ####                ####  ####      ####        ####
        ####            #######          ####      ####              #############     #########          ####    ####     ####        ####
        ####              #########       ####    ####               ###########       #########         ##############    ####        ####
        ####                   ######      ####  ####                ####     ####     ####             ################   ####        ####
        ####         ####        ####       ########                 ####      ####    ####            ####          ####  ####      #####
        ####          ####      ####         ######                  ####       ####   ##############  ####          ####  ##############
        ####            ##########            ####                   ####        ####  ##############  ####          ####  ###########

def tsvUnion(outputDir, samples, counting, minTotal=10, productiveOnly=True, source='repertoire'):
    '''
    ARGUMENTS: file name prefix convention
               list of sample names and locations
               counting (nucleotide or aminoAcid) for what sequences are being tracked
               minTotal parameter for the minimum number of counts to exclude
               productiveOnly (True/False) to include only productive sequences
    RETURNS: totalsfull: the full total counts regardless of minTotal
             totalsminTotal: the totals after the minTotal filter
             seqsfull: a list of sequence dictionaries regardless of minTotal
             seqsminTotal: a list of sequence dictionaries after minTotal filtering
             minTotEx: a list of sequence dictionaries excluded by minTotal
             exclude: a dictionary of sequences that had a 0 or less in the inputTemplateEstimate column that were not included in either sequence list
    PURPOSE: Reads all the data, and returns it in useable formats
    '''
    seqs = {}
    nucleatedCells=dict((key[0], 0) for key in samples)
    for sample in samples: # samples is a list of data files to read from (must be len 2)
        linenum=0
        #for ease of selection for batch mode, the user selects all samples and
        #run-analysis.sh does not unzip all of them because of space restrictions
        #they must be read in a way based on their zip status
        if '.gz' in sample[1][-3:]: #sample[1] is full file name and location
            openFile = gzip.open(sample[1], 'rb') #special way to read zipped files
        else:
            openFile = open(sample[1], 'r') #if file has already been unzipped
        for line in openFile:
            linenum+=1   #allows better error reporting
            if '#' in line and source =='repertoire': continue # lines at beginning with unnecessary info for this tool
            elif '#' in line and source =='nucleated':
                try:
                    if 'estTotalNucleatedCells' in line:
                        nucleatedCells[sample[0]] = (int(round(float(line.strip().split('=')[1]))))
                    if 'productionPCRAmountofTemplate' in line:
                        inputDNA = line.strip().split('=')[1]
                        if inputDNA != '':
                            nucleatedCells[sample[0]]=(int(round(float(inputDNA) * 153.846)))
                except ValueError:
                    print('nucleatedCell count for %s was not retreived due to a value error:\n\tline:'%sample[0] + line)
                    error_log(outputDir, 'nucleatedCell count for %s was not retreived due to a value error:\n\tline:'%sample[0] + line)
            elif ('nucleotide' in line) or ('rearrangement' in line): #this is the header line
                headers = line.split('\t') # each header is found by splitting the line
                try: seqi = headers.index(counting) # index of value matching counting (nucleotide or amino acid)
                except:
                    if counting == 'nucleotide':
                        seqi = headers.index('rearrangement')
                    elif counting == 'aminoAcid':
                        seqi = headers.index('amino_acid')
                try: cti=headers.index('inputTemplateEstimate') # find the index of the template counts column
                except:
                    try: cti=headers.index('templates')
                    except ValueError:
                        error_log(outputDir, '%s is incompatible with template estimates. It was omitted from all comparisons. Check the versioning of your files.\n' %(sample[0]))
                        return
                try: statusi = headers.index('sequenceStatus') # index of 'sequenceStatus' (column only includes 'In', 'Out', or 'Stop')
                except: statusi = headers.index('frame_type')
            else: #if it is not the header line (it is a data line)      ### else only occurs after the header line, so variables defined above hold
                fields = line.split('\t') #split values in same order as headers
                if fields[statusi] != 'In' and productiveOnly: continue #if it is Out or Stop and productiveOnly==True, continues to next loop iteration (so it only includes 'In')
                seq = fields[seqi] # data in either nucleotide or amino acid column
                ct  = int(fields[cti]) # data in inputTemplateEstimate column
                if ct<=0: 
                    ct=-1
                    error_log(outputDir, "Error in %s sample file. Clone %s had a processing issue with TemplateEstimate = %s. The clone has been excluded from analysis." %(sample[0], seq, fields[cti]))
                    try: seqs[seq][sample[0]] = ct
                    except KeyError:
                        seqs[seq] = dict((key[0], 0) for key in samples)
                        seqs[seq][sample[0]] = ct
                else:
                    try:
                        if seqs[seq][sample[0]] != -1:
                            seqs[seq][sample[0]] += ct #int(value) = count preference value + any other matching sequence values
                    except KeyError:
                        seqs[seq] = dict((key[0], 0) for key in samples) # so it is added to the seqs dictionary
                        seqs[seq][sample[0]] = ct #if seq not already in the dictionary, add it, making the value ct
        openFile.close()
    if len(seqs)<=0:
        error_log(outputDir, 'Due to productivityOnly settings, all sequences are excluded from the comparison of %s vs %s. The comparison was omitted.\n' %(samples[0][0], samples[1][0]))
        return
    seqsfull = [dict({'sequence':sequ}, **seqs[sequ]) for sequ in seqs] # reformat the dictionary to have 'sequence' as the key, and the sequence as its value
    totalsfull = dict((key[0], 0) for key in samples) #set start value for total counts

    ######## Generate Totals Dictionary
    exnums=[]
    for seq in xrange(len(seqsfull)):
        if -1 not in list(seqsfull[seq].values()):
            totalsfull[samples[0][0]] += seqsfull[seq][samples[0][0]] #add counts each loop to find total counts
            totalsfull[samples[1][0]] += seqsfull[seq][samples[1][0]] #add counts each loop to find total counts
        else:
            exnums.append(seq) #exclude samples that had erros
    for indnum in reversed(exnums): #reversed so that the indices stay the same
        del seqsfull[indnum] #delete the sequence from analysis

    if len(seqsfull)<=0:
        error_log(outputDir, 'The files %s and %s are incompatible with this tool. Their comparison was omitted. Check file versioning and inputTemplateEstimate count values.' %(samples[0][0], samples[1][0]))
        return

    ####### Min Totals Filtering
    minTotEx=[]
    totalsminTotal = dict((key[0], 0) for key in samples) #set start value for total counts
    if minTotal==0:
        seqsminTotal=seqsfull
        totalsminTotal=totalsfull
    else:
        minTotExapp=minTotEx.append
        seqsminTotal=[]
        seqsminTotalapp=seqsminTotal.append
        for seq in seqsfull:
            if (seq[samples[0][0]]+seq[samples[1][0]])>=minTotal:
                seqsminTotalapp(dict(seq))
                totalsminTotal[samples[0][0]] += seq[samples[0][0]]
                totalsminTotal[samples[1][0]] += seq[samples[1][0]]
            else:
                minTotExapp(dict(seq))

    if source == 'nucleated' and nucleatedCells.values() != [0,0]:
        totalsfull=nucleatedCells
    elif source == 'nucleated' and nucleatedCells.values() == [0,0]:
        print('nucleated cell counts encountered errors. nucleatedCells = '+nucleatedCells+(' the pair: %s vs %s was ignored.' %(samples[0][0], samples[1][0])))
        error_log(outputDir, 'nucleated cell counts encountered errors. nucleatedCells = '+nucleatedCells+(' the pair: %s vs %s was ignored.' %(samples[0][0], samples[1][0])))
        return

    return totalsfull, totalsminTotal, seqsfull, seqsminTotal, minTotEx



 #############     ##############  ##################         ####         #############     ####  ####      ####      ########      #####      #####  ####         ####         ####
 ###############   ##############  ##################        ######        ###############   ####  #####     ####    ############    #######  #######  ####        ######        ####
 ####        ####  ####                   ####              ########       ####        ####  ####  ######    ####   ####      ####   ################  ####       ########       ####
 ####       ####   ####                   ####             ####  ####      ####       ####   ####  #######   ####  ####        ####  #### ###### ####  ####      ####  ####      ####
 #############     #########              ####            ####    ####     #############     ####  ########  ####  ####        ####  ####  ####  ####  ####     ####    ####     ####
 ##############    #########              ####           ##############    ##############    ####  #### #### ####  ####        ####  ####        ####  ####    ##############    ####
 ####       ####   ####                   ####          ################   ####       ####   ####  ####  ########  ####        ####  ####        ####  ####   ################   ####
 ####        ####  ####                   ####         ####          ####  ####        ####  ####  ####   #######   ####      ####   ####        ####  ####  ####          ####  ####
 ###############   ##############         ####         ####          ####  ###############   ####  ####    ######    ############    ####        ####  ####  ####          ####  ##############
 #############     ##############         ####         ####          ####  #############     ####  ####     #####      ########      ####        ####  ####  ####          ####  ##############

def find(str, ch):
    return [i for i, ltr in enumerate(str) if ltr == ch]


#### unique to betabinomial ####

# beta binomial log likelihood
def bbll(k, n, mu, v):

    alpha = mu*v
    beta = (1-mu)*v
    result = -scipy.log(n+1) - betaln(n-k+1, k+1) + betaln(k+alpha, n-k+beta) - betaln(alpha, beta)
    return result

# calculate bb dispersion parameter for a given frequency and variance
def bbnDisp(mu,v):
    alpha = mu*v
    beta = (1-mu)*v
    hyperparameter = 1/(alpha + beta + 1)
    return hyperparameter

#variance function for bbn
def v(n, coefs):
    # note: kappa > 2 is required by unimodality (alpha, beta > 1)
    return 2 + scipy.poly1d(coefs)(scipy.log(n))**2

# joint log likelihood
def l(dat, mu, coefs):
    #print coefs
    result = sum(ct*bbll(n1, n1+n2, mu, v(n1+n2, coefs)) for (n1, n2, ct) in dat)
    #result = sum(ct*term(n1, n2, coefs) for (n1, n2, ct) in dat)
    #print coefs, result
    assert result <= 0.
    return result

def objective(x, mu, dat):
    return -l(dat, mu, x)

def MLE(dat):
    #for line in dat:
    #    print '\t'.join(map(str, line))
    #print dat
    mu = sum(n[0]*n[2] for n in dat)/sum(float((n[0]+n[1])*n[2]) for n in dat)
    x0 = (1., 1.)
    result = minimize(objective, x0, args=(mu, dat), method='Nelder-Mead')#, method='L-BFGS-B', bounds=[(0., None)]*len(coefs))#, options={'maxfev':1e10}, tol=1e-3)
    if not result.success:
        sys.exit(result)
    return result.x

def trainCoefs(replicateSamplesFile, outbase, items):
    '''train coefficients for betabinomial correction'''
    from scipy.special import betaln, gammaln, digamma, polygamma, logit, expit
    from scipy.optimize import minimize
    assert replicateSamplesFile is not None
    indx = find(replicateSamplesFile, '/')
    sampleloc = replicateSamplesFile[:indx[-1]+1]
    repList = [x.rstrip().split(',') for x in open(replicateSamplesFile, 'rU').readlines()]
    replicates = [[x[0], sampleloc + x[1]] for x in repList]
    assert len(replicates) == 2

    repTotals, totals, repSeqs, datoriginal, minTotEx = tsvUnion(outbase, replicates, #totals, dat = totals, seqs from tsvUnion
                                                                       items['count'],
                                                                       items['minTotal'],
                                                                       items['productiveOnly'])

    datDict = {}
    for (n1, n2) in [tuple(seq[sample[0]] for sample in replicates) for seq in repSeqs]:
        if (n1, n2) not in datDict:
            datDict[(n1, n2)] = 1
        else:
            datDict[(n1, n2)] += 1
    coefs = MLE([(x[0], x[1], datDict[x]) for x in datDict]) ## coefs are the beta1 and beta2 parameter of the beta distribution which describes a density function for a given probability
    print ('betabinomial coeficients = ' +str(coefs))
    return(coefs)


    ##########     ####    ##########     ####      ####  ####  ##############  ####    ##########            ####         ####      ####    ##########     ##############
  ####      ####   ####  ##############   #####     ####  ####  ##############  ####  ##############         ######        #####     ####  ##############   ##############
 ####        ####  ####  ####       ####  ######    ####  ####  ####            ####  ####       ####       ########       ######    ####  ####       ####  ####
 ######            ####  ####             #######   ####  ####  ####            ####  ####                 ####  ####      #######   ####  ####             ####
    #######        ####  ####             ########  ####  ####  #########       ####  ####                ####    ####     ########  ####  ####             #########
      #########    ####  ####    #######  #### #### ####  ####  #########       ####  ####               ##############    #### #### ####  ####             #########
           ######  ####  ####    #######  ####  ########  ####  ####            ####  ####              ################   ####  ########  ####             ####
 ####        ####  ####  ####       ####  ####   #######  ####  ####            ####  ####       ####  ####          ####  ####   #######  ####       ####  ####
  ####      ####   ####  ##############   ####    ######  ####  ####            ####  ##############   ####          ####  ####    ######  ##############   ##############
    ##########     ####    ##########     ####     #####  ####  ####            ####    ##########     ####          ####  ####     #####    ##########     ##############

def Pval(n2, n1, N2, N1, method='binomial', alternative='two-sided', coefs = None):
    #order switched from n1, n2, N1, N2 so that sample2>sample1 is displayed when the 'greater' alternative is selected
    '''
    ARGUMENTS: sample1 count, sample2 count, sample1 full total, sample2 full total, method, alternative
    RETURNS: pValue for those two sample counts
    PURPOSE: calculate pValue
    ***Adjusted in analyst version for betabinomial
    '''
    if method == 'binomial':
        p = N1/(1.*(N1+N2))
        if alternative == 'less':
            return binom.cdf(n1, n1+n2, p)
        elif alternative == 'greater':
            return binom.sf(n1-1, n1+n2, p) #n1-1 is to include the sequence itself
        elif alternative == 'two-sided':
            return binom_test(n1, n1+n2, p)
    if method == 'fisher':
        return fisher_exact([[n1, n2], [N1-n1, N2-n2]], alternative)[1]
    if method == 'betabinomial':
        assert coefs is not None
        mu = N1/(1.*(N1+N2))
        if alternative == 'greater':
            return sum(scipy.exp(bbll(x, n1+n2, mu, v(n1+n2, coefs))) for x in range(n1+1))
        elif alternative == 'less':
            return sum(scipy.exp(bbll(x, n1+n2, mu, v(n1+n2, coefs))) for x in range(n1,n1+n2+1))
        elif alternative == 'two-sided':
            thisProb = scipy.exp(bbll(n1, n1+n2, mu, v(n1+n2, coefs)))
            tot = 0
            for x in range(n1+n2+1):
                term = scipy.exp(bbll(x, n1+n2, mu, v(n1+n2, coefs)))
                if term <= thisProb:
                    tot += term
            return tot

def significance(pVals, correction='none', alpha=.01):
    '''
    ARGUMENTS: pVals list of sequence dictionaries with pValue included, correction, alpha (significance threshold)
    RETURNS: dictionary of lists of number rejected, significant sequences, and insignificant sequences
    PURPOSE: find which sequences are differentially abundant based on pValues (and correct as specified)
    ***Adjusted in analyst version for benjamini yekutieli
    '''
    pVals = sorted(pVals, key=lambda x: x['pvalue'])
    m = len(pVals)
    k = -1
    try: #allows alpha = 1
        if correction == 'bh':
            while pVals[k+1]['pvalue'] <= alpha*(k+2)/(1.*m): k += 1
        elif correction == 'by':
            Cm = sum([1./i for i in xrange(1,m+1)])                         #added 2/23/17
            while pVals[k+1]['pvalue'] <= alpha*((k+2)/(m*Cm)): k += 1
        elif correction == 'bonferroni':
            while pVals[k+1]['pvalue'] <= alpha/(1.*m): k += 1
        elif correction == 'none':
            while pVals[k+1]['pvalue'] <= alpha: k += 1
        return {'rejections':k+1, 'significant':pVals[:(k+1)], 'insignificant':pVals[(k+1):]}
    except IndexError: # this allows alpha = 1 (or alpha greater than expressed pVals)
        return {'rejections':m, 'significant':pVals, 'insignificant':[]}


 ###########       ####  ####          ####  ##############  #############        ##########     ####  ##################  ####        ####
 ##############    ####  ####          ####  ##############  ###############    ####      ####   ####  ##################   ####      ####
 ####      #####   ####  ####          ####  ####            ####        ####  ####        ####  ####         ####           ####    ####
 ####        ####  ####   ####        ####   ####            ####       ####   ######            ####         ####            ####  ####
 ####        ####  ####    ####      ####    #########       #############        #######        ####         ####             ########
 ####        ####  ####     ####    ####     #########       ###########            #########    ####         ####               ####
 ####        ####  ####      ####  ####      ####            ####     ####               ######  ####         ####               ####
 ####      #####   ####       ########       ####            ####      ####    ####        ####  ####         ####               ####
 ##############    ####        ######        ##############  ####       ####    ####      ####   ####         ####               ####
 ###########       ####         ####         ##############  ####        ####     ##########     ####         ####               ####

def diversity_calcs(totals, seqs, samples): #these are seqsfull and totalfull, so not impacted by minTotal
    '''
    ARGUMENTS: full totals, full sequences, samples list (for sample names)
    RETURNS: Morisita Index, Jaccard Index, tcr overlap, and Bhattacharyya's Coefficient of the set
    '''
    totaltcr=0
    overlaptcr=0
    totaljac=0
    overlapjac=0

    D1=0 # D1: Simpson's index of sample 1: (sum of (count/totalcount)^2)
    D2=0
    sumproduct=0
    BC = 0 #Bhattacharyya Coefficient

    for seq in seqs:
        s1count = seq[samples[0][0]]
        s2count = seq[samples[1][0]]
        s1freq = s1count/(1.*totals[samples[0][0]])
        s2freq = s2count/(1.*totals[samples[1][0]])

        #tcr and jaccard
        totaltcr+=(s1count+ s2count)
        totaljac+=1
        if s1count!=0 and s2count!=0:
            overlaptcr+=(s1count+ s2count)
            overlapjac+=1

        #morisita
        D1+=(s1freq)**2 #frequency = count s1/total counts s1
        D2+=(s2freq)**2 #frequency = count s2/total counts s2
        sumproduct+=(s1count)*(s2count) #count s1 * count s2

        #Bhattacharyya #sum of the square root of the product of frequencies
        BC+=((s1freq)*(s2freq))**(.5)

    mor = (2.*sumproduct)/(1.*(D1+D2)*((totals[samples[0][0]])*totals[samples[1][0]])) #m=(2.*sum(xi*yi))/((D1+D2)*X*Y)
    jac = (overlapjac/(1.*totaljac))
    tcr = (overlaptcr/(1.*totaltcr))

    return mor, jac, tcr, BC


 ##############  ####    ##########     ####       ####  #############     ##############     ##########
 ##############  ####  ##############   ####       ####  ###############   ##############   ####      ####
 ####            ####  ####       ####  ####       ####  ####        ####  ####            ####        ####
 ####            ####  ####             ####       ####  ####       ####   ####            ######
 #########       ####  ####             ####       ####  #############     #########          #######
 #########       ####  ####    #######  ####       ####  ###########       #########            #########
 ####            ####  ####    #######  ####       ####  ####     ####     ####                      ######
 ####            ####  ####       ####  ####       ####  ####      ####    ####            ####        ####
 ####            ####  ##############   ##############   ####       ####   ##############   ####      ####
 ####            ####    ##########       ##########     ####        ####  ##############     ##########

def addbarlabels(heights):
    '''
    ARGUMENTS: artist rectangles of the bar graph
    RETURNS: NONE
    PURPOSE: Add count labels to each bar of the graph
    '''
    horiz=[.05, .15, .25, .35, .45, .55, .65, .75, .85, .95]
    for i in xrange(len(heights)):
        plt.gca().text(horiz[i], heights[i],
                '%s' % ("{:,}".format(int(heights[i]))),
                horizontalalignment='center', verticalalignment='bottom', zorder=10, size=6)


def histogram(dat, outbase, base='linear', labels='on', plotformat='pdf'):
    '''
    ARGUMENTS: list of sequence dictionaries, output location, y-axis scale (linear or symlog)
    RETURNS: NONE
    PURPOSE: save a pdf figure of a histogram sorted by pValues
    '''
    plt.figure(1)
    plt.gca().xaxis.grid('off')
    plt.gca().spines['bottom'].set_color('#7d7d7d')
    barheights=plt.hist([seq['pvalue'] for seq in dat], bins=scipy.linspace(0, 1, num=11), color='#f27a63', edgecolor='white', zorder=10)[0]
    plt.tick_params(axis='both', which='both', bottom='off', top='off', right='off', left='off')
    plt.xlabel('p value'.upper()),  plt.ylabel('number of clones'.upper())
    plt.xlim(xmin=-.02, xmax=1.02)
    if base=='linear':
        plt.ylim(ymin=.75, ymax=1.02*max(barheights))#, ymax=1.1*maxall)
        plt.yscale('linear')
    elif base=='log':
        plt.ylim(ymin=.75, ymax=1.2*max(barheights))#, ymax=1.1*maxall)
        plt.yscale('log')
    if labels=='on':
        addbarlabels(barheights)
    plt.savefig(outbase+'.pvalueHistogram.'+plotformat)
    plt.close()

def minTotalsPlot(seqs, samples, outbase, plotformat='pdf'):
    '''
    ARGUMENTS: list of sequence dictionaries, list of sample names and locations, output location
    RETURNS: dictionary of how many sequences are included for each minTotal value
    PURPOSE: save a pdf figure of how minTotal effects sequence inclusion
             save a tsv reflecting the same information
    '''
    newdic={}
    for seq in seqs:
        sumcounts=seq[samples[0][0]]+seq[samples[1][0]]
        while sumcounts>=0:
            try: #faster than the if in / else
                newdic[sumcounts] += 1
            except KeyError:
                newdic[sumcounts] = 1
            sumcounts-=1
    plt.figure(5)
    plt.plot(list(newdic.keys())[1:], list(newdic.values())[1:], color='#f27a63', linewidth='3', zorder=5)
    plt.xlim(xmin=1, xmax=(max(list(newdic.keys()))+1)),    plt.ylim(ymin=.75)
    plt.yscale('log'),    plt.xscale('log')
    plt.gca().spines['bottom'].set_color('#7d7d7d'),    plt.gca().spines['bottom'].set_zorder(4)
    plt.tick_params(axis='both', which='both', bottom='off', top='off', right='off', left='off')
    plt.xlabel('minTotal'.upper()),    plt.ylabel('number of sequences'.upper())
    plt.savefig(outbase+'.minTotalPlot.'+plotformat)
    plt.close()

    totalseqs= newdic[0]
    with open(outbase+'.minTotal.tsv', 'w') as f: # write the minTotal data to a tsv
        f.write('minTotal value\t# of included sequences\t% of included sequences\n')
        for key in newdic:
            f.write('%s\t%s\t%.5f\n' %(key, "{:,}".format(int(newdic[key])), newdic[key]/(1.*totalseqs)*100))
    return newdic #for testing purposes only

def splitTwoSides(dat, minTotEx, tots1, tots2, samples, verbose):
    '''
    ARGUMENTS: list of sequence dictionaries, list of minTotal-excluded data, sample1 full total, sample2 full total, list of sample names and locations, verbosity of scatter
    RETURNS: dictionaries of abundance coordinates for sequences that are greater in s1 than 2, greater in s2 than 1, not differentially abundant, excluded by minTotal, or equal in the samples
    PURPOSE: sort the sequences for coloring in the scatter plot
    '''
    greater12={}
    greater21={}
    insig={}
    minEx={}
    equal={}
    for seq in dat['significant']:
        if seq[samples[0][0]]/(1.*tots1)<seq[samples[1][0]]/(1.*tots2): #orange
            try: greater21[seq[samples[0][0]], seq[samples[1][0]]]+=1
            except KeyError: greater21[seq[samples[0][0]], seq[samples[1][0]]]=1
        elif seq[samples[0][0]]/(1.*tots1)>seq[samples[1][0]]/(1.*tots2): #blue
            try: greater12[seq[samples[0][0]], seq[samples[1][0]]]+=1
            except KeyError: greater12[seq[samples[0][0]], seq[samples[1][0]]]=1
        else:
            try: equal[seq[samples[0][0]], seq[samples[1][0]]]+=1
            except KeyError: equal[seq[samples[0][0]], seq[samples[1][0]]]=1
    if verbose or verbose=='extra':
        for seq in dat['insignificant']:
            try: insig[seq[samples[0][0]], seq[samples[1][0]]]+=1
            except KeyError: insig[seq[samples[0][0]], seq[samples[1][0]]]=1
        if verbose=='extra':
            for seq in minTotEx:
                try: minEx[seq[samples[0][0]], seq[samples[1][0]]]+=1
                except KeyError: minEx[seq[samples[0][0]], seq[samples[1][0]]]=1
    return greater12, greater21, minEx, insig, equal

def scatterplot(dat, totals, outbase, samples, minTotEx, items):
    '''
    ARGUMENTS: list of sequence dictionaries, full totals, output location, list of sample names and location, list of sequence dictionaries excluded by minTotal, scatter plot settings
    RETURNS: list of dictionaries of plotted points for testing
    PURPOSE: save a pdf figure of a scatter plot with user settings
    '''
    [scale, base, verbose, line, minTotLine, points, size, lverbose] = items['scatter']
    minTotal=items['minTotal']
    #orange=(0.9450980392156862, 0.47843137254901963, 0.4)
    #blue=(0.34901960784313724, 0.7215686274509804, 0.807843137254902)
    #purple=(0.5607843137254902, 0.3686274509803922, 0.6039215686274509)
    tots=[totals[samples[0][0]], totals[samples[1][0]]]
    greater12, greater21, minEx, insig, equal = splitTwoSides(dat, minTotEx, tots[0], tots[1], samples, verbose)

    alldata=[greater12, greater21, equal, insig, minEx]
    datasets=[]
    for i in alldata:
        if len(i)>0:            datasets.append(i)
    if len(datasets)==0:
        items['scatter'][2]='extra'
        error_log(outputDir, 'scatter plot verbosity was changed to extra because nothing was plotted under current settings.')
        return scatterplot(dat, totals, outbase, samples, minTotEx, items)
    maxx=0
    maxy=0
    plt.figure(2, figsize=(6*1.25, 6)) #scatterplot
    cts=list(np.zeros((len(datasets)), dtype=np.int))
    for num in xrange(len(datasets)):
        cts[num]=sum(list(datasets[num].values()))
        outline='none'
        width=size/40
        if points =='spwo' or points=='sp':
            if len(datasets[num])!=4: # this is in case of a very, very weird glitch, where if the array is the same length as the number of points, the points are all weird colors
                rgba_colors = np.zeros((1, 4)) #in order to combat this, check if it is the same length and if it isn't, the array may be made
            else:
                rgba_colors = np.zeros((2,4)) #if it is the same length, the array must be a different size
            alphas=1
            if points=='spwo': #'wo' = with outline
                outline='black'
        else: #transparency, transparent outlines, white transparent, white width
            rgba_colors = np.zeros((len(datasets[num]),4))
            add=2.0
            for i in datasets[num]:
                datasets[num][i]=(datasets[num][i]/(datasets[num][i]+add))
            if points=='wpwwo':
                alphas = 1
                width=list(datasets[num].values())
                width=[(size/4.0)*w for w in width]
            else: #white, transparency, transparent outlines
                alphas = list(datasets[num].values())
        if datasets[num]==greater12:    rgba_colors[:,0:3] =list(matplotlib.colors.colorConverter.to_rgb('#4a99ad')) #Plot Dark Blue
        elif datasets[num]==greater21:    rgba_colors[:,0:3] = list(matplotlib.colors.colorConverter.to_rgb('#f27a63')) #Plot Red
        elif datasets[num]==insig:    rgba_colors[:,0:3] = list(matplotlib.colors.colorConverter.to_rgb('#b2adaa')) #Plot Stone
        elif datasets[num]==minEx:    rgba_colors[:,0:3] = list(matplotlib.colors.colorConverter.to_rgb('#ececeb')) #Plot Stone Lighter
        elif datasets[num]==equal:    rgba_colors[:,0:3] = list(matplotlib.colors.colorConverter.to_rgb('#8f5e9a')) #Purple .. Both
        rgba_colors[:, 3] = alphas
        maxx1 = max(v[0] for v in list(datasets[num].keys()))
        maxy1 = max(v[1] for v in list(datasets[num].keys()))
        if scale=='frequency':
            maxx1=1.*maxx1/tots[0]
            maxy1=1.*maxy1/tots[1]
        maxx=max(maxx, maxx1)
        maxy=max(maxy, maxy1)
        if points=='tpwo':
            outline='black'
        elif points=='wpwto' or points=='wpwwo':
            outline=rgba_colors
            rgba_colors='white'
        if scale == 'frequency':
            datasets[num]['plot']=plt.scatter([1.0*v[0]/tots[0] for v in list(datasets[num].keys())], [1.0*v[1]/tots[1] for v in list(datasets[num].keys())],
                     color=rgba_colors, edgecolor=outline, linewidth=width, clip_on=False, zorder=5, s=size**2)
        else:
            datasets[num]['plot']=plt.scatter([v[0] for v in list(datasets[num].keys())], [v[1] for v in list(datasets[num].keys())],
                     color=rgba_colors, edgecolor=outline, linewidth=width, clip_on=False, zorder=5, s=size**2)
        if datasets[num]==greater12:      greater12['title']='A > B:                '+ "{:,}".format(cts[num]) #this adds the counts (with commas) to the titles
        elif datasets[num]==greater21:    greater21['title']='B > A:                '+ "{:,}".format(cts[num]) #this adds the counts (with commas) to the titles
        elif datasets[num]==insig:        insig['title'] ='Not Significant:  '+ "{:,}".format(cts[num]) #this adds the counts (with commas) to the titles
        elif datasets[num]==minEx:        minEx['title']='Excluded:           '+ "{:,}".format(cts[num]) #this adds the counts (with commas) to the titles
        elif datasets[num]==equal:        equal['title']='A = B:                 '+ "{:,}".format(cts[num]) #this adds the counts (with commas) to the titles
    maxall=max(maxx, maxy)
    a=plt.scatter(-5, -5, color='white', edgecolor='white')
    lns=[a]
    lnnames=['']
    if line=='frequency':
        lnnames.append('Frequency Equality')
    elif line == 'abundance':
        lnnames.append('Abundance Equality')
    if line=='frequency' and scale=='abundance':
        x=np.linspace(0, max(1.1*maxall*tots[0], 1.1*maxall*tots[1]), maxall)
        ln1,=(plt.plot(x/(1.*tots[1]), x/(1.*tots[0]),  linestyle='--', color='#7d7d7d'))
    elif line=='abundance' and scale == 'abundance':
        x=np.linspace(0, maxall, maxall)
        ln1,=(plt.plot(x, x, linestyle='--', color='#7d7d7d'))
    elif line=='frequency' and scale == 'frequency':
        x=np.logspace(np.log(1.0/max(tots)), maxall, int(maxall*max(tots)+5)) #logscale base = 10, 1st and 2nd entries are exponent bounds, third entry is number of points.
        ln1,=(plt.plot(x, x, linestyle='--', color='#7d7d7d'))
    elif line == 'abundance' and scale=='frequency':
        x=np.linspace(0, 1)
        ln1,=(plt.plot(x*(1.*tots[1]), x*(1.*tots[0]),  linestyle='--', color='#7d7d7d'))
    if minTotLine=='on' and minTotal != 0:
        if scale=='abundance':
            cx=np.linspace(0,minTotal+1)
            cy=np.linspace(minTotal,-1)
        elif scale=='frequency':
            cx=np.linspace(0,1.*(minTotal)/tots[0])
            cy=np.linspace(1.*minTotal/tots[1],0)
        ln2,=(plt.plot(cx, cy, '--', color='#c2141a', alpha=.5, zorder=4))
        lnnames.append('Threshold for\nStatistical Comparison')
    try: lns.append(ln1)
    except NameError: pass
    try: lns.append(ln2)
    except NameError: pass
    plt.gca().spines['bottom'].set_color('#7d7d7d')
    plt.gca().spines['bottom'].set_zorder(4)
    plt.tick_params(axis='both', which='both', bottom='off', top='off', right='off', left='off')
    if base=='symlog':
        if scale=='frequency':
            linthreshex=1.0/tots[0]
            linthreshey=1.0/tots[1]
        elif scale=='abundance':
            linthreshex=1
            linthreshey=1
        plt.xscale('symlog', linthreshx=linthreshex, linscalex=.33)
        plt.yscale('symlog', linthreshy=linthreshey, linscaley=.33)
    elif base=='linear': plt.xscale('linear'), plt.yscale('linear')
    plt.xlim(xmin=0, xmax=1.1*maxall),   plt.ylim(ymin=0, ymax=1.1*maxall)
    if scale == 'frequency':
        plt.xlabel(('A: '+samples[0][0]+' '+scale+' (%)').upper()),    plt.ylabel(('B: '+samples[1][0]+' '+scale+' (%)').upper())
    else:
        plt.xlabel(('A: '+samples[0][0]+' '+scale).upper()),    plt.ylabel(('B: '+samples[1][0]+' '+scale).upper())
    if base =='symlog':
        xticks = plt.gca().xaxis.get_major_ticks()
        xticks[0].label1.set_visible(False)
        yticks = plt.gca().yaxis.get_major_ticks()
        yticks[0].label1.set_visible(False)
        if scale == 'frequency': plt.text(-linthreshex/2.5, -linthreshey/2.7, '0', fontsize=8)
        elif scale == 'abundance': plt.text(-.38, -.34, '0', fontsize=8)

    #legend:
    box = plt.gca().get_position()
    plt.gca().set_position([.1*box.x0+box.width*.1, box.y0, box.width*.83, box.height]) #rescale the plot to be square, so that the legend is outside (to the right)
    pts=[v['plot'] for v in datasets]
    ptnames=[v['title'] for v in datasets]
    if not lverbose or len(lns)==1: # ==1 because they are defined as filled lists for the extra spacing
        plt.legend([p for p in pts], [p for p in ptnames],
                   loc='center', bbox_to_anchor=(1.205, .5),
                   handletextpad=.2, borderpad=.8)
    elif lverbose:
        lnspts=pts+lns
        lnsptsnames=ptnames+lnnames
        plt.legend([l for l in lnspts], [l for l in lnsptsnames],
                   loc='center', bbox_to_anchor=(1.205, .5),
                   handletextpad=.2, borderpad=.8)

    plt.savefig(outbase+'.scatter.'+items['format'])
    plt.close()
    return datasets # for testing purposes

def fixrcplotdefaults():
    '''
    ARGUMENTS: NONE
    RETURNS: NONE
    PURPOSE: set default matplotlib settings so that they don't have to be set in each plotting function
    '''
    matplotlib.style.use('classic')
    matplotlib.rc('font', family='sans-serif', size=10)
    matplotlib.rcParams['text.color']='#7d7d7d'
    matplotlib.rc('axes', grid=True, labelsize=10, labelcolor='#7d7d7d', edgecolor='#d9dad8')
    matplotlib.rc('xtick', labelsize=8, color='#7d7d7d')
    matplotlib.rc('ytick', labelsize=8, color='#7d7d7d')
    matplotlib.rc('grid', color='#d9dad8', linestyle='-')
    matplotlib.rc('legend', fancybox=True, numpoints=1, scatterpoints=1, fontsize=8)
    matplotlib.rcParams['figure.figsize']= 6, 6
    matplotlib.rcParams['savefig.format']='pdf'
    matplotlib.rcParams['savefig.transparent']=True
    #edit font print type so they export as words in illustrator and other formats
    matplotlib.rcParams['ps.fonttype'] = 42
    matplotlib.rcParams['pdf.fonttype'] = 42
    matplotlib.rcParams['svg.fonttype'] = 'none'
    #'eps', 'jpg', 'jpeg', 'pdf', 'png', 'ps', 'raw', 'svg', 'svgz'



 #####      #####         ####         ####  ####      ####
 #######  #######        ######        ####  #####     ####
 ################       ########       ####  ######    ####
 #### ###### ####      ####  ####      ####  #######   ####
 ####  ####  ####     ####    ####     ####  ########  ####
 ####        ####    ##############    ####  #### #### ####
 ####        ####   ################   ####  ####  ########
 ####        ####  ####          ####  ####  ####   #######
 ####        ####  ####          ####  ####  ####    ######
 ####        ####  ####          ####  ####  ####     #####

def main(outputDir, ID, samples, items, coefs = None):
    outputDir += '/'
    fixrcplotdefaults()

    if not items['batch']:
        samples.sort()    # remove this for batchmode
        if items['reverse']=='yes':
            samples.reverse()
    outbase = outputDir +'/'+ ID # ID is subjectID_s1name_VS_s2name, same as directory # outputDir is a folder, so this creates files inside it
    assert len(samples) == 2, 'somehow a number other than two samples were selected' #make sure user only selected two samples
    try: totalfull, totals, seqs, datoriginal, minTotEx = tsvUnion(outputDir, samples, #totals, dat = totals, seqs from tsvUnion
                                                                       items['count'],
                                                                       items['minTotal'],
                                                                       items['productiveOnly'],
                                                                       items['source'])
    except TypeError:
        if items['batch']:
            if os.path.isfile(outputDir+'error.log'):
                error = ('Errors in pairing of %s vs. %s:\n'%(samples[0][0], samples[1][0])) + ('').join(open(outputDir+'error.log', 'r').readlines())
                os.remove(outputDir+'error.log')
                os.rmdir(outputDir) #this will be an empty file because the error file was already deleted
            return error
        else: sys.exit()
    nclones = len(datoriginal) #number of nucleotide or amino acid sequences
    if nclones == 0: # if there was nothing left after count filtering, reset count filter to 0
        if items['batch']:
            if os.path.isfile(outputDir+'error.log'):
                os.remove(outputDir+'error.log')
            return ('No rearrangments met the minTotal threshold in %s vs %s. The comparison was omitted.\n' %(samples[0][0], samples[1][0]))    #####for batch mode only
        error_log(outputDir, ' WARNING: minTotal reset to 0 because no rearrangments met the desired threshold of %s' %(items['minTotal']))
        datoriginal = seqs
        totals = totalfull
        nclones = len(datoriginal)
        items['minTotal'] = 0
        minTotEx=[]
    tot1=totalfull[samples[0][0]]
    tot2=totalfull[samples[1][0]]
    pHash = {}
    for seq in datoriginal: #aka for seq in seqs
        n1, n2 = [seq[sample[0]] for sample in samples] #counts of the samples respectively
        if (n1, n2) not in pHash: #if the comparison is not in the dictionary already, put it there
            pHash[(n1, n2)] = Pval(n1, n2,
                                   tot1, tot2,
                                   items['method'], items['alternative'], coefs)
        seq['pvalue'] = pHash[(n1, n2)]
    dat = significance(deepcopy(datoriginal),
                       items['correction'],
                       items['alpha'])
    mor, jac, tcr, bhattC=diversity_calcs(totalfull, seqs, samples)
    if bhattC == 0:
        bhattD= "inf"
    else:
        bhattD= -math.log(bhattC)

    if items['batch']:
        if items['batch_count_minTot_seqs']:
            batchdic=batch(dat['significant'], dat['insignificant'], minTotEx, samples, tot1, tot2)
        else:
            batchdic=batch(dat['significant'], dat['insignificant'], [], samples, tot1, tot2)
        batchdic["morisita_index"]=mor
        batchdic["tcr_overlap"]=tcr
        batchdic["jaccard_index"]=jac
        batchdic["bhattacharyya_distance"]= bhattD

    if items['verbose']=='scatter':
        scatterplot(deepcopy(dat), totalfull, outbase, samples, minTotEx, items)
        if items['batch']:
            return batchdic
        return

    minTotalsPlot(seqs, samples, outbase, items['format'])
    histogram(datoriginal, outbase, items['hist'][0], items['hist'][1], items['format'])
    scatterplot(deepcopy(dat), totalfull, outbase, samples, minTotEx, items)
    with open(outbase+'.differentialAbundance.tsv', 'w') as f: #write all the information to a tsv
        f.write('\n'.join(['# sample A: %s' %samples[0][0],
                           '# sample B: %s' %samples[1][0],
                           '# count: %-100s' %(items['count']),
                           '# minimum total filter: %-60s' %(items['minTotal']),
                           '# productive only: %-65s'%(items['productiveOnly']),
                           '# method: %-74s' %(items['method']),
                           '# alternative: %-69s'%(items['alternative']),
                           '# multiple test correction: %-56s'%(items['correction']),
                           '# alpha (significance threshold): %-45s' %(items['alpha']),
                           '# number of differentially abundant clones: %s'%(dat['rejections']),
                           '# number of clones: %s' %("{:,}".format(nclones)),
                           #'# sample1 total count (post minTotal): %s' %(totals[samples[0][0]]),
                           #'# sample2 total count (post minTotal): %s' %(totals[samples[1][0]]),
                           '# sample A total count: %s' %("{:,}".format(tot1)),
                           '# sample B total count: %s' %("{:,}".format(tot2)),
                           '# % of sample A excluded by minTotal value: ' +str(round(100-(totals[samples[0][0]]/(1.*tot1))*100, 2)),
                           '# % of sample B excluded by minTotal value: ' +str(round(100-(totals[samples[1][0]]/(1.*tot2))*100, 2)),
                           '# TCR Overlap: %s'%(tcr),
                           "# Morisita index: %s"%(mor),
                           '# Jaccard index: %s'%(jac),
                           '# Bhattacharyya Distance: %s'%(bhattD)
                           ])+'\n')
        f.write('\nsequence\tsampleA_abundance\tsampleB_abundance\tsampleA_frequency\tsampleB_frequency\tpValue\tsignificance\tfold_change\n')
        for seq in dat['significant']:
            f.write(getLine(samples, seq, tot1, tot2, 'significant'))
        if items['verbose'] or items['verbose']=='extra':
            for seq in dat['insignificant']:
                f.write(getLine(samples, seq, tot1, tot2, 'insignificant'))
            if items['verbose'] =='extra':
                for seq in minTotEx:
                    f.write(getLine(samples, seq, tot1, tot2, 'none'))
    matplotlib.rcdefaults() #return matplotlib defaults on computer
    if items['batch']:
        return batchdic

def getLine(samples, seq, tot1, tot2, significance):
    '''
    ARGUMENTS: list of samples name and path, dictionary of sequences, totals of sample 1 and 2, whether the seq set is significant, insignificant, or excluded
    RETURNS: line to be printed in tsv output file
    PURPOSE: create the line to be included in the ouptu file
    '''
    s1ct = seq[samples[0][0]]
    s2ct = seq[samples[1][0]]
    s1freq = s1ct/(1.*tot1)
    s2freq = s2ct/(1.*tot2)
    try:
        fold = s2freq/(1.*s1freq)
        if fold<1:
            fold = -1/(1.*fold) #this is equivalent to -s1ct/s2ct
    except ZeroDivisionError:
        if s2freq != 0 and s1freq == 0:
            fold = 'inf'
        elif s2freq == 0 and s1freq != 0:
            fold = '-inf'
        else:
            fold = 'N/A'
    if significance == 'significant':
        if (s2ct/(1.*tot2))>(s1ct/(1.*tot1)):       sig='B > A'
        elif (s2ct/(1.*tot2))<(s1ct/(1.*tot1)):     sig='A > B'
        else:                                       sig='A = B'
        pval = seq['pvalue']
    elif significance == 'insignificant':
        sig = 'False'
        pval = seq['pvalue']
    else:
        sig = 'N/A'
        pval = 'N/A'
    return ('\t'.join([seq['sequence'], str(s1ct), str(s2ct), str(s1freq), str(s2freq), str(pval), str(sig), str(fold) + '\n']))


 #############            ####         ##################    ##########     ####      ####
 ###############         ######        ##################  ##############   ####      ####
 ####        ####       ########              ####         ####       ####  ####      ####
 ####       ####       ####  ####             ####         ####             ####      ####
 #############        ####    ####            ####         ####             ##############
 ##############      ##############           ####         ####             ##############
 ####       ####    ################          ####         ####             ####      ####
 ####        ####  ####          ####         ####         ####       ####  ####      ####
 ###############   ####          ####         ####         ##############   ####      ####
 #############     ####          ####         ####           ##########     ####      ####

def batch(sig, insig, minTotEx, samples, tot1, tot2):
    '''
    ARGUMENTS: list of significant sequences, list of insignificant sequences, list of sequences excluded by minTotal, list of sample names, total of samples 1 and 2
    RETURNS: dictionary of sample metrics for batch output
    PURPOSE: calculate batch statistics

    '''
    share=0
    in1=0
    in2=0
    contract=0
    expand=0
    maxexpand=0
    maxcontract=0
    maxfold = 0

    for seq in sig:
        cnt1=seq[samples[0][0]]
        cnt2=seq[samples[1][0]]
        if cnt1!=0: in1+=1
        if cnt2!=0: in2+=1
        if cnt1!=0 and cnt2!=0: share+=1
        freq1=cnt1/(1.*tot1)
        freq2=cnt2/(1.*tot2)

        ###contraction and expansion
        diff=freq2-freq1
        if diff>0:
            expand+=1
            if diff>maxexpand: maxexpand=diff
        elif diff<0:
            contract+=1
            if diff<maxcontract: maxcontract=diff
        ###fold change
        try:
            fold = freq2/(1.*freq1)
            if fold<1:
                fold = -1/(1.*fold) #this is equivalent to -s1ct/s2ct
            if math.fabs(fold)>math.fabs(maxfold):
                maxfold = fold
        except ZeroDivisionError:
            pass

    for seq in insig:
        cnt1=seq[samples[0][0]]
        cnt2=seq[samples[1][0]]
        if cnt1!=0: in1+=1
        if cnt2!=0: in2+=1
        if cnt1!=0 and cnt2!=0: share+=1
    for seq in minTotEx:                            ## Added by JR on 11/29/16; includes clones below minTotal for calculating shared clones, percent shared, etc, in batch function
        #minTotEx is empty if config "batch ignore minTotal = False"
        #minTotEx is full if config "batch ignore minTotal = True"
        cnt1=seq[samples[0][0]]
        cnt2=seq[samples[1][0]]
        if cnt1!=0: in1+=1
        if cnt2!=0: in2+=1
        if cnt1!=0 and cnt2!=0: share+=1
    if share != 0: #if it does = 0 then in1 or in2 will be zero resulting in an error
        fracshare=share/(1.*(in1+in2-share)) # "-share" removes overlap
        fracshare1=share/(1.*in1)
        fracshare2=share/(1.*in2)
    else:
        fracshare = 0
        fracshare1 = 0
        fracshare2 = 0
    if maxfold == 0:
        maxfold = 'undefined fold change from or to 0'
    return {'clones_shared':share, 'percent_shared': 100.*fracshare, 'percent_shared_sample_a': 100.*fracshare1, 'percent_shared_sample_b':100.*fracshare2, #calculated from all: sig, insig, minTotal-Excluded (if requested in config)
             'number_expanded':expand, 'maximum_expansion(%)': 100.*maxexpand, 'number_contracted': contract, 'maximum_contraction(%)': 100.*maxcontract, 'maximum_fold_change': maxfold} #calculated only from significant clones


 ##############  #############     #############         ########      #############     ####                ########        ##########
 ##############  ###############   ###############     ############    ###############   ####              ############    ##############
 ####            ####        ####  ####        ####   ####      ####   ####        ####  ####             ####      ####   ####       ####
 ####            ####       ####   ####       ####   ####        ####  ####       ####   ####            ####        ####  ####
 #########       #############     #############     ####        ####  #############     ####            ####        ####  ####
 #########       ###########       ###########       ####        ####  ###########       ####            ####        ####  ####    #######
 ####            ####     ####     ####     ####     ####        ####  ####     ####     ####            ####        ####  ####    #######
 ####            ####      ####    ####      ####     ####      ####   ####      ####    ####             ####      ####   ####       ####
 ##############  ####       ####   ####       ####     ############    ####       ####   ##############    ############    ##############
 ##############  ####        ####  ####        ####      ########      ####        ####  ##############      ########        ##########

def error_log(outputDir, errornote):
    '''
    ARGUMENTS: output directory name, detailed error string
    RETURNS: NULL
    PURPOSE: adds error to errorlog file
    '''
    if errornote == '':
        return
    errornote+='\n'
    # appends to the error log file when called ('a' opens and adds to end of file)
    with open(outputDir+"error.log", "a") as error_log:
        error_log.write(errornote)


if __name__ == "__main__":
    #main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6]) #analyst version
    sys.exit("must be called from rundiffabBatch script")
