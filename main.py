#!/usr/bin/env python3
"""Runs blast on subject and query files, creates an alignment object, creates a distance matrix
and creates heatmaps/line plot"""

'''Python 3 hates me, please run with python 2.X.'''

'''TypeError: a bytes-like object is required, not str on line 59 when I run with python 3.
I've tried using .encode('UTF-8) but it still throws me an error. Please use python2. It runs, 
I promise.
'''
from collections import defaultdict
import sys
import subprocess
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

__author__ = "James Pinney"
__date__ = "2017-4-25"
__course__ = "BIO-184"
__assignment__ = "Final project - Visualizing sequence homology"

print("\n".join([__author__, __date__, __course__, __assignment__]))

class Main(object):
    """ Runs BLASTN with 3 different parameter sets given two files specificed via command line argument"""
    def __init__(self):
        self.query = (sys.argv[1])
        self.subject = (sys.argv[2])
        self.parameters = [[self.query, self.subject, ' 12', ' 100000', ' 0.001', ' 7', ' 3', '-5', ' 1'], [self.query, self.subject, '10', '1000', '0.0001', '5', '2', '-3', '2']]
        self.BLAST = ["blastn", "-query", '{0}', "-subject", '{1}', "-word_size", '{2}', "-max_hsps", '{3}', "-evalue", '{4}', "-gapopen", '{5}', "-gapextend", '{6}', "-penalty", '{7}', "-reward", '{8}', '-outfmt','10']

         
    def doBlast(self):
        """ Runs BLASTN as part of Main """
        #Converts BLAST argument list into a string
        list_to_string = ' '.join(self.BLAST)
        #Adds in first set of parameter values to string
        Add_Values = (list_to_string.format(*self.parameters[0]))
        #Turns string back into a list, so its usable with subprocess
        string_to_list = Add_Values.split()

        #Converts BLAST argument list into a string once again
        list_to_string2 = ' '.join(self.BLAST)
        #Adds in second set of parameter values to string
        Add_Values2 = (list_to_string2.format(*self.parameters[1]))
        #Turns string back into a list, so its usable with subprocess
        string_to_list2 = Add_Values2.split()

        #Runs BLASTN with default values, pipes STDout into variable. shell = false is default for subprocess, no need to explicitly specify
        BLAST_default = subprocess.Popen(["blastn", "-query", self.query, "-subject", self.subject, '-outfmt','10'], stdout=subprocess.PIPE)
        #Runs BLASTN with first parameter set
        BLAST2 = subprocess.Popen(string_to_list, stdout=subprocess.PIPE)
        #Runs BLASTN with second parameter set
        BLAST3 = subprocess.Popen(string_to_list2, stdout=subprocess.PIPE)
       
        #Initializing lists to store BLASTN output
        r_BLAST_default = []
        r_BLAST2 = []
        r_BLAST3 = []
        
        #Turns stdout from BLASTN into list of lists
        for line in iter(BLAST_default.stdout.readline,''):
            r_BLAST_default.append(line.strip().split(','))
        #Check if the list has anything in it as a sanity check
        if len(r_BLAST_default) == 0:
            print('Run failed. :( oops')
        else:
            print ('First run successful!')
        #Turns stdout from BLASTN into list of lists
        for line in iter(BLAST2.stdout.readline,''):
            r_BLAST2.append(line.strip().split(','))
        #Check if the list has anything in it as a sanity check
        if len(r_BLAST2) == 0:
            print('Run failed. :( oops')
        else:
            print ('Second run successful!')
        #Turns stdout from BLASTN into list of lists
        for line in iter(BLAST3.stdout.readline,''):
            r_BLAST3.append(line.strip().split(','))
        #Check if the list has anything in it as a sanity check
        if len(r_BLAST3) == 0:
            print('Run failed. :( oops')
        else:
            print('Third run successful!')
        #Return list of lists
        return [r_BLAST_default, r_BLAST2, r_BLAST3]


SET = Main().doBlast()
#Calls blast, gets respective list. Probably a better way to implement this. 
blast1 = SET[0]
blast2 = SET[1]
blast3 = SET[2]
      

'''TASK 2'''
class Alignment(object):
    '''Creates alignment objects from blast output'''
    def __init__(self, Align):
        self.Align = Align
        #Add extra element for distance calculation to be added later
        self.Align.append(' ')
        self.qseqid = self.Align[0]
        self.sseqid = self.Align[1]
        self.pident = self.Align[2]
        self.length = self.Align[3]
        self.mismatch = self.Align[4]
        self.gapopen = self.Align[5]
        self.qstart = self.Align[6]
        self.qend = self.Align[7]
        self.sstart = self.Align[8]
        self.send = self.Align[9]
        self.evalue = self.Align[10]
        self.bitscore = self.Align[11]
        self.distance = self.Align[12]
     
    def __iter__(self):
        #Makes Alignment object iterable
        return iter(self.Align)


class CALL(object):
    #Pipes BLAST output into Alignment Class, creates list of alignment objects where NT length is not less than 500
    def __init__(self, BLAST):
        self.BLAST = BLAST
    def convert(self):
        filtered = []
        for x in self.BLAST:
            Y = Alignment(x)
            if int(Y.length) >= 500:
                filtered.append(Y)
        return filtered
            

#Calls CALL for each BLAST output which calls Alignment
filter1 = CALL(blast1).convert()
filter2 = CALL(blast2).convert()
filter3 = CALL(blast3).convert()

        

        
class distance(object):
    #Calculates distance from ID% and length for every alignment object
    def __init__(self, filterlist):
        self.filterlist = filterlist
    def stats (self):
        distance = []
        for x in self.filterlist:
        
            DIFF = (100 - float(x.pident)) / 100 
            LENGTH = float(x.length)
            DISTANCE = (DIFF*LENGTH)
            x.distance = str(DISTANCE)   
            distance.append(x)
        return distance
#Calls distance for each of the filtered list
dist1 = distance(filter1).stats()
dist2 = distance(filter2).stats()
dist3 = distance(filter3).stats()



class dict_create(object):
    #Create dict of SUBJECT:DISTANCE pairs. Groups by subject ID and sums distance of same ID
    def __init__(self, dist):
        self.dist = dist
    def sortsum(self):
        totals = defaultdict(int)
        for c in self.dist:
            totals[c.sseqid] += float(c.distance)
            tuples = totals.items()
            SDPAIR = dict(tuples) 
        return SDPAIR
#Calls dict create for each distance list
SDPAIR1 = dict_create(dist1).sortsum()
SDPAIR2 = dict_create(dist2).sortsum()
SDPAIR3 = dict_create(dist3).sortsum()


#Sanity check
if len(SDPAIR1) == 0 or len(SDPAIR2) == 0 or len(SDPAIR3) == 0:
        print('Unsuccessful generation of distance matrix. IDK what happened.')
else:
        print ('Distance matrix for all runs generated!')



class plot (object):
    '''Turns alignment objects into plots '''
    def __init__(self, pair, filename):
        self.pair = pair
        self.filename = filename
       
    def heatplot(self):
        #Creates a heatmap for SUBJECT:DISTANCE pairs.
        #Does not generate axis label. Was not specified in directions, so I didn't do it.
        #Adding axis labels made my heatmaps weird. 
        width = [1]
        vals = np.array(self.pair.values())[:,None]
        axi = plt.matshow(vals, interpolation='none', cmap='hot')
        plt.colorbar(axi, shrink=0.8) # shrink size of colorbar #plt.colorbar(axi)
        plt.title(self.filename, y=1.2) # raise title higher #plt.title('Example of heatmap with colorbar')
        plt.savefig(self.filename + '.png', bbox_inches='tight')

    def lineplot(self):
        #Creates lineplot but only for 1 thing so far. Doesn't save it. 
        x = []
        y1 = []
        
        for c in filter1:
            bit = c.bitscore
            y1.append(bit)
        x = len(filter1)
        plot (x,y1)
        









#Creates a heatmap for each BLAST run, creates a lineplot 
plot(SDPAIR1, 'Heatmap1').heatplot()
plot(SDPAIR2, 'Heatmap2').heatplot()
plot(SDPAIR3, 'Heatmap3').heatplot()

plot(filter1, 'lineplot').lineplot()
print ('Check working directory for plots')

