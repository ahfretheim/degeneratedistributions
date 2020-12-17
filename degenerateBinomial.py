# -*- coding: utf-8 -*-
"""
degenerateBinomial package

Executes the Degenerate Binomial methodology for studying political polarization

The Degenerate Binomial Methodology calculates the difference between the sample
of actual groups and the prediction of the binomial distribution, which represents
the frequency of all possible groups rather than self-selected political reality.

To achieve this, the following object model is used:
    
    predictionModel class - a class, generated from a csv file of polling inputs
    with group ID's (should be same format as the sample file included with this
    package), an optional tuple input of the true probabilities (say, from
    election results), and an optional yes text for greater adaptability. 
    Computes overall polarization and provides graphing methods, which are built 
    upon ggplot2 using Rpy2 (requires R and the ggplot2 package on local machine).
    Might also include information theory implementations. When used to analyze 
    an entire population, esp. one where the final population percentages are not 
    known, the optional input should be left blank, whereas when analyzing the 
    effect of an institution's political biases, you should use the true 
    probabilities input.
    
    group class - represents each group in the predictionModel. Provides analysis
    of the internal polarization and likely aggressive or alienating tendencies
    within a group with respect to politics
    
Note that all binomial distributions are calculated with respect to the yes vote.

NOTE: The Binomial describes all possible groups, so the degeneracy is clearly
created not by an excess of polarized groups, but by a lack of unpolarized
groups. Therefore the most sensible metric may be an estimate of how many
unpolarized groups should exist but don't:
    
    A[i] = Actual frequency at point i
    T[i] = Theoretical probability/likelihood at point i
    A[i]/T[i] = Degree of denominator inflamation at point i
    D[i,j] = Degeneracy at point j measured from inflamed point i = (A[i]/T[i])*T[j] - A[j]

NEXT STEPS:
    
    1. Add functionality for fleshing out missing groups (?)
    2. Add chi-squared test functionality.
    
    We'll handle graphing capabilities and external file creation LAST.

Created on Sat Nov 14 18:47:02 2020

@author: Alexander Fretheim
"""

import scipy.stats as ss;
import csv as c;
from math import floor;

class predictionModel:
    
    def __init__(self, csv, ptruth = -1.0, yes = "Yes", grp = "Group ID", vtext = "Vote", chicut = 5):
        #opening the csv file:
        file = open(csv, newline="", errors="ignore");
        data = c.DictReader(file);
        sampleyes = 0.0;
        samplesize = 0.0;
        
        #initializing object dictionaries:
        self.groups = {}; #all groups in sample, group objects.
        self.ilks = {}; #types of group, each a list of groups therein.
        self.modes = [];
        self.degeneracy = 0.0;
        self.vacated = [];
        self.deglist = {};
        chiobs = [];
        chiexp = [];
        
        #reading the csv:
        for row in data:
            gid = row[grp];
            if gid in self.groups:
                self.groups[gid].add(row);
            else:
                self.groups[gid] = group(row, yes, grp, vtext);
            if row[vtext] == yes:
                sampleyes += 1.0;
            samplesize += 1.0;
        
        #housekeeping:
        file.close();
        
        #Temporary testing code:
        print("CSV read complete. Groups all read and established. Building table.")
        print("Total sample size read: " + str(samplesize))
        print("Total yes votes read: " + str(sampleyes))
        
        #running our table builder
        self.table = self.buildECDFTable();
        
        #self.B = ss.binom(); #creating our binomial class
        
        #creating the binom field - the underlying distribution of the non-polarized ideal
        if ptruth == -1.0:
            self.p = sampleyes/samplesize
            self.populationset = True
        else:
            self.p = ptruth
            self.populationset = False
        
        self.samplesize = samplesize;
        
        #finding which parts of ECDF are over and under predicted, and categorizing accordingly:
        #for reference, keys structured: [#yes]-[#total]. Numbers are quantity - multiply probabilities by samplesize to match
        
        ratiolist = []; #list of degeneracy ratios to calculate supremum from.
        
        for key in self.table:
            print(key) #temporary testing code
            a = key.find("-");
            k = int(key[0:a])
            n = int(key[a+1:]);
            B = ss.binom(p = self.p, n=n);
            
            #sorts ECDF entries by whether they are above or below the predicted distribution, with below first:
            print(B.pmf(k))
            chiobs.append(self.table[key])
            chiexp.append(B.pmf(k)*self.samplesize)
            if B.pmf(k)*self.samplesize < self.table[key]:
                self.vacated.append(key)
            else:
                self.modes.append(key)
                ratiolist.append(self.table[key]/(B.pmf(k)*self.samplesize))
                
        
        M = self.calcSupremum(ratiolist);
        degeneracy = 0.0;
        

        #Calculate total set degeneracy here, from vacated list
        for v in self.vacated:
            a = v.find('-')
            k = float(v[0:a])
            n = float(v[a+1:]);
            B = ss.binom(p = self.p, n=n);
            d = B.pmf(k)*M - self.table[v]
            self.deglist[v] = d
            degeneracy += abs(d)
        
        self.degeneracy = degeneracy/M
        self.supremum = M
        
        #Running chi-squared test:
        i = 0;
        elim = 0;
        chiof = [];
        chief = [];
        for C in chiobs:
            if C > chicut:
                if chiexp[i] > chicut:
                    chiof.append(C)
                    chief.append(chiexp[i])
                else:
                    elim += 1;
            else:
                elim += 1;
            i += 1;
        self.chi = ss.chisquare(chiof, f_exp=chief)
        self.groupstested = len(self.groups)-elim
        
    #A method, run in the constructor, to build the basic empirical distribution:
    #NOTE: This will perform better if we return empirical quantities instead of probabilities:
    def buildECDFTable(self):
        table = {};
        i = 0.0
        for group in self.groups:
            #the basic format of a ecdf entry:
            E = str(floor(self.groups[group].yes)) + "-" + str(floor(self.groups[group].size))
            i += 1.0
            
            if E in table:
                table[E] += 1.0
            else:
                table[E] = 1.0
        
        #Temporary testing code:

        print("Constructed table:")
        print(table.items())
        print("Total tabulated groups: " + str(i))
        
        return table;
    
    #A method for calculating the supremum of the degeneracy ratios, takes a list:
    
    def calcSupremum(self, lst):
        lst.sort(reverse=True)
        M = [];
        ret = []; #not really what it returns, but the primary mechanism of return
        for i in range(len(lst)-1):
            M.append(lst[i+1] - lst[i])
        V = ss.zscore(M)
        i = 0

        #this for loop removes outliers using the z-score of the difference
        for v in V:
            if abs(v) <= 2.0: #2 standard deviation cut-off for differences in the sorted list
                ret.append(lst[i])
            i += 1;
        
        ret.append(lst[i]) #always append the last item in the sorted list
        
        return max(ret); #now that the list has been cleaned of outliers, the supremum can be returned

class group:
    def __init__(self, row, yes, grp, vtext):
        self.affirmative = yes;
        self.vote = vtext;
        if row[vtext] == self.affirmative:
            self.yes = 1.0;
        else:
            self.yes = 0.0;
        self.name = row[grp]
        self.size = 1.0;
    
    def add(self, row):
        if row[self.vote] == self.affirmative:
            self.yes += 1.0
        self.size += 1.0
    
    def getYes(self):
        return self.yes/self.size

    def getNo(self):
        return (self.size - self.yes)/self.size