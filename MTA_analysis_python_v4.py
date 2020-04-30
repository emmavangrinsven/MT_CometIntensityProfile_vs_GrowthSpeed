# -*- coding: utf-8 -*-
"""

@authors: Emma van Grinsven, Eugene Katrukha
"""
#MTA-analysis python v3

import numpy 
import itertools
import math
import matplotlib.pyplot as plt

def findmaxmin(xyArr):
    p = numpy.polyfit(xyArr[:,0], xyArr[:,1], 1)
    arrlength = numpy.size(xyArr[:,1])
    linearfit = numpy.polyval(p, xyArr[:,0])
    diff = xyArr[:, 1] - linearfit
    #find min and max
    indmin = numpy.argmin(diff)
    indmax = numpy.argmax(diff)
    epoch = []
    #check if min and max is not a borderpoint
    if indmin > 1 and indmin < arrlength-2:
        epoch = [indmin]
    if indmax > 1 and indmax < arrlength-2:
        epoch_temp = [indmax]
        epoch = epoch + epoch_temp
    
    #special case
    if len(epoch) == 0:
        #let's approximate line with points at the ends
        #line parameters
        slope = (xyArr[arrlength-1,1]-xyArr[0,1])/(xyArr[arrlength-1,0]-xyArr[0,0])
        b = xyArr[0,1]-slope*xyArr[0,0]
        p= [slope, b]
        #repeat
        linearfit = numpy.polyval(p, xyArr[:,0])
        diff = xyArr[:, 1] - linearfit
        #find min and max
        indmin = numpy.argmin(diff)
        indmax = numpy.argmax(diff)
        epoch = []
        #check if min and max is not a borderpoint
        if indmin > 1 and indmin < arrlength-2:
            epoch = [indmin]
        if indmax > 1 and indmax < arrlength-2:
            #epoch = numpy.concatenate((epoch, indmax), axis=0)
            epoch_temp = [indmax]
            epoch = epoch + epoch_temp
        #okay, very special case. reached the line
        #gonna return nothing!
        #but display warning
        #if numpy.shape(epoch) ==0:
            #print ('Error! Full line approximation is reached')
    return (epoch)
         
   
def findepochs(xyArr, max_num_of_intervals):
     arrlength = numpy.size(xyArr[:, 0])
     epoches = [0, arrlength-1]
     nCount = 0
     #find splitting
     while numpy.size(epoches) < max_num_of_intervals+1:
         currEpochesNum = numpy.size(epoches)
         newepoches = []
         i=0
         while i<=(currEpochesNum - 2):
             newepoches = findmaxmin(xyArr[(epoches[i]):(epoches[i+1]+1),:])
             if i > 0:
                 newepoches = newepoches + (epoches[i])
             if numpy.size(newepoches)>0:
                epoches = numpy.concatenate([epoches,newepoches])  
             i=i+1
         epoches = sorted(epoches)
         nCount = nCount+1
         if nCount > 10:
             #print('reached')
             break
     return (epoches)
           

def get_norm_error(xyArr, epoches):
    num_epoches = len(epoches)
    p = numpy.polyfit(xyArr[:,0], xyArr[:,1], 1)
    linearfit = numpy.polyval(p, xyArr[:,0])
    error_zero = math.sqrt(sum((xyArr[:,1]-linearfit)**2))
    
    error_fit = 0
    for i in range(0, num_epoches-1):
        curr_interval = numpy.arange(epoches[i], epoches[i+1]+1)
        p = numpy.polyfit(xyArr[curr_interval,0],xyArr[curr_interval,1],1)
        linearfit = numpy.polyval(p,xyArr[curr_interval,0])
        error_fit = error_fit + sum((xyArr[curr_interval,1]-linearfit)**2)
    error_fit = math.sqrt(error_fit)
    error_fit = -math.log(error_fit/error_zero)/(num_epoches-2)
    return (error_fit)

def get_rms(xyArr, epoches):
    num_epoches = len(epoches)
    
    error_fit = 0
    for i in range(0, num_epoches-1):
        curr_interval = numpy.arange(epoches[i], epoches[i+1]+1)
        if curr_interval[-1] > xyArr[-1,0]:
            curr_interval[-1] = curr_interval[-1]-1
        p = numpy.polyfit(xyArr[curr_interval,0], xyArr[curr_interval,1], 1)
        linearfit = numpy.polyval(p, xyArr[curr_interval,0])
        error_fit = error_fit + sum((xyArr[curr_interval,1]-linearfit)**2)
    
    error_fit = math.sqrt(error_fit)
    return (error_fit)

def findapproximation(xyArr, max_num_of_intervals):
    epochs = findepochs(xyArr, max_num_of_intervals)
    numepoch = len(epochs)-2
    if numepoch>0:
        # get all possible combinations
        combin_list = [epochs]    
        for i in range(1, numepoch+1):
            combin_list = combin_list + [[epochs[0], epochs[i], epochs[numepoch+1]]]
        nCount = numepoch + 1
        for i in range(2, numepoch):
            combis = list(itertools.combinations(epochs[1:numepoch+1], i))
            sizer = len(combis)
            for j in range(0, sizer):
                to_add =  [[epochs[0]], list(combis[j]), [epochs[numepoch+1]]]
                to_add = sum(to_add, [])
                combin_list = combin_list + [to_add]
                nCount = nCount + 1
        numepoch = len(combin_list)
        error_list = numpy.zeros((numepoch, 1))
        for i in range(0, numepoch):
            error_list[i] = get_norm_error(xyArr, combin_list[i])
        indMax = numpy.argmax(error_list)
        finalepoches = combin_list[indMax]        
    else:
        finalepoches=epochs
    rms_error = get_rms(xyArr, finalepoches)        
    return ([finalepoches,rms_error])

def getapproximation(xyArr, epoches):   
    arrlength = numpy.size(xyArr[:,0])
    currEpochesNum = len(epoches)
    xyApprox = numpy.zeros((arrlength+currEpochesNum-2, 2))
    slopes = numpy.zeros((currEpochesNum-1,1))
    b_coeff = numpy.zeros((currEpochesNum-1,1))
    
    for i in range(0, currEpochesNum-1):
        curr_interval = range(epoches[i],(epoches[i+1]+1))
        p = numpy.polyfit(xyArr[curr_interval,0],xyArr[curr_interval,1],1)
        slopes[i]=p[0]
        b_coeff[i]=p[1]
        linearfit = numpy.polyval(p,xyArr[curr_interval,0])
        curr_approx_interval = range((epoches[i]+i),(epoches[i+1]+i+1))
        xyApprox[curr_approx_interval,0]=xyArr[curr_interval,0]
        xyApprox[curr_approx_interval,1]=linearfit       
        
    return (xyApprox,slopes,b_coeff)
                                                       
def mta_analysis(xyseries, rel_rms_improvement = 0.05, max_depth_of_tree = 10, max_num_of_intervals = 5):
    xyArr = xyseries
    if (numpy.size(xyArr,1) > numpy.size(xyArr,0)):
        xyArr.transpose()
    if (numpy.size(xyArr,1) != 2):
        print('Two dimentional array is required!')
        return(0,0,0)
   
    approx_tree = [None]*max_depth_of_tree
    rms_tree = numpy.zeros((max_depth_of_tree,1))
    segments_number = numpy.zeros((max_depth_of_tree,1))
    arrlength = numpy.size(xyArr[:, 0])
    epoches = [0, arrlength-1]
    approx_tree[0] = epoches
    #number of segments
    segments_number[0] = 1;
    rms_tree[0] = get_rms(xyArr,epoches)
    
    #diving down
    #building a tree of more detailed splitting of original data
    for i in range(1,max_depth_of_tree):
        #previous approximation
        previous_level_epoch = approx_tree[i-1]
        currEpochesNum = numpy.size(previous_level_epoch)
        rms_segments = numpy.zeros((currEpochesNum-1,1))
        new_epoches= [None]*(currEpochesNum-1)

        #go and split every new segment     
        j=0;
        while j<=(currEpochesNum-2):
            if (previous_level_epoch[j+1] - previous_level_epoch[j]) > 2:
                epoches, rms_segments[j] = findapproximation(xyArr[previous_level_epoch[j]:(previous_level_epoch[j+1]+1),:], max_num_of_intervals)
                epoches = [x+previous_level_epoch[j] for x in epoches]
                new_num = numpy.size(epoches)    
                epoches = sorted(epoches)
                epoches= numpy.concatenate((previous_level_epoch,epoches[1:(new_num-1)]));
                epoches = sorted(epoches)
                epoches=[int(el) for el in epoches];
                rms_segments[j] = get_rms(xyArr,epoches)
                new_epoches[j]= epoches
            #cannot split further, no improvement
            else:
                #print ('else')
                new_epoches[j] = previous_level_epoch
                rms_segments[j] = rms_tree[i-1]
            j=j+1
        indMin = numpy.argmin(rms_segments)
        approx_tree[i] = new_epoches[indMin] 
        #print(approx_tree[i]);
        segments_number[i]= numpy.size(new_epoches[indMin])-1  
        rms_tree[i] = rms_segments[indMin]
    #print('Tree building done')
    #now errors are in rms_tree
    #and multi-scale approximations are in approx_tree
    #let's find optimal approximation
    #relative improvement
    rmsini=rms_tree[0];
    rms_plot = [x/rmsini for x in rms_tree]
    opt_level = -1;
    i=0;
    while opt_level<0:
        if(rms_plot[i]-rms_plot[i+1]<rel_rms_improvement):
           opt_level = i;         
        i = i+1;
    optimal_epoches = approx_tree[opt_level];
    xyApprox,slopes, b_coeff = getapproximation(xyArr, optimal_epoches);
    #print('all done')
    return (optimal_epoches,slopes,xyApprox)

#Main BODY
#import data
xyArr = numpy.loadtxt("xyArr_file").reshape(367,2)
#run analysis
optimal_epoches,slopes,xyApprox = mta_analysis(xyArr, rel_rms_improvement = 0.05, max_depth_of_tree = 10, max_num_of_intervals = 5)
#plot input data
plt.plot(xyArr[:,0],xyArr[:,1])
#plot approximation
nIntervals = len(slopes)
for i in range(0, nIntervals):    
    rangeApprox = numpy.arange((optimal_epoches[i]+i),(optimal_epoches[i+1]+i+1))
    plt.plot(xyApprox[rangeApprox,0],xyApprox[rangeApprox,1])
plt.legend(['input data', 'approximation'],loc='upper right')
plt.show()


