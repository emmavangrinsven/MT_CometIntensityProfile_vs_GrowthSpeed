import PIL
import numpy
import pandas
import os
import math
import matplotlib.pyplot as plt
from read_roi import read_roi_file
import scipy.ndimage
from scipy.optimize import curve_fit
from scipy.optimize import OptimizeResult
from scipy import special
from MTA_analysis_python_v4 import *
import statistics 


# Select directory with files to analyse. 
folder_name = "D:/Chao_tip_analysis/code/analyzed/20191112 EB123KO Chao construct 1 EB3-GFP/01/kymos_analyze"
# Note: No nFilePairs, manually supply function with path_folder_name.

def load_all_filenames(folder_name):
	tif_files = []
	for file in os.listdir(folder_name):
		if file.endswith(".tif"):
			tif_files = tif_files + [file]

	nTotFiles = len(tif_files)
	marked_files = numpy.zeros(nTotFiles, dtype=int)

	ch1_filenames = []
	for nFile in range(0, (len(tif_files))):
		fn1 = tif_files[nFile]
		fl1 = len(fn1)
		fl1 = fl1-4
		ch1_filenames = ch1_filenames + [fn1[:fl1]]
	return (ch1_filenames)

# matlabs smooth() function
def smooth(a,WSZ):
    # a: NumPy 1-D array containing the data to be smoothed
    # WSZ: smoothing window size needs, which must be odd number,
    # as in the original MATLAB implementation
    out0 = numpy.convolve(a,numpy.ones(WSZ,dtype=int),'valid')/WSZ
    r = numpy.arange(1,WSZ-1,2)
    start = numpy.cumsum(a[:WSZ-1])[::2]/r
    stop = (numpy.cumsum(a[:-WSZ:-1])[::2]/r)[::-1]
    return numpy.concatenate((  start , out0, stop  ))

# fit_kymograph_EB_conv_with_lattice_v3

# parameters
halfwidthright = 15
halfwidthleft = 35
bMakeMovie = True

nUmPerPixel = 0.065
nMinPerFrame = 0.1/60

# filebase = "SLAIN2 100NM A_K002"
ch1_filenames = load_all_filenames(folder_name)

# going through all profiles
nTotFiles = len(ch1_filenames)
nSumOfFit = []
SumOfProfiles = []
nProfilesTotal = -1
for nFile in range(0,nTotFiles):
    filebase = folder_name + "/" + ch1_filenames[nFile]
    sWork = "working on (" + str(nFile+1) + "/" + str(nTotFiles) + ") " + ch1_filenames[nFile]
    print(sWork)
    # ADD WAITBAR
    
    # read tiff file
    filenameG = filebase + ".tif"
    
    I_G = plt.imread(filenameG)
    
    szG = I_G.shape
    
    imwidth = I_G.shape[1]
    imheight = I_G.shape[0]
    
    #read Kymo_ROI
    sROI = read_roi_file(filebase + ".roi")
    if sROI[ch1_filenames[nFile]]["type"] == "line":
        xcoord = list()
        ycoord = list()
        xcoord = int(sROI[ch1_filenames[nFile]]['x1']), (int(sROI[ch1_filenames[nFile]]['x2']))
        ycoord = int(sROI[ch1_filenames[nFile]]['y1']), (int(sROI[ch1_filenames[nFile]]['y2']))
        sROI[ch1_filenames[nFile]]["x"] = xcoord
        sROI[ch1_filenames[nFile]]["y"] = ycoord
        del sROI[ch1_filenames[nFile]]['x1']
        del sROI[ch1_filenames[nFile]]['x2']
        del sROI[ch1_filenames[nFile]]['y1']
        del sROI[ch1_filenames[nFile]]['y2']
        
    nPoints=len(sROI[ch1_filenames[nFile]]['x'])
    
    nTotProfiles = sROI[ch1_filenames[nFile]]['y'][-1] - sROI[ch1_filenames[nFile]]['y'][0]
    
    print (sROI[ch1_filenames[nFile]]['y'][-1], sROI[ch1_filenames[nFile]]['y'][0])
    
    #fitting functions: exponential convolved with gaussian of specified width CHECK IF OK!!
    def erfampG(x, bkg, amp1, amp2, mu, lambdA): 
        return (bkg+0.5*amp1*numpy.exp((1.72*1.72*lambdA*lambdA*0.5)+lambdA*x-lambdA*mu)*(1-special.erf((lambdA*1.72*1.72+x-mu)/(math.sqrt(2)*1.72)))+0.5*amp2*(special.erfc((x-mu)/(1.72*math.sqrt(2)))))
    
       
    # ADD MOVIE
    # clear 'MovieProfG'
    # MovieProfG(nTotProfiles) = struct('cdata',[],'colormap',[]);
    
    #arrays for storage of fitting results
    nFitParamNG = 5
    
    expfitresultG = numpy.zeros((nTotProfiles, nFitParamNG+4))
    xinterval=numpy.zeros((nTotProfiles,2))
    iterNG=numpy.zeros((nTotProfiles,1))
    fitoutG=numpy.empty((nTotProfiles,1), dtype=object)
    
    i=0
    print (i)
    for nPoint in range(1,(nPoints)): 
        #along comet
        x1=sROI[ch1_filenames[nFile]]['x'][nPoint-1]
        y1=sROI[ch1_filenames[nFile]]['y'][nPoint-1]
        x2=sROI[ch1_filenames[nFile]]['x'][nPoint]
        y2=sROI[ch1_filenames[nFile]]['y'][nPoint]
        
        invslope=(x2-x1)/(y2-y1)
        
        #figure
        
        # ADD WAITBAR
        

        
        for ycurr in range((y1+1),y2+1):
            xbeg=int(invslope*(ycurr-y1)+x1-halfwidthleft)
            xend=int(invslope*(ycurr-y1)+x1+halfwidthright)
            if xend > imwidth:
                xend = imwidth
            if xbeg < 0:
                xbeg = 0
        
            #store fitting range
            xinterval[i][0]=int(xbeg)
            xinterval[i][1]=int(xend)
            ## Green Channel Fitting
            fitx = range(xbeg,(xend+1))    
            fity = I_G[ycurr, xbeg:(xend+1)]
            miny=numpy.amin(fity)
            maxy=numpy.amax(fity)
            maxind = numpy.argmax(fity)
            maxamprange = (maxy) - (miny)
            #taking care of initial conditions
            ini = numpy.zeros((nFitParamNG,1))
            ini[0] = miny
            ini[1] = 2*maxamprange
            ini[2] = 0.2*maxamprange           
            ini[3] = fitx[maxind]
            ini[4] = 1.0
            
            #debug stuff
            if ycurr == 280:
                ycurr = ycurr+1
                ycurr = ycurr-1
        
        
            fitxflat = numpy.asarray(fitx).ravel()
            fityflat = numpy.asarray(fity).ravel()
            
            # ADD BOUNDARIES
            
            bFitSuccess = True
            
            try:
                popt, pcov = curve_fit(erfampG, fitxflat, fityflat, p0=ini, maxfev = 1000)
                message = 'Success, but fitting stopped because change in residuals less than tolerance (TolFun).'
            except RuntimeError:
                bFitSuccess = False
                print ('NaN warning')
                message = Exception
            
                                  
         
          #   iterNG[i]=output.iterations;
            fitoutG[i,0]= message
            if bFitSuccess:
                coeff = popt
                expfitresultG[i,:] = [coeff[0],coeff[1],coeff[2],coeff[3],coeff[4],ycurr,0,0,nFile]
                fittedy = erfampG(fitx, coeff[0],coeff[1],coeff[2],coeff[3],coeff[4])
                #store profile
                fittedMax=numpy.amax(fittedy)
                #normalize intensity
                normfity = numpy.zeros((fity.shape))
                for m in range(fity.shape[0]):
                    normfity[m] = ((fity[m]-coeff[0])/(fittedMax-coeff[0]))
                shiftedfitx = numpy.zeros((len(fitx)))
                for l in range(len(fitx)):
                    shiftedfitx[l] = fitx[l]-coeff[3]
                nProfilesTotal = nProfilesTotal + 1
                SumOfProfiles= SumOfProfiles + [[shiftedfitx,normfity]]
              
            else:
                coeff = ini
                #expfitresultG[i,:] = ["NaN","NaN","NaN","NaN","NaN",ycurr,0,0,nFile]
                expfitresultG[i,:] = [coeff[0],coeff[1],coeff[2],coeff[3],coeff[4],ycurr,0,0,nFile]
                
                
            if bMakeMovie:                
                plt.plot(fitx,fity)
                plt.plot(fitx,fittedy,'--')
                plt.show()
                plt.close()
        
        
            i = i + 1
   
    speedgrowth= numpy.diff(expfitresultG, axis = 0)
    speedgrowth = speedgrowth[:,3]
    for p in range(len(speedgrowth)):
        expfitresultG[p, nFitParamNG+2] = speedgrowth[p]
        expfitresultG[p, nFitParamNG+1] = nUmPerPixel*(1./expfitresultG[p,4])
    
    #outliers removal
    cometposition = expfitresultG[:,3]
    cometsliding = smooth(cometposition,5)
    stdval = numpy.median(numpy.absolute(cometposition-cometsliding))
    diff_x = numpy.absolute(cometposition-cometsliding)

    #MTA analysis
    xvals=expfitresultG[:,5]*nMinPerFrame
    yvals=smooth(expfitresultG[:,3],11)*nUmPerPixel
    xyArr= numpy.column_stack((xvals,yvals))
    optimal_epoches,slopes,xyApprox = mta_analysis(xyArr, rel_rms_improvement = 0.05, max_depth_of_tree = 10, max_num_of_intervals = 5)
    
    yvals = expfitresultG[:,3]*nUmPerPixel

    #fill average speed values
    MTASpeed = numpy.zeros((len(xvals),1));
    nEpN = len(slopes)
    for ind,val in enumerate(slopes):
        MTASpeed[optimal_epoches[ind]:optimal_epoches[ind+1],0]= slopes[ind,0]
    
    #average slope number in the transition (epoch) points
    if nEpN>2:
        for ind in range(1,nEpN):
            MTASpeed[optimal_epoches[ind],0] = statistics.mean([slopes[ind-1,0], slopes[ind,0]])
     
    
    #plot input data
    plt.plot(xyArr[:,0],xyArr[:,1])
    #plot approximation
    nIntervals = len(slopes)
    for i in range(0, nIntervals):   
        rangeApprox = numpy.arange((optimal_epoches[i]+i),(optimal_epoches[i+1]+i+1))
        plt.plot(xyApprox[rangeApprox,0],xyApprox[rangeApprox,1])
    plt.legend(['input data', 'approximation'],loc='upper right')
    plt.title(ch1_filenames[nFile])
    plt.show()
    expfitresultG[:,nFitParamNG+2]=MTASpeed[:,0]
    if numpy.size(nSumOfFit) == 0:
        nSumOfFit = expfitresultG
    else:
        nSumOfFit = numpy.vstack((nSumOfFit,expfitresultG))

nComLength_Speed = numpy.column_stack((nSumOfFit[:,6],nSumOfFit[:,7]))
nComLengthMSDN = numpy.column_stack((statistics.mean(nComLength_Speed[:,0]), numpy.std(nComLength_Speed[:,1]), numpy.size(nComLength_Speed[:,1])))
nSpeed = numpy.zeros((nFile, 2))
nCom = numpy.zeros((nFile, 2))
for i in range(0, nFile):
    filt = nSumOfFit[:, 8] == i
    nCom[i,0] = statistics.mean(nComLength_Speed[filt, 0])
    nCom[i,1] = numpy.std(nComLength_Speed[filt,0]/math.sqrt(numpy.size(nComLength_Speed[filt, 0])))
    nSpeed[i,0] = statistics.mean(nComLength_Speed[filt, 1])
    nSpeed[i,1] = numpy.std(nComLength_Speed[filt,1]/math.sqrt(numpy.size(nComLength_Speed[filt, 1])))

    
