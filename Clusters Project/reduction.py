# Procedure to fix known bad columns in CCD2 images.  2016 Oct 2 E. Gates

from astropy.io import fits
import numpy as np
import sys, getopt,os
import glob

# Create list of flat fielded data
def remove_bad_pixels(datain, dataout):

    # _bp in output file name stands for bad pixel corrected
    #dataout = [i[:-5]+ '_bp.fits' for i in datain]

    n=len(datain)
    # size of box for area around bad pixel to be averaged
    s=2

    # read in one image to get image size for bad pixel mask
    data,header=fits.getdata(datain[0],header=True)

    # make bad pixel mask
    mask=np.ma.make_mask(data,copy=True,shrink=True,dtype=np.bool_)
    mask[:,:]=False
    mask[:,255:257]=True
    mask[:,783:785]=True
    mask[:,1001:1003]=True

    # loop for all the data bad pixel correction 
    for k in range(0,n):
        data,header=fits.getdata(datain[k],header=True)
        mdata=np.ma.masked_array(data,mask=mask,fill_value=np.nan)
        dataFixed=data.copy()
        for i in range(0,mdata.shape[0]):
            for j in range(0,mdata.shape[1]):
                if np.math.isnan(mdata[i,j]):
                    x1=i-s
                    x2=i+s+1
                    y1=j-s
                    y2=j+s+1
                    if x1<0:
                        x1=0
                    if x2>mdata.shape[0]:
                        x2=mdata.shape[0]
                    if y1<0:
                        y1=0
                    if y2>mdata.shape[1]:
                        y2=mdata.shape[1]
                    dataFixed[i,j]=np.mean(mdata[x1:x2,y1:y2])
        header['HISTORY']='Bad columns replaced'
        fits.writeto(dataout[k],dataFixed,header)

def overscan_subtraction(f):

    # set fit = 'yes' to do legendre fit to overscan regions, 'no' to just use the median
    fit = 'yes'   
    
    # create input and output file name lists.  Output files will have unique names, leaving the input files unchanged.
    ifilelist = glob.glob(f)
    # _os stands for overscan subtracted in the output file names
    ofilelist = [i[:-5]+ '_os.fits' for i in ifilelist]
      
    # how many files
    numifiles = len(ifilelist)
    numofiles = len(ofilelist)
    if numifiles != numofiles:
        sys.exit('Input and output file lists have different numbers of files. Exiting.')

    # For each file in ifilelist, read in file, figure out overscan and data regions, fit 
    # overscan with desired function (if any), and subtract from data.  
    # Write data to ofilelist value.  

    for i in range(0,numifiles):
        ifile=ifilelist[i]
        ofile=ofilelist[i]
        data, header = fits.getdata(ifile,header=True)

        # change data to float
        data=data.astype('float32')

        # read necessary keywords from fits header

        #number of pixels in image
        xsize = header['NAXIS1']
        ysize = header['NAXIS2']
        #start column and row
        xorig = header['CRVAL1U']
        yorig = header['CRVAL2U']
        #binning and direction of reading pixels
        cdelt1 = header['CDELT1U']
        cdelt2 = header['CDELT2U']
        #number of overscan columns
        rover = header['ROVER']
        cover = header['COVER']
        #unbinned detector size
        detxsize = header['DNAXIS1']  
        detysize = header['DNAXIS2']
        #number of amplifiers
        ampsx = header['AMPSCOL']
        ampsy = header['AMPSROW']

        # determine number and sizes of overscan and data regions
        namps = ampsx*ampsy
        if rover > 0:
            over=rover
            sys.exit('Program does not yet deal with row overscans. Exiting.')
        else:
            over = cover
        if over == 0:
            sys.exit('No overscan region specified in FITS header. Exiting.')

        # single amplifier mode (assumes overscan is the righmost columns)
        if namps == 1:
            biassec = data[:,xsize-cover:xsize]
            datasec = data[0:,0:xsize-cover]

            # median overscan section
            bias=np.median(biassec, axis=1) 

            # legendre fit
            if fit == 'yes':
                # fit
                lfit = np.polynomial.legendre.legfit(range(0,len(bias)),bias,3)
                bias = np.polynomial.legendre.legval(range(0,len(bias)),lfit)

            # subtract overscan
            datanew = datasec
            for i in range(datasec.shape[1]):
                datanew[:,i] = datasec[:,i]-bias

        # two amplifier mode (assumes both amplifer overscans are at rightmost columns)
        if namps == 2:
            biasseca = data[:,xsize-cover*2:xsize-cover]
            biassecb = data[:,xsize-cover:xsize]

            # median overscan sections
            biasa=np.median(biasseca,axis=1)
            biasb=np.median(biassecb,axis=1)

            # legendre fit
            if fit == 'yes':
                lfita = np.polynomial.legendre.legfit(range(0,len(biasa)),biasa,3)
                lfitb = np.polynomial.legendre.legfit(range(0,len(biasb)),biasb,3)
                biasa = np.polynomial.legendre.legval(range(0,len(biasa)),lfita)
                biasb = np.polynomial.legendre.legval(range(0,len(biasb)),lfitb)

            # Extract data regions

            # determine boundary between amplifiers
            bd=detxsize/2/abs(cdelt1)

            # calculate x origin of readout in binned units if cdelt1 negative or positive
            if cdelt1 < 0:
                #if no binning x0=xorig-xsize-2*cover, with binning:
                x0=xorig/abs(cdelt1)- (xsize-2*cover) 
            else:
                x0=xorig/cdelt1
                
            xtest=x0+xsize-cover*2 # need to test if all data on one or two amplifiers

            # determine which columns are on which amplifier and subtract proper overscan region

            if xtest < bd: # all data on left amplifier
                datanew=data[:,0:xsize-cover*2]
                m=datanew.shape[1]
                for i in range(0,m):
                    datanew[:,i]=datanew[:,i]-biasa

            if x0 >= bd: # all data on right amplifier
                datanew=data[:,0:xsize-cover*2]
                m=datanew.shape[1]
                for i in range(0,m):
                    datanew[:,i]=datanew[:,i]-biasb

            if xtest >= bd and x0 < bd:  #data on both amplifiers
                x1=int(bd-x0)
                dataa=data[:,0:x1]
                datab=data[:,x1:-cover*2]
                ma=dataa.shape[1]
                mb=datab.shape[1]
                for i in range(0,ma):
                    dataa[:,i]=dataa[:,i]-biasa
                for i in range(0,mb):
                    datab[:,i]=datab[:,i]-biasb
                # merge dataa and datab into single image
                datanew=np.hstack([dataa,datab])
                
        if namps > 2:
            sys.exit('Program does not yet deal with more than two overscan regions. Exiting.')

        # add info to header
        header['HISTORY'] = 'Overscan subtracted'

        # write new fits file
        fits.writeto(ofile,datanew,header,overwrite=True)