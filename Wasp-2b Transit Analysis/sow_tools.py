# -*- coding: utf-8 -*-

import numpy as np

def mask_bad_pix(im, inst="NickelDIC", Nsig=7, thr_min=0.25):
    """
    Remove bad pixels by comparing their value to mean of neighbor pixels.
    im= image array to mark bad pixels in.
    inst= str with name of instrument im came from; only needed if instrument
            has known bad pixel locations.
    Nsig= # of standard deviations above the mean required to mark pixel as bad.
    thr_min= hard minimum value above the mean required to mark pixel as bad.
    """
    
    print ("\nMarking bad pixels (may take several minutes)...")
    
    # Create zero arrays same size as im but padded by 2 pixels on each side.
    mat = np.zeros((im.shape[0]+4, im.shape[1]+4))
    std_mat = np.zeros(mat.shape)
    # Shift and add im around edge of an 8-pixel square, with starting point
    # (0,0) Each pixel in mat is
    # the sum of the 8 nearest neighbors of the corresponding pixel in im.
    for yx in [(2,2), (1,2), (0,2), (0,1), (0,0), (1,0), (2,0), (2,1)]:
        mat[yx[0]:yx[0]+im.shape[0], yx[1]:yx[1]+im.shape[1]] += im.copy()
    # Calculate the mean of nearest neighbors.
    nebMean_mat = mat/8.
    # Shift and add squared difference from the mean of nearest neighbors.
    # Origin same as above. Needed to compute standard deviation of nearest
    # neighbors for each pixel in im.
    for yx in [(2,2), (1,2), (0,2), (0,1), (0,0), (1,0), (2,0), (2,1)]:
        std_mat[yx[0]:yx[0]+im.shape[0], yx[1]:yx[1]+im.shape[1]] += (im - nebMean_mat.copy()[yx[0]:yx[0]+im.shape[0], yx[1]:yx[1]+im.shape[1]])**2
    
    # Trim padding off of nebMean_mat so it matches shape of im.
    nebMeanArr = mat[1:im.shape[0]+1, 1:im.shape[1]+1]/8.
    # Finish calculating standard deviation of nearest neighbors and trim array.
    nebStdArr = np.sqrt(std_mat/8.)[1:im.shape[0]+1, 1:im.shape[1]+1]
    # Set threshold requirement for a "bad" pixel as Nsig standard deviations
    # above the mean or thr_min, whichever is smaller.
    thresh = Nsig*nebStdArr.copy()
    #thresh[thresh > thr_min] = thr_min
    # Find all pixels in im that exceed thresh relative to neighbors.
    diff = np.abs(nebMeanArr - im) - thresh
    whb = np.where(diff > 0.)
    # Replace all "bad" pixels with NaN.
    im[whb] = np.nan
    
    # Mask known bad pixels according to instrument specified.
    if inst=="NickelDIC":
        im[:, (256, 783, 784)] = np.nan
    
    print ("%d bad pixels (%.1f%%) marked as NaN (not counting known bad pixels)" % (whb[0].shape[0], 100*whb[0].shape[0]/float(im.size)))
    
    return im