

import numpy as np
import warnings
from scipy import optimize
from scipy import signal
from scipy import interpolate

def dextract(spectrum,traces,nspix,swext=2,npixavg=None):
    """Object spectral extraction routine
   
    Options:
    swext = 1: extract point source spectra by averaging over window
    swext = 2: extract full spatial orders with spatial interpolation
    swext = 3: extract full orders without interpolation 
    swext = 4: extract spectra by averaging over window

    Output:
    spectra: 2 or 3-d array of spectra of individual orders
    """
    
    # Case out the dimensionality of traces... 0 -> return
    if traces.ndim == 0:
        return 0
    elif traces.ndim == 1:
        nx = traces.size
        norders = 1
    else:
        norders, nx = traces.shape

    # Case out the option swext
    if swext == 0:
        return 0

    elif swext == 1:
        if npixavg is None:
            npixavg = 19
        spectra = np.empty((norders,nx), dtype=float)
        for io in range(norders):
            spectra[io,:] = dspecfit(spectrum, traces[io,:], bwidth=nspix, 
                                     extwidth=npixavg)
    
    elif swext == 2:
        spectra = np.empty((nspix, nx, norders), dtype=float)
        fnspix = np.arange(nspix, dtype=float) - nspix/2
        for io in range(norders):
            for ix in range(nx):
                # Interpolation:
                xt = traces[io,ix].astype(int) + fnspix
                ut = traces[io, ix] + fnspix
                vector = spectrum[xt, ix]
                ##IDL: return = INTERPOL(V, X, XOUT)
                #tvector = interpol(vector, xt, ut)
                f = interpolate.interp1d(xt, vector, bounds_error=False)
                tvector = f(ut)
                spectra[:,ix,io] = tvector

    elif swext == 3:
        spectra = np.empty((nspix, nx, norders), dtype=float)
        inspix = np.arange(nspix, dtype=int) - nspix/2
        for io in range(norders):
            for ix in range(ix):
                # W/O Interpolation:
                xt = traces[io,ix].astype(int) + inspix
                spectra[:,ix,io] = spectrum[xt, ix]

    elif swext == 4:
        if npixavg is None:
            npixavg = nspix
        spectra = np.empty((norders, nx), dtype=float)

        for io in range(norders):
            spectra[io,:] = specavg(spectrum, traces[io,:], npixavg)
    
    else:
        print("Silly user, you can't do that.")
        return 0

    return spectra



def gaussfit_func(x, a0, a1, a2, a3):
 
    # Silence RuntimeWarning for overflow, this function only
    warnings.simplefilter('ignore', RuntimeWarning)
    z = (x - a1) / a2
    y = a0 * np.exp(-z**2 / a2) + a3
    return y

def dflines(image, thresh=20., mark=False, title=''):
    """Automatically find and centroid lines in a 1-row image
 
    :image:
    :thresh: 20 DN above background
    :mark:
    :title:

    Returns:
    :centers: List of line centers (pixel #)
    :fwhm: The computed FWHM
    """

    warnings.simplefilter('ignore', optimize.OptimizeWarning)

    nx, ny = image.shape
    avgj = np.ndarray.flatten(image)

    print(f"Shapes of image: {image.shape}, and avgj: {avgj.shape}")

    peaks = []
    fwhm = []
 
    # Create background from median value of the image:
    bkgd = np.median(image)
    print(f'  Background level: {bkgd:.1f}')

    # Step through the cut and identify peaks:
    fmax = 11     # Minimum separations
    fwin = 15
    fhalfmax = int(np.floor(fmax/2))
    fhalfwin = int(np.floor(fwin/2))
    findmax = 50
    j0 = 0

    if mark:
        centers = [0]
    else:
        for j in range(ny):
            # print(f"avgj[j]: {avgj[j]}, {avgj[j].shape}")
            if j > (ny - fmax):
                continue
            if avgj[j] > (bkgd + thresh):
                j1 = j
                if np.abs(j1 - j0) < fmax:      
                    continue
                for jf in range(findmax):
                    itmp0 = avgj[jf + j]
                    itmp1 = avgj[jf + j + 1]
                    if itmp1 < itmp0:
                        icntr = jf + j
                        break
                if (icntr < fmax/2) or (icntr > (ny - fmax/2 - 1)):
                    continue
                xx = np.arange(fwin, dtype=float) + float(icntr - fhalfwin)
                temp = avgj[(icntr) - fhalfwin : (icntr) + fhalfwin + 1]
                temp = signal.medfilt(temp, kernel_size=3)
                
                # Run the fit, with error checking
                try:
                    p0 = [1000, np.mean(xx), 3, bkgd]
                    aa, cv = optimize.curve_fit(gaussfit_func, xx, temp, p0=p0)
                except RuntimeError:
                    continue  # Just skip this one
                # print(f"Gaussfit parameters: {aa}")
                tempfit = gaussfit_func(xx, *aa)
                center = aa[1]
                fw = aa[2] * 1.177 * 2.
                pmax = aa[0]
                if fw > 1.0:
                    peaks.append(center)
                    fwhm.append(fw)
                j0 = jf + j

        centers = np.asarray(peaks)
        #fwhm = fwhm[1:]
    
    cc = np.where(np.logical_and(centers >=0, centers <=2100))
    centers = centers[cc]

    szc = len(centers)
    print(f" Number of lines: {szc}")

    print(f"At this point the code makes some plots.  Yippee.")
 
    """ 
    plot,avgj,xra=[0,ny+2],xsty=1,yra=[0,max(avgj)+0.2*max(avgj)], $
    title=title, xtitle='CCD column', ytitle='I (DN)'
    for id=0,szc(1)-1 do begin
        plots,centers(id),avgj(centers(id))+20,psym=1
        xyouts,centers(id),avgj(centers(id))+30,strtrim(centers(id),2), $
        /data,orientation=0.
    endfor

    """
 
    return (centers, fwhm)

def specavg(spectrum, trace, wsize):
    """Extract an average spectrum along trace of size wsize
    :param spectrum: input spectrum
    :param trace: the trace along which to extract
    :param wsize: the size of the extraction (usually odd)
    :return:
    """
    # Case out the dimensionality of traces... 0 -> return
    if spectrum.ndim == 0:
        return 0
    elif spectrum.ndim == 1:
        nx = spectrum.size
    else:
        ny, nx = spectrum.shape
    speca = np.empty(nx, dtype=float)

    whalfsize = int(np.floor(wsize/2))

    # Because of python indexing, we need to "+1" the upper limit in order
    #   to get the full wsize elements for the average
    for i in range(nx):
        speca[i] = np.average(spectrum[int(trace[i]) - whalfsize : 
                                       int(trace[i]) + whalfsize + 1, i])
    
    #print(f"The shape of speca: {speca.shape}")
    return speca.reshape((1,nx))


def dspecfit(spec, trace, bwidth=101, extwidth=19):
    """Fit the background to a spectrum and subtract it

    Default values:
    :bwidth: 101
    :extwidth: 19
    """

    """     
    tval=fix(trace)

    nfit=bwidth & win=[bwidth/2-extwidth/2,bwidth/2+extwidth/2] ;  w=0 over win
    sz=size(spec)
    fits=fltarr(nfit,sz(1)) & datas=fits & subs=fits

    for i=0,sz(1)-1 do begin

        ww=fltarr(nfit) & ww(*)=1. & ww(win(0):win(1))=0.
        data=spec(i,tval(i)-nfit/2:tval(i)+nfit/2)
        coef=polyfitw(findgen(nfit),data,ww,1)	; go with linear fit

        fit=poly(findgen(nfit),coef)
        fits(*,i)=fit
        datas(*,i)=data
        subs(*,i)=data-fit

    endfor

    gplott=fltarr(sz(1))

    for i=0,sz(1)-1 do gplott(i)=total(subs(win(0):win(1),i))

    return,gplott

    end
    """
    return 0
