import numpy as np
from scipy.interpolate import interp1d
from scipy import constants

c_kms=constansts.c/10.**3
#convert from m/s to km/s

class spectrum:


    def __init__(self, flux, dispersion_axis=None, dispersion_axis_type=None,dispersion_axis_units=None,flux_units=None,vacuum=True):
        #if vacuum=False, then the dispersion_axis is after being defracted by air
        
        #type_dispersion_axis can be wavlength, wavenumber, or frequency
        #units is for wavelength and frequency type
        
        #dispers_axis_units can be 'A' or 'Angstroms' for Angstroms
        #'nm' or 'nanometers' for nm
        #'microns' for microns 
        #'mm' or 'millimeters' for millimeters
        #'Hz' for hertz
        #'1/cm' for wavenumber
        
        
        #flux_units can be something like erg/ s cm^2 nu or normalized or counts        


        self.flux = flux
        self.dispersion_axis=dispersion_axis

        self.dispersion_axis_type=dispersion_axis_type


        if isinstance(dispersion_axis_units,str):
            self.dispersion_axis_units=dispersion_axis_units

    def findbests2n(self,width=20,edge=10,p=100):
        """Finds the best signal to noise in a spectrum by dividing the mean of a region 
        by the standard deviation in that region

        not suitable for spectra with a continuum at 0 or no continuum

        p is percentile from 0 to 100
        for best use 100 or 99 or something like that
        100 would return the very best s2n but might be suseptible to outliers

        width is in pixels
        edge in pixels; edges can have odd edge effects, sometimes due to correcting for 
        the blaze function the best s2n will never be at ends"""

        l=len(spec)
        spec=self.flux

        s2nstemp=[]
        means=[]
        stds=[]
        i=edge

        while i<(l-width-edge):
            rangeend=i+width
            arr=spec[i:rangeend]
            avg=np.mean(arr)
            std=np.std(arr,ddof=1)
            s2n=np.mean(arr)/np.std(arr,ddof=1)
            s2nstemp.append(s2n)
            means.append(avg)
            stds.append(std)
            i=i+1  

        s2nstemp=np.nan_to_num(s2nstemp)
        return np.percentile(s2nstemp,p)

    def xcor(self, template, lb, ub, dispersion=0):
        #template must be an instance of spectrum
        #ub and lb are upper and lower bounds of pixel shift
        f=self.flux
        g=template.flux

        l=ub-lb
        arrl=len(f)
        #a1=np.zeros(arrl+l+1)
        #a2=a1
        shifts=np.arange(lb, ub+1)
        shifts=shifts.astype(np.int)

        corrs=[]
        for shift in shifts:
            if shift>0:
                a1=g[0:arrl-shift]
                a2=f[shift:]
            elif shift==0:
                a1=g
                a2=f
            else:
                a1=g[-shift:]
                a2=f[0:arrl+shift]

            corrs.append(np.corrcoef(a1,a2)[0,1])

        if dispersion!=0:
            rvs=np.arange(lb,ub+1,1)*dispersion
            corrs=corrs,rvs
        return corrs    


class spec_vs_wl(spectrum):
    
    dispersion_axis_type='wavelength'
    


    def rotbroad(self,vsini, eps=0.6,nr=10,ntheta=100, dif=0):
        """this don't work at the edges yet"""

        flux=self.flux
        wl=self.dispersion_axis

        # based on CMJ's program rotint.pro edited May 11 1994
        # ;  This routine reads in a spectrum, s, on a wavelength scale, w, and a vsini
        # ;  with which to rotationally broaden the spectrum.  The rotationally broadened
        # ;  spectrum is returned in ns.  Parameters that can be set are the coefficient
        # ;  of the limb darkening law, eps (0.6 default), the number of radial steps on
        # ;  the disk, nr (default = 10), and the maximum number of steps in angle around
        # ;  the disk, ntheta (default = 100).  Final optional parameter dif allows for 
        # ;  differential rotation according to the law Omeg(th)/Omeg(eq) = (1. - dif/2
        # ;  - (dif/2) cos(2 th)).  Dif = .675 nicely reproduces the law proposed by
        # ;  Smith, 1994, A&A, in press. to unify WTTS and CTTS.  Dif = .23 is similar to
        # ;  observed solar differential rotation.  Note: the th in the above expression
        # ;  is the stellar co-latitude, not the same as the integration variable used
        # ;  below.  This is a disk integration routine.

        #note from LF: using interpolate instead of spline

        final_spec=np.zeros(len(flux))

        ns=np.zeros(len(flux))
        tarea=0.0

        dr=1./nr

        for j in range(nr):
            r=dr/2.+j*dr
            area=((r+dr/2.)**2-(r-dr/2.)**2)/int(ntheta*r)*(1.-eps+eps*np.cos(np.arcsin(r)))
            for k in range(int(ntheta*r)-1):
                th=np.pi/int(ntheta*r)+k*2.*np.pi/int(ntheta*r)
                if dif!=0:
                    vl=vsini*r*np.sin(th)*(1.-dif/2.-dif/2.*np.cos(2.*np.arccos(r*np.cos(th))))
                    interped = interp1d(wl+wl*vl/(3*10**5), flux,fill_value=1,bounds_error=False)
                    ns=ns+area*interped(wl) #interpol(flux,wl+wl*vl/3.e5,wl) 
                    tarea=tarea+area
                else:
                    vl=r*vsini*np.sin(th)
                    interped = interp1d(wl+wl*vl/(3*10**5), flux,fill_value=1,bounds_error=False)
                    ns=ns+area*interped(wl)#     ns=ns+area*interpol(flux,wl+wl*vl/3.e5,wl) 
                    tarea=tarea+area

        ns=ns/tarea
        return spec_vs_wl(ns,wl)
    
    def doppler_shift(self,shift,return_wl=False):
        #return_wl=True returns only the new wavelength scale
        #shift in km/s
        #note, this is non-relativistic
        
        #if np.abs(shift)>0.1*c_kms:
            ##if shift is high enough to be relativistic
            #raise Exception("The absolute value of the Doppler shift should be less than 10% of 
                            #the speed of light since this program does need use a relativistic 
                            #treatement.")
        
        new_wl = self.dispersion_axis * (1.0 + shift / c_kms)
        if return_wl:
            return new_wl
        else:
            return spec_vs_wl(self.flux,new_wl)



class spec_vs_wn(spectrum):
    
    dispersion_axis_type='wavenumber'
    


    def wn2wl(self,units='microns',return_wl=False):
        #if return_wl is True then just the wl is returned
        #otherwise, a new object is returned with wl, flux

        #units can be 'A' or 'Angstroms' for Angstroms
        #'nm' or 'nanometers' for nm
        #'microns' for microns (defaut)
        #'mm' or 'millimeters' for millimeters
        #assume wave number is in 1/cm

        conversion_factor=10.0**4
        if units=='nm' or units=='nanometers':
            conversion_factor=10.0**7
        elif units=='Angstroms'  or units=='A':
            conversion_factor=10.0**8
        elif units=='mm' or units=='millimeters':
            conversion_factor=10.0

        wl=conversion_factor/self.dispersion_axis
        if return_wl:
            return wl
        else:
            return spec_vs_wl(self.flux,wl,dispersion_axis_units=units)

