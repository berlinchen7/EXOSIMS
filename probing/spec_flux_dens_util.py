from astropy.io import fits
import astropy.constants as const
import astropy.units as u
import numpy as np
import re, os, pickle, pkg_resources


def starMag(Vmag, BV, lam):
    lam_um = lam.to('um').value
    if lam_um < .550:
        b = 2.20
    else:
        b = 1.54
    mV = Vmag + b*BV*(1./lam_um - 1.818)
    return mV

def calc_F0_Traubetal(Vmag, BV, lam):
    return (1e4*10**(4.01 - (lam/u.nm - 550)/770)*u.ph/u.s/u.m**2/u.nm) *\
            10**(-0.4 * starMag(Vmag, BV, lam))

def calc_F0_Riemann(sdat, lam, BW):
    # Reimann integration of spectrum within bandwidth, converted from
    # erg/s/cm**2/angstrom to ph/s/m**2/nm, where dlam in nm is the
    # variable of integration.
    lmin = lam*(1-BW/2)
    lmax = lam*(1+BW/2)

    #midpoint Reimann sum
    band = (sdat.WAVELENGTH >= lmin.to(u.Angstrom).value) & (sdat.WAVELENGTH <= lmax.to(u.Angstrom).value)
    ls = sdat.WAVELENGTH[band]*u.Angstrom
    Fs = (sdat.FLUX[band]*u.erg/u.s/u.cm**2/u.AA)*(ls/const.h/const.c)        
    F0 = (np.sum((Fs[1:]+Fs[:-1])*np.diff(ls)/2.)/(lmax-lmin)*u.ph).to(u.ph/u.s/u.m**2/u.nm)
    return F0
  
def scale_spec_by_Traubetal(sdat, Vmag, BV):
    '''Crude scaling of given spectrum to fit Traub et al., which claims to be fairly accurate
        for 400 nm < \lambda < 1000 nm.
    '''
    i_450 = np.argmax(sdat.WAVELENGTH > (450*u.nm).to(u.Angstrom).value)
    i_700 = np.argmax(sdat.WAVELENGTH > (700*u.nm).to(u.Angstrom).value)
    i_900 = np.argmax(sdat.WAVELENGTH > (900*u.nm).to(u.Angstrom).value)
    
    bw_450, bw_700, bw_900 = 0.2, 0.15, 0.1 # \Delta\lambda \approx 100 nm

    # Compute correction factors:
    C_450 = calc_F0_Traubetal(Vmag, BV, 450*u.nm) / calc_F0_Riemann(sdat, 450*u.nm, bw_450)
    C_700 = calc_F0_Traubetal(Vmag, BV, 700*u.nm) / calc_F0_Riemann(sdat, 700*u.nm, bw_450)
    C_900 = calc_F0_Traubetal(Vmag, BV, 900*u.nm) / calc_F0_Riemann(sdat, 900*u.nm, bw_450)
    
    # Dilate FLUX by the average of the 3 factors:
    C = np.mean([C_450, C_700, C_900])
    sdat.FLUX = sdat.FLUX * C
    
    return sdat
    
def calc_F0(lam, BW, Vmag, BV, spec = None):

    indexf =  pkg_resources.resource_filename('EXOSIMS.TargetList','pickles_index.pkl')
    assert os.path.exists(indexf), "Pickles catalog index file not found in TargetList directory."
    datapath = pkg_resources.resource_filename('EXOSIMS.TargetList','dat_uvk')
    assert os.path.isdir(datapath), 'Could not locate %s in TargetList directory.' %(datapath)

    # grab Pickles Atlas index
    with open(indexf, 'rb') as handle:
        specindex = pickle.load(handle)

    speclist = sorted(specindex.keys())
    specdatapath = datapath
    
    #spectral type decomposition
    #default string: Letter|number|roman numeral
    #number is either x, x.x, x/x
    #roman numeral is either 
    #either number of numeral can be wrapped in ()
    specregex1 = re.compile(r'([OBAFGKMLTY])\s*\(*(\d*\.\d+|\d+|\d+\/\d+)\)*\s*\(*([IV]+\/{0,1}[IV]*)')
    #next option is that you have something like 'G8/K0IV'
    specregex2 = re.compile(r'([OBAFGKMLTY])\s*(\d+)\/[OBAFGKMLTY]\s*\d+\s*\(*([IV]+\/{0,1}[IV]*)')
    #next down the list, just try to match leading vals and assume it's a dwarf
    specregex3 = re.compile(r'([OBAFGKMLTY])\s*(\d*\.\d+|\d+|\d+\/\d+)')
    #last resort is just match spec type
    specregex4 = re.compile(r'([OBAFGKMLTY])')

    romandict = {'I':1,'II':2,'III':3,'IV':4,'V':5}
    specdict = {'O':0,'B':1,'A':2,'F':3,'G':4,'K':5,'M':6}
    
    #everything in speclist is correct, so only need first regexp
    specliste = []
    for s in speclist:
        specliste.append(specregex1.match(s).groups())
    specliste = np.vstack(specliste)
    spectypenum = np.array([specdict[l] for l in specliste[:,0]])*10+ np.array(specliste[:,1]).astype(float) 


    if spec is not None:
        # Try to decmompose the input spectral type
        tmp = specregex1.match(spec)
        if not(tmp):
            tmp = specregex2.match(spec)
        if tmp:
            spece = [tmp.groups()[0], \
                    float(tmp.groups()[1].split('/')[0]), \
                    tmp.groups()[2].split('/')[0]]
        else:
            tmp = specregex3.match(spec) 
            if tmp:
                spece = [tmp.groups()[0], \
                         float(tmp.groups()[1].split('/')[0]),\
                         'V']
            else:
                tmp = specregex4.match(spec) 
                if tmp:
                    spece = [tmp.groups()[0], 0, 'V']
                else:
                    spece = None

        #now match to the atlas
        if spece is not None:
            lumclass = specliste[:,2] == spece[2]
            ind = np.argmin( np.abs(spectypenum[lumclass] - \
                                    (specdict[spece[0]]*10+spece[1]) ))
            specmatch = ''.join(specliste[lumclass][ind])
        else:
            specmatch = None
    else:
        specmatch = None

    if specmatch == None:
        F0 = calc_F0_Traubetal(Vmag, BV, lam)
    else:
        # Open corresponding spectrum
        with fits.open(os.path.join(specdatapath, specindex[specmatch])) as hdulist:
            sdat = hdulist[1].data

        sdat = scale_spec_by_Traubetal(sdat, Vmag, BV)          
        F0 = calc_F0_Riemann(sdat, lam, BW)
    
    return F0
