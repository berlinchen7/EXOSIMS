import pickle, re, pkg_resources, os, csv
from astropy.io import fits
import astropy.units as u
import astropy.constants as const
import numpy as np
import pysynphot

#set up stuff for spectral type conversion
# Paths
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
rev_romandict = dict([reversed(i) for i in romandict.items()])
specdict = {'O':0,'B':1,'A':2,'F':3,'G':4,'K':5,'M':6}

#everything in speclist is correct, so only need first regexp
specliste = []
for spec in speclist:
    specliste.append(specregex1.match(spec).groups())
specliste = np.vstack(specliste)
spectypenum = np.array([specdict[l] for l in specliste[:,0]])*10+ np.array(specliste[:,1]).astype(float) 


#Set up stuffs for Kurutz-Cateslli atlas:
ck_indexf = '/Users/berlinc/Documents/exosims/EXOSIMS_fork/probing/spec_lib_summary/ck_spec_types.csv' 
with open(ck_indexf, 'r') as infile: 
    reader = csv.reader(infile) 
    reader.__next__() # Get rid of the first line which contains col labels 
    ck_specindex = dict((rows[0],rows[1]) for rows in reader) 
ck_speclist = sorted(ck_specindex.keys())
ck_specliste = []
for ck_spec in ck_speclist:
    ck_specliste.append(specregex1.match(ck_spec).groups())
ck_specliste = np.vstack(ck_specliste)
ck_spectypenum = np.array([specdict[l] for l in ck_specliste[:,0]])*10+ np.array(ck_specliste[:,1]).astype(float)

def get_ck_spec(BW=0.16, lam=551*u.nm, spec = None, teff = None, logg = None, m = None):
    """
    This function calculates the spectral flux density for a given 
    spectral type. Assumes the Pickles Atlas is saved to TargetList:
        ftp://ftp.stsci.edu/cdbs/grid/pickles/dat_uvk/
    If spectral type is provided, tries to match based on luminosity class,
    then spectral type. If no type, or not match, defaults to fit based on 
    Traub et al. 2016 (JATIS), which gives spectral flux density of
    ~9.5e7 [ph/s/m2/nm] @ 500nm
    
    Args:
        BW (float):
            Bandwidth fraction
        lam (astropy Quantity):
            Central wavelength in units of nm
        Spec (spectral type string):
            Should be something like G0V
            
    Returns:
        astropy Quantity:
            Spectral flux density in units of ph/m**2/s/nm.
    """
    
    # Begin by matching lam with nearest bandpass, then assign central
    # wavelength of ref band (UBVRIJK)
    bandeff = [365, 445, 551, 658, 806, 1220, 2190]*u.nm
    refLam = bandeff[(np.abs(bandeff-lam).argmin())]
    band_dict = {365: 'u', 445: 'b', 551: 'v', 658: 'r', 806: 'i',
                        1220: 'j', 1630: 'h', 2190: 'k', 3450: 'l',
                        4750: 'm'}
    refband = band_dict[refLam.value]

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
            lmax = lam*(1+BW/2)                                                    
            lmin = lam*(1-BW/2)
            # Otherwise, resort to the Castellie and Kurutz 2004 Atlas:
            if True:
                # Match by surface gravity, teff, and optionally metallicity if they exist:
                if logg is not None and teff is not None:
                    specmatch = True
                    # Default metallicity to 0 (i.e., same as the sun)
                    m = 0 if m is None else m
                    teff = teff.to(u.K).value
                # Otherwise match by spectral type:
                else:
                    lumclass = spece[2]
                    # CK's README only indicated suggested templates for odd lum classes
                    lumclass = romandict[lumclass] if romandict[lumclass]%2 == 1\
                                                   else romandict[lumclass] - 1
                    lumclass = rev_romandict[lumclass]
                    assert lumclass in ['I', 'III', 'V']
                    spece[2] = lumclass
                    lumclass_ind = ck_specliste[:,2] == spece[2]
                    ind = np.argmin( np.abs(ck_spectypenum[lumclass_ind] - (specdict[spece[0]]*10+spece[1]) ))
                    specmatch = ''.join(ck_specliste[lumclass_ind][ind])
                    spec_name = ck_specindex[specmatch]
                    # Parse the filename to get temp and logg:
                    parsed = re.split('\[|_|\]', spec_name)
                    teff = float(parsed[1])
                    logg = float(parsed[2][1:])/10
                    m = 0 # metallicity is always 0 in the suggested templates
                sp = pysynphot.Icat('ck04models', teff, m, logg)
                # renormalize the spectra:
                print(refband)
                sp = sp.renorm(0, 'vegamag', pysynphot.ObsBandpass('johnson,' + refband ))
                # some templates have all zeros which indicate invalid templates
                if len(sp.flux[sp.flux == 0]) == len(sp.flux):
                    print('Oh no! Failed to find matching templates')
                    specmatch = None
                WAVELENGTH = sp.wave * u.AA
                FLUX = sp.flux * u.erg/u.s/u.cm**2/u.AA
                return WAVELENGTH, FLUX
        else:
            specmatch = None
    else:
        specmatch = None

    if specmatch == None:
        print('Oops! Failed to find matching templates.')
        return None


def get_pickles_spec(BW=0.16, lam=551*u.nm, spec = None, teff = None, logg = None, m = None):
    """
    This function calculates the spectral flux density for a given 
    spectral type. Assumes the Pickles Atlas is saved to TargetList:
        ftp://ftp.stsci.edu/cdbs/grid/pickles/dat_uvk/
    If spectral type is provided, tries to match based on luminosity class,
    then spectral type. If no type, or not match, defaults to fit based on 
    Traub et al. 2016 (JATIS), which gives spectral flux density of
    ~9.5e7 [ph/s/m2/nm] @ 500nm
    
    Args:
        BW (float):
            Bandwidth fraction
        lam (astropy Quantity):
            Central wavelength in units of nm
        Spec (spectral type string):
            Should be something like G0V
            
    Returns:
        astropy Quantity:
            Spectral flux density in units of ph/m**2/s/nm.
    """
    
    # Begin by matching lam with nearest bandpass, then assign central
    # wavelength of ref band (UBVRIJK)
    bandeff = [365, 445, 551, 658, 806, 1220, 2190]*u.nm
    refLam = bandeff[(np.abs(bandeff-lam).argmin())]
    band_dict = {365: 'u', 445: 'b', 551: 'v', 658: 'r', 806: 'i',
                        1220: 'j', 1630: 'h', 2190: 'k', 3450: 'l',
                        4750: 'm'}
    refband = band_dict[refLam.value]

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
            lmax = lam*(1+BW/2)                                                    
            lmin = lam*(1-BW/2)
            # The case when the band falls within wavelength range supported by the Pickles Atlas:
            if lmax.to(u.Angstrom).value < 25000 and lmin.to(u.Angstrom).value > 1150:
                lumclass = specliste[:,2] == spece[2]
                ind = np.argmin( np.abs(spectypenum[lumclass] - (specdict[spece[0]]*10+spece[1]) ))
                specmatch = ''.join(specliste[lumclass][ind])
                # Open corresponding spectrum
                with fits.open(os.path.join(specdatapath,specindex[specmatch])) as hdulist:
                    sdat = hdulist[1].data
                WAVELENGTH = sdat.WAVELENGTH * u.AA
                FLUX = sdat.FLUX * u.erg/u.s/u.cm**2/u.AA
                sp = pysynphot.spectrum.ArraySourceSpectrum(wave=WAVELENGTH.value, flux=FLUX.value,waveunits='angstrom',fluxunits='flam')
                sp = sp.renorm(0, 'vegamag', pysynphot.ObsBandpass('johnson,' + refband ))
                return sp.wave * u.AA, sp.flux * u.erg/u.s/u.cm**2/u.AA
            # Otherwise, resort to the Castellie and Kurutz 2004 Atlas:
            else:
                specmatch = None
        else:
            specmatch = None
    else:
        specmatch = None

    if specmatch == None:
        return None
