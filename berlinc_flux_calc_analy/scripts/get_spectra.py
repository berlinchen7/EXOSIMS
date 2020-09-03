import pickle, re, pkg_resources, os, csv, warnings, sys
import matplotlib.pyplot as plt
from astropy.io import fits
import astropy.units as u
import astropy.constants as const
import numpy as np
# set up environment for and then import pysynphot
pysynphot_datapath = pkg_resources.resource_filename('EXOSIMS.TargetList','trds')
assert os.path.isdir(pysynphot_datapath), 'Could not locate %s in TargetList directory.' %(pysynphot_datapath)
os.environ['PYSYN_CDBS'] = pysynphot_datapath
warnings.filterwarnings("ignore", message="Extinction files not found") # we don't use extinction files 
import pysynphot

#set up stuff for spectral type conversion
# Paths
indexf =  pkg_resources.resource_filename('EXOSIMS.TargetList','pickles_index.pkl')
assert os.path.exists(indexf), "Pickles catalog index file not found in TargetList directory."
ck_indexf = pkg_resources.resource_filename('EXOSIMS.TargetList','ck_index.pkl')
assert os.path.exists(ck_indexf), "Castelli and Kurucz 2004 atlas index file not found in TargetList directory."

datapath = pkg_resources.resource_filename('EXOSIMS.TargetList','dat_uvk')
assert os.path.isdir(datapath), 'Could not locate %s in TargetList directory.' %(datapath)
ck_datapath = pkg_resources.resource_filename('EXOSIMS.TargetList','ckp00')
assert os.path.isdir(ck_datapath), 'Could not locate %s in TargetList directory.' %(ck_datapath)

# grab Pickles Atlas index
with open(indexf, 'rb') as handle:
    specindex = pickle.load(handle)
# grab Castelli and Kurucz 2004 Atlas index
with open(ck_indexf, 'rb') as handle:
    ck_specindex = pickle.load(handle) 

speclist = sorted(specindex.keys())
specdatapath = datapath
ck_speclist = sorted(ck_specindex.keys())
ck_specdatapath = ck_datapath

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

#everything in speclist and ck_speclist is correct, so only need first regexp
specliste = []
for spec in speclist:
    specliste.append(specregex1.match(spec).groups())
specliste = np.vstack(specliste)
spectypenum = np.array([specdict[l] for l in specliste[:,0]])*10+ np.array(specliste[:,1]).astype(float)

ck_specliste = []
for ck_spec in ck_speclist:
    ck_specliste.append(specregex1.match(ck_spec).groups())
ck_specliste = np.vstack(ck_specliste)
ck_spectypenum = np.array([specdict[l] for l in ck_specliste[:,0]])*10+ np.array(ck_specliste[:,1]).astype(float)

def get_ck_spec(spec):    
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
            # CK's README only indicated suggested templates for odd lum classes;
            #  hence, assuming subgiants are giants and bright giants are supergiants
            lumclasse = romandict[spece[2]] if romandict[spece[2]]%2 == 1\
                                           else romandict[spece[2]] - 1
            lumclasse = rev_romandict[lumclasse]
            assert lumclasse in ['I', 'III', 'V']
            lumclass = ck_specliste[:,2] == lumclasse
            ind = np.argmin(np.abs(ck_spectypenum[lumclass] 
                                   - (specdict[spece[0]]*10+spece[1])))
            specmatch = ''.join(ck_specliste[lumclass][ind])

            if specmatch is None:
                return None
            # values of ck_specindex are of the form ckp00_yyyyy[gxx],
            #  and we want the field gxx in ckp00_yyyyy.fits
            ck_info = ck_specindex[specmatch]
            ck_info = re.split('\[|\]', ck_info)
            ck_filename = ck_info[0] + '.fits'
            ck_field = ck_info[1]
            sdat = fits.getdata(os.path.join(ck_specdatapath, ck_filename))
            wave = sdat['WAVELENGTH'].copy()
            flux = sdat[ck_field].copy()
            sp = pysynphot.spectrum.ArraySourceSpectrum(wave=wave, flux=flux,\
                    waveunits='angstrom',fluxunits='flam')
            sp = sp.renorm(0, 'vegamag', pysynphot.ObsBandpass('johnson,v'))
            WAVELENGTH = sp.wave * u.AA
            FLUX = sp.flux * u.erg/u.s/u.cm**2/u.AA
            return WAVELENGTH, FLUX
    return None

def get_pickles_spec(spec):
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
            ind = np.argmin( np.abs(spectypenum[lumclass] - (specdict[spece[0]]*10+spece[1]) ))
            specmatch = ''.join(specliste[lumclass][ind])
            # Open corresponding spectrum
            with fits.open(os.path.join(specdatapath,specindex[specmatch])) as hdulist:
                sdat = hdulist[1].data
            WAVELENGTH = sdat.WAVELENGTH * u.AA
            FLUX = sdat.FLUX * u.erg/u.s/u.cm**2/u.AA
            sp = pysynphot.spectrum.ArraySourceSpectrum(wave=WAVELENGTH.value, flux=FLUX.value,waveunits='angstrom',fluxunits='flam')
            return sp.wave * u.AA, sp.flux * u.erg/u.s/u.cm**2/u.AA
        return None

def plot_spectra(spec, MAX_L=50000*u.AA):
    plt.figure(figsize=(15, 7.5))

    p_wavelength, p_flux = get_pickles_spec(spec)
    ck_wavelength, ck_flux = get_ck_spec(spec)

    plt.plot(p_wavelength[p_wavelength < MAX_L], p_flux[p_wavelength < MAX_L], label='Pickles Library')
    plt.plot(ck_wavelength[ck_wavelength < MAX_L], ck_flux[ck_wavelength < MAX_L], label='Castelli-Kurucz Atlas')
    plt.xlabel('Wavelength (Angstrom)')
    plt.ylabel('Irradiance (erg/s/cm**2/AA )')
    plt.title(spec)
    plt.legend()
    filename = '../plots/spectra_plots_' + spec + '.png'
    print('[get_spectra] saving spectra plot to ' + filename)
    plt.savefig(filename, dpi=300)

if __name__ == '__main__':                                                     
    spec = sys.argv[1]
    plot_spectra(spec)
    
# Example usage:                                                               
# python get_spectra.py G2III
