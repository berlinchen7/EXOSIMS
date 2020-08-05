from retrieve_spectral_coverage_utils import *
from astroquery.simbad import Simbad
import numpy as np
import pandas as pd
import os
from astropy.io import fits
import pickle, os


def update_spec_type_dict_from_csv(csv_name, SUP_DICT):
    '''Helper for load_spec_type_map_from_csv.'''
    df = pd.read_csv(csv_name)
    tmp = df.set_index('Name').T.to_dict('list')
    # tmp is now {'star name': ['spec type with weird "\xao" characters'], ...}
    for key in tmp:
        tmp[key] = tmp[key][0]
        tmp[key] = tmp[key].replace(u'\\xa0', u' ').strip()
    SUP_DICT_tmp = {}
    for key in tmp:
        SUP_DICT_tmp[key.upper()] = tmp[key]
    SUP_DICT.update(SUP_DICT_tmp.copy())
    hist = []
    names = df['Name']
    for n in names:
        if n in hist:
            print(n)
        else:
            hist.append(n)
    return SUP_DICT

def load_spec_type_map_from_csv():
    ''' Loading cached data as a dictionary object (SUP_DICT)'''

    SUP_DICT = {}


    # data of curr_calspec_sup.csv pulled from:
    # - https://www.stsci.edu/hst/instrumentation/reference-data-for-calibration-and-tools/astronomical-catalogs/calspec
    SUP_DICT = update_spec_type_dict_from_csv('curr_calspec_sup.csv', SUP_DICT)

    # data of cohen_cw_spectra_spec_sup.csv scrapped from tem/dvl files in 'cohen' and 'cw_spectra'.
    # - see retrieve_spectral_coverage_utils.make_spec_type_csv_cohen_cw_spectra()
    SUP_DICT = update_spec_type_dict_from_csv('cohen_cw_spectra_spec_sup.csv', SUP_DICT)

    # data of curr_calspec_sup.csv pulled from:
    # - https://www.stsci.edu/hst/instrumentation/reference-data-for-calibration-and-tools/astronomical-catalogs/calspec
    SUP_DICT = update_spec_type_dict_from_csv('curr_calspec_sup.csv', SUP_DICT)

    # data of curr_calspec_sup.csv pulled from:
    # - https://www.stsci.edu/hst/instrumentation/reference-data-for-calibration-and-tools/astronomical-catalogs/castelli-and-kurucz-atlas
    # Note: ckp00_39000[g40] and ckp00_34000[g40] correspond to multiple associated spectral types (O6.5V/O6I and 
    #  O8.5V/O8I, respectively); for convenience, I deleted the O6V and O8.5V ones as other models have
    #  the same spectral type after standardization.
    SUP_DICT = update_spec_type_dict_from_csv('ck_spec_types.csv', SUP_DICT)
 
    # Load Pickle's index data:
    pickles_indexf= '../../EXOSIMS/TargetList/pickles_index.pkl'
    with open(pickles_indexf, 'rb') as handle:
        specindex = pickle.load(handle)
    for spec in specindex.keys():
        if spec == 'M10III': # Somehow pickles has 'M10III', which is not needed.
            continue
        SUP_DICT[specindex[spec]] = spec

    return SUP_DICT
    
    
def find_spec_type_given_name(name, SUP_DICT):
    
    name = name.replace(' ', '').upper()
    if name in list(SUP_DICT.keys()):
        return SUP_DICT[name]
    
    
    return None # Cannot find 'name' in cached csv files. Returning None.

def find_spec_types_given_names(names, SUP_DICT):
    
    names = np.array(names)
    specs = np.empty(names.shape, dtype='|U20')

    cache_mask = np.isin(names, list(SUP_DICT.keys()))
    cache_indices = np.where(cache_mask)[0]
    for i in cache_indices:
        specs[i] = find_spec_type_given_name(names[i], SUP_DICT)

    if len(names[~cache_mask]) != 0:
        print('WARNING: There exists names whose corresponding spectral types cannot be found.')
        specs[~cache_indices] = None

    return specs

def getSpCov_engelke():
    # Based on Table 1 of https://arxiv.org/pdf/astro-ph/0507139.pdf
    return {'K2III': (1000, 34933.3), 'K1III': (1009, 34933.3), 'K0III': (1000, 34933.3)}

def getSpCov_irs_cal():
    # Based on https://arxiv.org/pdf/1408.5922.pdf and engelke_zeromag_standard.tbl.txt
    # Explanation of SL1/2 - https://irsa.ipac.caltech.edu/data/SPITZER/docs/files/spitzer/irstr04002.pdf
    return {'A0V': (5126.52, 15348.30), 'G9III': (5126.52, 15348.30), 'K1III': (5126.52, 15348.30)}

def getSpCov_calspec(dir_path):
    # Based on https://www.stsci.edu/hst/instrumentation/reference-data-for-calibration-and-tools/astronomical-catalogs/calspec
    spec_cov = {}
    spec_map = load_spec_type_map_from_csv()
    for file in get_all_files_given_ext(dir_path, 'fits'):
 
        name = os.path.basename(file)
        name = name[:-5]
        # Surgery on 'name' so it matches with ones in the .csv files:
        extraneous = ['_stis', '_mod', '_nic', '_fos', '_00']
        for e in extraneous:
            if name.find(e) != -1:
                name = name[:name.find(e)]
                break # break here necessary

        spec = find_spec_type_given_name(name.upper(), spec_map)
        spec = standerize_spec_type(spec)
        if spec == None:
            print('Spec is none for the file ' + file)
            continue

        with fits.open(file) as hdulist:
                sdat = hdulist[1].data
        min_wv, max_wv = sdat.WAVELENGTH[0], sdat.WAVELENGTH[-1]
        min_wv, max_wv = min_wv/10, max_wv/10 # Convert from Angstrom to nm
        if spec in spec_cov:
            if max_wv > spec_cov[spec][1]:
                spec_cov[spec] = (min_wv, max_wv)
        else:
            spec_cov[spec] = (min_wv, max_wv)
    return spec_cov

def getSpCov_cohen(dir_path):
    spec_cov = {}
    files = get_all_files_given_ext(dir_path, 'tem') + get_all_files_given_ext(dir_path, 'dlv')
    names = []
    for file in files:
        name = os.path.basename(file)
        name = name[:-4]
        names.append(name.upper())
        
    spec_map = load_spec_type_map_from_csv()
    specs = find_spec_types_given_names(names, spec_map)

    for i, file in enumerate(files):    
        spec = specs[i]
        spec = standerize_spec_type(spec)
        if spec == None:
            print('Spec is none for the file ' + file)
            continue
        min_wv, max_wv = get_min_max_wavelength_from_tem_dlv(file)
        min_wv, max_wv = min_wv*1000, max_wv*1000 # Convert from microns to nm
        
        if spec in spec_cov:
            if max_wv > spec_cov[spec][1]:
                spec_cov[spec] = (min_wv, max_wv)
        else:
            spec_cov[spec] = (min_wv, max_wv)
    return spec_cov

def getSpCov_cw_spectra(dir_path):
    spec_cov = {}
    spec_map = load_spec_type_map_from_csv()
    files = get_all_files_given_ext(dir_path, 'tem')
    names = []
    for file in files:
        name = os.path.basename(file)
        name = name[:-4]
        names.append(name.upper())
        
    specs = find_spec_types_given_names(names, spec_map)

    for i, file in enumerate(files):    
        spec = specs[i]
        spec = standerize_spec_type(spec)
        if spec == None:
            print('Spec is none for the file ' + file)
            continue
        min_wv, max_wv = get_min_max_wavelength_from_tem_dlv(file)
        min_wv, max_wv = min_wv*1000, max_wv*1000 # Convert from microns to nm
        
        if spec in spec_cov:
            if max_wv > spec_cov[spec][1]:
                spec_cov[spec] = (min_wv, max_wv)
        else:
            spec_cov[spec] = (min_wv, max_wv)
    return spec_cov

def getSpCov_ck04models(dir_path):
    spec_cov = {}
    spec_map = load_spec_type_map_from_csv()
    files = get_all_files_given_ext(dir_path, 'fits')
    for file in files:
        name = os.path.basename(file)
        name = name[:-5] # Get rid of '.fits'
        
        g, spec = None, None
        for i in range(0, 51, 5):
            curr_spec = find_spec_type_given_name(name.upper() + "[G{:02d}]".format(i), spec_map)
            if curr_spec != None:
                g = "g{:02d}".format(i)
                spec = curr_spec
                break
                
        # In case when the model does not correspond to a spectral type:
        if spec == None:
            continue
            
        spec = standerize_spec_type(spec)

        with fits.open(file) as hdulist:
                sdat = hdulist[1].data
        min_wv, max_wv = sdat.WAVELENGTH[0], sdat.WAVELENGTH[-1]
        min_wv, max_wv = min_wv/10, max_wv/10 # Convert from Angstrom to nm
        if spec in spec_cov:
            if max_wv > spec_cov[spec][1]:
                spec_cov[spec] = (min_wv, max_wv)
        else:
            spec_cov[spec] = (min_wv, max_wv)
    return spec_cov

def getSpCov_pickles(indexf='../../EXOSIMS/TargetList/pickles_index.pkl'):
    with open(indexf, 'rb') as handle: 
        specindex = pickle.load(handle)
    ret = {}
    for spec_type in specindex.keys(): # Note spec_type is already in standardized format, so no need to standardize
        if spec_type == 'M10III': # Somehow pickles has 'M10III', which is not needed
            continue
        fitsfile = specindex[spec_type]
        fitsfile = os.path.dirname(indexf) + '/dat_uvk/' + fitsfile
        with fits.open(fitsfile) as hdulist:
                sdat = hdulist[1].data
        min_wv, max_wv = sdat.WAVELENGTH[0], sdat.WAVELENGTH[-1]
        min_wv, max_wv = min_wv/10, max_wv/10 # Convert from Angstrom to nm
        ret[spec_type] = (min_wv, max_wv)
    return ret


def getSpCov_bpgs(fdir='./bpgs'):
    ret = {}
    files = get_all_files_given_ext(fdir, 'fits')
    for f in files:
        with fits.open(f) as hdulist:
            header = hdulist[0].header
        sptype = header['MKTYPE']
        if sptype == '' or sptype == 'N': # Not sure what spectral type 'N' is, so skip
            continue
        ret[standerize_spec_type(sptype)] = (229, 2560)
    return ret
        
def retrieve_files_from_spec_type(spec, calstar_path):
    '''Exhaustively search through calstar_templates for template files corresponding
        to the given spectral type.'''

    ret = []
    spec_map = load_spec_type_map_from_csv()
    
    # engelke:
    engelke_map = {'K2III': [calstar_path + '/engelke/SAO17718.dat',
                             calstar_path + '/engelke/BD+681022.dat'], 
                   'K0III': [calstar_path + '/engelke/KF09T1.dat', 
                             calstar_path + '/engelke/KF08T3.dat'], 
                   'K1III': [calstar_path + '/engelke/KF06T2.dat', 
                             calstar_path + '/engelke/KF06T1.dat']}
    
    # irs_cal:
    irs_cal_map = {'A0V': [calstar_path + '/irs_cal/HR2194.txt'], 
                   'G9III': [calstar_path + '/irs_cal/HR6606.txt'],
                   'K1III': [calstar_path + '/irs_cal/HR7341.txt']}
    
    # cohen and cw_spectra:
    co_cws_map = {}
    files = get_all_files_given_ext(calstar_path + '/cw_spectra', 'tem')
    files += get_all_files_given_ext(calstar_path + '/cohen', 'dlv')
    files += get_all_files_given_ext(calstar_path + '/cohen', 'tem')
    for file in files:
        with open(file, 'r') as fobj:
            lines = fobj.readlines()   
        file_basename = os.path.basename(file)
        
        # Observation: for HD stars: spec type is the last info on the first line where there is the letter 'HD';
        #              for other stars: spec type is the last info on the first line.
        if file_basename.upper()[:2] == 'HD':
            for line in lines:
                if 'HD'in line:
                    curr_spec = line.split()[-1]
                    break
        else:
            curr_spec = lines[0].split()[-1]
        if curr_spec not in co_cws_map:
            co_cws_map[curr_spec] = [file]
        else:
            co_cws_map[curr_spec].append(file)
    
    # current_calspec:
    curr_calspec_map = {}
    files = get_all_files_given_ext(calstar_path + '/current_calspec/', 'fits')
    for file in files:
        name = os.path.basename(file)
        name = name[:-5]
        # Surgery on 'name' so it matches with ones in the .csv files:
        extraneous = ['_stis', '_mod', '_nic', '_fos', '_00']
        for e in extraneous:
            if name.find(e) != -1:
                name = name[:name.find(e)]
                break # break here necessary

        curr_spec = find_spec_type_given_name(name.upper(), spec_map)
        if curr_spec in curr_calspec_map:
            curr_calspec_map[curr_spec].append(file)
        else:
            curr_calspec_map[curr_spec] = [file]
    
    
    for key in engelke_map.keys():
        if standerize_spec_type(key) == spec:
            ret += engelke_map[key]
    for key in irs_cal_map.keys():
        if standerize_spec_type(key) == spec:
            ret += irs_cal_map[key]
    for key in co_cws_map.keys():
        if standerize_spec_type(key) == spec:
            ret += co_cws_map[key]
    for key in curr_calspec_map.keys():
        if standerize_spec_type(key) == spec:
            ret += curr_calspec_map[key]
            
    return ret
