from retrieve_spectral_coverage_utils import *

import pandas as pd
from astroquery.simbad import Simbad

SUP_DICT = {}


# data of stellar_spec_sup.csv pulled from:
# - https://arxiv.org/pdf/astro-ph/0507139.pdf
df = pd.read_csv('stellar_spec_sup.csv')
tmp = df.set_index('Name').T.to_dict('list')
for key in tmp:
    tmp[key] = tmp[key][0]
SUP_DICT_tmp = {}
for key in tmp:
    SUP_DICT_tmp[key.upper()] = tmp[key].strip()
SUP_DICT = SUP_DICT_tmp.copy()


# data of curr_calspec_sup.csv pulled from:
# - https://www.stsci.edu/hst/instrumentation/reference-data-for-calibration-and-tools/astronomical-catalogs/calspec
df = pd.read_csv('curr_calspec_sup.csv')
tmp = df.set_index('Name').T.to_dict('list')
for key in tmp:
    tmp[key] = tmp[key][0]
    tmp[key] = tmp[key].replace(u'\\xa0', u' ').strip()
SUP_DICT_tmp = {}
for key in tmp:
    SUP_DICT_tmp[key.upper()] = tmp[key]
SUP_DICT.update(SUP_DICT_tmp.copy())


# data of cohen_cw_spectra_spec_sup.csv scrapped from tem/dvl files in 'cohen' and 'cw_spectra'.
# - see retrieve_spectral_coverage_utils.make_spec_type_csv_cohen_cw_spectra()
df = pd.read_csv('cohen_cw_spectra_spec_sup.csv')
tmp = df.set_index('Name').T.to_dict('list')

for key in tmp:
    tmp[key] = tmp[key][0]
    tmp[key] = tmp[key].replace(u'\\xa0', u' ').strip()
SUP_DICT_tmp = {}
for key in tmp:
    SUP_DICT_tmp[key.upper()] = tmp[key]
SUP_DICT.update(SUP_DICT_tmp.copy())
    

# import os, warnings, astropy
# os.path.exists('EXOSIMS/StarCatalog/mission_exocat_2019.08.22_11.37.24.votable')
# from astropy.io.votable import parse
# catalogpath = 'EXOSIMS/StarCatalog/mission_exocat_2019.08.22_11.37.24.votable'
# with warnings.catch_warnings():
#     # warnings for IPAC votables are out of control 
#     #   they are not moderated by pedantic=False
#     #   they all have to do with units, which we handle independently anyway
#     warnings.simplefilter('ignore', 
#             astropy.io.votable.exceptions.VOTableSpecWarning)
#     warnings.simplefilter('ignore', 
#             astropy.io.votable.exceptions.VOTableChangeWarning)
#     votable = parse(catalogpath)
# table = votable.get_first_table()
# data = table.array
# names, specs = data['hd_name'], data['st_spttype'].astype(str)

# for n, s in zip(names, specs):
#     SUP_DICT[(n.decode("utf-8")).replace(" ", "").upper()] = s
    
# print(SUP_DICT.keys())

def find_spec_type_given_name(name):
    name = name.replace(' ', '').upper()
    if name in list(SUP_DICT.keys()):
        return SUP_DICT[name]
    
    print('Cannot find ' + name + ' in cached csv files. Resorting to Simbad.')
    customSimbad = Simbad()
    customSimbad.add_votable_fields('sptype')
    try:
        result = customSimbad.query_object(name)
        return result['SP_TYPE'][0].decode("utf-8")
    except:
        print(name + ' not found.')
        return None

import numpy as np
def find_spec_types_given_names(names):
    names = np.array(names)
    specs = np.empty(names.shape, dtype='|U20')

    cache_mask = np.isin(names, list(SUP_DICT.keys()))
    cache_indices = np.where(cache_mask)[0]
    for i in cache_indices:
        specs[i] = find_spec_type_given_name(names[i])

    if len(names[~cache_mask]) != 0:
        print('Some spectral types cannot be found in the cached csv files. Resorting to Simbad for those.')
        customSimbad = Simbad()
        customSimbad.add_votable_fields('sptype')
        result = customSimbad.query_objects(names[~cache_mask])
        result = [s.decode("utf-8") for s in result['SP_TYPE']]
        result = np.array(result)
        if result.shape[0] != names[~cache_mask].shape[0]:
            print('There exists name that cannot be found even with Simbad. Returning None.')
            return None
        specs[~cache_mask] = result
    return specs
    

# print(find_spec_types_given_names(['hd216386', 'hd218452', 'KF03T4']))



import os
from astropy.io import fits

def getSpCov_engelke():
    # Based on Table 1 of https://arxiv.org/pdf/astro-ph/0507139.pdf
    return {'K2III': (100, 3493.33), 'K1III': (100, 3493.33), 'K0III': (100, 3493.33)}

def getSpCov_irs_cal():
    # Based on https://arxiv.org/pdf/1408.5922.pdf
    # Explanation of SL1/2 - https://irsa.ipac.caltech.edu/data/SPITZER/docs/files/spitzer/irstr04002.pdf
    return {'A0V': (512.652, 1534.830), 'G9III': (512.652, 1534.830), 'K1III': (512.652, 1534.830)}

def getSpCov_calspec(dir_path):
    # Based on https://www.stsci.edu/hst/instrumentation/reference-data-for-calibration-and-tools/astronomical-catalogs/calspec
    spec_cov = {}
    for file in get_all_files_given_ext(dir_path, 'fits'):
 
        name = os.path.basename(file)
        name = name[:-5]
        # Surgery on 'name' so it matches with ones in the .csv files:
        extraneous = ['_stis', '_mod', '_nic', '_fos', '_00']
        for e in extraneous:
            if name.find(e) != -1:
                name = name[:name.find(e)]
                break # break here necessary
                
        spec = find_spec_type_given_name(name.upper())
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
        
    specs = find_spec_types_given_names(names)

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
    files = get_all_files_given_ext(dir_path, 'tem')
    names = []
    for file in files:
        name = os.path.basename(file)
        name = name[:-4]
        names.append(name.upper())
        
    specs = find_spec_types_given_names(names)

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

# print(getSpCov_calspec('/Users/berlinc/Downloads/calstar_templates/current_calspec'))
# print(getSpCov_cohen('/Users/berlinc/Downloads/calstar_templates/cohen'))
# print(getSpCov_cw_spectra('/Users/berlinc/Downloads/calstar_templates/cw_spectra'))


