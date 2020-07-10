import re
import os
import pandas as pd

def get_all_files_given_ext(dir_path, ext):
    import glob
    return glob.glob(dir_path +'/**/*.' + ext, recursive=True)

def get_min_max_wavelength_from_tem_dlv(file):
    f = open(file, 'r')
    # Get rid of empty lines:
    lines = [line for line in f.readlines() if line.strip()]
    f.close()
    data = lines
    for i, l in enumerate(lines):
        # Spacing differs from file to file, as well as capitalization of W
        if l.replace(' ', '').upper() == '(microns)(W/cm2/um)(W/cm2/um)(%)(%)\n'.upper():
            data = lines[i+2:]
            break
    if len(data) == len(lines):
        print('Error with parsing ' + file)
        return None

    min_wv = data[0].strip().split()
    min_wv = float(min_wv[0])
    max_wv = data[-1].strip().split()
    max_wv = float(max_wv[0])
    
    return min_wv, max_wv

def get_spectra_from_tem_dlv(file):
    f = open(file, 'r')
    # Get rid of empty lines:
    lines = [line for line in f.readlines() if line.strip()]
    f.close()
    data = lines
    for i, l in enumerate(lines):
        # Spacing differs from file to file, as well as capitalization of W
        if l.replace(' ', '').upper() == '(microns)(W/cm2/um)(W/cm2/um)(%)(%)\n'.upper():
            data = lines[i+2:]
            break
    if len(data) == len(lines):
        print('Error with parsing ' + file)
        return None
    
    wavelength, irrad = [], []
    for d in data:
        wavelength.append(float(d.split()[0]))
        irrad.append(float(d.split()[1]))
    return wavelength, irrad
    
def get_spectra_from_ascii(file):
    f = open(file, 'r')
    # Get rid of empty lines:
    lines = [line for line in f.readlines() if line.strip()]
    f.close()
    
    wavelength, irrad = [], []
    for l in lines:
        entry = l.split()
        if len(entry) == 3: # for some reason the file engelke/KF08T3.dat 
                            #  has weird values such as '-' at the beginning of a line or
                            #  a missing decimal point, e.g., '1 48391' (supposed to be '1.48391')
            try:
                int(entry[0])
                entry[0] = entry[0] + '.' + entry[1]
                entry[1] = entry[2]
                entry = entry[:2]
            except:
                entry = entry[1:]                
        
        try: # for some reason the file engelke/KF08T3.dat has messed up wavelengths such as 2.322.7
            w = float(entry[0])
        except:
            print('Faulty wavelength: ' + entry[0] + ' in ' + file + '.')
            continue
        wavelength.append(w)
        irrad.append(float(entry[1]))

    return wavelength, irrad
    
def make_spec_type_csv_cohen_cw_spectra(calstar_tem_path):
    files = get_all_files_given_ext(calstar_tem_path + '/cohen', 'tem')
    files += get_all_files_given_ext(calstar_tem_path + '/cohen', 'dlv')
    files += get_all_files_given_ext(calstar_tem_path + '/cw_spectra', 'tem')

    Name, Spec = [], []
    for file in files:
        with open(file, 'r') as fobj:
            lines = fobj.readlines()
            
        file_basename = os.path.basename(file)
        Name.append(file_basename[:-4].upper())
        
        # Observation: for HD stars: spec type is the last info on the first line where there is the letter 'HD';
        #              for other stars: spec type is the last info on the first line.
        if file_basename.upper()[:2] == 'HD':
            for line in lines:
                if 'HD'in line:
                    Spec.append(line.split()[-1])
                    break
            assert len(Spec) == len(Name)
        else:
            Spec.append(lines[0].split()[-1])
    
    # Remove redundencies:
    unique_Name, unique_Spec = [], []
    for i, n in enumerate(Name):
        if n in Name[:i]:
            j = Name[:i].index(n)
            if Spec[i] != Spec[j]:
                print('Conflicting spec types for ' + n + ' (either ' + Spec[j] + ' or ' + Spec[i] + ' ). Opt for '\
                       + Spec[j] + '.')
            continue
        unique_Name.append(n)
        unique_Spec.append(Spec[i])
    
    csv_input = {'Name': unique_Name, 'Spec': unique_Spec}
    df = pd.DataFrame(csv_input)
    df.to_csv('cohen_cw_spectra_spec_sup.csv', index=False)
            
            
def standerize_spec_type(spec):
    if len(spec) == 0:
#         print('spectral type ' + spec + ' cannot be standerized.')
        return None
    
    if spec[:2] == 'sd': # subdwarfs
        spec = spec[2:] + '.0VI' if '.' not in spec else spec[2:] + '0VI'
    elif spec[0] == 'D': # white dwarfs
        spec = spec[1:] + '.0VII'if '.' not in spec else spec[1:] + '0VII'
    #spectral type decomposition
    #default string: Letter|number|roman numeral
    #number is either x, x.x, x/x
    #roman numeral is either 
    #either number of numeral can be wrapped in ()
    specregex1 = re.compile(r'([OBAFGKMLTY])\s*\(*(\d*\.\d+|\d+|\d+\/\d+)\)*\s*\(*([IV]+\/{0,1}[IV]*)')
    #next option is that you have something like 'G8/K0IV'
    specregex2 = re.compile(r'([OBAFGKMLTY])\s*(\d+)\/[OBAFGKMLTY]\s*\d+\s*\(*([IV]+\/{0,1}[IV]*)')
    #next down the list, just try to match leading vals and assume it's a dwarf
    specregex3 = re.compile(r'([OBAFGKMLTY])\s*(\d*\.\d+|\d+\/\d+|\d+)')
    #last resort is just match spec type
    specregex4 = re.compile(r'([OBAFGKMLTY])')

    # Try to decompose the input spectral type
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
#                 print('spectral type ' + spec + ' cannot be standerized.')
                return None  
    return spece[0] + str(int(spece[1])) + spece[2]

def test():
    assert standerize_spec_type('sdF8.5') == 'F8VI'
    assert standerize_spec_type('F8.5') == 'F8V'
    assert standerize_spec_type('sdO') == 'O0VI'
    assert standerize_spec_type('...') == None
    assert standerize_spec_type('M2.5IIIFe-0.5') == 'M2III'
#test()
