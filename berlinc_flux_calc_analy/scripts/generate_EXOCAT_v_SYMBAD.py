import sys
import EXOSIMS,EXOSIMS.MissionSim,os.path
import astropy.units as u
import pandas as pd
from astroquery.simbad import Simbad
import numpy as np


def generate_EXOCAT_v_SYMBAD(scriptfile):
    print('[generate_EXOCAT_v_SYMBAD] Instantiating EXOSIMS object.')
    sim = EXOSIMS.MissionSim.MissionSim(scriptfile) 
    
    targetListNames = sim.TargetList.Name.copy()
    targetListSpecs = sim.TargetList.Spec.copy()

    print('[generate_EXOCAT_v_SYMBAD] Retrieving corresponding SIMBAD spec types.') 
    SB_in = []                                                                 
    for name in targetListNames:
        name = name.split()                                                    
        SB_in.append(name[0] + ' ' + name[1])                                  
    customSimbad = Simbad()
    customSimbad.add_votable_fields('sptype')                                  
    result = customSimbad.query_objects(SB_in)                                 
    result = [s.decode("utf-8") for s in result['SP_TYPE']]                    
    SIMBADSpecs = result

    df_in = {'name': targetListNames, 'spec_EXO': targetListSpecs, 'spec_SIMBAD': SIMBADSpecs}
    df = pd.DataFrame(df_in)
    filename = '../plots/spec_EXO_vs_SIMBAD.csv'
    print('[generate_EXOCAT_v_SYMBAD] Saving to ' + filename)
    df.to_csv(filename, index=False)

if __name__ == '__main__':
    scriptfile = sys.argv[1]
    generate_EXOCAT_v_SYMBAD(scriptfile)

# Example usage:
# python generate_EXOCAT_v_SYMBAD.py EXOSIMS/Scripts/HabEx_4m_coroOnly_DulzE_promo_EPRVmwants_lucky_charAsnr5_binMenn_20190929_mod.json
