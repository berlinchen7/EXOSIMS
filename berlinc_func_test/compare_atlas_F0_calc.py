import sys
import EXOSIMS,EXOSIMS.MissionSim,os.path
import astropy.units as u
import pandas as pd
from astroquery.simbad import Simbad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns; sns.set()
import numpy as np
import plotly.express as px


def compare_atlas_F0_calc(scriptfile):
    # pickles_band and CK_band are very much the same band, except that
    # CK_band exceeds 2.5 um just a little to trick TargetList.F0()
    # into using CK atlas for zero flux calculation.
    Pickles_band = (2190*u.nm, 0.28)
    CK_band = (2190*u.nm, 0.29)

    print('[compare_atlas_F0_calc] Instantiating EXOSIMS object.')
    sim = EXOSIMS.MissionSim.MissionSim(scriptfile) 
    
    targetListNames = sim.TargetList.Name.copy()
    targetListSpecs = sim.TargetList.Spec.copy()

    print('[compare_atlas_F0_calc] Checking TargetList.Spec against SYMBAD.')                                                                             
    SB_in = []                                                                 
    for name in targetListNames:
        name = name.split()                                                    
        SB_in.append(name[0] + ' ' + name[1])                                  
    customSimbad = Simbad()
    customSimbad.add_votable_fields('sptype')                                  
    result = customSimbad.query_objects(SB_in)                                 
    result = [s.decode("utf-8") for s in result['SP_TYPE']]                    
    targetListSpecs = result
    
    print('[compare_atlas_F0_calc] Calling TargetList.F0() to compute stellar zero flux.')
    star_name, specs, CK_F0, Pickles_F0 = [], [], [], []                                              
    for stInd, spec in enumerate(targetListSpecs):
        curr_CK_F0, _ = sim.TargetList.F0(CK_band[1], CK_band[0], stInd, spec)   
        curr_Pickles_F0, _ = sim.TargetList.F0(Pickles_band[1], Pickles_band[0], stInd, spec)

        star_name.append(targetListNames[stInd])
        specs.append(spec)
        CK_F0.append(curr_CK_F0.value)
        Pickles_F0.append(curr_Pickles_F0.value)

    F0_residual = np.array([Pickles_F0[i] - CK_F0[i] for i in range(len(Pickles_F0))])
    print('\n Statistics of distributions of Pickles F0 - CK F0 at central wavelength = 2190 nm; BW ~= 0.3:')
    print('  Mean: ' + str(np.mean(F0_residual)) + ' ph/m2/nm/s')
    print('  Median: ' + str(np.median(F0_residual)) + ' ph/m2/nm/s')
    print('  Standard deviation: ' + str(np.std(F0_residual)) + ' ph/m2/nm/s')
    print('\n')

    print('[compare_atlas_F0_calc] Generating scatterplot comparing Pickles F0 to CK F0.')
    df_in = {'Name': star_name, 'Spec': specs, 'CK F0': CK_F0, 'Pickles F0': Pickles_F0}
    df = pd.DataFrame(df_in)

    fig = px.scatter(df, x="CK F0", y="Pickles F0",
                     hover_data=['Name', 'Spec'],
                     title='Comparing zero-flux calculation at central wavelength = 2190 nm; BW ~= 0.3')
    # Add identity:
    minVal = min(min(df['CK F0']), min(df['Pickles F0']))
    maxVal = max(max(df['CK F0']), max(df['Pickles F0']))
    
    fig.update_layout(shapes=[
                dict(
                    type="line",
                    yref='y1',
                    y0=minVal,
                    y1=maxVal,
                    xref='x1',
                    x0=minVal,
                    x1=maxVal,
                    line=dict(
                        color="Red",
                        width=1,
                    )
                )
            ])
    
    filename = "F0_compare_scatterplot.html"
    fig.write_html(filename)
    print('[compare_atlas_F0_calc] Scatterplot saved as ' + filename + '.')

if __name__ == '__main__':
    scriptfile = sys.argv[1]
    compare_atlas_F0_calc(scriptfile)

# Example usage:
# python compare_atlas_F0_calc.py EXOSIMS/Scripts/HabEx_4m_coroOnly_DulzE_promo_EPRVmwants_lucky_charAsnr5_binMenn_20190929_mod.json
