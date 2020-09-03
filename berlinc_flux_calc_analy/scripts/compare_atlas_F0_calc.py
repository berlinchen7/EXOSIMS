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


def compare_atlas_F0_calc(scriptfile, useSYMBAD=True, verbose=True):
    bands = [(500*u.nm, 0.2), (625*u.nm, 0.2), (910*u.nm, 0.2), (1200*u.nm, 0.2), (2000*u.nm, 0.2)]
    CK_band = (2190*u.nm, 0.2)

    print('[compare_atlas_F0_calc] Instantiating EXOSIMS object.')
    sim = EXOSIMS.MissionSim.MissionSim(scriptfile) 
    
    targetListNames = sim.TargetList.Name.copy()
    targetListSpecs = sim.TargetList.Spec.copy()

    if useSYMBAD:
        print('[compare_atlas_F0_calc] using SYMBAD instead of TargetList.Spec.') 
        SB_in = []                                                                 
        for name in targetListNames:
            name = name.split()                                                    
            SB_in.append(name[0] + ' ' + name[1])                                  
        customSimbad = Simbad()
        customSimbad.add_votable_fields('sptype')                                  
        result = customSimbad.query_objects(SB_in)                                 
        result = [s.decode("utf-8") for s in result['SP_TYPE']]                    
        targetListSpecs = result

    for lamb, bw in bands:        
        if verbose:
            print('[compare_atlas_F0_calc] Calling TargetList.F0() to compute stellar zero flux.')
        star_name, specs, CK_F0, Pickles_F0 = [], [], [], []                                              
        for stInd, spec in enumerate(targetListSpecs):
            curr_CK_F0, _ = sim.TargetList.F0(bw, lamb, stInd, spec, usePickles=False)   
            curr_Pickles_F0, _ = sim.TargetList.F0(bw, lamb, stInd, spec, usePickles=True)
    
            star_name.append(targetListNames[stInd])
            specs.append(spec)
            CK_F0.append(curr_CK_F0.value)
            Pickles_F0.append(curr_Pickles_F0.value)
    
        F0_residual = np.array([Pickles_F0[i] - CK_F0[i] for i in range(len(Pickles_F0))])
        print('\n Statistics of distributions of Pickles F0 - CK F0 at central wavelength = ' + str(lamb.value) + ' nm;' +
                ' BW = ' + str(bw) + ':')
        print('  Mean: ' + str(np.mean(F0_residual)) + ' ph/m2/nm/s')
        print('  Standard deviation: ' + str(np.std(F0_residual)) + ' ph/m2/nm/s')
        print('  Q1: ' + str(np.quantile(F0_residual, .25)) + ' ph/m2/nm/s')
        print('  Median: ' + str(np.median(F0_residual)) + ' ph/m2/nm/s')
        print('  Q3: ' + str(np.quantile(F0_residual, .75)) + ' ph/m2/nm/s')
        print('\n')
    
        if verbose:
            print('[compare_atlas_F0_calc] Generating scatterplot comparing Pickles F0 to CK F0.')
        df_in = {'Name': star_name, 'Spec': specs, 'CK F0': CK_F0, 'Pickles F0': Pickles_F0}
        df = pd.DataFrame(df_in)
    
        fig = px.scatter(df, x="CK F0", y="Pickles F0",
                         hover_data=['Name', 'Spec'],
                         title='Comparing zero-flux calculation at central wavelength = ' + str(lamb.value) + ' nm;' +
                                ' BW = ' + str(bw) + ':')
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
        
        filename = '../plots/F0_compare_scatterplot_' + str(int(lamb.value)) + '_' + str(bw).replace('.', 'p')
        if useSYMBAD:
            filename = filename + '_SYMBAD.html'
        else:
            filename = filename + '.html'
        fig.write_html(filename)
        if verbose:
            print('[compare_atlas_F0_calc] Scatterplot saved as ' + filename + '.')

if __name__ == '__main__':
    scriptfile = sys.argv[1]
    compare_atlas_F0_calc(scriptfile)

# Example usage:
# python compare_atlas_F0_calc.py EXOSIMS/Scripts/HabEx_4m_coroOnly_DulzE_promo_EPRVmwants_lucky_charAsnr5_binMenn_20190929_mod.json
