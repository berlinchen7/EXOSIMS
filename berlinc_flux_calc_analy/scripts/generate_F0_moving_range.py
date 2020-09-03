import sys
import EXOSIMS,EXOSIMS.MissionSim,os.path
import astropy.units as u 
import pandas as pd
from astroquery.simbad import Simbad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns; sns.set()

def generate_F0_moving_range(scriptfile, useSYMBAD=True):
    # Bands to test:
    bandeff = [365, 445, 551, 658, 806, 1220, 1630, 2190, 3450, 4750]*u.nm
    band_BW_dict = {365*u.nm: 0.18, 445*u.nm: 0.21, 551*u.nm: 0.16, 658*u.nm: 0.21, 806*u.nm: 0.18,
                        1220*u.nm: 0.17, 1630*u.nm: 0.18, 2190*u.nm: 0.18, 3450*u.nm: 0.14,
                        4750*u.nm: 0.1}
    band_dict = {365*u.nm: 'U', 445*u.nm: 'B',
                     551*u.nm: 'V', 658*u.nm: 'R',
                     806*u.nm: 'I', 1220*u.nm: 'J',
                     1630*u.nm: 'H', 2190*u.nm: 'K',
                     3450*u.nm: 'L', 4750*u.nm: 'M'}
    
    # A priori knowledge on which bands exceeds 2.5 um, in which case the CK atlas is used:
    exceeding_bands = ['L', 'M']

#    sys.path += ['/home/berlinc/anaconda3/envs/es/lib/python3.7/site-packages']
    print('[generate_F0_moving_range] Instantiating EXOSIMS object.')
    sim = EXOSIMS.MissionSim.MissionSim(scriptfile)
    
    targetListNames = sim.TargetList.Name.copy()
    targetListSpecs = sim.TargetList.Spec.copy()
    
    if useSYMBAD:
        print('[generate_F0_moving_range] Using SYMBAD outputs instead of TargetList.Spec.')
        SB_in = []
        for name in targetListNames:
            name = name.split()
            SB_in.append(name[0] + ' ' + name[1])
        customSimbad = Simbad()
        customSimbad.add_votable_fields('sptype')
        result = customSimbad.query_objects(SB_in)
        result = [s.decode("utf-8") for s in result['SP_TYPE']]
        targetListSpecs = result
    
    print('[generate_F0_moving_range] Calling TargetList.F0() to compute stellar zero flux.')
    band, F0 = [], []
    for LAM in bandeff:
        for stInd, spec in enumerate(targetListSpecs):
            f0, mag = sim.TargetList.F0(band_BW_dict[LAM], LAM, stInd, spec)
            F0.append(f0.value)
            band.append(band_dict[LAM])
    
    df_in = {'Band': band, 'F0 (ph/m2/nm/s)': F0}
    df = pd.DataFrame(df_in)

    print('[generate_F0_moving_range] Generating lineplot of F0 data.')
    plt.figure(figsize=(15, 7))
    ax = sns.lineplot(x="Band", y="F0 (ph/m2/nm/s)", data=df, ci="sd", 
                  estimator="mean", sort=False) 
    plt.title('F0 calc of target list - mean with 1 sigma moving range')
    filename = '../plots/F0_calc_lineplot'
    if useSYMBAD:
        filename = filename + '_SYMBAD.png'
    else:
        filename = filename + '.png'
    plt.savefig(filename, dpi=300)
    print('[generate_F0_moving_range] Lineplot saved as ' + filename + '.')

if __name__ == '__main__':
    scriptfile = sys.argv[1]
    generate_F0_moving_range(scriptfile)

# Example usage:
# python generate_F0_moving_range.py EXOSIMS/Scripts/HabEx_4m_coroOnly_DulzE_promo_EPRVmwants_lucky_charAsnr5_binMenn_20190929_mod.json
