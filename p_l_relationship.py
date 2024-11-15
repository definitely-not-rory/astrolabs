from imports import *



def get_periods():
    
    df = pd.read_csv("periodmagdata.txt", delimiter=" ", header=None)
    
    mags = df[0].to_numpy()
    mags_errs = df[1].to_numpy()
    periods = df[2].to_numpy()
    periods_errs = df[3].to_numpy()
    
    
    return mags, mags_errs, periods, periods_errs


def get_parallaxes():
    
    objs=next(os.walk('.'))[1]
    objs = objs[1:-1]
    
    gaiadr3_ids = []
    
    Simbad.add_votable_fields('ids')
    
    for obj in objs:
        
        result = Simbad.query_object(obj)
        
        for id in result['ids'][0].split('|'):
            
            if 'Gaia DR3' in id:
                
                gaia_id = id.split(' ')[-1]
            
                break
            
        gaiadr3_ids = np.append(gaiadr3_ids,gaia_id)
            
    
    
    
        
    query = f"SELECT parallax \
    FROM gaiadr3.gaia_source \
    WHERE source_id = " + gaiadr3_ids[0] + "\
    OR source_id = " + gaiadr3_ids[1] + "\
    OR source_id = " + gaiadr3_ids[2] + "\
    OR source_id = " + gaiadr3_ids[3] + "\
    OR source_id = " + gaiadr3_ids[4] + "\
    OR source_id = " + gaiadr3_ids[5] + "\
    OR source_id = " + gaiadr3_ids[6] + "\
    OR source_id = " + gaiadr3_ids[7] + "\
    OR source_id = " + gaiadr3_ids[8] + "\
    OR source_id = " + gaiadr3_ids[9] + "\
    OR source_id = " + gaiadr3_ids[10] + "\
    OR source_id = " + gaiadr3_ids[11] + "\
    OR source_id = " + gaiadr3_ids[12] + "\
    OR source_id = " + gaiadr3_ids[13] + "\
    OR source_id = " + gaiadr3_ids[14] + "\
    OR source_id = " + gaiadr3_ids[15] + "\
    OR source_id = " + gaiadr3_ids[16] + "\
    OR source_id = " + gaiadr3_ids[17] + "\
    OR source_id = " + gaiadr3_ids[18] + "\
    OR source_id = " + gaiadr3_ids[19] + "\
    OR source_id = " + gaiadr3_ids[20] + "\
    OR source_id = " + gaiadr3_ids[21] + "\
    OR source_id = " + gaiadr3_ids[22] + "\
    OR source_id = " + gaiadr3_ids[23]


    job     = Gaia.launch_job_async(query)
    results = job.get_results()
    
    parallaxes = []
    
    for i in range(24):
        
        parallaxes = np.append(parallaxes,results['parallax'][i])
        
    dist_pc = 1/(0.001 * parallaxes)                                        #GAIA parallaxes are in mas, convert to arcseconds here
    
    return dist_pc[0:8]



        

def plot_pl(mags, mags_errs, periods, periods_errs, dist_pc):
    
    abs_mags = mags - 5 * np.log10(dist_pc) + 5
    
    plt.errorbar(periods, abs_mags, xerr=periods_errs, yerr=(0.05*np.abs(abs_mags) + 0.05), linestyle=" ")
    
    plt.show()
    
    
    
mags, mags_errs, periods, periods_errs = get_periods()

print(mags)
print(mags_errs)

dist_pc = get_parallaxes()

plot_pl(mags, mags_errs, periods, periods_errs, dist_pc)





