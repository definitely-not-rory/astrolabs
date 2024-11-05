from imports import *



def get_periods():
    
    df = pd.read_csv('starsdata.csv',delimiter='\t',header=None)
    
    arrays = df.to_numpy()[:,:]
    
    print(arrays)
    
    mags = np.array(arrays[:,0])
    mags_errs  = np.array(arrays[:,1])
    periods = np.array(arrays[:,2])
    periods_errs = np.array(arrays[:,3])
    
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
            
    
    
    #retrieval_type = 'ALL'          # Options are: 'EPOCH_PHOTOMETRY', 'MCMC_GSPPHOT', 'MCMC_MSC', 'XP_SAMPLED', 'XP_CONTINUOUS', 'RVS', 'ALL'
    #data_structure = 'RAW'     # Options are: 'INDIVIDUAL' or 'RAW'
    #data_release   = 'Gaia DR3'     # Options are: 'Gaia DR3' (default), 'Gaia DR2'
    datalink = Gaia.load_data(ids=gaiadr3_ids, data_release='Gaia DR3', retrieval_type='ALL', data_structure='RAW')
    
    '''
    dl_keys  = [inp for inp in datalink.keys()]
    dl_keys.sort()
    print(f'The following Datalink products have been downloaded:')
    for dl_key in dl_keys:
        print(f' * {dl_key}')
    '''
    
    print(vars(datalink["EPOCH_PHOTOMETRY_RAW.xml"]))
    
    
    
    
    
#mags, mags_errs, periods, periods_errs = get_periods()

get_parallaxes()



