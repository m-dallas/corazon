#%%
import numpy as np
import s3fs
import lightkurve as lk
import ast
import time
import matplotlib.pyplot as plt
import fsspec
from fsspec.implementations.local import LocalFileSystem
import s3fs

from corazon import run_pipeline

import coiled


cp_s3_df = pd.read_excel('/Users/smullally/Science/tess_false_alarms/search2024/coil_bls/catalog_table_s3_paths.xlsx')
#cp_s3_df = pd.read_excel('/Users/smullally/Science/tess_false_alarms/search2024/coil_bls/cp_s3_locs.xlsx')
cp_s3_df = cp_s3_df.drop('Unnamed: 0', axis=1)

#This creates 3 array of information based on the s3 path names
s3_locs = []
sectors = []
tics = []
for i, row in enumerate(cp_s3_df['s3_paths']):
    tic = cp_s3_df['TIC'][i]
    for s3_path in ast.literal_eval(row):
        s3_locs.append(s3_path)
        sectors.append(int(s3_path[33:35]))
        tics.append(tic)

print(s3_locs[0:3], tics[0:3])


#%% ------------
#Arrange information in to a list of lists
targetinfo = []

for i, row in enumerate(cp_s3_df['s3_paths']):
    tic = cp_s3_df['TIC'][i]
    for s3_path in ast.literal_eval(row):   #don't do this, better to remove the string

        info = [tic, s3_path, int(s3_path[33:35])]
        targetinfo.append(info)
print(len(targetinfo))
#%%
import corazon.pipeline as pipeline
from datetime import datetime
import os
from exovetter import vetters
import matplotlib.pyplot as plt
import corazon.gen_lightcurve as genlc
import s3fs
import json
from fsspec.implementations.local import LocalFileSystem
from astropy.utils.misc import JsonCustomEncoder
import matplotlib
def run_one_froms3_coiled(ticid, s3_location, sector, write_config,  
                config = None, plot=False, vetter_list=None, local=True):
    """
    Run the full bls search on a list of ticids stored in a file.
    Read from s3 and store output files in s3.  write_config contains the info.

    Parameters
    ----------
    ticid : int
       tess input catalog number
    s3_location : str
        string of MAST s3 location of TGLC lightcurve
    sector : int
       tess sector to search
    write_config : dictionary
        A dictionary containing 'fs', A fsspec filesystem. Either s3fs.S3FileSystem 
         or a LocalFileSystem, and 'outdir' for either the 
         directory location or the s3 bucket name, and 'run_tag' a string to put in the filename
         set to None to get the current day and lc_author as a string for keeping track
    config : dictionary
        A dictionary to describe how to detrend and run the BLS
    vetter_list: list
        List of vetters from exovetter, e.g. [vetters.LeoTransitEvents()]

    Returns
    -------
    None.

    """
  # Set up where the output goes by reading in the 
    fs = write_config['fs']
    lc_author = write_config['lc_author']
    out_dir = write_config['outdir']

    #if not os.path.exists(out_dir):
    #    os.mkdir(out_dir)
    if local:
        if not fs.exists(out_dir):
            fs.mkdir(out_dir)

    if write_config['run_tag'] is None:
        now = datetime.now()
        run_tag = now.strftime("crz%m%d%Y")
    else:
        run_tag = write_config['run_tag']    

    
    if plot:
        matplotlib.use('Agg')
    
    # Default BLS parameters
    config = {
        "det_window" : 95,  #window used for detrending
        "noise_window" : 19, #window used for running outlier rejection
        "n_sigma" : 3.0,  #noise/outlier reject sigma
        "max_period_days" : 11,
        "min_period_days" : 0.8,
        "bls_durs_hrs" : [1,2,4,6,8,10,12],
        "minSnr" : [1],
        "maxTces" : 20,
        "fracRemain" : 0.7
        }
    
    # load exovetter vetters to get some statistics on the generated tces
    vetter_list = [vetters.LeoTransitEvents(), vetters.Sweet(), 
                   vetters.TransitPhaseCoverage()]

  

    # Set up where the output goes
    target_dir = "/tic%09is%02i/" % (int(ticid), sector)
    #if not os.path.exists(out_dir+target_dir):
    #    os.mkdir(out_dir+target_dir)
    if local:
        if not fs.exists(out_dir + target_dir):
            fs.mkdir(out_dir + target_dir)

    log_name = out_dir + target_dir + "tic%09i-%s.log" % (ticid, run_tag)
    output_file = out_dir + target_dir + "tic%09i-%s-tcesum.csv" % (ticid, run_tag)
    
    try:
        
        lcdata = genlc.from_S3(s3_location)
        
        if lc_author == 'qlp':
            lcdata['quality'] = lcdata['quality'].value & 2237
         
        tce_list, result_strings, metrics_list = pipeline.search_and_vet_one(ticid, 
                                sector, lcdata, config, 
                                vetter_list, plot=plot)
        
        if plot:
        
            plotfilename = "tic%09i-%s-plot.png" % (ticid, 
                                                    run_tag)
            #plotfilename = "s3://" + out_dir + target_dir + plotfilename
            plotfilename = out_dir + target_dir + plotfilename 
            print(f"Writing to {plotfilename}")
            with fs.open(plotfilename, 'wb') as fp:
                plt.savefig(fp, bbox_inches='tight')
        
        # new way to save csv with vetters
        with fs.open(output_file, 'w') as fp: 
            for i,r in enumerate(result_strings):
                newstr = ", %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f \n" % (metrics_list[i]['MES'],
                                           metrics_list[i]['SHP'],
                                           metrics_list[i]['CHI'],
                                           metrics_list[i]['med_chases'],
                                           metrics_list[i]['mean_chases'],
                                           metrics_list[i]['mean_chases'],
                                           metrics_list[i]['max_SES'],
                                           metrics_list[i]['DMM'],
                                           metrics_list[i]['amp'][2][0], # last array in Sweet (amplitude to uncertainty ratio): half-period
                                           metrics_list[i]['amp'][2][1], # period
                                           metrics_list[i]['amp'][2][2], # twice the period
                                           metrics_list[i]['transit_phase_coverage'],
                                           metrics_list[i]['snr'])
                newr = r[:-1] + newstr
                fp.write(newr)
    
        
        #Write TCEs
        for tce in tce_list:
            tcefilename = "tic%09i-%02i-%s.json" % (ticid, 
                                                    int(tce['event']), 
                                                    run_tag)
    
            full_filename = out_dir + target_dir + tcefilename
            tce['lc_author'] = lc_author
            with fs.open(full_filename, 'w') as fp:
                json.dump(tce, fp, cls=JsonCustomEncoder)

        with fs.open(log_name, 'w') as fp:
            fp.write("Success!")

        return tce_list

    except Exception as e:
        with fs.open(log_name,'w') as fp:
            fp.write("Failed to create TCEs for TIC %i for Sector %i \n" % (ticid, sector))
            fp.write(str(e))
            #return ["failed", output_file, str(e)]
            print(output_file)
            print(str(e))
            raise



#%% -----------------
#Run it locally first.
write_config={
        'fs': LocalFileSystem(),
        'outdir': '/Users/smullally/Science/tess_false_alarms/search2024/test',
        'run_tag': 'tmp_20241205b',
        'lc_author': 'TGLC'
        }
tinfo = targetinfo[13]
result = run_one_froms3_coiled(tinfo[0],tinfo[1],tinfo[2], 
                        write_config, plot=True, local=True)
print(result)

#%%
#Run Just one.
import s3fs
write_config={
        'fs': s3fs.S3FileSystem(anon=False, profile="default"),
        'outdir': 'corazon-search-2024/search-test2024',
        'run_tag': 'tmp_20241205c',
        'lc_author': 'TGLC'
        }
tinfo = targetinfo[9]
result = run_one_froms3_coiled(tinfo[0],tinfo[1],tinfo[2], 
                        write_config, plot=True, local=False)
print(result)

#%%
@coiled.function()
def wrapped_run_s3_search(tinfo):
    try:
        write_config={
        'fs': s3fs.S3FileSystem(anon=False),
        'outdir': 'corazon-search-2024/search-tglc-dec2024_full',
        'run_tag': '20241213b',
        'lc_author': 'TGLC'
        }
        result = run_one_froms3_coiled(tinfo[0], tinfo[1], tinfo[2], 
                               write_config,  plot=False, local=False)

    except Exception as e:
        return str(e)
    
    return result

#%%
start = time.time()
num=50000
search_result = list(wrapped_run_s3_search.map(targetinfo[10000:num]))

m = (time.time()-start)/60
print(f'finished in {m} m for {num} in my list')
# %%

search_result = wrapped_run_s3_search(targetinfo[5])
# %%
