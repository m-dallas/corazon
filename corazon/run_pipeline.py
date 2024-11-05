import corazon.pipeline as pipeline
from datetime import datetime
import os
from exovetter import vetters
import matplotlib.pyplot as plt
import corazon.gen_lightcurve as genlc
#sys.path[2] = '/Users/smullally/Python_Code/lightkurve/lightkurve'


def run_write_one(ticid, sector, out_dir, lc_author = 'qlp',local_dir = None,
               run_tag = None, config_file = None, plot=False):
    """
    Run the full bls search on a list of ticids stored in a file.

    Parameters
    ----------
    ticid : int
       tess input catalog number
    sector : int
       tess sector to search
    out_dir : string
        directory to store all the results. One dir per ticid will be created.
    lc_author : string
        'qlp' or 'tess-spoc'
    local_dir : string
        defaul is None and then pulls data from MAST API. Otherwise contains
        directory name for the location of the data files.
    run_tag : string, optional
        directory name and string to attach to output file names. 

    Returns
    -------
    None.

    """
    
    if run_tag is None:
        now = datetime.now()
        run_tag = now.strftime("crz%m%d%Y") + "_"+lc_author
    
    # Default BLS parameters
    config = {
        "det_window" : 95,  #window used for detrending
        "noise_window" : 19, #window used for running outlier rejection
        "n_sigma" : 4.5,  #noise/outlier reject sigma
        "max_period_days" : 11,
        "min_period_days" : 0.8,
        "bls_durs_hrs" : [1,2,4,6,8,10,12],
        "minSnr" : [1],
        "maxTces" : 20,
        "fracRemain" : 0.1, #0.7 # was 0.7
        }
    
    # load exovetter vetters to get some statistics on the generated tces
    vetter_list = [vetters.LeoTransitEvents(), vetters.Sweet(), vetters.TransitPhaseCoverage()]
    
    # Set up where the output goes
    target_dir = "/tic%09is%02i/" % (int(ticid), sector)
    log_name = out_dir + target_dir + "tic%09i-%s.log" % (ticid, run_tag)
    output_file = out_dir + target_dir + "tic%09i-%s-tcesum.csv" % (ticid, run_tag)
    
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)

    try:
        os.mkdir(out_dir+target_dir)   
    except FileExistsError:
        pass
    except PermissionError as e:
        log_obj = open(log_name,'w+')
        log_obj.write("Permission Error on Target Directory ")
        log_obj.write(e)
        log_obj.close()
        
    try:
        
        # This will change to genlc.get_lc_from_S3()
        lcdata = genlc.hlsp(ticid, sector, author=lc_author,local_dir = local_dir)
        
        print(lcdata)
        
        if lc_author == 'qlp':
            lcdata['quality'] = lcdata['quality'].value & 2237
         
        tce_list, result_strings, metrics_list = pipeline.search_and_vet_one(ticid, 
                                sector, lcdata, config, 
                                vetter_list, plot=plot)

        # metrics_list is currently not being used!
        
        if plot:
            plotfilename = "tic%09i-%s-plot.png" % (ticid, 
                                                    run_tag)
            plt.savefig(out_dir + target_dir + plotfilename, bbox_inches='tight')
            plt.close()
        
        output_obj = open(output_file, 'w')
        for r in result_strings:
            output_obj.write(r)
        
    
        output_obj.close()
        
        #Write TCEs
        for tce in tce_list:
            tcefilename = "tic%09i-%02i-%s.json" % (ticid, 
                                                    int(tce['event']), 
                                                    run_tag)
    
            full_filename = out_dir + target_dir + tcefilename
            tce['lc_author'] = lc_author
            tce.to_json(full_filename)

        log_obj = open(log_name, 'w+')
        log_obj.write("Success.")
        log_obj.close()

    except Exception as e:
        log_obj = open(log_name,'w+')
        log_obj.write("Failed to create TCEs for TIC %i for Sector %i \n" % (ticid, sector))
        log_obj.write(str(e))
        log_obj.close() 

def run_write_one_from_s3(ticid, s3_location, sector, out_dir, lc_author = 'TGLC',local_dir = None,
               run_tag = None, config_file = None, plot=False):
    """
    Run the full bls search on a list of ticids stored in a file.

    Parameters
    ----------
    ticid : int
       tess input catalog number
    sector : int
       tess sector to search
    out_dir : string
        directory to store all the results. One dir per ticid will be created.
    lc_author : string
        'qlp' or 'tess-spoc'
    local_dir : string
        defaul is None and then pulls data from MAST API. Otherwise contains
        directory name for the location of the data files.
    run_tag : string, optional
        directory name and string to attach to output file names. 

    Returns
    -------
    None.

    """
    
    if run_tag is None:
        now = datetime.now()
        run_tag = now.strftime("crz%m%d%Y") + "_"+lc_author
    
    # Default BLS parameters
    config = {
        "det_window" : 95,  #window used for detrending
        "noise_window" : 19, #window used for running outlier rejection
        "n_sigma" : 4.5,  #noise/outlier reject sigma
        "max_period_days" : 11,
        "min_period_days" : 0.8,
        "bls_durs_hrs" : [1,2,4,6,8,10,12],
        "minSnr" : [1],
        "maxTces" : 20,
        "fracRemain" : 0.7
        }
    
    # load exovetter vetters to get some statistics on the generated tces
    vetter_list = [vetters.LeoTransitEvents(), vetters.Sweet(), vetters.TransitPhaseCoverage()]
    
    # Set up where the output goes
    target_dir = "/tic%09is%02i/" % (int(ticid), sector)
    log_name = out_dir + target_dir + "tic%09i-%s.log" % (ticid, run_tag)
    output_file = out_dir + target_dir + "tic%09i-%s-tcesum.csv" % (ticid, run_tag)
    
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)

    try:
        os.mkdir(out_dir+target_dir)   
    except FileExistsError:
        pass
    except PermissionError as e:
        log_obj = open(log_name,'w+')
        log_obj.write("Permission Error on Target Directory ")
        log_obj.write(e)
        log_obj.close()
        
    try:
        
        # This will change to genlc.get_lc_from_S3()
        #lcdata = genlc.hlsp(ticid, sector, author=lc_author,local_dir = local_dir)

        # so I changed it:
        lcdata = genlc.from_S3(s3_location)
        
        
        if lc_author == 'qlp':
            lcdata['quality'] = lcdata['quality'].value & 2237
         
        tce_list, result_strings, metrics_list = pipeline.search_and_vet_one(ticid, 
                                sector, lcdata, config, 
                                vetter_list, plot=plot)

        # metrics_list is currently not being used!
        
        if plot:
            plotfilename = "tic%09i-%s-plot.png" % (ticid, 
                                                    run_tag)
            plt.savefig(out_dir + target_dir + plotfilename, bbox_inches='tight')
            plt.close()
        
        output_obj = open(output_file, 'w')
        for r in result_strings:
            output_obj.write(r)
        
    
        output_obj.close()
        
        #Write TCEs
        for tce in tce_list:
            tcefilename = "tic%09i-%02i-%s.json" % (ticid, 
                                                    int(tce['event']), 
                                                    run_tag)
    
            full_filename = out_dir + target_dir + tcefilename
            tce['lc_author'] = lc_author
            tce.to_json(full_filename)

        log_obj = open(log_name, 'w+')
        log_obj.write("Success.")
        log_obj.close()

    except Exception as e:
        log_obj = open(log_name,'w+')
        log_obj.write("Failed to create TCEs for TIC %i for Sector %i \n" % (ticid, sector))
        log_obj.write(str(e))
        log_obj.close() 