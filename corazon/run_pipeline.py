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
               run_tag = None, local=True, plot=False):
    """
    Run the full bls search on a list of ticids stored in a file.

    Parameters
    ----------
    ticid : int
       tess input catalog number
    s3_location : str
        string of MAST s3 location of TGLC lightcurve
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
    local : bool
        Specifies the file system to use

    Returns
    -------
    None.

    """
    
    #print('running!')
    if run_tag is None:
        now = datetime.now()
        run_tag = now.strftime("crz%m%d%Y") + "_"+lc_author

    if not local:
        fs = s3fs.S3FileSystem(anon=False, profile="default")
    else:
        fs = LocalFileSystem()
    
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
        
        #print(lcdata)
        
        if lc_author == 'qlp':
            lcdata['quality'] = lcdata['quality'].value & 2237
         
        tce_list, result_strings, metrics_list = pipeline.search_and_vet_one(ticid, 
                                sector, lcdata, config, 
                                vetter_list, plot=plot)

        # metrics_list is currently not being used!
        # print(metrics_list) # see which ones we want to use
        
        if plot:
            plotfilename = "tic%09i-%s-plot.png" % (ticid, 
                                                    run_tag)
            plt.savefig(out_dir + target_dir + plotfilename, bbox_inches='tight')
            plt.close()
        
        # change to fs.open(output_file, 'w') as fp:
        # output_obj = open(output_file, 'w')
        # for r in result_strings:
        #     output_obj.write(r)
        #output_obj.close()
        
        # new way to save csv with vetters
        with fs.open(output_file, 'w') as fp: 
            for i,r in enumerate(result_strings):
                newstr = ", %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f " % (metrics_list[i]['MES'],
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
                newr = r[:-1]+newstr
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
            #tce.to_json(full_filename)

        with fs.open(log_name, 'w') as fp:
            fp.write("Success!")
        # log_obj = open(log_name, 'w+')
        # log_obj.write("Success.")
        # log_obj.close()

    except Exception as e:
        with fs.open(log_name,'w') as fp:
            fp.write("Failed to create TCEs for TIC %i for Sector %i \n" % (ticid, sector))
            fp.write(str(e))
            return ["failed", output_file, str(e)]
        # log_obj = open(log_name,'w+')
        # log_obj.write("Failed to create TCEs for TIC %i for Sector %i \n" % (ticid, sector))
        # log_obj.write(str(e))
        # log_obj.close() 

def run_one_froms3_coiled(ticid, s3_location, sector, write_config,  
                config = None, plot=False, vetter_list=None):
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
    if write_config['run_tag'] is None:
        now = datetime.now()
        run_tag = now.strftime("crz%m%d%Y")
    else:
        run_tag = write_config['run_tag']    

    
    if plot and (fs[1:2] != "Lo"):
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
            with fs.open("s3://" + out_dir + target_dir + plotfilename) as fp:
                plt.savefig(fp, bbox_inches='tight')
        
        # new way to save csv with vetters
        with fs.open(output_file, 'w') as fp: 
            for i,r in enumerate(result_strings):
                newstr = ", %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f " % (metrics_list[i]['MES'],
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
            return ["failed", output_file, str(e)]

    