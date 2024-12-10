__all__ = ['search_and_vet_one', 'vet_tce','vet_all_tces','get_disposition']

import corazon.planetSearch as ps
import corazon.gen_lightcurve as genlc
import matplotlib.pyplot as plt
import exovetter.tce as TCE
import astropy.units as u
import exovetter.const as const
import lightkurve as lk
import numpy as np
from exovetter import vetters


def search_and_vet_one(ticid, sector, lcdata, config, vetter_list, plot=True):
    """
    Search and vet one ticid using config and vetter list
    
    Parameters
    ----------
    ticid : int
         TIC Identification number
    sector : int
        Sector of the TESS data to use for analysis
    lcdata : lightkkurve obect
        time and flux and quality flags populated
    config : dict
        configuration dictionary
    vetter_list : list
        list of vetters from exovetter to run

    Returns
    -------
    tce_tces : list
        list of exovetter TCEs for this target
    result_strings : str
       string version of tce and decision
      
    metrics_list : list
        all metrics, one per tce

    
    """
    
    time = lcdata['time'].value
    flux = lcdata['flux'].value
    # flags = lcdata['quality'].value

    # Take a lightcurve and clean it (This is susan's method, changing it to the basic lightkurve package detrending)
    # good_time, meddet_flux = ps.clean_timeseries(time, flux, flags,
    #                                       config["det_window"], 
    #                                       config["noise_window"], 
    #                                       config["n_sigma"], 
    #                                       sector)

    # basic lightkurve based cleaning
    lcdata = lcdata[~np.isnan(lcdata["flux"])] # remove nans
    lcdata = lcdata.remove_outliers(sigma=config["n_sigma"]) # sigma clipping
    lcdata = lcdata.flatten() # Savitzky-Golay filter

    good_time = lcdata['time'].value
    meddet_flux = lcdata['flux'].value # "meddet_flux" might not realllly be the output anymore but I'll pass it to identify Tces etc rather than changing the name everywhere 
    flags = lcdata['quality'].value # I think we're removing points not just flagging them so might not actually be catching the flags- this is just used in plotting though 

    # Run BLS and get TCEs from it using the config dictionary
    tce_list, stats = ps.identifyTces(good_time, meddet_flux, 
                                      bls_durs_hrs=config["bls_durs_hrs"],
                                      minSnr=config["minSnr"], 
                                      fracRemain=config["fracRemain"], 
                                      maxTces=config["maxTces"], 
                                      minP=config["min_period_days"], 
                                      maxP=config["max_period_days"])
    
    if plot:
        plot_lc_tce(ticid, tce_list, time, flux, flags, good_time, 
                    meddet_flux, stats, sector)
    
    lcformat = lcdata['time'].format
    tce_lc = lk.LightCurve(time=good_time, flux=meddet_flux+1,
                        time_format=lcformat, meta={'sector':sector})
    
    result_strings, metrics_list, tce_tces = vet_all_tces(tce_lc, 
                                                    tce_list, ticid, 
                                                    vetter_list,
                                                    plot=False)
    
    return tce_tces, result_strings, metrics_list
    
def vet_all_tces(lc, tce_dict_list, ticid, vetter_list, plot=False):
    lcformat = lc['time'].format
    result_list = []
    metrics_list = []
    tce_list = []
    pn = 1
    for item in tce_dict_list:
        tce = TCE.Tce(period = item[0]*u.day, epoch=item[1]*u.day, 
                      depth=item[2] * const.frac_amp,
                      duration=item[3]*u.day, 
                      epoch_offset=const.string_to_offset[lcformat],
                      snr=item[4],
                      target = f"TIC {ticid}",
                      sector = lc.sector,
                      event = f"{pn}")

        metrics = vet_tce(tce, lc, vetter_list, plot=plot) # dictionary of all vetting metrics

        metrics['snr'] = tce['snr'] # Potentially could steal this from LEO? Although LEO will give you a nan if it doesn't work well on the given tce

        result_string = make_result_string(tce)
        tce_list.append(tce)
        result_list.append(result_string)
        metrics_list.append(metrics)
        pn=pn+1
        
    return result_list, metrics_list, tce_list

def vet_tce(tce, tce_lc, vetter_list, plot=False):
    """Create a dictionary with all vetting parameters returned by each vetter in vetter_list"""
    metrics = dict()
    for v in vetter_list:
        vetter = v
        
        try:
            vetter_dict = vetter.run(tce, tce_lc)
        except ValueError:
            print('Did not Vet')
            pass
        if plot:
            vetter.plot()
        metrics.update(vetter_dict)
        
    return metrics

def make_result_string(tce):
    """
   Create a string that summarizes the TCE and its disposition
    Parameters
    ----------
    tce : TYPE
        DESCRIPTION.
    disposition : string
        DESCRIPTION.
    reason : string
        DESCRIPTION.

    Returns
    -------
    None.

    """
    st = "%s, %s, %i, %8.4f, %9.4f, %8.3f, %5.3f, %5.2f\n" % \
                                        (tce['target'], tce['event'],
                                               tce['sector'],
                                               tce['period'].value, 
                                               tce['epoch'].value,
                                               tce['depth'].value*1e6,
                                               tce['duration'].value*24.0,
                                               tce['snr'])
                                               
    return st

def plot_lc_tce(ticid, tce_list, time, flux, flags, good_time, 
                good_flux, stats, sector):
    col = ['tab:orange','tab:green','tab:purple','tab:brown',
               'gold','magenta','lightpink']
    plt.figure(figsize=(10,6))
    plt.subplot(211)
    plt.plot(good_time, good_flux,'.')
    plt.title("Lightcurve for TIC %i in S%i" % (int(ticid), int(sector)))
   
    axes = plt.gca()
    y_min, y_max = axes.get_ylim()
    x_min, x_max = axes.get_xlim()
    for n,s in enumerate(stats):
        plt.vlines(stats[n]['transit_times'], y_min, y_max, 
                   colors=col[n], zorder=1, label=str(n+1))
    plt.legend()
    plt.subplot(212)
    plt.plot(time, flux,'.', label="original lc", color='red')
    plt.plot(good_time, good_flux,'.', label="detrended lc", color='green')
    #plt.plot(time[flags!=0], flux[flags!=0],'o', ms=3, label='flagged')
    plt.legend()
    plt.xlim(x_min, x_max)

def open_output_file(filename, headerlist):
    fobj = open(filename, 'a')
    
    
    header = headerlist[0]
    for h in headerlist[1:]:
        header = header + ", " + h
        
    fobj.write(header + '\n')
    
    return fobj

