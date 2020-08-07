"""
Created on Fri Jun 26 11:12:24 2020

@author: Anna Murphree (annammurphree@gmail.com)
For Summer 2020 Internship with Dr. Yaping Zhou

This code opens output files from findgrids.py and makes pairplots of their variables.
"""

def open_event_info(file, stats, ysn):
    from netCDF4 import Dataset
    import numpy as np
    import matplotlib.pyplot as plt
    
    # open the event stats from IMERG:
    stats = Dataset(stats, mode='r')
    params = stats.variables['event_parameters'][:]
    
    # create holders for event stats & averages:
    duras = []
    vols = []
    areas = []
    vmax_lats = []
    vmax_lons = []
    vavg_lats = []
    vavg_lons = []
    qv2ms_avg = []
    qv2mp_avg = []
    qv2me_avg = []
    q250s_avg = []
    q250p_avg = []
    q250e_avg = []
    q500s_avg = []
    q500p_avg = []
    q500e_avg = []
    q850s_avg = []
    q850p_avg = []
    q850e_avg = []
    u500s_avg = []
    u500p_avg = []
    u500e_avg = []
    v500s_avg = []
    v500p_avg = []
    v500e_avg = []
    
    for event in params:
        dstart = event[0]    # start time (days) of the year
        dend = event[1]      # end time (days) of the year
        year = event[22]     # year of event, which I add 1 to because the stats files are a year off
        year = int(year + 1)
        #print('Dates: ',dstart,'-', dend, ',',year)
        long = event[2]      # longitude (deg East)
        #print('IMERG long: ',long)
        lati = event[3]      # latitude (deg North)
        #print('IMERG lat: ',lati)
        dura = event[4]      # Extreme event duration [hours]
        vol = event[6]       # Extreme event rain volume - Heavy rain only [m^3]
        area = event[7]      # Extreme event total area - Heavy rain only [km^2]
        vmax_lat = event[18] # Maximum latitude component of system velocity [km hr^-1]
        vmax_lon = event[19] # Maximum longitude component of system velocity [km hr^-1]
        vavg_lat = event[20] # Mean latitude component of system velocity [km hr^-1]
        vavg_lon = event[21] # Mean longitude component of system velocity [km hr^-1]
        
        tmax = event[34]    # time into event of max rain volume (hours) 
        index = int(event[28])   # event index #
        #print('Index: ',index)
        tstart = int((dstart - int(dstart))*24)
        #print('Start hour:',tstart)
        tmax_days = dstart + (tmax/24)
        tmax_hrs = int((tmax_days - int(tmax_days))*24)
        #print('Peak hour: ',tmax_hrs)
        tend = int((dend - int(dend))*24)
        
        # small: dura < 6hrs
        # med  : 6 <= dura < 24hrs
        # big  : dura > 24hrs
        data = Dataset(file, mode='r')
        indices = len(data.variables['index'])
        
        if index < indices:
            qtest0 = data.variables['QV2M'][index,tstart,lati,long]
            qtest1 = data.variables['QV2M'][index,tmax_hrs,lati,long]
            qtest2 = data.variables['QV2M'][index,tend,lati,long]
            if qtest0 != np.nan and qtest1 != np.nan and qtest2 != np.nan:
                duras.append(dura)
                vols.append(vol)
                areas.append(area)
                vmax_lats.append(vmax_lat)
                vmax_lons.append(vmax_lon)
                vavg_lats.append(vavg_lat)
                vavg_lons.append(vavg_lon)
                
                qv2ms_avg.append(data.variables['QV2M'][index,tstart,lati,long])
                qv2mp_avg.append(data.variables['QV2M'][index,tmax_hrs,lati,long])
                qv2me_avg.append(data.variables['QV2M'][index,tend,lati,long])
                q250s_avg.append(data.variables['Q250'][index,tstart,lati,long])
                q250p_avg.append(data.variables['Q250'][index,tmax_hrs,lati,long])
                q250e_avg.append(data.variables['Q250'][index,tend,lati,long])
                q500s_avg.append(data.variables['Q500'][index,tstart,lati,long])
                q500p_avg.append(data.variables['Q500'][index,tmax_hrs,lati,long])
                q500e_avg.append(data.variables['Q500'][index,tend,lati,long])
                q850s_avg.append(data.variables['Q850'][index,tstart,lati,long])
                q850p_avg.append(data.variables['Q850'][index,tmax_hrs,lati,long])
                q850e_avg.append(data.variables['Q850'][index,tend,lati,long])
                u500s_avg.append(data.variables['U500'][index,tstart,lati,long])
                u500p_avg.append(data.variables['U500'][index,tmax_hrs,lati,long])
                u500e_avg.append(data.variables['U500'][index,tend,lati,long])
                v500s_avg.append(data.variables['V500'][index,tstart,lati,long])
                v500p_avg.append(data.variables['V500'][index,tmax_hrs,lati,long])
                v500e_avg.append(data.variables['V500'][index,tend,lati,long])
                
                data.close()
 
    qv2ms_avg = [np.round(num, 10) for num in qv2ms_avg]
    qv2mp_avg = [np.round(num, 10) for num in qv2mp_avg]
    qv2me_avg = [np.round(num, 10) for num in qv2me_avg]
    q250s_avg = [np.round(num, 10) for num in q250s_avg]
    q250p_avg = [np.round(num, 10) for num in q250p_avg]
    q250e_avg = [np.round(num, 10) for num in q250e_avg]
    q500s_avg = [np.round(num, 10) for num in q500s_avg]
    q500p_avg = [np.round(num, 10) for num in q500p_avg]
    q500e_avg = [np.round(num, 10) for num in q500e_avg]
    q850s_avg = [np.round(num, 10) for num in q850s_avg]
    q850p_avg = [np.round(num, 10) for num in q850p_avg]
    q850e_avg = [np.round(num, 10) for num in q850e_avg]
    u500s_avg = [np.round(num, 10) for num in u500s_avg]
    u500p_avg = [np.round(num, 10) for num in u500p_avg]
    u500e_avg = [np.round(num, 10) for num in u500e_avg]
    v500s_avg = [np.round(num, 10) for num in v500s_avg]
    v500p_avg = [np.round(num, 10) for num in v500p_avg]
    v500e_avg = [np.round(num, 10) for num in v500e_avg]
    
    import pandas as pd
    import seaborn as sns
    sns.set(style="ticks", color_codes=True)
    ''
    dfs = pd.DataFrame({'Duration': duras, 'QV2Ms_Avg': qv2ms_avg,'Q250s_Avg': q250s_avg, 
                        'Q500s_Avg': q500s_avg, 'U500s_Avg': u500s_avg,
                        'V500s_Avg': v500s_avg, 'Area': areas})
    
    pairs = sns.pairplot(dfs, height=3)
    pairs.axes[2,0].set_ylim([np.nanmin(q250s_avg)-np.nanstd(q250s_avg),np.nanmax(q250s_avg)+np.nanstd(q250s_avg)])
    pairs.axes[2,1].set_ylim([np.nanmin(q250s_avg)-np.nanstd(q250s_avg),np.nanmax(q250s_avg)+np.nanstd(q250s_avg)])
    pairs.axes[2,3].set_ylim([np.nanmin(q250s_avg)-np.nanstd(q250s_avg),np.nanmax(q250s_avg)+np.nanstd(q250s_avg)])
    pairs.axes[2,4].set_ylim([np.nanmin(q250s_avg)-np.nanstd(q250s_avg),np.nanmax(q250s_avg)+np.nanstd(q250s_avg)])
    pairs.axes[2,5].set_ylim([np.nanmin(q250s_avg)-np.nanstd(q250s_avg),np.nanmax(q250s_avg)+np.nanstd(q250s_avg)])
    pairs.axes[2,6].set_ylim([np.nanmin(q250s_avg)-np.nanstd(q250s_avg),np.nanmax(q250s_avg)+np.nanstd(q250s_avg)])
    pairs.axes[3,0].set_ylim([np.nanmin(q500s_avg)-np.nanstd(q500s_avg),np.nanmax(q500s_avg)+np.nanstd(q500s_avg)])
    pairs.axes[3,1].set_ylim([np.nanmin(q500s_avg)-np.nanstd(q500s_avg),np.nanmax(q500s_avg)+np.nanstd(q500s_avg)])
    pairs.axes[3,2].set_ylim([np.nanmin(q500s_avg)-np.nanstd(q500s_avg),np.nanmax(q500s_avg)+np.nanstd(q500s_avg)])
    pairs.axes[3,4].set_ylim([np.nanmin(q500s_avg)-np.nanstd(q500s_avg),np.nanmax(q500s_avg)+np.nanstd(q500s_avg)])
    pairs.axes[3,5].set_ylim([np.nanmin(q500s_avg)-np.nanstd(q500s_avg),np.nanmax(q500s_avg)+np.nanstd(q500s_avg)])
    pairs.axes[3,6].set_ylim([np.nanmin(q500s_avg)-np.nanstd(q500s_avg),np.nanmax(q500s_avg)+np.nanstd(q500s_avg)])
    plt.show()
    pairs.savefig(f'{ysn}_pairplots_s.png')
    ''
    dfp = pd.DataFrame({'Duration': duras, 'QV2Mp_Avg': qv2mp_avg,'Q250p_Avg': q250p_avg, 
                        'Q500s_Avg': q500p_avg, 'U500p_Avg': u500p_avg,
                        'V500s_Avg': v500p_avg, 'Area': areas})
    #print(dfp)
    pairp = sns.pairplot(dfp, height=3)
    pairp.axes[2,0].set_ylim([np.nanmin(q250s_avg)-np.nanstd(q250s_avg),np.nanmax(q250s_avg)+np.nanstd(q250s_avg)])
    pairp.axes[2,1].set_ylim([np.nanmin(q250s_avg)-np.nanstd(q250s_avg),np.nanmax(q250s_avg)+np.nanstd(q250s_avg)])
    pairp.axes[2,3].set_ylim([np.nanmin(q250s_avg)-np.nanstd(q250s_avg),np.nanmax(q250s_avg)+np.nanstd(q250s_avg)])
    pairp.axes[2,4].set_ylim([np.nanmin(q250s_avg)-np.nanstd(q250s_avg),np.nanmax(q250s_avg)+np.nanstd(q250s_avg)])
    pairp.axes[2,5].set_ylim([np.nanmin(q250s_avg)-np.nanstd(q250s_avg),np.nanmax(q250s_avg)+np.nanstd(q250s_avg)])
    pairp.axes[2,6].set_ylim([np.nanmin(q250s_avg)-np.nanstd(q250s_avg),np.nanmax(q250s_avg)+np.nanstd(q250s_avg)])
    pairp.axes[3,0].set_ylim([np.nanmin(q500s_avg)-np.nanstd(q500s_avg),np.nanmax(q500s_avg)+np.nanstd(q500s_avg)])
    pairp.axes[3,1].set_ylim([np.nanmin(q500s_avg)-np.nanstd(q500s_avg),np.nanmax(q500s_avg)+np.nanstd(q500s_avg)])
    pairp.axes[3,2].set_ylim([np.nanmin(q500s_avg)-np.nanstd(q500s_avg),np.nanmax(q500s_avg)+np.nanstd(q500s_avg)])
    pairp.axes[3,4].set_ylim([np.nanmin(q500s_avg)-np.nanstd(q500s_avg),np.nanmax(q500s_avg)+np.nanstd(q500s_avg)])
    pairp.axes[3,5].set_ylim([np.nanmin(q500s_avg)-np.nanstd(q500s_avg),np.nanmax(q500s_avg)+np.nanstd(q500s_avg)])
    pairp.axes[3,6].set_ylim([np.nanmin(q500s_avg)-np.nanstd(q500s_avg),np.nanmax(q500s_avg)+np.nanstd(q500s_avg)])
    plt.show()
    pairp.savefig(f'{ysn}_pairplots_p.png')
    
    dfe = pd.DataFrame({'Duration': duras, 'QV2Me_Avg': qv2me_avg,'Q250e_Avg': q250e_avg, 
                        'Q500e_Avg': q500e_avg, 'U500e_Avg': u500e_avg,
                        'V500s_Avg': v500e_avg, 'Area': areas})
    paire = sns.pairplot(dfe, height=3)
    paire.axes[2,0].set_ylim([np.nanmin(q250s_avg)-np.nanstd(q250s_avg),np.nanmax(q250s_avg)+np.nanstd(q250s_avg)])
    paire.axes[2,1].set_ylim([np.nanmin(q250s_avg)-np.nanstd(q250s_avg),np.nanmax(q250s_avg)+np.nanstd(q250s_avg)])
    paire.axes[2,3].set_ylim([np.nanmin(q250s_avg)-np.nanstd(q250s_avg),np.nanmax(q250s_avg)+np.nanstd(q250s_avg)])
    paire.axes[2,4].set_ylim([np.nanmin(q250s_avg)-np.nanstd(q250s_avg),np.nanmax(q250s_avg)+np.nanstd(q250s_avg)])
    paire.axes[2,5].set_ylim([np.nanmin(q250s_avg)-np.nanstd(q250s_avg),np.nanmax(q250s_avg)+np.nanstd(q250s_avg)])
    paire.axes[2,6].set_ylim([np.nanmin(q250s_avg)-np.nanstd(q250s_avg),np.nanmax(q250s_avg)+np.nanstd(q250s_avg)])
    paire.axes[3,0].set_ylim([np.nanmin(q500s_avg)-np.nanstd(q500s_avg),np.nanmax(q500s_avg)+np.nanstd(q500s_avg)])
    paire.axes[3,1].set_ylim([np.nanmin(q500s_avg)-np.nanstd(q500s_avg),np.nanmax(q500s_avg)+np.nanstd(q500s_avg)])
    paire.axes[3,2].set_ylim([np.nanmin(q500s_avg)-np.nanstd(q500s_avg),np.nanmax(q500s_avg)+np.nanstd(q500s_avg)])
    paire.axes[3,4].set_ylim([np.nanmin(q500s_avg)-np.nanstd(q500s_avg),np.nanmax(q500s_avg)+np.nanstd(q500s_avg)])
    paire.axes[3,5].set_ylim([np.nanmin(q500s_avg)-np.nanstd(q500s_avg),np.nanmax(q500s_avg)+np.nanstd(q500s_avg)])
    paire.axes[3,6].set_ylim([np.nanmin(q500s_avg)-np.nanstd(q500s_avg),np.nanmax(q500s_avg)+np.nanstd(q500s_avg)])
    plt.show()
    paire.savefig(f'{ysn}_pairplots_e.png')
    
    dfv = pd.DataFrame({'V_Avg Lats': vavg_lats, 'V_Avg Lons': vavg_lons,
                        'U500s_Avg': u500s_avg,'V500s_Avg': v500s_avg, 
                        'U500p_Avg': u500p_avg,'V500p_Avg': v500p_avg, 
                        'U500e_Avg': u500e_avg,'V500e_Avg': v500e_avg})
    
    pairv = sns.pairplot(dfv, height=3)
    pairv.savefig(f'{ysn}_pairplots_v.png')
    
#open_event_info('2016_JJA_var_avgs.nc4', 'extreme_event_stats_2016_JJA.nc')