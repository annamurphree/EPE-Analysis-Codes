"""
Created on Fri Jun 26 15:04:15 2020

@author: Anna Murphree (annammurphree@gmail.com)
For Summer 2020 Internship with Dr. Yaping Zhou


"""
def plot_stats(axs, n1, n2, v1, v2i, v2s, v2p, v2e, v2t, var1, scale=10):
    import numpy as np
    axs[n1,n2].scatter(v2i, v1, color='black', alpha=0.5, label=r't$_{start}-12$', s=scale)
    axs[n1,n2].scatter(v2s, v1, color='red', alpha=0.5, label=r't$_{start}$', s=scale)
    axs[n1,n2].scatter(v2p, v1, color='blue', alpha=0.5, label=r't$_{peak}$', s=scale)
    axs[n1,n2].scatter(v2e, v1, color='green', alpha=0.5, label=r't$_{end}$', s=scale)
    axs[n1,n2].scatter(v2t, v1, color='orange', alpha=0.5, label=r't$_{end}+12$', s=scale)
    axs[n1,n2].set_ylim([np.nanmin(v1)-np.nanstd(v1),np.nanmax(v1)+np.nanstd(v1)])
    axs[n1,n2].set_xlim([np.nanmin(v2p)-np.nanstd(v2p),np.nanmax(v2p)+np.nanstd(v2p)])
    #axs[0,0].set(xlabel='Event Duration (hours)', ylabel=r'1000 * QV2M Average (kg kg$^{-1}$)')
    if n1 == 0:
        axs[n1,n2].set_title(r'%s'%(var1), fontsize=12)

def percent_plot(arr, params, data, var, ysn, dlimb, dlima, nd, lev=0):
    import numpy as np
    from netCDF4 import Dataset
    import matplotlib.pyplot as plt
    # find percentiles of rain volume (input arr):
    #pers = [90, 80, 70, 60, 50, 40, 30, 20, 10, 1, 0.1, 0.01, 0.001]
    e90 = arr[(np.abs(arr - np.percentile(arr, 90))).argmin()]
    e80 = arr[(np.abs(arr - np.percentile(arr, 80))).argmin()]
    e70 = arr[(np.abs(arr - np.percentile(arr, 70))).argmin()]
    e60 = arr[(np.abs(arr - np.percentile(arr, 60))).argmin()]
    e50 = arr[(np.abs(arr - np.percentile(arr, 50))).argmin()]
    e40 = arr[(np.abs(arr - np.percentile(arr, 40))).argmin()]
    e30 = arr[(np.abs(arr - np.percentile(arr, 30))).argmin()]
    e20 = arr[(np.abs(arr - np.percentile(arr, 20))).argmin()]
    e10 = arr[(np.abs(arr - np.percentile(arr, 10))).argmin()]
    e1 = arr[(np.abs(arr - np.percentile(arr, 1))).argmin()]
    e01 = arr[(np.abs(arr - np.percentile(arr, 0.1))).argmin()]
    e001 = arr[(np.abs(arr - np.percentile(arr, 0.01))).argmin()]
    e0001 = arr[(np.abs(arr - np.percentile(arr, 0.001))).argmin()]
    
    # find the index of the events where the rain volume is those percentiles:
    i90 = [int(event[28]) for event in params if event[6] == e90][0]
    i80 = [int(event[28]) for event in params if event[6] == e80][0]
    i70 = [int(event[28]) for event in params if event[6] == e70][0]
    i60 = [int(event[28]) for event in params if event[6] == e60][0]
    i50 = [int(event[28]) for event in params if event[6] == e50][0]
    i40 = [int(event[28]) for event in params if event[6] == e40][0]
    i30 = [int(event[28]) for event in params if event[6] == e30][0]
    i20 = [int(event[28]) for event in params if event[6] == e20][0]
    i10 = [int(event[28]) for event in params if event[6] == e10][0]
    i1 = [int(event[28]) for event in params if event[6] == e1][0]
    i01 = [int(event[28]) for event in params if event[6] == e01][0]
    i001 = [int(event[28]) for event in params if event[6] == e001][0]
    i0001 = [int(event[28]) for event in params if event[6] == e0001][0]
    
    inds = [i90, i80, i70, i60, i50, i40, i30, i20, i10, i1, i01, i001, i0001]
    #print(percents)
    dat = Dataset(data, mode='r')
    indices = len(dat.variables['index'])
    #dat.variables[f'{var}']
        
    vi = []
    vs = []
    vp = []
    ve = []
    vt = []
    vit = []
    vst = []
    vpt = []
    vet = []
    vtt = []
    # for each percentile, find the var averages/stddevs (matched by index):
    gen = (i for i in inds if i < indices)
    if nd == 2:
        for i in gen:
            '''
            ds = [event[0] for event in params if event[28] == i][0]
            de = [event[1] for event in params if event[28] == i][0]
            tp = [event[34] for event in params if event[28] == i][0]
            #print(ds, de, tp)
            tstart = int((ds - int(ds))*24)
            tmax_days = ds + (tp/24)
            tmax_hrs = int((tmax_days - int(tmax_days))*24)
            tend = int((de - int(de))*24)
            '''
            vi.append(dat.variables[f'{var}'][i,0,0,0])
            vs.append(dat.variables[f'{var}'][i,1,0,0])
            vp.append(dat.variables[f'{var}'][i,2,0,0])
            ve.append(dat.variables[f'{var}'][i,3,0,0])
            vt.append(dat.variables[f'{var}'][i,4,0,0])
            
            vit.append(dat.variables[f'{var}'][i,0,0,1])
            vst.append(dat.variables[f'{var}'][i,1,0,1])
            vpt.append(dat.variables[f'{var}'][i,2,0,1])
            vet.append(dat.variables[f'{var}'][i,3,0,1])
            vtt.append(dat.variables[f'{var}'][i,4,0,1])
    elif nd == 3:
        for i in gen:
            vi.append(dat.variables[f'{var}'][i,0,lev,0])
            vs.append(dat.variables[f'{var}'][i,1,lev,0])
            vp.append(dat.variables[f'{var}'][i,2,lev,0])
            ve.append(dat.variables[f'{var}'][i,3,lev,0])
            vt.append(dat.variables[f'{var}'][i,4,lev,0])
            
            vit.append(dat.variables[f'{var}'][i,0,lev,1])
            vst.append(dat.variables[f'{var}'][i,1,lev,1])
            vpt.append(dat.variables[f'{var}'][i,2,lev,1])
            vet.append(dat.variables[f'{var}'][i,3,lev,1])
            vtt.append(dat.variables[f'{var}'][i,4,lev,1])
    '''    
    vi = [v for v in vi if v < 1e30]
    vs = [v for v in vs if v < 1e30]
    vp = [v for v in vp if v < 1e30]
    ve = [v for v in ve if v < 1e30]
    vt = [v for v in vt if v < 1e30]
    vit = [v for v in vit if v < 1e30]
    vst = [v for v in vst if v < 1e30]
    vpt = [v for v in vpt if v < 1e30]
    vet = [v for v in vet if v < 1e30]
    vtt = [v for v in vtt if v < 1e30]
    '''
    xlabs = ['90', '80', '70', '60', '50', '40', '30', '20', '10', '1', '0.1', '0.01', '0.001']
    xs = np.arange(13)
    plt.xticks(xs, xlabs)
    plt.plot(xs, vi, color='yellow', alpha=0.5, label=r't$_{start}$-12')
    plt.errorbar(xs, vi, yerr=vit)
    plt.plot(xs, vs, color='r', alpha=0.5, label=r't$_{start}$')
    plt.errorbar(xs, vs, yerr=vst)
    plt.plot(xs, vp, color='b', alpha=0.5, label=r't$_{peak}$')
    plt.errorbar(xs, vp, yerr=vpt)
    plt.plot(xs, ve, color='g', alpha=0.5, label=r't$_{end}$')
    plt.errorbar(xs, ve, yerr=vet)
    plt.plot(xs, vt, color='black', alpha=0.5, label=r't$_{end}$+12')
    plt.errorbar(xs, vt, yerr=vtt)
    plt.ylim(top=np.nanmax(vi))
    plt.ylabel(f'{var}')
    plt.xlabel('Rain Volume Percentile')
    plt.legend()
    plt.tight_layout(pad=0.1)
    plt.savefig(f'{ysn}_percentplot_{var}_dlimb{dlimb}_dlima{dlima}_lev{lev}_new.png')
    plt.show()
''
#def stat_list2d(var, name):
#    var_0 = var[0,0,0]

def open_event_info(file, stats, ysn, dlimb=0, dlima=1000):
    from netCDF4 import Dataset
    import numpy as np
    import matplotlib.pyplot as plt
    
    # open the event stats from IMERG:
    stat = Dataset(stats, mode='r')
    params = stat.variables['event_parameters'][:]
    
    # create holders for event stats:
    duras = []
    vols = []
    areas = []
    #vmax_lats = []
    #vmax_lons = []
    vavg_lats = []
    vavg_lons = []
    ''
    # variable averages (2D):
    qv2mi_avg = []
    qv2ms_avg = []
    qv2mp_avg = []
    qv2me_avg = []
    qv2mt_avg = []
    
    q250i_avg = []
    q250s_avg = []
    q250p_avg = []
    q250e_avg = []
    q250t_avg = []
    
    q500i_avg = []
    q500s_avg = []
    q500p_avg = []
    q500e_avg = []
    q500t_avg = []
    
    q850i_avg = []
    q850s_avg = []
    q850p_avg = []
    q850e_avg = []
    q850t_avg = []
    
    u500i_avg = []
    u500s_avg = []
    u500p_avg = []
    u500e_avg = []
    u500t_avg = []
    
    v500i_avg = []
    v500s_avg = []
    v500p_avg = []
    v500e_avg = []
    v500t_avg = []
    
    tqvi_avg = []
    tqvs_avg = []
    tqvp_avg = []
    tqve_avg = []
    tqvt_avg = []
    
    t2mi_avg = []
    t2ms_avg = []
    t2mp_avg = []
    t2me_avg = []
    t2mt_avg = []
    ''
    
    data = Dataset(file, mode='r')
    indices = len(data.variables['index'])
    #print(data.variables)
    
    for event in params:
        '''
        dstart = event[0]    # start time (days) of the year
        dend = event[1]      # end time (days) of the year
        year = event[22]     # year of event, which I add 1 to because the stats files are a year off
        year = int(year + 1)
        #print('Dates: ',dstart,'-', dend, ',',year)
        #long = event[2]      # longitude (deg East)
        #print('IMERG long: ',long)
        #lati = event[3]      # latitude (deg North)
        #print('IMERG lat: ',lati)
        '''
        dura = event[4]      # Extreme event duration [hours]
        vol = event[6]       # Extreme event rain volume - Heavy rain only [m^3]
        area = event[7]      # Extreme event total area - Heavy rain only [km^2]
        #vmax_lat = event[18] # Maximum latitude component of system velocity [km hr^-1]
        #vmax_lon = event[19] # Maximum longitude component of system velocity [km hr^-1]
        vavg_lat = event[20] # Mean latitude component of system velocity [km hr^-1]
        vavg_lon = event[21] # Mean longitude component of system velocity [km hr^-1]
        
        index = int(event[28])   # event index #
        #print('Index: ',index)
        '''
        tmax = event[34]    # time into event of max rain volume (hours) 
        tstart = int((dstart - int(dstart))*24)
        #print('Start hour:',tstart)
        tmax_days = dstart + (tmax/24)
        tmax_hrs = int((tmax_days - int(tmax_days))*24)
        #print('Peak hour: ',tmax_hrs)
        tend = int((dend - int(dend))*24)
        '''
        # small: dura < 6hrs       (dlimb=0,  dlima=6)
        # med  : 6 <= dura < 24hrs (dlimb=6,  dlima=24)
        # big  : dura > 24hrs      (dlimb=24, dlima=inf)
        if index < indices and (dlimb <= dura <= dlima):
        #if index < indices:
            qtest0 = data.variables['QV2M'][index,1,0,0]
            #print(qtest0)
            qtest1 = data.variables['QV2M'][index,2,0,0]
            #print(qtest1)
            qtest2 = data.variables['QV2M'][index,3,0,0]
            #print(qtest2)
            if qtest0 != np.nan and qtest1 != np.nan and qtest2 != np.nan:
                duras.append(dura)
                vols.append(vol)
                areas.append(area)
                #vmax_lats.append(vmax_lat)
                #vmax_lons.append(vmax_lon)
                vavg_lats.append(vavg_lat)
                vavg_lons.append(vavg_lon)
                '''
                # 2D:
                qv2m = data.variables['QV2M'][index]
                q250 = data.variables['Q250'][index]
                q500 = data.variables['Q500'][index]
                q850 = data.variables['Q850'][index]
                u500 = data.variables['U500'][index]
                v500 = data.variables['V500'][index]
                tqv = data.variables['TQV'][index]    # total precipitable water vapor
                t2m = data.variables['T2M'][index]    # air temperature 2m
                # 3D:
                vv = data.variables['OMEGA'][index]   # vertical pressure velocity
                rh = data.variables['RH'][index]      # relative humidity after moist
                qv = data.variables['QV'][index]      # specific humidity
                u = data.variables['U'][index]        # eastward wind
                v = data.variables['V'][index]        # northward wind
                temp = data.variables['T'][index]     # temperature
                
                '''
                qv2mi_avg.append(data.variables['QV2M'][index,0,0,0])
                qv2ms_avg.append(data.variables['QV2M'][index,1,0,0])
                qv2mp_avg.append(data.variables['QV2M'][index,2,0,0])
                qv2me_avg.append(data.variables['QV2M'][index,3,0,0])
                qv2mt_avg.append(data.variables['QV2M'][index,4,0,0])
                
                q250i_avg.append(data.variables['Q250'][index,0,0,0])
                q250s_avg.append(data.variables['Q250'][index,1,0,0])
                q250p_avg.append(data.variables['Q250'][index,2,0,0])
                q250e_avg.append(data.variables['Q250'][index,3,0,0])
                q250t_avg.append(data.variables['Q250'][index,4,0,0])
                
                q500i_avg.append(data.variables['Q500'][index,0,0,0])
                q500s_avg.append(data.variables['Q500'][index,1,0,0])
                q500p_avg.append(data.variables['Q500'][index,2,0,0])
                q500e_avg.append(data.variables['Q500'][index,3,0,0])
                q500t_avg.append(data.variables['Q500'][index,4,0,0])
                
                q850i_avg.append(data.variables['Q850'][index,0,0,0])
                q850s_avg.append(data.variables['Q850'][index,1,0,0])
                q850p_avg.append(data.variables['Q850'][index,2,0,0])
                q850e_avg.append(data.variables['Q850'][index,3,0,0])
                q850t_avg.append(data.variables['Q850'][index,4,0,0])
                
                u500i_avg.append(data.variables['U500'][index,0,0,0])
                u500s_avg.append(data.variables['U500'][index,1,0,0])
                u500p_avg.append(data.variables['U500'][index,2,0,0])
                u500e_avg.append(data.variables['U500'][index,3,0,0])
                u500t_avg.append(data.variables['U500'][index,4,0,0])
                
                v500i_avg.append(data.variables['V500'][index,0,0,0])
                v500s_avg.append(data.variables['V500'][index,1,0,0])
                v500p_avg.append(data.variables['V500'][index,2,0,0])
                v500e_avg.append(data.variables['V500'][index,3,0,0])
                v500t_avg.append(data.variables['V500'][index,4,0,0])
                
                tqvi_avg.append(data.variables['TQV'][index,0,0,0])
                tqvs_avg.append(data.variables['TQV'][index,1,0,0])
                tqvp_avg.append(data.variables['TQV'][index,2,0,0])
                tqve_avg.append(data.variables['TQV'][index,3,0,0])
                tqvt_avg.append(data.variables['TQV'][index,4,0,0])
                
                t2mi_avg.append(data.variables['T2M'][index,0,0,0])
                t2ms_avg.append(data.variables['T2M'][index,1,0,0])
                t2mp_avg.append(data.variables['T2M'][index,2,0,0])
                t2me_avg.append(data.variables['T2M'][index,3,0,0])
                t2mt_avg.append(data.variables['T2M'][index,4,0,0])
                ''
                
                #data.close()
    '''
    percent_plot(vols, params, file, 'QV2M', ysn, dlimb, dlima, 2)
    percent_plot(vols, params, file, 'Q250', ysn, dlimb, dlima, 2)
    percent_plot(vols, params, file, 'Q500', ysn, dlimb, dlima, 2)
    percent_plot(vols, params, file, 'Q850', ysn, dlimb, dlima, 2)
    percent_plot(vols, params, file, 'U500', ysn, dlimb, dlima, 2)
    percent_plot(vols, params, file, 'V500', ysn, dlimb, dlima, 2)
    percent_plot(vols, params, file, 'TQV', ysn, dlimb, dlima, 2)
    percent_plot(vols, params, file, 'T2M', ysn, dlimb, dlima, 2)
    ''
    percent_plot(vols, params, file, 'OMEGA', ysn, dlimb, dlima, 3, 1)
    percent_plot(vols, params, file, 'RH', ysn, dlimb, dlima, 3, 1)
    percent_plot(vols, params, file, 'QV', ysn, dlimb, dlima, 3, 1)
    percent_plot(vols, params, file, 'U', ysn, dlimb, dlima, 3, 1)
    percent_plot(vols, params, file, 'V', ysn, dlimb, dlima, 3, 1)
    percent_plot(vols, params, file, 'T', ysn, dlimb, dlima, 3, 1)
    
    percent_plot(vols, params, file, 'OMEGA', ysn, dlimb, dlima, 3, 2)
    percent_plot(vols, params, file, 'RH', ysn, dlimb, dlima, 3, 2)
    percent_plot(vols, params, file, 'QV', ysn, dlimb, dlima, 3, 2)
    percent_plot(vols, params, file, 'U', ysn, dlimb, dlima, 3, 2)
    percent_plot(vols, params, file, 'V', ysn, dlimb, dlima, 3, 2)
    percent_plot(vols, params, file, 'T', ysn, dlimb, dlima, 3, 2)
    
    percent_plot(vols, params, file, 'OMEGA', ysn, dlimb, dlima, 3, 3)
    percent_plot(vols, params, file, 'RH', ysn, dlimb, dlima, 3, 3)
    percent_plot(vols, params, file, 'QV', ysn, dlimb, dlima, 3, 3)
    percent_plot(vols, params, file, 'U', ysn, dlimb, dlima, 3, 3)
    percent_plot(vols, params, file, 'V', ysn, dlimb, dlima, 3, 3)
    percent_plot(vols, params, file, 'T', ysn, dlimb, dlima, 3, 3)
    ''
    fig, axs = plt.subplots(2, 6, figsize=(30,10), sharey='row')
    #areas[areas > 1e8] = np.nan
    #areas = [np.nan for num in areas if num > 1e8]
    plot_stats(axs, 0, 0, vavg_lats, qv2ms_avg, qv2mp_avg, qv2me_avg, 'QV2M Avg (kg kg$^{-1}$)')
    plot_stats(axs, 0, 1, vavg_lats, q250s_avg, q250p_avg, q250e_avg, 'Q250 Avg (kg kg$^{-1}$)')
    plot_stats(axs, 0, 2, vavg_lats, q500s_avg, q500p_avg, q500e_avg, 'Q500 Avg (kg kg$^{-1}$)')
    plot_stats(axs, 0, 3, vavg_lats, q850s_avg, q850p_avg, q850e_avg, 'Q850 Avg (kg kg$^{-1}$)')
    plot_stats(axs, 0, 4, vavg_lats, u500s_avg, u500p_avg, u500e_avg, 'U500 Avg (kg kg$^{-1}$)')
    plot_stats(axs, 0, 5, vavg_lats, v500s_avg, v500p_avg, v500e_avg, 'V500 Avg (kg kg$^{-1}$)')
    axs[0,0].legend(fontsize=12)
    axs[0,0].set_ylabel('V$_{avg,lats}$ (km s$^{-1}$)')
    plot_stats(axs, 1, 0, vavg_lons, qv2ms_avg, qv2mp_avg, qv2me_avg, 'QV2M Avg (kg kg$^{-1}$)')
    plot_stats(axs, 1, 1, vavg_lons, q250s_avg, q250p_avg, q250e_avg, 'Q250 Avg (kg kg$^{-1}$)')
    plot_stats(axs, 1, 2, vavg_lons, q500s_avg, q500p_avg, q500e_avg, 'Q500 Avg (kg kg$^{-1}$)')
    plot_stats(axs, 1, 3, vavg_lons, q850s_avg, q850p_avg, q850e_avg, 'Q850 Avg (kg kg$^{-1}$)')
    plot_stats(axs, 1, 4, vavg_lons, u500s_avg, u500p_avg, u500e_avg, 'U500 Avg (kg kg$^{-1}$)')
    plot_stats(axs, 1, 5, vavg_lons, v500s_avg, v500p_avg, v500e_avg, 'V500 Avg (kg kg$^{-1}$)')
    axs[1,0].set_ylabel('V$_{avg,lons}$ (km s$^{-1}$)')
    plt.show()
    fig.savefig(f'{ysn}_scatplots_dlimb{dlimb}_dlima{dlima}_vels.png')
    '''
    fig1, axs1 = plt.subplots(2, 8, figsize=(50,20), sharey='row')
    #areas[areas > 1e8] = np.nan
    #areas = [np.nan for num in areas if num > 1e8]
    plot_stats(axs1, 0, 0, areas, qv2mi_avg, qv2ms_avg, qv2mp_avg, qv2me_avg, qv2mt_avg, 'QV2M Avg (kg kg$^{-1}$)')
    plot_stats(axs1, 0, 1, areas, q250i_avg, q250s_avg, q250p_avg, q250e_avg, q250t_avg, 'Q250 Avg (kg kg$^{-1}$)')
    plot_stats(axs1, 0, 2, areas, q500i_avg, q500s_avg, q500p_avg, q500e_avg, q500t_avg, 'Q500 Avg (kg kg$^{-1}$)')
    plot_stats(axs1, 0, 3, areas, q850i_avg, q850s_avg, q850p_avg, q850e_avg, q850t_avg, 'Q850 Avg (kg kg$^{-1}$)')
    plot_stats(axs1, 0, 4, areas, u500i_avg, u500s_avg, u500p_avg, u500e_avg, u500t_avg, 'U500 Avg (kg kg$^{-1}$)')
    plot_stats(axs1, 0, 5, areas, v500i_avg, v500s_avg, v500p_avg, v500e_avg, v500t_avg, 'V500 Avg (kg kg$^{-1}$)')
    plot_stats(axs1, 0, 6, areas, tqvi_avg, tqvs_avg, tqvp_avg, tqve_avg, tqvt_avg, 'TQV Avg')
    plot_stats(axs1, 0, 7, areas, t2mi_avg, t2ms_avg, t2mp_avg, t2me_avg, t2mt_avg, 'T2M Avg')
    axs1[0,0].legend(fontsize=12)
    axs1[0,0].set_ylabel('Area (km$^2$)')
    plot_stats(axs1, 1, 0, vols, qv2mi_avg, qv2ms_avg, qv2mp_avg, qv2me_avg, qv2mt_avg, 'QV2M Avg (kg kg$^{-1}$)')
    plot_stats(axs1, 1, 1, vols, q250i_avg, q250s_avg, q250p_avg, q250e_avg, q250t_avg, 'Q250 Avg (kg kg$^{-1}$)')
    plot_stats(axs1, 1, 2, vols, q500i_avg, q500s_avg, q500p_avg, q500e_avg, q500t_avg, 'Q500 Avg (kg kg$^{-1}$)')
    plot_stats(axs1, 1, 3, vols, q850i_avg, q850s_avg, q850p_avg, q850e_avg, q850t_avg, 'Q850 Avg (kg kg$^{-1}$)')
    plot_stats(axs1, 1, 4, vols, u500i_avg, u500s_avg, u500p_avg, u500e_avg, u500t_avg, 'U500 Avg (kg kg$^{-1}$)')
    plot_stats(axs1, 1, 5, vols, v500i_avg, v500s_avg, v500p_avg, v500e_avg, v500t_avg, 'V500 Avg (kg kg$^{-1}$)')
    plot_stats(axs1, 1, 6, vols, tqvi_avg, tqvs_avg, tqvp_avg, tqve_avg, tqvt_avg, 'TQV Avg')
    plot_stats(axs1, 1, 7, vols, t2mi_avg, t2ms_avg, t2mp_avg, t2me_avg, t2mt_avg, 'T2M Avg')
    axs1[1,0].set_ylabel('Volume (m$^3$)')
    plt.show()
    fig1.savefig(f'{ysn}_scatplots_dlimb{dlimb}_dlima{dlima}_areavol_new.png')
    '''
    sns_plot(ysn, dlimb, dlima, qv2mi_avg, q250i_avg, q500i_avg, u500i_avg, v500i_avg, tqvi_avg, t2mi_avg, q850i_avg, duras, areas, vols,'12hrs Before Start')
    sns_plot(ysn, dlimb, dlima, qv2ms_avg, q250s_avg, q500s_avg, u500s_avg, v500s_avg, tqvs_avg, t2ms_avg, q850s_avg, duras, areas, vols,'Start')
    sns_plot(ysn, dlimb, dlima, qv2mp_avg, q250p_avg, q500p_avg, u500p_avg, v500p_avg, tqvp_avg, t2mp_avg, q850p_avg, duras, areas, vols,'Peak')
    sns_plot(ysn, dlimb, dlima, qv2me_avg, q250e_avg, q500e_avg, u500e_avg, v500e_avg, tqve_avg, t2me_avg, q850e_avg, duras, areas, vols,'End')
    sns_plot(ysn, dlimb, dlima, qv2mt_avg, q250t_avg, q500t_avg, u500t_avg, v500t_avg, tqvt_avg, t2mt_avg, q850t_avg, duras, areas, vols,'12hrs After End')
    '''
    sns_plot(ysn, dlimb, dlima, qv2ms_avg, q250s_avg, q500s_avg, u500s_avg, v500s_avg, tqvs_avg, t2ms_avg, q850s_avg, 
                                qv2mp_avg, q250p_avg, q500p_avg, u500p_avg, v500p_avg, tqvp_avg, t2mp_avg, q850p_avg, 
                                qv2me_avg, q250e_avg, q500e_avg, u500e_avg, v500e_avg, tqve_avg, t2me_avg, q850e_avg, 
                                duras, areas, vols,'SPE')
    
#open_event_info('2016_JJA_var_avgs.nc4', 'extreme_event_stats_2016_JJA.nc')
    
    
def sns_plot(ysn, dlimb, dlima, qv2ms, q250s, q500s, u500s, v500s, tqvs, t2ms, q850s, 
                                qv2mp, q250p, q500p, u500p, v500p, tqvp, t2mp, q850p,
                                qv2me, q250e, q500e, u500e, v500e, tqve, t2me, q850e,
             dura, area, vol, t):
    import pandas as pd
    import seaborn as sns
    import matplotlib.pyplot as plt
    import numpy as np
    
    sns.set(style="ticks")
    #sns.set(style="ticks", color_codes=True)
    
    dfs = pd.DataFrame({'Duration': dura, 'QV2M': qv2ms,'Q250': q250s,'Q500': q500s, 
                       'U500': u500s,'V500': v500s, 'Area': area, 'Volume': vol, 
                       'TQV': tqvs, 'T2M': t2ms, 'Q850': q850s})
    
    dfp = pd.DataFrame({'Duration': dura, 'QV2M': qv2mp,'Q250': q250p,'Q500': q500p, 
                       'U500': u500p,'V500': v500p, 'Area': area, 'Volume': vol, 
                       'TQV': tqvp, 'T2M': t2mp, 'Q850': q850p})
    
    dfe = pd.DataFrame({'Duration': dura, 'QV2M': qv2me,'Q250': q250e,'Q500': q500e, 
                       'U500': u500e,'V500': v500e, 'Area': area, 'Volume': vol, 
                       'TQV': tqve, 'T2M': t2me, 'Q850': q850e})
    concat = pd.concat([dfs.assign(dataset='start'),dfp.assign(dataset='peak'),dfe.assign(dataset='end')])
    #concat = pd.concat([dfs,dfp,dfe], keys=['start','peak','end'],names=['time'])
    print(concat)
    '''
    pairs = sns.pairplot(df, height=3, x_vars=["Duration", "Volume", "Area"],
                 y_vars=["QV2M", "Q250", "Q500", "Q850", "U500", "V500", "TQV", "T2M"])
    
    pairs.axes[1,0].set_ylim([np.nanmin(q250)-np.nanstd(q250),np.nanmax(q250)+np.nanstd(q250)])
    pairs.axes[1,1].set_ylim([np.nanmin(q250)-np.nanstd(q250),np.nanmax(q250)+np.nanstd(q250)])
    pairs.axes[1,2].set_ylim([np.nanmin(q250)-np.nanstd(q250),np.nanmax(q250)+np.nanstd(q250)])
    
    plt.suptitle(f'{t} Time', y=0.998)
    plt.show()
    pairs.savefig(f'{ysn}_pairplots_{t}_dlimb{dlimb}_dlima{dlima}.png')
    '''
    pairc = sns.pairplot(concat, hue='dataset')
    plt.show()
    
    pairs = sns.pairplot(concat, hue='dataset',height=3, x_vars=["Duration", "Volume", "Area"],
                 y_vars=["QV2M", "Q250", "Q500", "Q850", "U500", "V500", "TQV", "T2M"])
    pairs.axes[1,0].set_ylim([np.nanmin(q250s)-np.nanstd(q250s),np.nanmax(q250s)+np.nanstd(q250s)])
    pairs.axes[1,1].set_ylim([np.nanmin(q250p)-np.nanstd(q250p),np.nanmax(q250p)+np.nanstd(q250p)])
    pairs.axes[1,2].set_ylim([np.nanmin(q250e)-np.nanstd(q250e),np.nanmax(q250e)+np.nanstd(q250e)])
    plt.suptitle(f'{t} Time', y=0.998)
    plt.show()
    
    pairc.savefig(f'{ysn}_pairplots_{t}_dlimb{dlimb}_dlima{dlima}.png')
    pairs.savefig(f'{ysn}_pairplot_{t}_dlimb{dlimb}_dlima{dlima}.png')