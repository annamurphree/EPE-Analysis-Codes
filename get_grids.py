"""
Created on Fri Jun 26 14:35:02 2020

@author: Anna Murphree (annammurphree@gmail.com)
For Summer 2020 Internship with Dr. Yaping Zhou

This code reads in EPE statistics from events created from IMERG data and outputs 
a new .nc4 file with grids of meteorological variables from MERRA-2 data files. 
The size of the grids depends on the event duration:
    0-6hrs:  5x5 degrees
    6-24hrs: 10x10 degrees
    24+hrs:  15x15 degrees
"""
# this makes a list of all the MERRA-2 data files for a season:
def files_list(core, yr, m0, m1, m2):
    #open the MERRA-2 data files in their directories:
    import os
    # the first month of the season:
    directory0 = os.fsencode(f'/discover/nobackup/projects/gmao/merra2/data/products/MERRA2_{core}/Y{yr}/M{m0}')
    # the second month of the season:
    directory1 = os.fsencode(f'/discover/nobackup/projects/gmao/merra2/data/products/MERRA2_{core}/Y{yr}/M{m0}')
    # the third month of the season:
    directory2 = os.fsencode(f'/discover/nobackup/projects/gmao/merra2/data/products/MERRA2_{core}/Y{yr}/M{m0}')
    # the 2D averaged files:
    files = [os.path.join(os.fsdecode(directory0), os.fsdecode(fi)) 
             for fi in os.listdir(directory0) if f'{core}.tavg1_2d_slv_Nx' in os.fsdecode(fi)]
    files1 = [os.path.join(os.fsdecode(directory1), os.fsdecode(fi)) 
             for fi in os.listdir(directory1) if f'{core}.tavg1_2d_slv_Nx' in os.fsdecode(fi)]
    files2 = [os.path.join(os.fsdecode(directory2), os.fsdecode(fi)) 
             for fi in os.listdir(directory2) if f'{core}.tavg1_2d_slv_Nx' in os.fsdecode(fi)]
    # the 3D instantaneous files:
    files3 = [os.path.join(os.fsdecode(directory0), os.fsdecode(fi)) 
             for fi in os.listdir(directory0) if f'{core}.inst3_3d_asm_Np' in os.fsdecode(fi)]
    files4 = [os.path.join(os.fsdecode(directory1), os.fsdecode(fi)) 
             for fi in os.listdir(directory1) if f'{core}.inst3_3d_asm_Np' in os.fsdecode(fi)]
    files5 = [os.path.join(os.fsdecode(directory2), os.fsdecode(fi)) 
             for fi in os.listdir(directory2) if f'{core}.inst3_3d_asm_Np' in os.fsdecode(fi)]
    # add these all into 1 list:
    files.extend(files1)
    files.extend(files2)
    files.extend(files3)
    files.extend(files4)
    files.extend(files5)
    
    # write this data into a new text file:
    #with open(f'{ysn}_files.txt', 'w') as filename:
    #    filename.writelines('%s\n' % file for file in files)
    
    return files

# this picks out the date from the MERRA-2 data file:
def pick_day(fi):
    from netCDF4 import Dataset
    import datetime
    fh = Dataset(fi, mode='r')
    date = fh.RangeBeginningDate
    dy = int(date[0:4])  # year
    dm = int(date[5:7])  # month
    dd = int(date[8:10]) # day
    day = datetime.datetime(dy, dm, dd).date()
    return day

# this fills in the meteorological data for a time step of the event:
def ft(v0, v1, v2, v3, v4, v5, v6, v7, lat_range, lon_range, index, t,
       nc_var0, nc_var1, nc_var2, nc_var3, nc_var4, nc_var5, nc_var6, nc_var7, plotsave='plot'):
    import numpy as np
    
    if plotsave == 'save':
        # add the averages to the holders in the .nc file:
        for la in lat_range:
            for lo in lon_range:
                nc_var0[index,t,la,lo] = v0[t,la,lo]
                nc_var1[index,t,la,lo] = v1[t,la,lo]
                nc_var2[index,t,la,lo] = v2[t,la,lo]
                nc_var3[index,t,la,lo] = v3[t,la,lo]
                nc_var4[index,t,la,lo] = v4[t,la,lo]
                nc_var5[index,t,la,lo] = v5[t,la,lo]
                nc_var6[index,t,la,lo] = v6[t,la,lo]
                nc_var7[index,t,la,lo] = v7[t,la,lo]
        
        return nc_var0, nc_var1, nc_var2, nc_var3, nc_var4, nc_var5, nc_var6, nc_var7
    
    elif plotsave =='plot':
        dv0 = []
        dv1 = []
        dv2 = []
        dv3 = []
        dv4 = []
        dv5 = []
        dv6 = []
        dv7 = []
        
        for la in lat_range:
            for lo in lon_range:
                dv0.append(v0[t,la,lo])
                dv1.append(v1[t,la,lo])
                dv2.append(v2[t,la,lo])
                dv3.append(v3[t,la,lo])
                dv4.append(v4[t,la,lo])
                dv5.append(v5[t,la,lo])
                dv6.append(v6[t,la,lo])
                dv7.append(v7[t,la,lo])
        ''
        dv0 = np.array(dv0)
        d0 = np.reshape(dv0, (len(lat_range), len(lon_range)))
        dv1 = np.array(dv1)
        d1 = np.reshape(dv1, (len(lat_range), len(lon_range)))
        dv2 = np.array(dv2)
        d2 = np.reshape(dv2, (len(lat_range), len(lon_range)))
        dv3 = np.array(dv3)
        d3 = np.reshape(dv3, (len(lat_range), len(lon_range)))
        dv4 = np.array(dv4)
        d4 = np.reshape(dv4, (len(lat_range), len(lon_range)))
        dv5 = np.array(dv5)
        d5 = np.reshape(dv5, (len(lat_range), len(lon_range)))
        dv6 = np.array(dv6)
        d6 = np.reshape(dv6, (len(lat_range), len(lon_range)))
        dv7 = np.array(dv7)
        d7 = np.reshape(dv7, (len(lat_range), len(lon_range)))
    
        return d0, d1, d2, d3, d4, d5, d6, d7

def ft3(v8, v9, v10, v11, v12, v13, lat_range, lon_range, index, t, l,
    nc_var0, nc_var1, nc_var2, nc_var3, nc_var4, nc_var5, plotsave='plot'):
    import numpy as np
    # levels: 21=250hPa, 16=500hPa, 6=850hPa
    
    if plotsave == 'save':
        # add the averages to the holders in the .nc file:
        for la in lat_range:
            for lo in lon_range:
                nc_var0[index,t,la,lo] = v8[t,la,lo]
                nc_var1[index,t,la,lo] = v9[t,la,lo]
                nc_var2[index,t,la,lo] = v10[t,la,lo]
                nc_var3[index,t,la,lo] = v11[t,la,lo]
                nc_var4[index,t,la,lo] = v12[t,la,lo]
                nc_var5[index,t,la,lo] = v13[t,la,lo]
        
        return nc_var0, nc_var1, nc_var2, nc_var3, nc_var4, nc_var5
    
    elif plotsave == 'plot':
        dv8 = []
        dv9 = []
        dv10 = []
        dv11 = []
        dv12 = []
        dv13 = []
        
        for la in lat_range:
            for lo in lon_range:
                dv8.append(v8[t,l,la,lo])
                dv9.append(v9[t,l,la,lo])
                dv10.append(v10[t,l,la,lo])
                dv11.append(v11[t,l,la,lo])
                dv12.append(v12[t,l,la,lo])
                dv13.append(v13[t,l,la,lo])
        ''
        dv8 = np.array(dv8)
        d8 = np.reshape(dv8, (len(lat_range), len(lon_range)))
        dv9 = np.array(dv9)
        d9 = np.reshape(dv9, (len(lat_range), len(lon_range)))
        dv10 = np.array(dv10)
        d10 = np.reshape(dv10, (len(lat_range), len(lon_range)))
        dv11 = np.array(dv11)
        d11 = np.reshape(dv11, (len(lat_range), len(lon_range)))
        dv12 = np.array(dv12)
        d12 = np.reshape(dv12, (len(lat_range), len(lon_range)))
        dv13 = np.array(dv13)
        d13 = np.reshape(dv13, (len(lat_range), len(lon_range)))
        
        return d8, d9, d10, d11, d12, d13

def map_plot(d0_i, d0_s, d0_p, d0_e, d0_t, m, fig, axs, col, v0, pick='n'):
    import numpy as np
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    #d0 = [d0_i, d0_s, d0_p, d0_e, d0_t]
    d0 = [d0_s, d0_p, d0_e]
    d0_min = np.min(d0) 
    d0_max = np.max(d0)
    if pick == 'n':
        m.ax = axs[0,col]
        axs[0,col].set_title(f'{v0}')
    elif pick == 'y':
        m.ax = axs[0]
        axs[0].set_title(f'{v0}')
    m.drawparallels(np.arange(-80., 81., 10.), labels=[1,0,0,0], fontsize=10)
    m.drawmeridians(np.arange(-180., 181., 10.), labels=[0,0,0,1], fontsize=10)
    m.drawcoastlines()
    m.drawstates()
    m.drawcountries()
    m.imshow(d0_i, origin='lower', vmin=d0_min, vmax=d0_max)
    if pick == 'n':
        m.ax = axs[1,col]
    elif pick == 'y':
        m.ax = axs[1]
    m.drawparallels(np.arange(-80., 81., 10.), labels=[1,0,0,0], fontsize=10)
    m.drawmeridians(np.arange(-180., 181., 10.), labels=[0,0,0,1], fontsize=10)
    m.drawcoastlines()
    m.drawstates()
    m.drawcountries()
    m.imshow(d0_s, origin='lower', vmin=d0_min, vmax=d0_max)
    if pick == 'n':
        m.ax = axs[2,col]
    elif pick == 'y':
        m.ax = axs[2]
    m.drawparallels(np.arange(-80., 81., 10.), labels=[1,0,0,0], fontsize=10)
    m.drawmeridians(np.arange(-180., 181., 10.), labels=[0,0,0,1], fontsize=10)
    m.drawcoastlines()
    m.drawstates()
    m.drawcountries()
    m.imshow(d0_p, origin='lower', vmin=d0_min, vmax=d0_max)
    if pick == 'n':
        m.ax = axs[3,col]
    elif pick == 'y':
        m.ax = axs[3]
    m.drawparallels(np.arange(-80., 81., 10.), labels=[1,0,0,0], fontsize=10)
    m.drawmeridians(np.arange(-180., 181., 10.), labels=[0,0,0,1], fontsize=10)
    m.drawcoastlines()
    m.drawstates()
    m.drawcountries()
    m.imshow(d0_e, origin='lower', vmin=d0_min, vmax=d0_max)
    if pick == 'n':
        m.ax = axs[4,col]
    elif pick == 'y':
        m.ax = axs[4]
    m.drawparallels(np.arange(-80., 81., 10.), labels=[1,0,0,0], fontsize=10)
    m.drawmeridians(np.arange(-180., 181., 10.), labels=[0,0,0,1], fontsize=10)
    m.drawcoastlines()
    m.drawstates()
    m.drawcountries()
    pic = m.imshow(d0_t, origin='lower', vmin=d0_min, vmax=d0_max)
    #cax = fig.add_axes([.13, 0.07, .08, .01])
    if pick == 'n':
        div = make_axes_locatable(axs[4,col])
    elif pick == 'y':
        div = make_axes_locatable(axs[4])
    
    cax = div.append_axes("bottom", size="5%", pad="15%")
    cbar = fig.colorbar(pic, cax=cax, orientation="horizontal", ticks=[d0_min, 0, d0_max])
    # stagger the colorbars so they don't overlap:
    if col % 2 == 0:
        cbar.ax.tick_params(pad=10, labelsize=9) 
    else:
        cbar.ax.tick_params(pad=15, labelsize=9) 

# this huge function plots 5 meteorological variables for each event
# each variable is mapped for the start, peak, and end times of the event
def plot_vars(v0, d0_i, d0_s, d0_p, d0_e, d0_t, v1, d1_i, d1_s, d1_p, d1_e, d1_t, v2, d2_i, d2_s, d2_p, d2_e, d2_t,
              v3, d3_i, d3_s, d3_p, d3_e, d3_t, v4, d4_i, d4_s, d4_p, d4_e, d4_t, v5, d5_i, d5_s, d5_p, d5_e, d5_t,
              lon, lat, lo, la, ind, ysn, ts, tp, te, lev=0, pickvar='all'):
    import matplotlib.pyplot as plt
    #from mpl_toolkits.axes_grid1 import make_axes_locatable
    from mpl_toolkits.basemap import Basemap
    #import numpy as np
    
    # create figure and axes:
    if pickvar == 'all':
        fig, axs = plt.subplots(5, 6, figsize=(15, 15), sharey='row', sharex='col')
    elif pickvar != 'all':
        fig, axs = plt.subplots(5, 1, figsize=(15, 15), sharey='row', sharex='col')
    plt.subplots_adjust(wspace=0.25, hspace=0.05)
    
    # setup stereographic basemap:
    # lon_0,lat_0 is the center point
    m = Basemap(llcrnrlon=(lo-len(lon)/2),llcrnrlat=(la-len(lat)/2),urcrnrlon=(lo+len(lon)/2),urcrnrlat=(la+len(lat)/2), 
                resolution='l',projection='stere',lat_0=la,lon_0=lo)
    if pickvar == 'all':
        map_plot(d0_i, d0_s, d0_p, d0_e, d0_t, m, fig, axs, 0, v0)
        map_plot(d1_i, d1_s, d1_p, d1_e, d1_t, m, fig, axs, 1, v1)
        map_plot(d2_i, d2_s, d2_p, d2_e, d2_t, m, fig, axs, 2, v2)
        map_plot(d3_i, d3_s, d3_p, d3_e, d3_t, m, fig, axs, 3, v3)
        map_plot(d4_i, d4_s, d4_p, d4_e, d4_t, m, fig, axs, 4, v4)
        map_plot(d5_i, d5_s, d5_p, d5_e, d5_t, m, fig, axs, 5, v5)
        
    elif pickvar == 0:
        map_plot(d0_i, d0_s, d0_p, d0_e, d0_t, m, fig, axs, 0, v0, 'y')
    elif pickvar == 1:
        map_plot(d1_i, d1_s, d1_p, d1_e, d1_t, m, fig, axs, 0, v1, 'y')
    elif pickvar == 2:
        map_plot(d2_i, d2_s, d2_p, d2_e, d2_t, m, fig, axs, 0, v2, 'y')
    elif pickvar == 3:
        map_plot(d3_i, d3_s, d3_p, d3_e, d3_t, m, fig, axs, 0, v3, 'y')
    elif pickvar == 4:
        map_plot(d4_i, d4_s, d4_p, d4_e, d4_t, m, fig, axs, 0, v4, 'y')
    elif pickvar == 5:
        map_plot(d5_i, d5_s, d5_p, d5_e, d5_t, m, fig, axs, 0, v5, 'y')
        
    #fig.text(0.09, 0.5, 'Latitude (degrees)', ha='center', va='center', rotation='vertical', fontsize=14)
    #fig.text(0.5, 0.07, 'Longitude (degrees)', ha='center', va='center', fontsize=14)
    fig.text(0.09, 0.83, '12hrs Before: ', ha='center', va='center', rotation='vertical', fontsize=12)
    fig.text(0.09, 0.65, 'Start Time: ', ha='center', va='center', rotation='vertical', fontsize=12)
    #fig.text(0.09, 0.77, f'{ts}', ha='center', va='center', rotation='vertical', fontsize=10)
    fig.text(0.09, 0.5, 'Peak Time: ', ha='center', va='center', rotation='vertical', fontsize=12)
    #fig.text(0.09, 0.5, f'{tp}', ha='center', va='center', rotation='vertical', fontsize=10)
    fig.text(0.09, 0.34, 'End Time: ', ha='center', va='center', rotation='vertical', fontsize=12)
    #fig.text(0.09, 0.23, f'{te}', ha='center', va='center', rotation='vertical', fontsize=10)
    fig.text(0.09, 0.18, '12hrs After: ', ha='center', va='center', rotation='vertical', fontsize=12)
    fig.text(0.5, 0.9, f'Event {ind}', ha='center', va='center', fontsize=12)
    
    #fig.tight_layout(h_pad=0.05, w_pad=0.35)
    #fig.tight_layout()
    fig.savefig(f'varmaps_{ysn}_event{ind}_level{lev}_{pickvar}.png', bbox_inches = "tight")
    plt.show()

def find_IMERG_in_MERRA2(imerg, core, yr, sn, output, ysn, pick, plotsave='plot', box=0, lev=0, levin=0, pickvar='all'):
    # imerg = the stats found from events in IMERG data
    # filelist = list of MERRA-2 data files created by files_list
    # output = new .nc4 file with each event's index, times, and averages
    # ysn = 'YEAR_SEA', like '2017_MAM'; just for labeling saved files
    # pick = index of event you want to look at
    # plotsave = either 'plot' or 'save' the event's data
    # box = size of box around the event's center lat/lon
    # lev = 3D level; in this code, 0=2D files, 1=250hPa, 2=500hPa, 3=850hPa
    # levin = 3D level index in MERRA-2 data; 21=250, 16=500, 6=850
    
    from netCDF4 import Dataset
    import numpy as np
    import datetime
    import os
    
    # open the event stats from IMERG:
    stats = Dataset(imerg, mode='r')
    params = stats.variables['event_parameters'][:] 
    
    sms = {0: ['12', '01', '02'], 1: ['03', '04', '05'], 2: ['06', '07', '08'], 3: ['09', '10', '11']}
    files = files_list(core, yr, sms[sn][0], sms[sn][1], sms[sn][2],)
    
    if plotsave == 'save':
        # open a new NetCDF4 file to write the data to:
        nc_out = Dataset(f'{output}', 'w', format='NETCDF4')   
        nc_out.createDimension('index', None)                 # create # of indices
        ind = nc_out.createVariable('index', 'i4', ('index',))
        
        # Using our previous dimension information, we can create the new dimensions
        nc_in = Dataset(files[0], mode='r')
        nc_dims = [dim for dim in nc_in.dimensions]
        data = {}
        # this sets the same dimensions as the MERRA-2 data files:
        for dim in nc_dims:
            nc_out.createDimension(dim, nc_in.variables[dim].size)
            data[dim] = nc_out.createVariable(dim, nc_in.variables[dim].dtype,\
                (dim,))
        
        # open a new NetCDF4 file to write the data to:
        nc_out = Dataset(f'{output}', 'w', format='NETCDF4')    
        nc_out.createDimension('time', 5)       # create 5 timesteps (prets, ts, tp, te, poste)  
        nc_out.createDimension('index', None)   # create # of indices
        ind = nc_out.createVariable('index', 'i4', ('index',))
        nc_out.createDimension('lev', 4)        # levels; 0=avg2d, 1=250hPa, 2=500hPa, 3=850hPa
        nc_out.createDimenstion('lat', None)
        nc_out.createDimenstion('lon', None)
        
        # 2D:
        nc_var0 = nc_out.createVariable('QV2M', 'f8', ('index', 'time', 'lev', 'lat', 'lon',))
        nc_var1 = nc_out.createVariable('Q250', 'f8', ('index', 'time', 'lev',  'lat', 'lon',))
        nc_var2 = nc_out.createVariable('Q500', 'f8', ('index', 'time',  'lev', 'lat', 'lon',))
        nc_var3 = nc_out.createVariable('Q850', 'f8', ('index', 'time',  'lev', 'lat', 'lon',))
        nc_var4 = nc_out.createVariable('U500', 'f8', ('index', 'time',  'lev', 'lat', 'lon',))
        nc_var5 = nc_out.createVariable('V500', 'f8', ('index', 'time',  'lev', 'lat', 'lon',))
        nc_var6 = nc_out.createVariable('TQV', 'f8', ('index', 'time',  'lev', 'lat', 'lon',))
        nc_var7 = nc_out.createVariable('T2M', 'f8', ('index', 'time',  'lev', 'lat', 'lon',))
        # 3D:
        nc_var8 = nc_out.createVariable('OMEGA', 'f8', ('index', 'time',  'lev', 'lat', 'lon',))
        nc_var9 = nc_out.createVariable('RH', 'f8', ('index', 'time',  'lev', 'lat', 'lon',))
        nc_var10 = nc_out.createVariable('QV', 'f8', ('index', 'time',  'lev', 'lat', 'lon',))
        nc_var11 = nc_out.createVariable('U', 'f8', ('index', 'time',  'lev', 'lat', 'lon',))
        nc_var12 = nc_out.createVariable('V', 'f8', ('index', 'time',  'lev', 'lat', 'lon',))
        nc_var13 = nc_out.createVariable('T', 'f8', ('index', 'time',  'lev', 'lat', 'lon',))
    
    # cycle through each event from IMERG:
    for event in params:
        dstart = event[0]       # start time (days) of the year
        dend = event[1]         # end time (days) of the year
        year = event[22]        # year of event
        year = int(year + 1)    # add 1 because the stats files are a year off
        long = event[2]         # IMERG longitude (deg East)
        lati = event[3]         # IMERG latitude (deg North)
        tmax = event[34]        # time into event of max rain volume (hours) 
        index = int(event[28])  # event index #
        
        if plotsave == 'save':
            ind[index] = index      # write this index # into the new .nc4 file
            
        dura = event[4]         # event duration (hrs)
        
        # convert start/end # days/year to datetime format:
        # datetime.datetime: (year, month, day, hour, minute, second, timezone info)
        # this format gives the full date and time:
        ds = datetime.datetime(year, 1, 1, 0, 0, 0) + datetime.timedelta(dstart - 1)
        #print('Start time: ',ds)
        dp = datetime.datetime(year, 1, 1, 0, 0, 0) + datetime.timedelta(dstart - 1, tmax*3600)
        #print('Peak time: ',dp)
        de = datetime.datetime(year, 1, 1, 0, 0, 0) + datetime.timedelta(dend - 1)
        #print('End time: ' ,de)
        
        dpre = ds - datetime.timedelta(hours=12)
        dpost = de + datetime.timedelta(hours=12)
        
        # this just gives the hours as integers, for indexing the MERRA-2 data:
        tstart = int((dstart - int(dstart))*24)  # start time
        dmax = dstart + (tmax/24)                
        tmax_hrs = int((dmax - int(dmax))*24)    # peak time
        tend = int((dend - int(dend))*24)        # end time
        
        prets = int(dpre.hour)         # 12hrs before tstart
        #print('Ts: ', ds, tstart)
        #print('Pre-Ts: ', dpre, prets)
        poste = int(dpost.hour)        # 12hrs after tend
        #print('Te: ', de, tend)
        #print('Post-Te: ', dpost, poste)
        
        # put each event day's file into a list:
        event_days = [fi for fi in files if (dpre.date() <= pick_day(fi) <= dpost.date())]
        #print(event_days)
        
        # can search for events here, by adding conditions:
        #if len(event_days) != 0 and (-96 <= long <= -90) and (36 <= lati <= 46):
        if len(event_days) != 0 and index == pick:
            #ind[index] = index      # write this index # into the new .nc4 file
            # if the MERRA-2 date is within the date range of the event:
            for dat in event_days:
                # 2D files:
                if '400.tavg1_2d_slv_Nx' in os.fsdecode(dat):
                    # print event info here, so it's only for the events matched:
                    print('2D Index: ',index)
                    
                    fh = Dataset(dat, mode='r')
                    #variables = fh.variables.keys()
                    #print(variables)
                    date = fh.RangeBeginningDate
                    dy = int(date[0:4])  # year
                    dm = int(date[5:7])  # month
                    dd = int(date[8:10]) # day
                    dy = datetime.datetime(dy, dm, dd).date()
                    #print('Event date:', dy, ds.date())
                    lon = fh.variables['lon'][:]    # longitude, 576 dimension: 0.625 deg
                    lat = fh.variables['lat'][:]    # latitude, 361 dimension: 0.5 deg
                    #time = fh.variables['time'][:] # in minutes, 0-23hrs
                    qv2m = fh.variables['QV2M'][:]  # 2m specific humidity
                    q250 = fh.variables['Q250'][:]  # 250m humidity
                    q500 = fh.variables['Q500'][:]  # 500m humidity
                    q850 = fh.variables['Q850'][:]  # 850m humidity 
                    u500 = fh.variables['U500'][:]  # eastward wind 500m
                    v500 = fh.variables['V500'][:]  # northward wind 500m
                    tqv = fh.variables['TQV'][:]    # total precipitable water vapor
                    t2m = fh.variables['T2M'][:]    # air temperature 2m
                    
                    # match up IMERG coordinates w/ MERRA-2 coordinates:
                    lon = np.asarray(lon)
                    close_lon = (np.abs(lon - long)).argmin()
                    elon = lon[close_lon]         # event longitude in MERRA-2
                    lat = np.asarray(lat)
                    close_lat = (np.abs(lat - lati)).argmin()
                    elat = lat[close_lat]         # event latitude in MERRA-2
                    
                    # if you don't specify how big the box is, default to these:
                    if box == 0:
                        if dura < 6:            # small events (<6hrs): box = 5x5deg
                            off = 2.5
                        elif 6 <= dura < 24:    # med events (6-24hrs): box = 10x10deg
                            off = 5
                        elif dura > 24:         # big events (>24hrs): box = 15x15deg
                            off = 7.5
                    # if you want to pick a specific box size:
                    elif box != 0:
                        off = box/2
                    # (box is the total side length, off is the offset around the center point)
                    lon0 = elon-off         # make a box around the event lat/lon
                    lon1 = elon+off
                    lat0 = elat-off
                    lat1 = elat+off
                    la_range = np.arange(int(lat0), int(lat1)) 
                    lo_range = np.arange(int(lon0), int(lon1))
                    #print('lons: ', lon0, elon, lon1)
                    #print('lats: ', lat0, elat, lat1)
                    lo0 = (8/5)*(lon0 + 180) + 1    # solve for the index where the lat/lon of the box corners are
                    la0 = 2*(lat0 + 90) + 1
                    lo1 = (8/5)*(lon1 + 180) + 1
                    la1 = 2*(lat1 + 90) + 1
                    #print('lon range: ', lo0, lo1)
                    #print('lat range: ', la0, la1)
                    lat_range = np.arange(int(la0), int(la1)) # now, the lat/lon ranges are indexing the lat/lon of the variables correctly
                    lon_range = np.arange(int(lo0), int(lo1))
                    
                    if plotsave == 'plot': 
                        # if it's the first day of the event, add the start time info
                        if dy == dpre.date():
                            (dv0i, dv1i, dv2i, dv3i, dv4i, dv5i, dv6i, dv7i) = ft(qv2m, q250, q500, q850, u500, v500, tqv, t2m, lat_range, lon_range, index, prets)
                            if dy == ds.date():
                                (dv0s, dv1s, dv2s, dv3s, dv4s, dv5s, dv6s, dv7s) = ft(qv2m, q250, q500, q850, u500, v500, tqv, t2m, lat_range, lon_range, index, tstart)
                                if dy == dp.date():
                                    (dv0p, dv1p, dv2p, dv3p, dv4p, dv5p, dv6p, dv7p) = ft(qv2m, q250, q500, q850, u500, v500, tqv, t2m, lat_range, lon_range, index, tmax_hrs)
                                    if dy == de.date():
                                        (dv0e, dv1e, dv2e, dv3e, dv4e, dv5e, dv6e, dv7e) = ft(qv2m, q250, q500, q850, u500, v500, tqv, t2m, lat_range, lon_range, index, tend)
                        elif dy == ds.date():
                            (dv0s, dv1s, dv2s, dv3s, dv4s, dv5s, dv6s, dv7s) = ft(qv2m, q250, q500, q850, u500, v500, tqv, t2m, lat_range, lon_range, index, tstart)
                            if dy == dp.date():
                                (dv0p, dv1p, dv2p, dv3p, dv4p, dv5p, dv6p, dv7p) = ft(qv2m, q250, q500, q850, u500, v500, tqv, t2m, lat_range, lon_range, index, tmax_hrs)
                                if dy == de.date():
                                    (dv0e, dv1e, dv2e, dv3e, dv4e, dv5e, dv6e, dv7e) = ft(qv2m, q250, q500, q850, u500, v500, tqv, t2m, lat_range, lon_range, index, tend)
                                    if dy == dpost.date():
                                        (dv0t, dv1t, dv2t, dv3t, dv4t, dv5t, dv6t, dv7t) = ft(qv2m, q250, q500, q850, u500, v500, tqv, t2m, lat_range, lon_range, index, poste)
                        # if it's the peak day of the event, add the peak time info
                        elif dy == dp.date():
                            (dv0p, dv1p, dv2p, dv3p, dv4p, dv5p, dv6p, dv7p) = ft(qv2m, q250, q500, q850, u500, v500, tqv, t2m, lat_range, lon_range, index, tmax_hrs)
                            if dy == de.date():
                                (dv0e, dv1e, dv2e, dv3e, dv4e, dv5e, dv6e, dv7e) = ft(qv2m, q250, q500, q850, u500, v500, tqv, t2m, lat_range, lon_range, index, tend)
                                if dy == dpost.date():
                                    (dv0t, dv1t, dv2t, dv3t, dv4t, dv5t, dv6t, dv7t) = ft(qv2m, q250, q500, q850, u500, v500, tqv, t2m, lat_range, lon_range, index, poste)
                        # if it's the last day of the event, add the end time info
                        elif dy == de.date():
                            (dv0e, dv1e, dv2e, dv3e, dv4e, dv5e, dv6e, dv7e) = ft(qv2m, q250, q500, q850, u500, v500, tqv, t2m, lat_range, lon_range, index, tend)
                            if dy == dpost.date():
                                (dv0t, dv1t, dv2t, dv3t, dv4t, dv5t, dv6t, dv7t) = ft(qv2m, q250, q500, q850, u500, v500, tqv, t2m, lat_range, lon_range, index, poste)
                        elif dy == dpost.date():
                            (dv0t, dv1t, dv2t, dv3t, dv4t, dv5t, dv6t, dv7t) = ft(qv2m, q250, q500, q850, u500, v500, tqv, t2m, lat_range, lon_range, index, poste)
                        # if it isn't the start, peak, or end day of the event, keep moving
                        else:
                            continue
                    elif plotsave == 'save':
                        # if it's the first day of the event, add the start time info
                        if dy == dpre.date():
                            (nc_var0[ind,0,0], nc_var1[ind,0,0], nc_var2[ind,0,0], nc_var3[ind,0,0], nc_var4[ind,0,0], nc_var5[ind,0,0], nc_var6[ind,0,0], nc_var7[ind,0,0]) = ft(qv2m, q250, q500, q850, u500, v500, tqv, t2m, lat_range, lon_range, index, prets)
                            if dy == ds.date():
                                (nc_var0[ind,1,0], nc_var1[ind,1,0], nc_var2[ind,1,0], nc_var3[ind,1,0], nc_var4[ind,1,0], nc_var5[ind,1,0], nc_var6[ind,1,0], nc_var7[ind,1,0]) = ft(qv2m, q250, q500, q850, u500, v500, tqv, t2m, lat_range, lon_range, index, tstart)
                                if dy == dp.date():
                                    (nc_var0[ind,2,0], nc_var1[ind,2,0], nc_var2[ind,2,0], nc_var3[ind,2,0], nc_var4[ind,2,0], nc_var5[ind,2,0], nc_var6[ind,2,0], nc_var7[ind,2,0]) = ft(qv2m, q250, q500, q850, u500, v500, tqv, t2m, lat_range, lon_range, index, tmax_hrs)
                                    if dy == de.date():
                                        (nc_var0[ind,3,0], nc_var1[ind,3,0], nc_var2[ind,3,0], nc_var3[ind,3,0], nc_var4[ind,3,0], nc_var5[ind,3,0], nc_var6[ind,3,0], nc_var7[ind,3,0]) = ft(qv2m, q250, q500, q850, u500, v500, tqv, t2m, lat_range, lon_range, index, tend)
                        elif dy == ds.date():
                            (nc_var0[ind,1,0], nc_var1[ind,1,0], nc_var2[ind,1,0], nc_var3[ind,1,0], nc_var4[ind,1,0], nc_var5[ind,1,0], nc_var6[ind,1,0], nc_var7[ind,1,0]) = ft(qv2m, q250, q500, q850, u500, v500, tqv, t2m, lat_range, lon_range, index, tstart)
                            if dy == dp.date():
                                (nc_var0[ind,2,0], nc_var1[ind,2,0], nc_var2[ind,2,0], nc_var3[ind,2,0], nc_var4[ind,2,0], nc_var5[ind,2,0], nc_var6[ind,2,0], nc_var7[ind,2,0]) = ft(qv2m, q250, q500, q850, u500, v500, tqv, t2m, lat_range, lon_range, index, tmax_hrs)
                                if dy == de.date():
                                    (nc_var0[ind,3,0], nc_var1[ind,3,0], nc_var2[ind,3,0], nc_var3[ind,3,0], nc_var4[ind,3,0], nc_var5[ind,3,0], nc_var6[ind,3,0], nc_var7[ind,3,0]) = ft(qv2m, q250, q500, q850, u500, v500, tqv, t2m, lat_range, lon_range, index, tend)
                                    if dy == dpost.date():
                                        (nc_var0[ind,4,0], nc_var1[ind,4,0], nc_var2[ind,4,0], nc_var3[ind,4,0], nc_var4[ind,4,0], nc_var5[ind,4,0], nc_var6[ind,4,0], nc_var7[ind,4,0]) = ft(qv2m, q250, q500, q850, u500, v500, tqv, t2m, lat_range, lon_range, index, poste)
                        # if it's the peak day of the event, add the peak time info
                        elif dy == dp.date():
                            (nc_var0[ind,2,0], nc_var1[ind,2,0], nc_var2[ind,2,0], nc_var3[ind,2,0], nc_var4[ind,2,0], nc_var5[ind,2,0], nc_var6[ind,2,0], nc_var7[ind,2,0]) = ft(qv2m, q250, q500, q850, u500, v500, tqv, t2m, lat_range, lon_range, index, tmax_hrs)
                            if dy == de.date():
                                (nc_var0[ind,3,0], nc_var1[ind,3,0], nc_var2[ind,3,0], nc_var3[ind,3,0], nc_var4[ind,3,0], nc_var5[ind,3,0], nc_var6[ind,3,0], nc_var7[ind,3,0]) = ft(qv2m, q250, q500, q850, u500, v500, tqv, t2m, lat_range, lon_range, index, tend)
                                if dy == dpost.date():
                                    (nc_var0[ind,4,0], nc_var1[ind,4,0], nc_var2[ind,4,0], nc_var3[ind,4,0], nc_var4[ind,4,0], nc_var5[ind,4,0], nc_var6[ind,4,0], nc_var7[ind,4,0]) = ft(qv2m, q250, q500, q850, u500, v500, tqv, t2m, lat_range, lon_range, index, poste)
                        # if it's the last day of the event, add the end time info
                        elif dy == de.date():
                            (nc_var0[ind,3,0], nc_var1[ind,3,0], nc_var2[ind,3,0], nc_var3[ind,3,0], nc_var4[ind,3,0], nc_var5[ind,3,0], nc_var6[ind,3,0], nc_var7[ind,3,0]) = ft(qv2m, q250, q500, q850, u500, v500, tqv, t2m, lat_range, lon_range, index, tend)
                            if dy == dpost.date():
                                (nc_var0[ind,4,0], nc_var1[ind,4,0], nc_var2[ind,4,0], nc_var3[ind,4,0], nc_var4[ind,4,0], nc_var5[ind,4,0], nc_var6[ind,4,0], nc_var7[ind,4,0]) = ft(qv2m, q250, q500, q850, u500, v500, tqv, t2m, lat_range, lon_range, index, poste)
                        elif dy == dpost.date():
                            (nc_var0[ind,4,0], nc_var1[ind,4,0], nc_var2[ind,4,0], nc_var3[ind,4,0], nc_var4[ind,4,0], nc_var5[ind,4,0], nc_var6[ind,4,0], nc_var7[ind,4,0]) = ft(qv2m, q250, q500, q850, u500, v500, tqv, t2m, lat_range, lon_range, index, poste)
                        # if it isn't the start, peak, or end day of the event, keep moving
                        else:
                            continue
                # 3D files:
                elif '400.inst3_3d_asm_Np' in os.fsdecode(dat):
                    # print event info here, so it's only for the events matched:
                    print('3D Index: ',index)
                    
                    fh = Dataset(dat, mode='r')
                    #variables = fh.variables.keys()
                    #print(variables)
                    date = fh.RangeBeginningDate
                    dy = int(date[0:4])  # year
                    dm = int(date[5:7])  # month
                    dd = int(date[8:10]) # day
                    dy = datetime.datetime(dy, dm, dd).date()
                    #print('Event date:', dy, ds.date())
                    lon = fh.variables['lon'][:]    # longitude, 576 dimension: 0.625 deg
                    lat = fh.variables['lat'][:]    # latitude, 361 dimension: 0.5 deg
                    #time = fh.variables['time'][:] # in minutes, 0-23hrs
                    #lev = fh.variables['lev'][:]    # levels, 0-42; 21=250hPa, 16=500hPa, 6=850hPa
                    # variables: (time, lev, lat, lon)
                    vv = fh.variables['OMEGA'][:]   # vertical pressure velocity
                    rh = fh.variables['RH'][:]      # relative humidity after moist
                    qv = fh.variables['QV'][:]      # specific humidity
                    u = fh.variables['U'][:]        # eastward wind
                    v = fh.variables['V'][:]        # northward wind
                    temp = fh.variables['T'][:]     # temperature
                    
                    # match up IMERG coordinates w/ MERRA-2 coordinates:
                    lon = np.asarray(lon)
                    close_lon = (np.abs(lon - long)).argmin()
                    elon = lon[close_lon]         # event longitude in MERRA-2
                    lat = np.asarray(lat)
                    close_lat = (np.abs(lat - lati)).argmin()
                    elat = lat[close_lat]         # event latitude in MERRA-2
                    
                    #print('times: ',prets, tstart, tmax_hrs, tend, poste)
                    # need to match up times: 2D/3D are different 
                    # 2D times: 0-23hrs
                    # 3D times: 0, 3, 6, 9, 12, 15, 18, 21 hrs
                    #  indices: 0, 1, 2, 3, 4,  5,  6,  7
                    # I'll divide the matched t3d time by 3 to get the right index:
                    t3d = [0,3,6,9,12,15,18,21]
                    t3d = np.asarray(t3d)
                    close_pts = (np.abs(t3d - prets)).argmin()
                    prets = int(t3d[close_pts]/3)
                    #print('PreTs: ',prets)
                    close_ts = (np.abs(t3d - tstart)).argmin()
                    tstart = int(t3d[close_ts]/3)
                    #print('Ts: ',tstart)
                    close_tm = (np.abs(t3d - tmax_hrs)).argmin()
                    tmax_hrs = int(t3d[close_tm]/3)
                    #print('Tm: ',tmax_hrs)
                    close_te = (np.abs(t3d - tend)).argmin()
                    tend = int(t3d[close_te]/3)
                    #print('Te: ',tend)
                    close_tp = (np.abs(t3d - poste)).argmin()
                    poste = int(t3d[close_tp]/3)
                    #print('PosTe: ',poste)
                    
                    # if you don't specify how big the box is, default to these:
                    if box == 0:
                        if dura < 6:            # small events (<6hrs): box = 5x5deg
                            off = 2.5
                        elif 6 <= dura < 24:    # med events (6-24hrs): box = 10x10deg
                            off = 5
                        elif dura > 24:         # big events (>24hrs): box = 15x15deg
                            off = 7.5
                    # if you want to pick a specific box size:
                    elif box != 0:
                        off = box/2
                    # (box is the total side length, off is the offset around the center point)
                        
                    lon0 = elon-off         # make a box around the event lat/lon
                    lon1 = elon+off
                    lat0 = elat-off
                    lat1 = elat+off
                    la_range = np.arange(int(lat0), int(lat1))
                    lo_range = np.arange(int(lon0), int(lon1))
                    #print('lons: ', lon0, elon, lon1)
                    #print('lats: ', lat0, elat, lat1)
                    lo0 = (8/5)*(lon0 + 180) + 1    # solve for the index where the lat/lon of the box corners are
                    la0 = 2*(lat0 + 90) + 1
                    lo1 = (8/5)*(lon1 + 180) + 1
                    la1 = 2*(lat1 + 90) + 1
                    
                    #print('lon range: ', lo0, lo1)
                    #print('lat range: ', la0, la1)
                    lat_range = np.arange(int(la0), int(la1)) # now, the lat/lon ranges are indexing the lat/lon of the variables correctly
                    lon_range = np.arange(int(lo0), int(lo1))
                    
                    # levels: 21=250hPa, 16=500hPa, 6=850hPa
                    
                    if plotsave == 'plot':
                        # if it's the first day of the event, add the start time info
                        if dy == dpre.date():
                            (d3v0i, d3v1i, d3v2i, d3v3i, d3v4i, d3v5i) = ft3(vv, rh, qv, u, v, temp, lat_range, lon_range, index, prets, levin)
                            if dy == ds.date():
                                (d3v0s, d3v1s, d3v2s, d3v3s, d3v4s, d3v5s) = ft3(vv, rh, qv, u, v, temp, lat_range, lon_range, index, tstart, levin)
                                if dy == dp.date():
                                    (d3v0p, d3v1p, d3v2p, d3v3p, d3v4p, d3v5p) = ft3(vv, rh, qv, u, v, temp, lat_range, lon_range, index, tmax_hrs, levin)
                                    if dy == de.date():
                                        (d3v0e, d3v1e, d3v2e, d3v3e, d3v4e, d3v5e) = ft3(vv, rh, qv, u, v, temp, lat_range, lon_range, index, tend, levin)
                        elif dy == ds.date():
                            (d3v0s, d3v1s, d3v2s, d3v3s, d3v4s, d3v5s) = ft3(vv, rh, qv, u, v, temp, lat_range, lon_range, index, tstart, levin)
                            if dy == dp.date():
                                (d3v0p, d3v1p, d3v2p, d3v3p, d3v4p, d3v5p) = ft3(vv, rh, qv, u, v, temp, lat_range, lon_range, index, tmax_hrs, levin)
                                if dy == de.date():
                                    (d3v0e, d3v1e, d3v2e, d3v3e, d3v4e, d3v5e) = ft3(vv, rh, qv, u, v, temp, lat_range, lon_range, index, tend, levin)
                                    if dy == dpost.date():
                                        (d3v0t, d3v1t, d3v2t, d3v3t, d3v4t, d3v5t) = ft3(vv, rh, qv, u, v, temp, lat_range, lon_range, index, poste, levin)
                        # if it's the peak day of the event, add the peak time info
                        elif dy == dp.date():
                            (d3v0p, d3v1p, d3v2p, d3v3p, d3v4p, d3v5p) = ft3(vv, rh, qv, u, v, temp, lat_range, lon_range, index, tmax_hrs, levin)
                            if dy == de.date():
                                (d3v0e, d3v1e, d3v2e, d3v3e, d3v4e, d3v5e) = ft3(vv, rh, qv, u, v, temp, lat_range, lon_range, index, tend, levin)
                                if dy == dpost.date():
                                    (d3v0t, d3v1t, d3v2t, d3v3t, d3v4t, d3v5t) = ft3(vv, rh, qv, u, v, temp, lat_range, lon_range, index, poste, levin)
                        # if it's the last day of the event, add the end time info
                        elif dy == de.date():
                            (d3v0e, d3v1e, d3v2e, d3v3e, d3v4e, d3v5e) = ft3(vv, rh, qv, u, v, temp, lat_range, lon_range, index, tend, levin)
                            if dy == dpost.date():
                                (d3v0t, d3v1t, d3v2t, d3v3t, d3v4t, d3v5t) = ft3(vv, rh, qv, u, v, temp, lat_range, lon_range, index, poste, levin)
                        elif dy == dpost.date():
                            (d3v0t, d3v1t, d3v2t, d3v3t, d3v4t, d3v5t) = ft3(vv, rh, qv, u, v, temp, lat_range, lon_range, index, poste, levin)
                        # if it isn't the start, peak, or end day of the event, keep moving
                        else:
                            continue
                    elif plotsave == 'save':
                        # if it's the first day of the event, add the start time info
                        if dy == dpre.date():
                            (nc_var0[ind,0,0], nc_var1[ind,0,0], nc_var2[ind,0,0], nc_var3[ind,0,0], nc_var4[ind,0,0], nc_var5[ind,0,0], nc_var6[ind,0,0], nc_var7[ind,0,0]) = ft3(qv2m, q250, q500, q850, u500, v500, tqv, t2m, lat_range, lon_range, index, prets)
                            if dy == ds.date():
                                (nc_var0[ind,1,0], nc_var1[ind,1,0], nc_var2[ind,1,0], nc_var3[ind,1,0], nc_var4[ind,1,0], nc_var5[ind,1,0], nc_var6[ind,1,0], nc_var7[ind,1,0]) = ft3(qv2m, q250, q500, q850, u500, v500, tqv, t2m, lat_range, lon_range, index, tstart)
                                if dy == dp.date():
                                    (nc_var0[ind,2,0], nc_var1[ind,2,0], nc_var2[ind,2,0], nc_var3[ind,2,0], nc_var4[ind,2,0], nc_var5[ind,2,0], nc_var6[ind,2,0], nc_var7[ind,2,0]) = ft3(qv2m, q250, q500, q850, u500, v500, tqv, t2m, lat_range, lon_range, index, tmax_hrs)
                                    if dy == de.date():
                                        (nc_var0[ind,3,0], nc_var1[ind,3,0], nc_var2[ind,3,0], nc_var3[ind,3,0], nc_var4[ind,3,0], nc_var5[ind,3,0], nc_var6[ind,3,0], nc_var7[ind,3,0]) = ft3(qv2m, q250, q500, q850, u500, v500, tqv, t2m, lat_range, lon_range, index, tend)
                        elif dy == ds.date():
                            (nc_var0[ind,1,0], nc_var1[ind,1,0], nc_var2[ind,1,0], nc_var3[ind,1,0], nc_var4[ind,1,0], nc_var5[ind,1,0], nc_var6[ind,1,0], nc_var7[ind,1,0]) = ft3(qv2m, q250, q500, q850, u500, v500, tqv, t2m, lat_range, lon_range, index, tstart)
                            if dy == dp.date():
                                (nc_var0[ind,2,0], nc_var1[ind,2,0], nc_var2[ind,2,0], nc_var3[ind,2,0], nc_var4[ind,2,0], nc_var5[ind,2,0], nc_var6[ind,2,0], nc_var7[ind,2,0]) = ft3(qv2m, q250, q500, q850, u500, v500, tqv, t2m, lat_range, lon_range, index, tmax_hrs)
                                if dy == de.date():
                                    (nc_var0[ind,3,0], nc_var1[ind,3,0], nc_var2[ind,3,0], nc_var3[ind,3,0], nc_var4[ind,3,0], nc_var5[ind,3,0], nc_var6[ind,3,0], nc_var7[ind,3,0]) = ft3(qv2m, q250, q500, q850, u500, v500, tqv, t2m, lat_range, lon_range, index, tend)
                                    if dy == dpost.date():
                                        (nc_var0[ind,4,0], nc_var1[ind,4,0], nc_var2[ind,4,0], nc_var3[ind,4,0], nc_var4[ind,4,0], nc_var5[ind,4,0], nc_var6[ind,4,0], nc_var7[ind,4,0]) = ft3(qv2m, q250, q500, q850, u500, v500, tqv, t2m, lat_range, lon_range, index, poste)
                        # if it's the peak day of the event, add the peak time info
                        elif dy == dp.date():
                            (nc_var0[ind,2,0], nc_var1[ind,2,0], nc_var2[ind,2,0], nc_var3[ind,2,0], nc_var4[ind,2,0], nc_var5[ind,2,0], nc_var6[ind,2,0], nc_var7[ind,2,0]) = ft3(qv2m, q250, q500, q850, u500, v500, tqv, t2m, lat_range, lon_range, index, tmax_hrs)
                            if dy == de.date():
                                (nc_var0[ind,3,0], nc_var1[ind,3,0], nc_var2[ind,3,0], nc_var3[ind,3,0], nc_var4[ind,3,0], nc_var5[ind,3,0], nc_var6[ind,3,0], nc_var7[ind,3,0]) = ft3(qv2m, q250, q500, q850, u500, v500, tqv, t2m, lat_range, lon_range, index, tend)
                                if dy == dpost.date():
                                    (nc_var0[ind,4,0], nc_var1[ind,4,0], nc_var2[ind,4,0], nc_var3[ind,4,0], nc_var4[ind,4,0], nc_var5[ind,4,0], nc_var6[ind,4,0], nc_var7[ind,4,0]) = ft3(qv2m, q250, q500, q850, u500, v500, tqv, t2m, lat_range, lon_range, index, poste)
                        # if it's the last day of the event, add the end time info
                        elif dy == de.date():
                            (nc_var0[ind,3,0], nc_var1[ind,3,0], nc_var2[ind,3,0], nc_var3[ind,3,0], nc_var4[ind,3,0], nc_var5[ind,3,0], nc_var6[ind,3,0], nc_var7[ind,3,0]) = ft3(qv2m, q250, q500, q850, u500, v500, tqv, t2m, lat_range, lon_range, index, tend)
                            if dy == dpost.date():
                                (nc_var0[ind,4,0], nc_var1[ind,4,0], nc_var2[ind,4,0], nc_var3[ind,4,0], nc_var4[ind,4,0], nc_var5[ind,4,0], nc_var6[ind,4,0], nc_var7[ind,4,0]) = ft3(qv2m, q250, q500, q850, u500, v500, tqv, t2m, lat_range, lon_range, index, poste)
                        elif dy == dpost.date():
                            (nc_var0[ind,4,0], nc_var1[ind,4,0], nc_var2[ind,4,0], nc_var3[ind,4,0], nc_var4[ind,4,0], nc_var5[ind,4,0], nc_var6[ind,4,0], nc_var7[ind,4,0]) = ft3(qv2m, q250, q500, q850, u500, v500, tqv, t2m, lat_range, lon_range, index, poste)
                        # if it isn't the start, peak, or end day of the event, keep moving
                        else:
                            continue
                ''
                fh.close()
                
        elif len(event_days) != 0 and index <= pick:
            continue
        else: 
            break
        
    if lev == 0:       #2D data 
        plot_vars('QV2M', dv0i, dv0s, dv0p, dv0e, dv0t, 'Q850', dv3i, dv3s, dv3p, dv3e, dv3t, 'Q500', dv2i, dv2s, dv2p, dv2e, dv2t, 
                  'Q250', dv1i, dv1s, dv1p, dv1e, dv1t, 'U500', dv4i, dv4s, dv4p, dv4e, dv4t, 'V500', dv5i, dv5s, dv5p, dv5e, dv5t,
                  lo_range, la_range, elon, elat, index-1, ysn, tstart, tmax_hrs, tend, pickvar=pickvar)
        # other variables you could substitute in:
        # 'T2M', dv7i, dv7s, dv7p, dv7e, dv7t
        # 'TQV', dv6i, dv6s, dv6p, dv6e, dv6t
        
    elif lev != 0:     #3D data
        plot_vars('OMEGA', d3v0i, d3v0s, d3v0p, d3v0e, d3v0t, 'RH', d3v1i, d3v1s, d3v1p, d3v1e, d3v1t, 'QV', d3v2i, d3v2s, d3v2p, d3v2e, d3v2t,
                  'U', d3v3i, d3v3s, d3v3p, d3v3e, d3v3t, 'V', d3v4i, d3v4s, d3v4p, d3v4e, d3v4t, 'T', d3v5i, d3v5s, d3v5p, d3v5e, d3v5t,
                  lo_range, la_range, elon, elat, index-1, ysn, ds, dp, de, lev, pickvar)
            
    #nc_out.close()  # close the new file
