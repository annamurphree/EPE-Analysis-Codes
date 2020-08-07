"""
Created on Fri Jun 26 14:35:02 2020

@author: Anna Murphree (annammurphree@gmail.com)
For Summer 2020 Internship with Dr. Yaping Zhou

This code reads in EPE statistics from events created from IMERG data and outputs 
a new .nc4 file with averages and standard deviations of grids of meteorological variables 
from MERRA-2 data files. 
The size of the collected grids depends on the event duration:
    0-6hrs:  5x5 degrees
    6-24hrs: 10x10 degrees
    24+hrs:  15x15 degrees
"""
# this makes a list of all the MERRA-2 data files for a season:
def files_list(ym0, ym1, ym2, ysn):
    # open the MERRA-2 data files in their directories:
    import os
    # the first month of the season:
    directory0 = os.fsencode('/discover/nobackup/projects/gmao/merra2/data/products/MERRA2_400/%s'%(ym0))
    # the second month of the season:
    directory1 = os.fsencode('/discover/nobackup/projects/gmao/merra2/data/products/MERRA2_400/%s'%(ym1))
    # the third month of the season:
    directory2 = os.fsencode('/discover/nobackup/projects/gmao/merra2/data/products/MERRA2_400/%s'%(ym2))
    # the 2D averaged files:
    files = [os.path.join(os.fsdecode(directory0), os.fsdecode(fi)) 
             for fi in os.listdir(directory0) if '400.tavg1_2d_slv_Nx' in os.fsdecode(fi)]
    files1 = [os.path.join(os.fsdecode(directory1), os.fsdecode(fi)) 
             for fi in os.listdir(directory1) if '400.tavg1_2d_slv_Nx' in os.fsdecode(fi)]
    files2 = [os.path.join(os.fsdecode(directory2), os.fsdecode(fi)) 
             for fi in os.listdir(directory2) if '400.tavg1_2d_slv_Nx' in os.fsdecode(fi)]
    # the 3D instantaneous files:
    files3 = [os.path.join(os.fsdecode(directory0), os.fsdecode(fi)) 
             for fi in os.listdir(directory0) if '400.inst3_3d_asm_Np' in os.fsdecode(fi)]
    files4 = [os.path.join(os.fsdecode(directory1), os.fsdecode(fi)) 
             for fi in os.listdir(directory1) if '400.inst3_3d_asm_Np' in os.fsdecode(fi)]
    files5 = [os.path.join(os.fsdecode(directory2), os.fsdecode(fi)) 
             for fi in os.listdir(directory2) if '400.inst3_3d_asm_Np' in os.fsdecode(fi)]
    # add these all into 1 list:
    files.extend(files1)
    files.extend(files2)
    files.extend(files3)
    files.extend(files4)
    files.extend(files5)
    # write this data into a new text file:
    with open(f'{ysn}_files.txt', 'w') as filename:
        filename.writelines('%s\n' % file for file in files)

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

# this fills in the meteorological data for each time of the event, for 2D variables:
def ft(v0, v1, v2, v3, v4, v5, v6, v7, lat_range, lon_range, index, ti, to, 
        nc_var0, nc_var1, nc_var2, nc_var3, nc_var4, nc_var5, nc_var6, nc_var7):
    import numpy as np
    
    v0s = [v0[ti,la,lo] for la in lat_range for lo in lon_range]
    v1s = [v1[ti,la,lo] for la in lat_range for lo in lon_range]
    v2s = [v2[ti,la,lo] for la in lat_range for lo in lon_range]
    v3s = [v3[ti,la,lo] for la in lat_range for lo in lon_range]
    v4s = [v4[ti,la,lo] for la in lat_range for lo in lon_range]
    v5s = [v5[ti,la,lo] for la in lat_range for lo in lon_range]
    v6s = [v6[ti,la,lo] for la in lat_range for lo in lon_range]
    v7s = [v7[ti,la,lo] for la in lat_range for lo in lon_range]
    
    # add the averages to the holders in the .nc file:
    nc_var0[index,to,0,0] = np.nanmean(v0s)
    nc_var1[index,to,0,0] = np.nanmean(v1s)
    nc_var2[index,to,0,0] = np.nanmean(v2s)
    nc_var3[index,to,0,0] = np.nanmean(v3s)
    nc_var4[index,to,0,0] = np.nanmean(v4s)
    nc_var5[index,to,0,0] = np.nanmean(v5s)
    nc_var6[index,to,0,0] = np.nanmean(v6s)
    nc_var7[index,to,0,0] = np.nanmean(v7s)
    
    nc_var0[index,to,0,1] = np.nanstd(v0s)
    nc_var1[index,to,0,1] = np.nanstd(v1s)
    nc_var2[index,to,0,1] = np.nanstd(v2s)
    nc_var3[index,to,0,1] = np.nanstd(v3s)
    nc_var4[index,to,0,1] = np.nanstd(v4s)
    nc_var5[index,to,0,1] = np.nanstd(v5s)
    nc_var6[index,to,0,1] = np.nanstd(v6s)
    nc_var7[index,to,0,1] = np.nanstd(v7s)
    
    #return v0_(tn), v1s, v2s, v3s, v4s, v5s
    
def ft3d(v8, v9, v10, v11, v12, v13, lat_range, lon_range, index, ti, to,
        nc_var8, nc_var9, nc_var10, nc_var11, nc_var12, nc_var13):
    import numpy as np
    # levels: 21=250hPa, 16=500hPa, 6=850hPa
    # level: 250hPa
    v82 = [v8[ti,21,la,lo] for la in lat_range for lo in lon_range]
    v92 = [v9[ti,21,la,lo] for la in lat_range for lo in lon_range]
    v102 = [v10[ti,21,la,lo] for la in lat_range for lo in lon_range]
    v112 = [v11[ti,21,la,lo] for la in lat_range for lo in lon_range]
    v122 = [v12[ti,21,la,lo] for la in lat_range for lo in lon_range]
    v132 = [v13[ti,21,la,lo] for la in lat_range for lo in lon_range]
    # add the averages & stddevs to the holders in the .nc file:
    nc_var8[index,to,1,0] = np.nanmean(v82)
    nc_var9[index,to,1,0] = np.nanmean(v92)
    nc_var10[index,to,1,0] = np.nanmean(v102)
    nc_var11[index,to,1,0] = np.nanmean(v112)
    nc_var12[index,to,1,0] = np.nanmean(v122)
    nc_var13[index,to,1,0] = np.nanmean(v132)
    nc_var8[index,to,1,1] = np.nanstd(v82)
    nc_var9[index,to,1,1] = np.nanstd(v92)
    nc_var10[index,to,1,1] = np.nanstd(v102)
    nc_var11[index,to,1,1] = np.nanstd(v112)
    nc_var12[index,to,1,1] = np.nanstd(v122)
    nc_var13[index,to,1,1] = np.nanstd(v132)
    # level: 500hPa
    v85 = [v8[ti,16,la,lo] for la in lat_range for lo in lon_range]
    v95 = [v9[ti,16,la,lo] for la in lat_range for lo in lon_range]
    v105 = [v10[ti,16,la,lo] for la in lat_range for lo in lon_range]
    v115 = [v11[ti,16,la,lo] for la in lat_range for lo in lon_range]
    v125 = [v12[ti,16,la,lo] for la in lat_range for lo in lon_range]
    v135 = [v13[ti,16,la,lo] for la in lat_range for lo in lon_range]
    # add the averages & stddevs to the holders in the .nc file:
    nc_var8[index,to,2,0] = np.nanmean(v85)
    nc_var9[index,to,2,0] = np.nanmean(v95)
    nc_var10[index,to,2,0] = np.nanmean(v105)
    nc_var11[index,to,2,0] = np.nanmean(v115)
    nc_var12[index,to,2,0] = np.nanmean(v125)
    nc_var13[index,to,2,0] = np.nanmean(v135)
    nc_var8[index,to,2,1] = np.nanstd(v85)
    nc_var9[index,to,2,1] = np.nanstd(v95)
    nc_var10[index,to,2,1] = np.nanstd(v105)
    nc_var11[index,to,2,1] = np.nanstd(v115)
    nc_var12[index,to,2,1] = np.nanstd(v125)
    nc_var13[index,to,2,1] = np.nanstd(v135)
    # level: 850hPa
    v88 = [v8[ti,6,la,lo] for la in lat_range for lo in lon_range]
    v98 = [v9[ti,6,la,lo] for la in lat_range for lo in lon_range]
    v108 = [v10[ti,6,la,lo] for la in lat_range for lo in lon_range]
    v118 = [v11[ti,6,la,lo] for la in lat_range for lo in lon_range]
    v128 = [v12[ti,6,la,lo] for la in lat_range for lo in lon_range]
    v138 = [v13[ti,6,la,lo] for la in lat_range for lo in lon_range]
    # add the averages & stddevs to the holders in the .nc file:
    nc_var8[index,to,3,0] = np.nanmean(v88)
    nc_var9[index,to,3,0] = np.nanmean(v98)
    nc_var10[index,to,3,0] = np.nanmean(v108)
    nc_var11[index,to,3,0] = np.nanmean(v118)
    nc_var12[index,to,3,0] = np.nanmean(v128)
    nc_var13[index,to,3,0] = np.nanmean(v138)
    nc_var8[index,to,3,1] = np.nanstd(v88)
    nc_var9[index,to,3,1] = np.nanstd(v98)
    nc_var10[index,to,3,1] = np.nanstd(v108)
    nc_var11[index,to,3,1] = np.nanstd(v118)
    nc_var12[index,to,3,1] = np.nanstd(v128)
    nc_var13[index,to,3,1] = np.nanstd(v138)
    
    #return v0_(tn), v1s, v2s, v3s, v4s, v5s

def find_IMERG_in_MERRA2(imerg, filelist, output):
    # imerg = the stats found from events in IMERG data
    # filelist = list of MERRA-2 data files created by files_list
    # output = new .nc4 file with each event's index, times, and averages
    
    from netCDF4 import Dataset
    import numpy as np
    import datetime
    import os
    
    # open the event stats from IMERG:
    stats = Dataset(imerg, mode='r')
    params = stats.variables['event_parameters'][:]
    
    # open the list of MERRA-2 data files created by files_list:
    f = open(filelist).readlines()
    files = [(i.split()[0]) for i in f][:]        
    
    # open a new NetCDF4 file to write the data to:
    nc_out = Dataset(f'{output}', 'w', format='NETCDF4')    
    nc_out.createDimension('time', 5)       # create 5 timesteps (prets, ts, tp, te, poste)  
    nc_out.createDimension('index', None)   # create # of indices
    ind = nc_out.createVariable('index', 'i4', ('index',))
    nc_out.createDimension('stats', 2)      # statistics; 0=avg, 1=stddev
    nc_out.createDimension('lev', 4)        # levels; 0=avg2d, 1=250hPa, 2=500hPa, 3=850hPa
    
    # create our departure variables:
    # 2D:
    nc_var0 = nc_out.createVariable('QV2M', 'f8', ('index', 'time', 'lev', 'stats',))
    nc_var1 = nc_out.createVariable('Q250', 'f8', ('index', 'time', 'lev',  'stats',))
    nc_var2 = nc_out.createVariable('Q500', 'f8', ('index', 'time',  'lev', 'stats',))
    nc_var3 = nc_out.createVariable('Q850', 'f8', ('index', 'time',  'lev', 'stats',))
    nc_var4 = nc_out.createVariable('U500', 'f8', ('index', 'time',  'lev', 'stats',))
    nc_var5 = nc_out.createVariable('V500', 'f8', ('index', 'time',  'lev', 'stats',))
    nc_var6 = nc_out.createVariable('TQV', 'f8', ('index', 'time',  'lev', 'stats',))
    nc_var7 = nc_out.createVariable('T2M', 'f8', ('index', 'time',  'lev', 'stats',))
    # 3D:
    nc_var8 = nc_out.createVariable('OMEGA', 'f8', ('index', 'time',  'lev', 'stats',))
    nc_var9 = nc_out.createVariable('RH', 'f8', ('index', 'time',  'lev', 'stats',))
    nc_var10 = nc_out.createVariable('QV', 'f8', ('index', 'time',  'lev', 'stats',))
    nc_var11 = nc_out.createVariable('U', 'f8', ('index', 'time',  'lev', 'stats',))
    nc_var12 = nc_out.createVariable('V', 'f8', ('index', 'time',  'lev', 'stats',))
    nc_var13 = nc_out.createVariable('T', 'f8', ('index', 'time',  'lev', 'stats',))
    
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
        if len(event_days) != 0:
            
            # if the MERRA-2 date is within the date range of the event:
            for dat in event_days:
                # 2D files:
                if '400.tavg1_2d_slv_Nx' in os.fsdecode(dat):
                    # print event info here, so it's only for the events matched:
                    print('2D Index: ',index)
                    fh = Dataset(dat, mode='r')
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
                    
                    '''
                    if dura < 6:            # small events (<6hrs): box = 5x5deg
                        lat_range = np.array(list(np.arange(int(np.round(elat-2.5)),int(np.round(elat+2.5)))))
                        lon_range = np.array(list(np.arange(int(np.round(elon-2.5)),int(np.round(elon+2.5)))))
                    elif 6 <= dura < 24:    # med events (6-24hrs): box = 10x10deg
                        lat_range = np.array(list(np.arange(int(np.round(elat-5)),int(np.round(elat+5)))))
                        lon_range = np.array(list(np.arange(int(np.round(elon-5)),int(np.round(elon+5)))))
                    elif dura > 24:         # big events (>24hrs): box = 15x15deg
                        lat_range = np.array(list(np.arange(int(np.round(elat-7.5)),int(np.round(elat+7.5)))))
                        lon_range = np.array(list(np.arange(int(np.round(elon-7.5)),int(np.round(elon+7.5)))))
                    ''
                    if dura < 6:            # small events (<6hrs): box = 5x5deg
                        lon0 = elon-2.5
                        lon1 = elon+2.5
                        lat0 = elat-2.5
                        lat1 = elat+2.5
                        lo0 = (8/5)*(lon0 + 180) + 1    # solve for the index where the lat/lon of the box corners are
                        la0 = 2*(lat0 + 90) + 1
                        lo1 = (8/5)*(lon1 + 180) + 1
                        la1 = 2*(lat1 + 90) + 1
                        lat_range = np.arange(la0, la1)
                        lon_range = np.arange(lo0, lo1)
                        #lat_range = np.array(list(np.arange(int(np.round(elat-2.5)),int(np.round(elat+2.5)))))
                        #lon_range = np.array(list(np.arange(int(np.round(elon-2.5)),int(np.round(elon+2.5)))))
                    elif 6 <= dura < 24:    # med events (6-24hrs): box = 10x10deg
                        lat_range = np.array(list(np.arange(int(np.round(elat-5)),int(np.round(elat+5)))))
                        lon_range = np.array(list(np.arange(int(np.round(elon-5)),int(np.round(elon+5)))))
                    elif dura > 24:         # big events (>24hrs): box = 15x15deg
                        lat_range = np.array(list(np.arange(int(np.round(elat-7.5)),int(np.round(elat+7.5)))))
                        lon_range = np.array(list(np.arange(int(np.round(elon-7.5)),int(np.round(elon+7.5)))))
                    '''
                    if dura < 6:            # small events (<6hrs): box = 5x5deg
                        off = 2.5
                    elif 6 <= dura < 24:    # med events (6-24hrs): box = 10x10deg
                        off = 5
                    elif dura > 24:         # big events (>24hrs): box = 15x15deg
                        off = 7.5
                    lon0 = elon-off         # make a box around the event lat/lon
                    lon1 = elon+off
                    lat0 = elat-off
                    lat1 = elat+off
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
                    
                    # if it's the first day of the event, add the start time info
                    if dy == dpre.date():
                        ft(qv2m, q250, q500, q850, u500, v500, tqv, t2m, lat_range, lon_range, index, prets, 0,
                            nc_var0, nc_var1, nc_var2, nc_var3, nc_var4, nc_var5, nc_var6, nc_var7)
                        if dy == ds.date():
                            ft(qv2m, q250, q500, q850, u500, v500,  tqv, t2m, lat_range, lon_range, index, tstart, 1,
                               nc_var0, nc_var1, nc_var2, nc_var3, nc_var4, nc_var5, nc_var6, nc_var7)
                            if dy == dp.date():
                                ft(qv2m, q250, q500, q850, u500, v500,  tqv, t2m, lat_range, lon_range, index, tmax_hrs, 2,
                                   nc_var0, nc_var1, nc_var2, nc_var3, nc_var4, nc_var5, nc_var6, nc_var7)
                                if dy == de.date():
                                    ft(qv2m, q250, q500, q850, u500, v500,  tqv, t2m, lat_range, lon_range, index, tend, 3,
                                       nc_var0, nc_var1, nc_var2, nc_var3, nc_var4, nc_var5, nc_var6, nc_var7)
                    elif dy == ds.date():
                        ft(qv2m, q250, q500, q850, u500, v500,  tqv, t2m, lat_range, lon_range, index, tstart, 1,
                            nc_var0, nc_var1, nc_var2, nc_var3, nc_var4, nc_var5, nc_var6, nc_var7)
                        if dy == dp.date():
                            ft(qv2m, q250, q500, q850, u500, v500,  tqv, t2m, lat_range, lon_range, index, tmax_hrs, 2,
                               nc_var0, nc_var1, nc_var2, nc_var3, nc_var4, nc_var5, nc_var6, nc_var7)
                            if dy == de.date():
                                ft(qv2m, q250, q500, q850, u500, v500,  tqv, t2m, lat_range, lon_range, index, tend, 3,
                                   nc_var0, nc_var1, nc_var2, nc_var3, nc_var4, nc_var5, nc_var6, nc_var7)
                                if dy == dpost.date():
                                    ft(qv2m, q250, q500, q850, u500, v500,  tqv, t2m, lat_range, lon_range, index, poste, 4,
                                       nc_var0, nc_var1, nc_var2, nc_var3, nc_var4, nc_var5, nc_var6, nc_var7)
                    # if it's the peak day of the event, add the peak time info
                    elif dy == dp.date():
                        ft(qv2m, q250, q500, q850, u500, v500,  tqv, t2m, lat_range, lon_range, index, tmax_hrs, 2,
                           nc_var0, nc_var1, nc_var2, nc_var3, nc_var4, nc_var5, nc_var6, nc_var7)
                        if dy == de.date():
                            ft(qv2m, q250, q500, q850, u500, v500,  tqv, t2m, lat_range, lon_range, index, tend, 3,
                               nc_var0, nc_var1, nc_var2, nc_var3, nc_var4, nc_var5, nc_var6, nc_var7)
                            if dy == dpost.date():
                                ft(qv2m, q250, q500, q850, u500, v500,  tqv, t2m, lat_range, lon_range, index, poste, 4,
                                   nc_var0, nc_var1, nc_var2, nc_var3, nc_var4, nc_var5, nc_var6, nc_var7)
                    # if it's the last day of the event, add the end time info
                    elif dy == de.date():
                        ft(qv2m, q250, q500, q850, u500, v500,  tqv, t2m, lat_range, lon_range, index, tend, 3,
                           nc_var0, nc_var1, nc_var2, nc_var3, nc_var4, nc_var5, nc_var6, nc_var7)
                        if dy == dpost.date():
                            ft(qv2m, q250, q500, q850, u500, v500,  tqv, t2m, lat_range, lon_range, index, poste, 4,
                               nc_var0, nc_var1, nc_var2, nc_var3, nc_var4, nc_var5, nc_var6, nc_var7)
                    elif dy == dpost.date():
                        ft(qv2m, q250, q500, q850, u500, v500,  tqv, t2m, lat_range, lon_range, index, poste, 4,
                            nc_var0, nc_var1, nc_var2, nc_var3, nc_var4, nc_var5, nc_var6, nc_var7)
                    # if it isn't the start, peak, or end day of the event, keep moving
                    else:
                        continue
                # 3D files:
                elif '400.inst3_3d_asm_Np' in os.fsdecode(dat):
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
                    #time = fh.variables['time'][:] # in minutes, 0-21hrs (every 3hrs: index 0-7)
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
                    '''
                    if dura < 6:            # small events (<6hrs): box = 5x5deg
                        lat_range = np.array(list(np.arange(int(np.round(elat-2.5)),int(np.round(elat+2.5)))))
                        lon_range = np.array(list(np.arange(int(np.round(elon-2.5)),int(np.round(elon+2.5)))))
                    elif 6 <= dura < 24:    # med events (6-24hrs): box = 10x10deg
                        lat_range = np.array(list(np.arange(int(np.round(elat-5)),int(np.round(elat+5)))))
                        lon_range = np.array(list(np.arange(int(np.round(elon-5)),int(np.round(elon+5)))))
                    elif dura > 24:         # big events (>24hrs): box = 15x15deg
                        lat_range = np.array(list(np.arange(int(np.round(elat-7.5)),int(np.round(elat+7.5)))))
                        lon_range = np.array(list(np.arange(int(np.round(elon-7.5)),int(np.round(elon+7.5)))))
                    '''
                    if dura < 6:            # small events (<6hrs): box = 5x5deg
                        off = 2.5
                    elif 6 <= dura < 24:    # med events (6-24hrs): box = 10x10deg
                        off = 5
                    elif dura > 24:         # big events (>24hrs): box = 15x15deg
                        off = 7.5
                    lon0 = elon-off         # make a box around the event lat/lon
                    lon1 = elon+off
                    lat0 = elat-off
                    lat1 = elat+off
                    lo0 = (8/5)*(lon0 + 180) + 1    # solve for the index where the lat/lon of the box corners are
                    la0 = 2*(lat0 + 90) + 1
                    lo1 = (8/5)*(lon1 + 180) + 1
                    la1 = 2*(lat1 + 90) + 1
                    lat_range = np.arange(int(la0), int(la1)) # now, the lat/lon ranges are indexing the lat/lon of the variables correctly
                    lon_range = np.arange(int(lo0), int(lo1))
                    
                    # if it's the first day of the event, add the start time info
                    if dy == dpre.date():
                        ft3d(vv, rh, qv, u, v, temp, lat_range, lon_range, index, prets, 0,
                             nc_var8, nc_var9, nc_var10, nc_var11, nc_var12, nc_var13)
                        if dy == ds.date():
                            ft3d(vv, rh, qv, u, v, temp, lat_range, lon_range, index, tstart, 1, 
                                 nc_var8, nc_var9, nc_var10, nc_var11, nc_var12, nc_var13)
                            if dy == dp.date():
                                ft3d(vv, rh, qv, u, v, temp, lat_range, lon_range, index, tmax_hrs, 2,
                                     nc_var8, nc_var9, nc_var10, nc_var11, nc_var12, nc_var13)
                                if dy == de.date():
                                    ft3d(vv, rh, qv, u, v, temp, lat_range, lon_range, index, tend, 3,
                                         nc_var8, nc_var9, nc_var10, nc_var11, nc_var12, nc_var13)
                    elif dy == ds.date():
                        ft3d(vv, rh, qv, u, v, temp, lat_range, lon_range, index, tstart, 1, 
                             nc_var8, nc_var9, nc_var10, nc_var11, nc_var12, nc_var13)
                        if dy == dp.date():
                            ft3d(vv, rh, qv, u, v, temp, lat_range, lon_range, index, tmax_hrs, 2, 
                                 nc_var8, nc_var9, nc_var10, nc_var11, nc_var12, nc_var13)
                            if dy == de.date():
                                ft3d(vv, rh, qv, u, v, temp, lat_range, lon_range, index, tend, 3, 
                                     nc_var8, nc_var9, nc_var10, nc_var11, nc_var12, nc_var13)
                                if dy == dpost.date():
                                    ft3d(vv, rh, qv, u, v, temp, lat_range, lon_range, index, poste, 4, 
                                         nc_var8, nc_var9, nc_var10, nc_var11, nc_var12, nc_var13)
                    # if it's the peak day of the event, add the peak time info
                    elif dy == dp.date():
                        ft3d(vv, rh, qv, u, v, temp, lat_range, lon_range, index, tmax_hrs, 2,
                             nc_var8, nc_var9, nc_var10, nc_var11, nc_var12, nc_var13)
                        if dy == de.date():
                            ft3d(vv, rh, qv, u, v, temp, lat_range, lon_range, index, tend, 3,
                                 nc_var8, nc_var9, nc_var10, nc_var11, nc_var12, nc_var13)
                            if dy == dpost.date():
                                ft3d(vv, rh, qv, u, v, temp, lat_range, lon_range, index, poste, 4,
                                     nc_var8, nc_var9, nc_var10, nc_var11, nc_var12, nc_var13)
                    # if it's the last day of the event, add the end time info
                    elif dy == de.date():
                        ft3d(vv, rh, qv, u, v, temp, lat_range, lon_range, index, tend, 3,
                             nc_var8, nc_var9, nc_var10, nc_var11, nc_var12, nc_var13)
                        if dy == dpost.date():
                            ft3d(vv, rh, qv, u, v, temp, lat_range, lon_range, index, poste, 4,
                                 nc_var8, nc_var9, nc_var10, nc_var11, nc_var12, nc_var13)
                    elif dy == dpost.date():
                        ft3d(vv, rh, qv, u, v, temp, lat_range, lon_range, index, poste, 4,
                             nc_var8, nc_var9, nc_var10, nc_var11, nc_var12, nc_var13)
                    # if it isn't the start, peak, or end day of the event, keep moving
                    else:
                        continue
                    
                fh.close()
            
    nc_out.close()  # close the new file
