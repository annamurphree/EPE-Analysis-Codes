"""
Created on Thu Jun 25 16:26:36 2020

@author: Anna Murphree (annammurphree@gmail.com)
For Summer 2020 Internship with Dr. Yaping Zhou

This code locates EPEs from a statistics file within the MERRA-2 dataset, collects 
variables in grids around the events' centers, and saves those grids into a new .nc4 file.
MERRA-2 files used: .tavg1_2d_slv_Nx (2D single-level, 1-hourly)
Areas collected: 5x5deg for 0-6hr events, 10x10deg for 6-24hr events, 15x15deg for >24hr events
Times collected: start, peak, and end times of the events. 
"""

# this fills in the meteorological data for the start time of the event:
def fts(qv2m, q250, q500, q850, u850, v850, 
        lat_range, lon_range, lat, lon, index, tstart, 
        nc_var0, nc_var1, nc_var2, nc_var3, nc_var4, nc_var5):
    import numpy as np
    print('start day')
    # make temporary holders
    v0s = []
    v1s = []
    v2s = []
    v3s = []
    v4s = []
    v5s = []
    for la in lat_range:
        for lo in lon_range:
            v0i = qv2m[tstart, la, lo]
            v0s.append(v0i)
            v1i = q250[tstart, la, lo]
            v1s.append(v1i)
            v2i = q500[tstart, la, lo]
            v2s.append(v2i)
            v3i = q850[tstart, la, lo]
            v3s.append(v3i)
            v4i = u850[tstart, la, lo]
            v4s.append(v4i)
            v5i = v850[tstart, la, lo]
            v5s.append(v5i)
    
    # add the averages to the holders in the .nc file:
    nc_var0[index,tstart,lat,lon] = np.nanmean(v0s)
    nc_var1[index,tstart,lat,lon] = np.nanmean(v1s)
    nc_var2[index,tstart,lat,lon] = np.nanmean(v2s)
    nc_var3[index,tstart,lat,lon] = np.nanmean(v3s)
    nc_var4[index,tstart,lat,lon] = np.nanmean(v4s)
    nc_var5[index,tstart,lat,lon] = np.nanmean(v5s)
    
    return v0s, v1s, v2s, v3s, v4s, v5s

# this fills in the meteorological data for the peak time of the event:
def ftp(qv2m, q250, q500, q850, u850, v850, 
        lat_range, lon_range, lat, lon, index, tmax_hrs, 
        nc_var0, nc_var1, nc_var2, nc_var3, nc_var4, nc_var5):
    import numpy as np
    print('peak day')
    # make temporary holders
    v0p = []
    v1p = []
    v2p = []
    v3p = []
    v4p = []
    v5p = []
    for la in lat_range:
        for lo in lon_range:
            v0i = qv2m[tmax_hrs, la, lo]
            v0p.append(v0i)
            v1i = q250[tmax_hrs, la, lo]
            v1p.append(v1i)
            v2i = q500[tmax_hrs, la, lo]
            v2p.append(v2i)
            v3i = q850[tmax_hrs, la, lo]
            v3p.append(v3i)
            v4i = u850[tmax_hrs, la, lo]
            v4p.append(v4i)
            v5i = v850[tmax_hrs, la, lo]
            v5p.append(v5i)
    
    # add the averages to the holders in the .nc file:
    nc_var0[index,tmax_hrs,lat,lon] = np.nanmean(v0p)
    nc_var1[index,tmax_hrs,lat,lon] = np.nanmean(v1p)
    nc_var2[index,tmax_hrs,lat,lon] = np.nanmean(v2p)
    nc_var3[index,tmax_hrs,lat,lon] = np.nanmean(v3p)
    nc_var4[index,tmax_hrs,lat,lon] = np.nanmean(v4p)
    nc_var5[index,tmax_hrs,lat,lon] = np.nanmean(v5p)
    
    return v0p, v1p, v2p, v3p, v4p, v5p

# this fills in the meteorological data for the end time of the event:
def fte(qv2m, q250, q500, q850, u850, v850, 
        lat_range, lon_range, lat, lon, index, tend, 
        nc_var0, nc_var1, nc_var2, nc_var3, nc_var4, nc_var5):            
    import numpy as np
    print('end day')
    # make temporary holders
    v0e = []
    v1e = []
    v2e = []
    v3e = []
    v4e = []
    v5e = []
    for la in lat_range:
        for lo in lon_range:
            v0i = qv2m[tend, la, lo]
            v0e.append(v0i)
            v1i = q250[tend, la, lo]
            v1e.append(v1i)
            v2i = q500[tend, la, lo]
            v2e.append(v2i)
            v3i = q850[tend, la, lo]
            v3e.append(v3i)
            v4i = u850[tend, la, lo]
            v4e.append(v4i)
            v5i = v850[tend, la, lo]
            v5e.append(v5i)
    
    # add the averages to the holders in the .nc file:
    nc_var0[index,tend,lat,lon] = np.nanmean(v0e)
    nc_var1[index,tend,lat,lon] = np.nanmean(v1e)
    nc_var2[index,tend,lat,lon] = np.nanmean(v2e)
    nc_var3[index,tend,lat,lon] = np.nanmean(v3e)
    nc_var4[index,tend,lat,lon] = np.nanmean(v4e)
    nc_var5[index,tend,lat,lon] = np.nanmean(v5e)
    
    return v0e, v1e, v2e, v3e, v4e, v5e

# this makes a list of all the MERRA-2 data files for a season:
def files_list(core, yr, m0, m1, m2):
    # open the MERRA-2 data files in their directories:
    import os
    # MERRA 2 cores: 100 = 1980-1990, 200 = 1991-1999, 300 = 2000-2009, 400 = 2010-2019
    # the first month of the season:
    directory0 = os.fsencode(f'/discover/nobackup/projects/gmao/merra2/data/products/MERRA2_{core}/Y{yr}/M{m0}')
    # the second month of the season:
    directory1 = os.fsencode(f'/discover/nobackup/projects/gmao/merra2/data/products/MERRA2_{core}/Y{yr}/M{m1}')
    # the third month of the season:
    directory2 = os.fsencode(f'/discover/nobackup/projects/gmao/merra2/data/products/MERRA2_{core}/Y{yr}/M{m2}')
    # the 2D averaged files:
    files = [os.path.join(os.fsdecode(directory0), os.fsdecode(fi)) 
             for fi in os.listdir(directory0) if f'{core}.tavg1_2d_slv_Nx' in os.fsdecode(fi)]
    files1 = [os.path.join(os.fsdecode(directory1), os.fsdecode(fi)) 
             for fi in os.listdir(directory1) if f'{core}.tavg1_2d_slv_Nx' in os.fsdecode(fi)]
    files2 = [os.path.join(os.fsdecode(directory2), os.fsdecode(fi)) 
             for fi in os.listdir(directory2) if f'{core}.tavg1_2d_slv_Nx' in os.fsdecode(fi)]
    # add these all into 1 list:
    files.extend(files1)
    files.extend(files2)
    
    return files

def find_IMERG_in_MERRA2(imerg, core, yr, m0, m1, m2, output):
    # imerg = the stats found from events in IMERG data
    # direct = directory to MERRA-2 data files
    
    from netCDF4 import Dataset
    import numpy as np
    import datetime
    
    # open the event stats from IMERG:
    stats = Dataset(imerg, mode='r')
    #variables = fh.variables.keys()
    #print(variables)
    #print(fh.dimensions)
    params = stats.variables['event_parameters'][:]
    #print(mask) # this is a masked array without a mask?
    #print(params.data)
    
    #open the MERRA-2 data files in their directories:
    files = files_list(core, yr, m0, m1, m2)
    
    # open a new .nc4 file to write the data to:
    nc_out = Dataset('%s'%(output), 'w', format='NETCDF4')
    
    # using our previous dimension information, we can create the new dimensions:
    nc_in = Dataset(files[0], mode='r')
    nc_dims = [dim for dim in nc_in.dimensions]
    data = {}
    # this sets the same dimensions as the MERRA-2 data files:
    for dim in nc_dims:
        nc_out.createDimension(dim, nc_in.variables[dim].size)
        data[dim] = nc_out.createVariable(dim, nc_in.variables[dim].dtype,\
            (dim,))
        #print(dim)
        # you can do this step yourself but someone else did the work for us:
        for ncattr in nc_in.variables[dim].ncattrs():
            data[dim].setncattr(ncattr, nc_in.variables[dim].getncattr(ncattr))
            #print(ncattr)
            
    nc_out.createDimension('index', None) # create # of indices
    ind = nc_out.createVariable('index', 'i4', ('index',))
    
    # create departure variables:
    # if you wanted to find other variables, change these to what you'd like:
    nc_var0 = nc_out.createVariable('QV2M', 'f8', ('index', 'time', 'lat', 'lon',))
    nc_var0.setncatts({'long_name': u"2mb specific humidity",\
                      'units': u"kg kg^-1", 'level_desc': u'2mb',\
                      'var_desc': u"2mb specific humidity",\
                      'statistic': u'instantaneous'})
    nc_var1 = nc_out.createVariable('Q250', 'f8', ('index', 'time', 'lat', 'lon',))
    nc_var1.setncatts({'long_name': u"250mb specific humidity at start",\
                      'units': u"kg kg^-1", 'level_desc': u'250mb',\
                      'var_desc': u"250mb specific humidity",\
                      'statistic': u'instantaneous'})
    nc_var2 = nc_out.createVariable('Q500', 'f8', ('index', 'time', 'lat', 'lon',))
    nc_var2.setncatts({'long_name': u"500mb specific humidity",\
                      'units': u"kg kg^-1", 'level_desc': u'500mb',\
                      'var_desc': u"500mb specific humidity",\
                      'statistic': u'instantaneous'})
    nc_var3 = nc_out.createVariable('Q850', 'f8', ('index', 'time', 'lat', 'lon',))
    nc_var3.setncatts({'long_name': u"850mb specific humidity",\
                      'units': u"kg kg^-1", 'level_desc': u'850mb',\
                      'var_desc': u"850mb specific humidity",\
                      'statistic': u'instantaneous'})
    nc_var4 = nc_out.createVariable('U500', 'f8', ('index', 'time', 'lat', 'lon',))
    nc_var4.setncatts({'long_name': u"500mb eastward wind",\
                      'units': u"m/s", 'level_desc': u'500mb',\
                      'var_desc': u"500mb eastward wind",\
                      'statistic': u'instantaneous'})
    nc_var5 = nc_out.createVariable('V500', 'f8', ('index', 'time', 'lat', 'lon',))
    nc_var5.setncatts({'long_name': u"500mb northward wind",\
                      'units': u"m/s", 'level_desc': u'500mb',\
                      'var_desc': u"500mb northward wind",\
                      'statistic': u'instantaneous'})
    ''
    # cycle through each event from IMERG:
    for event in params:
        dstart = event[0]   # start time (days) of the year
        dend = event[1]     # end time (days) of the year
        year = event[22]    # year of event, which I add 1 to because the stats files are a year off
        year = int(year + 1)
        #print('Dates: ',dstart,'-', dend, ',',year)
        long = event[2]     # longitude (deg East)
        #print('IMERG long: ',long)
        lati = event[3]     # latitude (deg North)
        dura = event[4]     # event duration (hrs)
        #print('IMERG lat: ',lati)
        tmax = event[34]    # time into event of max rain volume (hours) 
        index = int(event[28])   # event index #
        #print('Index: ',index)
        ind[index] = index
        # convert start/end # days/year to datetime format
        # datetime.datetime: (year, month, day, hour, minute, second, timezone info)
        # this format gives the full date and time:
        ds = datetime.datetime(year, 1, 1, 0, 0, 0) + datetime.timedelta(dstart - 1)
        #print('Start time: ',ds)
        dp = datetime.datetime(year, 1, 1, 0, 0, 0) + datetime.timedelta(dstart - 1, tmax*3600)
        #print('Peak time: ',dp)
        de = datetime.datetime(year, 1, 1, 0, 0, 0) + datetime.timedelta(dend - 1)
        #print('End time: ' ,de)
        # this just gives the hours as integers, for indexing the MERRA-2 data:
        tstart = int((dstart - int(dstart))*24)
        #print('Start hour:',tstart)
        tmax_days = dstart + (tmax/24)
        tmax_hrs = int((tmax_days - int(tmax_days))*24)
        #print('Peak hour: ',tmax_hrs)
        tend = int((dend - int(dend))*24)
        #print('End hour: ',tend)
        #tim_range = np.array(list([tstart, tmax_hrs, tend]))
        #print('Time range: ',tim_range)
        #day_count = (de.date() - ds.date()).days + 1
        #print('Day range: ',day_count)
        
        # find these events in the MERRA-2 data:
        # put each event day's file into a list:
        event_days = []
        for fi in files:
            fh = Dataset(fi, mode='r')
            #print(fh)
            #variables = fh.variables.keys()
            #print(variables)
            date = fh.RangeBeginningDate
            #print('MERRA-2 day: ',date)
            dy = int(date[0:4])  # year
            dm = int(date[5:7])  # month
            dd = int(date[8:10]) # day
            day = datetime.datetime(dy, dm, dd).date()
            
            #for day in (ds.date() + datetime.timedelta(n) for n in range(day_count)):
            if (ds.date() <= day <= de.date()):
                event_days.append(fi)
            fh.close()
        
        #print(event_days)
        # can search for events here, by adding conditions:
        #if len(event_days) != 0 and (-96 <= long <= -90) and (36 <= lati <= 46):
        if len(event_days) != 0:
            
            # if the MERRA-2 date is within the date range of the event:
            for day in event_days:
                # print event info here, so it's only for the events matched:
                #print('IMERG long: ',long)
                #print('IMERG lat: ',lati)
                print('Index: ',index)
                '''
                print('Start time: ',ds)
                print('Peak time: ',dp)
                print('End time: ' ,de)
                print('Day range: ',day_count)
                '''
                
                fh = Dataset(day, mode='r')
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
                #u850 = fh.variables['U850'][:]  # eastward wind 850m
                #v850 = fh.variables['V850'][:]  # northward wind 850m
                u500 = fh.variables['U500'][:]  # eastward wind 500m
                v500 = fh.variables['V500'][:]  # northward wind 500m
                
                # match up IMERG coordinates w/ MERRA-2 coordinates:
                lon = np.asarray(lon)
                close_lon = (np.abs(lon - long)).argmin()   # event lon's index in MERRA-2
                elon = lon[close_lon]         # event longitude in MERRA-2
                lat = np.asarray(lat)
                close_lat = (np.abs(lat - lati)).argmin()   # event lat's index in MERRA-2
                elat = lat[close_lat]         # event latitude in MERRA-2
                
                # default to these grid sizes, scaled by the event's duration:
                if dura < 6:            # small events (<6hrs): box = 5x5deg
                    off = 2.5
                elif 6 <= dura < 24:    # med events (6-24hrs): box = 10x10deg
                    off = 5
                elif dura > 24:         # big events (>24hrs) : box = 15x15deg
                    off = 7.5
                # (off is the offset around the center point)
                lon0 = elon-off         # make a box around the event lat/lon
                lon1 = elon+off
                lat0 = elat-off
                lat1 = elat+off
                #la_range = np.arange(int(lat0), int(lat1)) 
                #lo_range = np.arange(int(lon0), int(lon1))
                #print('lons: ', lon0, elon, lon1)
                #print('lats: ', lat0, elat, lat1)
                lo0 = (8/5)*(lon0 + 180) + 1    # solve for the index where the lon of the box corners are
                la0 = 2*(lat0 + 90) + 1
                lo1 = (8/5)*(lon1 + 180) + 1    # solve for the index where the lat of the box corners are
                la1 = 2*(lat1 + 90) + 1
                #print('lon range: ', lo0, lo1)
                #print('lat range: ', la0, la1)
                lat_range = np.arange(int(la0), int(la1)) # now, the lat/lon ranges are indexing the lat/lon of the variables correctly
                lon_range = np.arange(int(lo0), int(lo1))
                
                # if it's the first day of the event, add the start time info
                if dy == ds.date():
                    (v0s, v1s, v2s, v3s, v4s, v5s) = fts(qv2m, q250, q500, q850, u500, v500, lat_range, lon_range, elat, elon, index, tstart, 
                                                         nc_var0, nc_var1, nc_var2, nc_var3, nc_var4, nc_var5)
                    if dy == dp.date():
                        (v0p, v1p, v2p, v3p, v4p, v5p) = ftp(qv2m, q250, q500, q850, u500, v500, lat_range, lon_range,  elat, elon, index, tmax_hrs, 
                                                             nc_var0, nc_var1, nc_var2, nc_var3, nc_var4, nc_var5)
                        if dy == de.date():
                            (v0e, v1e, v2e, v3e, v4e, v5e) = fte(qv2m, q250, q500, q850, u500, v500, lat_range, lon_range,  elat, elon, index, tend, 
                                                                 nc_var0, nc_var1, nc_var2, nc_var3, nc_var4, nc_var5)
                # if it's the peak day of the event, add the peak time info
                elif dy == dp.date():
                    (v0p, v1p, v2p, v3p, v4p, v5p) = ftp(qv2m, q250, q500, q850, u500, v500, lat_range, lon_range,  elat, elon, index, tmax_hrs, 
                                                         nc_var0, nc_var1, nc_var2, nc_var3, nc_var4, nc_var5)
                    if dy == de.date():
                        (v0e, v1e, v2e, v3e, v4e, v5e) = fte(qv2m, q250, q500, q850, u500, v500, lat_range, lon_range,  elat, elon, index, tend, 
                                                             nc_var0, nc_var1, nc_var2, nc_var3, nc_var4, nc_var5)
                # if it's the last day of the event, add the end time info
                elif dy == de.date():
                    (v0e, v1e, v2e, v3e, v4e, v5e) = fte(qv2m, q250, q500, q850, u500, v500, lat_range, lon_range,  elat, elon, index, tend, 
                                                         nc_var0, nc_var1, nc_var2, nc_var3, nc_var4, nc_var5)
                # if it isn't the start, peak, or end day of the event, keep moving
                else:
                    continue
                
                fh.close()
            
    nc_out.close()  # close the new file

