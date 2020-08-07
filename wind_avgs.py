"""
Created on Thu Jul 16 10:16:01 2020

@author: Anna Murphree (annammurphree@gmail.com)
For Summer 2020 Internship with Dr. Yaping Zhou

This code computes the seasonal mean and trend of u and v from 1980 to 2019, or over a
chosen range of years. You can also specify a latitude/longitude range.
MERRA-2 files used:
    files_list: .tavg1_2d_slv_Nx (2D single-level, 1-hourly)
    files_list_monthly: .instM_3d_asm_Np (3D, monthly averages)
"""

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
    # write this data into a new text file:
    #with open(f'{ysn}_files.txt', 'w') as filename:
    #    filename.writelines('%s\n' % file for file in files)
    
    return files

# this makes a list of a season's MERRA-2 monthly average files:
def files_list_monthly(core, yr, m0, m1, m2):
    # open the MERRA-2 data files in their directories:
    import os
    # MERRA 2 cores: 100 = 1980-1990, 200 = 1991-1999, 300 = 2000-2009, 400 = 2010-2019
    # the first month of the season:
    directory0 = os.fsencode(f'/discover/nobackup/projects/gmao/merra2/data/products/MERRA2_{core}/Y{yr}/M{m0}')
    # the second month of the season:
    directory1 = os.fsencode(f'/discover/nobackup/projects/gmao/merra2/data/products/MERRA2_{core}/Y{yr}/M{m1}')
    # the third month of the season:
    directory2 = os.fsencode(f'/discover/nobackup/projects/gmao/merra2/data/products/MERRA2_{core}/Y{yr}/M{m2}')
    # the 3D monthly average files:
    files = [os.path.join(os.fsdecode(directory0), os.fsdecode(fi)) 
             for fi in os.listdir(directory0) if f'{core}.instM_3d_asm_Np' in os.fsdecode(fi)]
    files1 = [os.path.join(os.fsdecode(directory1), os.fsdecode(fi)) 
             for fi in os.listdir(directory1) if f'{core}.instM_3d_asm_Np' in os.fsdecode(fi)]
    files2 = [os.path.join(os.fsdecode(directory2), os.fsdecode(fi)) 
             for fi in os.listdir(directory2) if f'{core}.instM_3d_asm_Np' in os.fsdecode(fi)]
    # add these all into 1 list:
    files.extend(files1)
    files.extend(files2)
    # write this data into a new text file:
    #with open(f'{ysn}_files.txt', 'w') as filename:
    #    filename.writelines('%s\n' % file for file in files)
    
    return files

# This finds the right grids in MERRA-2 and averages them for a whole season or year
# yrsn is either 'yr' or 'sn', depending on what you want to average
# avgmap is either 'avg' or 'map', depending on if you want an overall average or a map of averages    
def get_avgs(core, yr, m0, m1, m2, nc_u, nc_v, yrsn, avgmap, lonmin=-65, lonmax=-125, latmin=25, latmax=50, sn=0, month='yes'):
    from netCDF4 import Dataset
    import numpy as np
    # if you don't want to use monthly averages (slower):
    if month == 'no':
        if yrsn == 'sn':
            # for just a season:
            files = files_list(core, yr, m0, m1, m2)
        elif yrsn == 'yr':
            # for a year:
            files = files_list(core, yr, f'0{m0}', f'0{m1}', f'0{m2}')
            files1 = files_list(core, yr, f'0{m0+3}', f'0{m1+3}', f'0{m2+3}')
            files2 = files_list(core, yr, f'0{m0+6}', f'0{m1+6}', f'0{m2+6}')
            files3 = files_list(core, yr, m0+9, m1+9, m2+9)
            files.extend(files1)
            files.extend(files2)
            files.extend(files3)
    # if you do want to use monthly averages (faster, default):
    elif month == 'yes':
        if yrsn == 'sn':
            # for just a season:
            files = files_list_monthly(core, yr, m0, m1, m2)
        elif yrsn == 'yr':
            # for a year:
            files = files_list_monthly(core, yr, f'0{m0}', f'0{m1}', f'0{m2}')
            files1 = files_list_monthly(core, yr, f'0{m0+3}', f'0{m1+3}', f'0{m2+3}')
            files2 = files_list_monthly(core, yr, f'0{m0+6}', f'0{m1+6}', f'0{m2+6}')
            files3 = files_list_monthly(core, yr, m0+9, m1+9, m2+9)
            files.extend(files1)
            files.extend(files2)
            files.extend(files3)
    
    #print(files)
    # create holders for the U and V values
    u = []
    v = []
    # loop through each data file in the season:
    for fi in files:
        fh = Dataset(fi, mode='r')
        lon = fh.variables['lon'][:]    # longitude, 576 dimension: 0.625 deg
        lat = fh.variables['lat'][:]    # latitude, 361 dimension: 0.5 deg
        # define the U and V variables in the file:
        # (if you wanted to run this on other variables, you'd just need to change what you're indexing here)
        if month == 'no':
            u500 = fh.variables['U500'][:]  # eastward wind 500m
            v500 = fh.variables['V500'][:]  # northward wind 500m
        elif month == 'yes':
            u500 = fh.variables['U'][0,16]  # eastward wind (0=time, 16=index for 500hPa)
            v500 = fh.variables['V'][0,16]  # northward wind (0=time, 16=index for 500hPa)
        
        # CONUS: 25->50 deg latitude, -65->-125 deg longitude
        # the variables are indexed as var(time->24,lat->361,lon->576)
        lon = np.asarray(lon)
        lat = np.asarray(lat)
        ind_lonmin = (np.abs(lon -lonmin)).argmin()   # index for minimum longitude (default is -65)
        ind_lonmax = (np.abs(lon -lonmax)).argmin()   # index for maximum longitude (default is -125)
        ind_latmin = (np.abs(lat - latmin)).argmin()  # index for minimum latitude (default is 25)
        ind_latmax = (np.abs(lat - latmax)).argmin()  # index for maximum latitude (default is 50)
        lat_range = np.arange(int(ind_latmin), int(ind_latmax))
        lon_range = np.arange(int(ind_lonmax), int(ind_lonmin))  # list them backwards (index for -125 is bigger than index for -65)
        #print(lon_range)
        '''
        # finds the MERRA-2 lat/lon values at those indices:
        lonmin = lon[ind_lonmin]
        lonmax = lon[ind_lonmax]
        latmin = lat[ind_latmin]
        latmax = lat[ind_latmax]
        
        # to find lat/lon of a point by its indices i & j:
        long = -180 + (5/8)*(i - 1)
        lati = -90 + (1/2)*(j - 1)
        '''
        
        if month == 'no':
            time_range = np.arange(0,24)
            # for an average over the whole time/lat/lon domain:
            if avgmap == 'avg':
                for t in time_range: 
                    for la in lat_range:
                        for lo in lon_range:
                            # add the U and V values to their holders:
                            u.append(u500[t,la,lo])
                            v.append(v500[t,la,lo])
            # for a map of each coordinate's average over the time domain:
            elif avgmap == 'map':
                for t in time_range: 
                    u5 = [u500[t,la,lo] for la in lat_range for lo in lon_range]
                    v5 = [v500[t,la,lo] for la in lat_range for lo in lon_range]
                    u5 = np.array(u5)
                    ui = np.reshape(u5, (len(lat_range), len(lon_range)))
                    #ui = np.abs(ui)   # try taking absval before averaging arrays
                    u.append(ui)
                    v5 = np.array(v5)
                    vi = np.reshape(v5, (len(lat_range), len(lon_range)))
                    #vi = np.abs(vi)   # try taking absval before averaging arrays
                    v.append(vi)
        elif month == 'yes':
            if avgmap == 'avg':
                for la in lat_range:
                    for lo in lon_range:
                        u.append(u500[la,lo])
                        v.append(v500[la,lo])
            elif avgmap == 'map':
                u5 = [u500[la,lo] for la in lat_range for lo in lon_range]
                v5 = [v500[la,lo] for la in lat_range for lo in lon_range]
                u5 = np.array(u5)
                ui = np.reshape(u5, (len(lat_range), len(lon_range)))
                #ui = np.abs(ui)   # try taking absval before averaging arrays
                u.append(ui)
                v5 = np.array(v5)
                vi = np.reshape(v5, (len(lat_range), len(lon_range)))
                #vi = np.abs(vi)   # try taking absval before averaging arrays
                v.append(vi)
    
    if avgmap == 'avg':
        print('U size: ', len(u))
        print('V size: ', len(v))
        u_avg = np.nanmean(u)
        print('U Avg: ', u_avg)
        v_avg = np.nanmean(v)
        print('V Avg: ', v_avg)
        yr_index = int(yr-1980)  # set index for outfile (0-39)
        print('Year: ', yr, yr_index)
        # add the averages to the new .nc4 file variables:
        if yrsn == 'sn':
            nc_u[yr_index,sn] = round(u_avg,8)
            nc_v[yr_index,sn] = round(v_avg,8)
        elif yrsn == 'yr':
            nc_u[yr_index] = round(u_avg,8)
            nc_v[yr_index] = round(v_avg,8)
            
        return round(u_avg,8), round(v_avg,8)
    
    elif avgmap == 'map':
        ugrids = np.array(u)
        #print(ugrids)
        vgrids = np.array(v)
        uavgd = np.mean(ugrids,axis=0)
        #print('U Avgd: ', uavgd)
        vavgd = np.mean(vgrids,axis=0)
        # to compute change map: 
        times = np.arange(ugrids.shape[0])          # of days x 24hrs
        uslope = []
        vslope = []
        for la in range(ugrids.shape[1]):           # 2D dimension (latitude)
            for lo in range(ugrids.shape[2]):       # 2D dimension (longitude)
                ut = []                             # holder for each pixel's time series
                for t in range(ugrids.shape[0]):    # 3rd dimension (# of 2D arrays)
                    ut.append(ugrids[t,la,lo])        
                #print(ut)
                mu, bu = np.polyfit(times, ut, 1)   # fits a line to the pixel's time series
                #print(mu)
                uslope.append(mu)
        for la in range(vgrids.shape[1]):
            for lo in range(vgrids.shape[2]):
                vt = []
                for t in range(vgrids.shape[0]):
                    vt.append(vgrids[t,la,lo])        
                #print(vt)
                mv, bv = np.polyfit(times, vt, 1)
                #print(mv)
                vslope.append(mv)
        # make sure the slope arrays are in the right shape:
        uslope = np.array(uslope)
        uslope = np.reshape(uslope, (len(lat_range), len(lon_range)))
        vslope = np.array(vslope)
        vslope = np.reshape(vslope, (len(lat_range), len(lon_range)))
        if yrsn == 'sn':
            print('U size: ',len(u), ugrids.shape, uavgd.shape)
            print('V size: ',len(v), vgrids.shape, vavgd.shape)
            # plot the 
            map_plot1(uavgd, 'U500', yr, sn, lonmin, lonmax, latmin, latmax, grid_change=uslope)
            map_plot1(vavgd, 'V500', yr, sn, lonmin, lonmax, latmin, latmax, grid_change=vslope)
            map_plot1(uslope, 'U500_Change', yr, sn, lonmin, lonmax, latmin, latmax, grid_change=uavgd)
            map_plot1(vslope, 'V500_Change', yr, sn, lonmin, lonmax, latmin, latmax, grid_change=vavgd)
            yr_index = int(yr-1980)  # set index for outfile (0-39)
            print('Year: ', yr, yr_index)
            ''
            for la in range(0,50):
                for lo in range(0,96):
                    nc_u[yr_index,sn,la,lo] = uavgd[la,lo]
                    nc_v[yr_index,sn,la,lo] = vavgd[la,lo]
                    ''
        elif yrsn == 'yr':
            print('U size: ',len(u), ugrids.shape, uavgd.shape)
            print('V size: ',len(v), vgrids.shape, vavgd.shape)
            map_plot1(uavgd, 'U500', yr, lonmin, lonmax, latmin, latmax, grid_change=uslope)
            map_plot1(vavgd, 'V500', yr, lonmin, lonmax, latmin, latmax, grid_change=vslope)
            map_plot1(uslope, 'U500_Change', yr, lonmin, lonmax, latmin, latmax, grid_change=uavgd)
            map_plot1(vslope, 'V500_Change', yr, lonmin, lonmax, latmin, latmax, grid_change=vavgd)
            yr_index = int(yr-1980)  # set index for outfile (0-39)
            print('Year: ', yr, yr_index)
            for la in range(0,50):
                for lo in range(0,96):
                    nc_u[yr_index,la,lo] = uavgd[la,lo]
                    nc_v[yr_index,la,lo] = vavgd[la,lo]
    

# yrsn is either 'yr' or 'sn', depending on what you want to average
# avgmap is either 'avg' or 'map', depending on if you want an overall average or a map of averages
def main(outfile, yrsn, avgmap, lonmin=-65, lonmax=-125, latmin=25, latmax=50, pick_c='no', pick_yr='no', pick_sn='no'):
    from netCDF4 import Dataset
    import numpy as np
    # open a new NetCDF4 file to write the data to:
    nc_out = Dataset(f'{outfile}', 'w', format='NETCDF4')    
    # MERRA-2 cores: 100 = 1980-1990, 200 = 1991-1999, 300 = 2000-2009, 400 = 2010-2019
    # years = np.arange(1980, 2019)
    cyrs = {100:np.arange(1980,1992), 200:np.arange(1992,2001), 300:np.arange(2001,2011), 400:np.arange(2011,2020)}
    # sns : 0 (DJF), 1 (MAM), 2 (JJA), 3 (SON)
    # months = ['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12']
    sms = {0: ['12', '01', '02'], 1: ['03', '04', '05'], 2: ['06', '07', '08'], 3: ['09', '10', '11']}
    # if you want to pick a specific season:
    if pick_yr != 'no':
        nc_out.createDimension('season', 4)     # 0=DJF, 1=MAM, 2=JJA, 3=SON
        nc_out.createDimension('year', 40)      # 0=1980, ..., 38=2019
        nc_out.createDimension('lat', None)  
        nc_out.createDimension('lon', None)  
        nc_u = nc_out.createVariable('U500', 'f8', ('year', 'season', 'lat', 'lon',))
        nc_v = nc_out.createVariable('V500', 'f8', ('year', 'season', 'lat', 'lon',))
        
        get_avgs(pick_c, pick_yr, sms[pick_sn][0], sms[pick_sn][1], sms[pick_sn][2], nc_u, nc_v, yrsn, avgmap, lonmin, lonmax, latmin, latmax, sn=pick_sn, month='yes')
    # if you want to run it for every yr/sn since 1980:
    elif pick_yr == 'no':
        # if you want seasonal averages:
        if yrsn == 'sn':
            nc_out.createDimension('season', 4)     # 0=DJF, 1=MAM, 2=JJA, 3=SON
            nc_out.createDimension('year', 40)      # 0=1980, ..., 38=2019
            # if you want the average over the whole CONUS:
            if avgmap == 'avg':
                nc_u = nc_out.createVariable('U500', 'f8', ('year', 'season',))
                nc_v = nc_out.createVariable('V500', 'f8', ('year', 'season',))
                for c in cyrs:              # this should loop through every core
                    for yr in cyrs[c]:      # and then through each year
                        for sn in sms:      # month 1    month 2     month 3
                            nc_u[int(yr-1980),sn], nc_v[int(yr-1980),sn] = get_avgs(c, yr, sms[sn][0], sms[sn][1], sms[sn][2], nc_u, nc_v, 'sn', 'avg', lonmin, lonmax, latmin, latmax, sn)
                            print('Data written: ',nc_u[int(yr-1980),sn],nc_v[int(yr-1980),sn])
            # if you want a map of the averages over CONUS:
            elif avgmap == 'map':
                nc_out.createDimension('lat', None)  
                nc_out.createDimension('lon', None)  
                nc_u = nc_out.createVariable('U500', 'f8', ('year', 'season', 'lat', 'lon',))
                nc_v = nc_out.createVariable('V500', 'f8', ('year', 'season', 'lat', 'lon',))
                for c in cyrs:              # this should loop through every core
                    for yr in cyrs[c]:      # and then through each year
                        for sn in sms:      # month 1    month 2     month 3
                            get_avgs(c, yr, sms[sn][0], sms[sn][1], sms[sn][2], nc_u, nc_v, 'sn', 'map', lonmin, lonmax, latmin, latmax, sn)
                                         
        # if you want yearly averages:
        elif yrsn == 'yr':
            nc_out.createDimension('year', 40)      # 0=1980, ..., 38=2019
            if avgmap == 'avg':
                nc_u = nc_out.createVariable('U500', 'f8', ('year',))
                nc_v = nc_out.createVariable('V500', 'f8', ('year',))
                for c in cyrs:                  # loop through each core
                    for yr in cyrs[c]:          # then through each year
                        nc_u[int(yr-1980)], nc_v[int(yr-1980)] = get_avgs(c, yr, 1, 2, 3, nc_u, nc_v, 'yr', 'avg', lonmin, lonmax, latmin, latmax)
                        print('Data written: ',nc_u[int(yr-1980)],nc_v[int(yr-1980)])
            elif avgmap == 'map':
                nc_out.createDimension('lat', None)  
                nc_out.createDimension('lon', None)  
                nc_u = nc_out.createVariable('U500', 'f8', ('year', 'lat', 'lon',))
                nc_v = nc_out.createVariable('V500', 'f8', ('year', 'lat', 'lon',))
                for c in cyrs:                  # loop through each core
                    for yr in cyrs[c]:          # then through each year
                        get_avgs(c, yr, 1, 2, 3, nc_u, nc_v, 'yr', 'map', lonmin, lonmax, latmin, latmax, month='yes')
                
def sn_plot(us, vs, yi, ye, sn, lat0=25, lat1=50):
    import matplotlib.pyplot as plt
    import numpy as np
    from scipy import stats
    
    u = us[(yi-1980):(ye-1980)]
    v = vs[(yi-1980):(ye-1980)]
    label_distu = 1.5*np.std(u)
    label_distv = 3*np.std(v)
    xlabs = np.arange(yi, ye, step=4)
    xlocs = np.arange(len(u), step=4)
    
    plt.figure(figsize=(6,6))
    xs = np.arange(len(u))
    plt.xticks(xlocs, xlabs)
    
    # fit a line to the datasets (linear regression):
    cu = np.polyfit(xs,u,1)
    mu, bu  = np.polyfit(xs,u,1)
    plt.text(np.mean(xs), np.mean(u)+label_distu, 'y = {:.4f}*x + {:.4f}'.format(mu,bu))
    polyu = np.poly1d(cu) # a function which takes in yrs and returns an estimate for u
    plt.plot(xs, u, 'b', label='U')
    plt.plot(xs, polyu(xs), '--m')
    cv = np.polyfit(xs,v,1)
    mv, bv  = np.polyfit(xs,v,1)
    plt.text(np.mean(xs), np.mean(v)+label_distv, 'y = {:.4f}*x + {:.4f}'.format(mv,bv))
    polyv = np.poly1d(cv) # a function which takes in yrs and returns an estimate for v
    plt.plot(xs, v, 'k', label='V')
    plt.plot(xs, polyv(xs), '--m')
    
    # do student's t-test for significance:
    tstatu, pvalu = stats.ttest_ind(xs,u)
    plt.text(np.mean(xs), np.mean(u)-2*label_distu, 'T stat={:.6f}, p={:.6f}'.format(tstatu,pvalu))
    tstatv, pvalv = stats.ttest_ind(xs,v)
    plt.text(np.mean(xs), np.mean(v)-label_distv, 'T stat={:.6f}, p={:.6f}'.format(tstatv,pvalv))
    
    plt.ylabel('Wind Speed')
    plt.xlabel('Year')
    plt.legend()
    #plt.tight_layout(pad=0.1)
    plt.title(f'{sn}')
    if lat0 != 25:
        plt.savefig(f'wind_avgs_{yi}_{ye}_{sn}_{lat0}_{lat1}_new.png')
    else:
        plt.savefig(f'wind_avgs_{yi}_{ye}_{sn}_new.png')
    plt.show()
    

def plot_avgs(file, yi=1980, ye=2019, yrsn='sn', lat0=25, lat1=50):
    from netCDF4 import Dataset
    import matplotlib.pyplot as plt
    import numpy as np
    from scipy import stats
    import pandas as pd
    
    data = Dataset(file, mode='r')
    
    u500 = data.variables['U500'][:]
    v500 = data.variables['V500'][:]
    
    if yrsn == 'sn':
        u = [sn for yr in u500 for sn in yr]
        v = [sn for yr in v500 for sn in yr]
        u0 = [yr[0] for yr in u500]  # DJF season
        v0 = [yr[0] for yr in v500]  
        u1 = [yr[1] for yr in u500]  # MAM season
        v1 = [yr[1] for yr in v500]
        u2 = [yr[2] for yr in u500]  # JJA season
        v2 = [yr[2] for yr in v500]
        u3 = [yr[3] for yr in u500]  # SON season
        v3 = [yr[3] for yr in v500]
        us = u[(yi-1980)*4:(ye-1980)*4]
        vs = v[(yi-1980)*4:(ye-1980)*4]
        label_distu = 1.5*np.std(us)
        label_distv = 3*np.std(vs)
        xlabs = np.arange(yi, ye, step=4)
        xlocs = np.arange(len(us), step=16)
        
    elif yrsn == 'yr':
        u = [yr for yr in u500]
        v = [yr for yr in v500]
        us = u[(yi-1980):(ye-1980)]
        vs = v[(yi-1980):(ye-1980)]
        label_distu = -5*np.std(us)
        label_distv = 5*np.std(vs)
        xlabs = np.arange(yi, ye, step=4)
        xlocs = np.arange(len(us), step=4)
    # turn the data lists into pandas dataframes to get rid of Nan values:
    datalists = [u0, u1, u2, u3, us, v0, v1, v2, v3, vs]
    for i in datalists:
        df = pd.DataFrame(i)
        notnan = df.dropna()
        i = notnan[0].values.tolist()
    
    sn_plot(u0, v0, yi, ye, 'DJF', lat0, lat1)
    sn_plot(u1, v1, yi, ye, 'MAM', lat0, lat1)
    sn_plot(u2, v2, yi, ye, 'JJA', lat0, lat1)
    sn_plot(u3, v3, yi, ye, 'SON', lat0, lat1)
    
    plt.figure(figsize=(6,6))
    xs = np.arange(len(us))
    plt.xticks(xlocs, xlabs)
    ''
    # fit a line to the datasets (linear regression):
    cu = np.polyfit(xs,us,1)
    mu, bu  = np.polyfit(xs,us,1)
    plt.text(np.mean(xs), np.mean(us)+label_distu, 'y = {:.4f}*x + {:.4f}'.format(mu,bu))
    polyu = np.poly1d(cu) # a function which takes in yrs and returns an estimate for u
    plt.plot(xs, us, 'k', label='U')
    plt.plot(xs, polyu(xs), '--m')
    cv = np.polyfit(xs,vs,1)
    mv, bv  = np.polyfit(xs,vs,1)
    plt.text(np.mean(xs), np.mean(vs)+label_distv, 'y = {:.4f}*x + {:.4f}'.format(mv,bv))
    polyv = np.poly1d(cv) # a function which takes in yrs and returns an estimate for v
    plt.plot(xs, vs, 'b', label='V')
    plt.plot(xs, polyv(xs), '--m')
    
    # do student's t-test for significance:
    tstatu, pvalu = stats.ttest_ind(xs,us)
    plt.text(np.mean(xs), np.mean(us)-label_distu, 'T stat={:.6f}, p={:.6f}'.format(tstatu,pvalu))
    tstatv, pvalv = stats.ttest_ind(xs,vs)
    plt.text(np.mean(xs), np.mean(vs)-label_distv, 'T stat={:.6f}, p={:.6f}'.format(tstatv,pvalv))
    
    ''' 
    # doesn't work for some reason?
    #from scipy.stats import linregress
    #linregress(xs,u) 
    from mlxtend.plotting import plot_linear_regression
    intu, slopeu, corr_coeffu = plot_linear_regression(xs, us)
    intv, slopev, corr_coeffv = plot_linear_regression(xs, vs)
    #plt.plot(xs, vs, color='b', alpha=0.5, label='V')
    #plt.plot(xs, us, color='r', alpha=0.5, label='U')
    ''
    # this one does work, but I think it's too complicated
    from sklearn.linear_model import LinearRegression
    # create a linear regression model
    model = LinearRegression()
    xs = np.array(xs).reshape((-1,1))  # reshape x values for model to work right
    model.fit(xs, us)
    # predict y from the data
    x_new = np.linspace(0, 156, 100)
    y_new = model.predict(x_new[:, np.newaxis])
    # plot the results
    ax = plt.axes()
    ax.plot(xs, us)
    ax.plot(x_new, y_new, '--k')
    '''
    #plt.plot(us, np.poly1d(np.polyfit(xs, us, 1))(us))
    
    plt.ylabel('Wind Speed')
    plt.xlabel('Year')
    plt.legend()
    plt.tight_layout(pad=0.1)
    if lat0 != 25:
        plt.savefig(f'wind_avgs_{yi}_{ye}_{yrsn}_{lat0}_{lat1}_new.png')
    else:
        plt.savefig(f'wind_avgs_{yi}_{ye}_{yrsn}_new.png')
    plt.show()
    
def map_plot(file, var, yr):
    import numpy as np
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    import matplotlib.pyplot as plt
    from netCDF4 import Dataset
    from mpl_toolkits.basemap import Basemap
    
    data = Dataset(file)
    u500 = data.variables[f'{var}'][yr-1980]
    djf = u500[0]
    mam = u500[1]
    jja = u500[2]
    son = u500[3]
    
    # create figure and axes:
    fig, axs = plt.subplots(2, 2, figsize=(20, 20), sharey='row', sharex='col')
    plt.subplots_adjust(wspace=0.15, hspace=0.15)
    
    # setup stereographic basemap:
    # make sure this is plotting the right lat/lon ranges
    # lat range : 25, 50
    # lon range : -65, -125 (west)
    m = Basemap(projection='merc', lat_0=37.5, lon_0=-95, llcrnrlon=-125, llcrnrlat=25,urcrnrlon=-65, urcrnrlat=50, resolution='l')
    #projection='stere',lat_0=37.5,lon_0=-95
    sns = [djf, mam, jja, son]
    sn_min = np.min(sns) 
    sn_max = np.max(sns)
    
    fig.suptitle(f'{yr} {var} (m/s)')
    
    m.ax = axs[0,0]
    m.drawparallels(np.arange(-80., 81., 10.), labels=[1,0,0,0], fontsize=10)
    m.drawmeridians(np.arange(-180., 181., 10.), labels=[0,0,0,1], fontsize=10)
    m.drawcoastlines()
    m.drawstates()
    m.drawcountries()
    pic = m.imshow(djf, origin='lower', vmin=sn_min, vmax=sn_max)
    '''
    nlats = 50; nlons = 96
    #delta = 2.*np.pi/(nlons-1)
    lats = (0.5*np.pi-0.5*np.indices((nlats,nlons))[0,:,:])
    lons = (0.625*np.indices((nlats,nlons))[1,:,:])
    x, y = m(lons*180./np.pi, lats*180./np.pi)
    pic = m.contour(x, y, djf, np.arange(-1,15), origin='lower', vmin=sn_min, vmax=sn_max)
    '''
    div = make_axes_locatable(axs[0,0])
    cax = div.append_axes("bottom", size="5%", pad="15%")
    cbar = fig.colorbar(pic, cax=cax, orientation="horizontal", ticks=[sn_min, sn_max])
    cbar.ax.tick_params(labelsize=9) 
    axs[0,0].set_title('DJF')
    
    m.ax = axs[0,1]
    m.drawparallels(np.arange(-80., 81., 10.), labels=[1,0,0,0], fontsize=10)
    m.drawmeridians(np.arange(-180., 181., 10.), labels=[0,0,0,1], fontsize=10)
    m.drawcoastlines()
    m.drawstates()
    m.drawcountries()
    pic = m.imshow(mam, origin='lower', vmin=sn_min, vmax=sn_max)
    '''
    nlats = 50; nlons = 96
    #delta = 2.*np.pi/(nlons-1)
    lats = (0.5*np.pi-0.5*np.indices((nlats,nlons))[0,:,:])
    lons = (0.625*np.indices((nlats,nlons))[1,:,:])
    x, y = m(lons*180./np.pi, lats*180./np.pi)
    pic = m.contour(x, y, mam, np.arange(-1,15), origin='lower', vmin=sn_min, vmax=sn_max)
    '''
    div = make_axes_locatable(axs[0,1])
    cax = div.append_axes("bottom", size="5%", pad="15%")
    cbar = fig.colorbar(pic, cax=cax, orientation="horizontal", ticks=[sn_min, sn_max])
    cbar.ax.tick_params(labelsize=9) 
    axs[0,1].set_title('MAM')
    
    m.ax = axs[1,0]
    m.drawparallels(np.arange(-80., 81., 10.), labels=[1,0,0,0], fontsize=10)
    m.drawmeridians(np.arange(-180., 181., 10.), labels=[0,0,0,1], fontsize=10)
    m.drawcoastlines()
    m.drawstates()
    m.drawcountries()
    pic = m.imshow(jja, origin='lower', vmin=sn_min, vmax=sn_max)
    '''
    nlats = 50; nlons = 96
    #delta = 2.*np.pi/(nlons-1)
    lats = (0.5*np.pi-0.5*np.indices((nlats,nlons))[0,:,:])
    lons = (0.625*np.indices((nlats,nlons))[1,:,:])
    x, y = m(lons*180./np.pi, lats*180./np.pi)
    pic = m.contour(x, y, jja, np.arange(-1,15), origin='lower', vmin=sn_min, vmax=sn_max)
    '''
    div = make_axes_locatable(axs[1,0])
    cax = div.append_axes("bottom", size="5%", pad="15%")
    cbar = fig.colorbar(pic, cax=cax, orientation="horizontal", ticks=[sn_min, sn_max])
    cbar.ax.tick_params(labelsize=9) 
    axs[1,0].set_title('JJA')
    
    m.ax = axs[1,1]
    m.drawparallels(np.arange(-80., 81., 10.), labels=[1,0,0,0], fontsize=10)
    m.drawmeridians(np.arange(-180., 181., 10.), labels=[0,0,0,1], fontsize=10)
    m.drawcoastlines()
    m.drawstates()
    m.drawcountries()
    pic = m.imshow(son, origin='lower', vmin=sn_min, vmax=sn_max)
    '''
    nlats = 50; nlons = 96
    #delta = 2.*np.pi/(nlons-1)
    lats = (0.5*np.pi-0.5*np.indices((nlats,nlons))[0,:,:])
    lons = (0.625*np.indices((nlats,nlons))[1,:,:])
    x, y = m(lons*180./np.pi, lats*180./np.pi)
    pic = m.contour(x, y, son, np.arange(-1,15), origin='lower', vmin=sn_min, vmax=sn_max)
    '''
    cax = fig.add_axes([.13, 0.07, .08, .01])
    div = make_axes_locatable(axs[1,1])
    cax = div.append_axes("bottom", size="5%", pad="15%")
    cbar = fig.colorbar(pic, cax=cax, orientation="horizontal", ticks=[sn_min, sn_max])
    cbar.ax.tick_params(labelsize=9) 
    axs[1,1].set_title('SON')
    
    #plt.tight_layout(pad=0.1)
    fig.savefig(f'{yr}_{var}_snavgmaps_change.png')
    plt.show()
    
def map_plot1(grid, var, yr, season=4, lonmin=-65, lonmax=-125, latmin=25, latmax=50, grid_change=0):
    import numpy as np
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    import matplotlib.pyplot as plt
    from mpl_toolkits.basemap import Basemap, cm
    
    if season == 0:
        sn = 'DJF'
    elif season == 1:
        sn = 'MAM'
    elif season == 2:
        sn = 'JJA'
    elif season == 3:
        sn = 'SON'
    else:
        sn = ''
    
    # create figure and axes:
    fig, axs = plt.subplots(1, 1, figsize=(15, 10))
    
    # setup stereographic basemap:
    if latmax < 90:
        m = Basemap(projection='merc', llcrnrlon=lonmax, llcrnrlat=latmin, 
                    urcrnrlon=lonmin, urcrnrlat=latmax, resolution='l')
        m.ax = axs
        sn_min = np.min(grid) 
        sn_max = np.max(grid)
        ch_min = np.min(grid_change) 
        ch_max = np.max(grid_change)
        pic = m.imshow(grid_change, origin='lower', vmin=ch_min, vmax=ch_max)
        ny = grid.shape[0]; nx = grid.shape[1]
        lons, lats = m.makegrid(nx,ny)
        x, y = m(lons,lats)
        # plot the contours for all of the grid data as dashed lines:
        con = m.contour(x,y,grid,cmap=cm.s3pcpn, vmin=sn_min, vmax=0, linestyles='dashdot')
        plt.clabel(con, inline=True, fmt='%1.0f', fontsize=10)
        # plot the positive contours as solid lines:
        con1 = m.contour(x,y,grid,cmap=cm.s3pcpn, vmin=0, vmax=sn_max, linestyles='solid')
        plt.clabel(con1, inline=True, fmt='%1.0f', fontsize=10)
        m.drawparallels(np.arange(-90., 91., 10.), labels=[1,0,0,0], fontsize=12)
        m.drawmeridians(np.arange(-180., 181., 10.), labels=[0,0,0,1], fontsize=12)
        #pic = m.imshow(grid, origin='lower', vmin=sn_min, vmax=sn_max)
        m.drawcoastlines()
        m.drawstates()
        m.drawcountries()
    elif latmax == 90:
        m = Basemap(projection='cyl')
        m.ax = axs
        sn_min = np.min(grid) 
        sn_max = np.max(grid)
        ch_min = np.min(grid_change) 
        ch_max = np.max(grid_change)
        pic = m.imshow(grid_change, origin='lower', vmin=ch_min, vmax=ch_max)
        ny = grid.shape[0]; nx = grid.shape[1]
        lons, lats = m.makegrid(nx,ny)
        x, y = m(lons,lats)
        # plot the contours for all of the grid data as dashed lines:
        con = m.contour(x,y,grid,cmap=cm.s3pcpn, vmin=sn_min, vmax=0, linestyles='dashdot')
        plt.clabel(con, inline=True, fmt='%1.0f', fontsize=10)
        # plot the positive contours as solid lines:
        con1 = m.contour(x,y,grid,cmap=cm.s3pcpn, vmin=0, vmax=sn_max, linestyles='solid')
        plt.clabel(con1, inline=True, fmt='%1.0f', fontsize=10)
        m.drawparallels(np.arange(-90., 91., 30.), labels=[1,0,0,0], fontsize=12)
        m.drawmeridians(np.arange(-180., 181., 30.), labels=[0,0,0,1], fontsize=12)
        m.drawcoastlines()
        m.drawcountries()
    
    div = make_axes_locatable(axs)
    cax = div.append_axes("bottom", size="5%", pad="15%")
    cbar = fig.colorbar(pic, cax=cax, orientation="horizontal", ticks=[ch_min, 0, ch_max])
    cbar.ax.tick_params(labelsize=12) 
    axs.set_title(f'{yr} {sn} {var} (m/s)', fontsize=14)
    #plt.tight_layout(pad=0.1)
    fig.savefig(f'{yr}_{var}_{sn}_avgmap_lat{latmin}_{latmax}.png')
    plt.show()
    
# this code gets the trend of two variables over a season, from 1980 to 2019 (or to a year you pick: yr_pick)
# you can specify the level you want with lev (16 is the index for 500hPa)
# you can plot an estimate of each coordinate's total change since 1980 (trend='total_change')
#   or the slopes of each coordinate's regression since 1980 (trend='slope')
def get_trend(v1, v2, v1_label, v2_label, lev, sn, lonmin=180, lonmax=-180, latmin=-90, latmax=90, yr_pick='no', trend='total_change'):
    # for CONUS: lonmin=-65, lonmax=-125, latmin=25, latmax=50 
    from netCDF4 import Dataset
    import numpy as np
    
    cyrs = {100:np.arange(1980,1992), 200:np.arange(1992,2001), 300:np.arange(2001,2011), 400:np.arange(2011,2019)}
    sms = {0: ['12', '01', '02'], 1: ['03', '04', '05'], 2: ['06', '07', '08'], 3: ['09', '10', '11']}
    u_avgs = []
    v_avgs = []
    for c in cyrs:                  # loop through each core
        for yr in cyrs[c]:          # then through each year
            files = files_list_monthly(c, yr, sms[sn][0], sms[sn][1], sms[sn][2])
            #print(files)
            u = []
            v = []
            
            for fi in files:
                fh = Dataset(fi, mode='r')
                lon = fh.variables['lon'][:]    # longitude, 576 dimension: 0.625 deg
                lat = fh.variables['lat'][:]    # latitude, 361 dimension: 0.5 deg
                u500 = fh.variables[f'{v1}'][0,lev]  # eastward wind (0=time, 16=index for 500hPa)
                v500 = fh.variables[f'{v2}'][0,lev]  # northward wind (0=time, 16=index for 500hPa)
                
                # CONUS: 25->50 deg latitude, -65->-125 deg longitude
                # the variables are indexed as var(time->24,lat->361,lon->576)
                lon = np.asarray(lon)
                lat = np.asarray(lat)
                ind_lonmin = (np.abs(lon -lonmin)).argmin()   # index for minimum longitude (default is -65)
                ind_lonmax = (np.abs(lon -lonmax)).argmin()   # index for maximum longitude (default is -125)
                ind_latmin = (np.abs(lat - latmin)).argmin()  # index for minimum latitude (default is 25)
                ind_latmax = (np.abs(lat - latmax)).argmin()  # index for maximum latitude (default is 50)
                lat_range = np.arange(int(ind_latmin), int(ind_latmax))
                lon_range = np.arange(int(ind_lonmax), int(ind_lonmin))  # list them backwards (index for -125 is bigger than index for -65)
            
                u5 = [u500[la,lo] for la in lat_range for lo in lon_range]
                v5 = [v500[la,lo] for la in lat_range for lo in lon_range]
                u5 = np.array(u5)
                ui = np.reshape(u5, (len(lat_range), len(lon_range)))
                #ui = np.abs(ui)   # try taking absval before averaging arrays
                u.append(ui)
                v5 = np.array(v5)
                vi = np.reshape(v5, (len(lat_range), len(lon_range)))
                #vi = np.abs(vi)   # try taking absval before averaging arrays
                v.append(vi)
            
            ugrids = np.array(u)
            #print(ugrids)
            vgrids = np.array(v)
            uavgd = np.mean(ugrids,axis=0)
            #print('U Avgd: ', uavgd)
            vavgd = np.mean(vgrids,axis=0)
            u_avgs.append(uavgd)
            v_avgs.append(vavgd)
            
    u_grids = np.array(u_avgs)
    #print('U_avgs: ', u_grids)
    v_grids = np.array(v_avgs)
    # these are the seasonal averages from 1980-2019:
    print('U Sn Avgs: ',u_grids.shape)
    print('V Sn Avgs: ',v_grids.shape)
    # these are the averages of the seasonal averages from 1980-2019:
    u_avg = np.mean(u_grids,axis=0)
    v_avg = np.mean(v_grids,axis=0)
    
    # to compute change map: 
    times = np.arange(u_grids.shape[0])          # of days x 24hrs
    uslope = []
    vslope = []
    for la in range(u_grids.shape[1]):           # 2D dimension (latitude)
        for lo in range(u_grids.shape[2]):       # 2D dimension (longitude)
            ut = []                             # holder for each pixel's time series
            for t in range(u_grids.shape[0]):    # 3rd dimension (# of 2D arrays)
                ut.append(u_grids[t,la,lo])        
            #print(ut)
            mu, bu = np.polyfit(times, ut, 1)   # fits a line to the pixel's time series
            #print(mu)
            uslope.append(mu)
    for la in range(v_grids.shape[1]):
        for lo in range(v_grids.shape[2]):
            vt = []
            for t in range(v_grids.shape[0]):
                vt.append(v_grids[t,la,lo])        
            #print(vt)
            mv, bv = np.polyfit(times, vt, 1)
            #print(mv)
            vslope.append(mv)
    uslope = np.array(uslope)
    uslope = np.reshape(uslope, (len(lat_range), len(lon_range)))
    vslope = np.array(vslope)
    vslope = np.reshape(vslope, (len(lat_range), len(lon_range)))
    
    if yr_pick == 'no':
        map_plot1(u_avg, f'{v1_label}', 'Trend', sn, lonmin, lonmax, latmin, latmax, grid_change=uslope)
        map_plot1(v_avg, f'{v2_label}', 'Trend', sn, lonmin, lonmax, latmin, latmax, grid_change=vslope)
        map_plot1(uslope, f'{v1_label} Change', 'Trend', sn, lonmin, lonmax, latmin, latmax, grid_change=u_avg)
        map_plot1(vslope, f'{v2_label} Change', 'Trend', sn, lonmin, lonmax, latmin, latmax, grid_change=v_avg)
    # if you've picked a year to end at:
    elif yr_pick != 'no':
        # if you want to plot an estimate of each coordinate's total change since 1980:
        if trend == 'total_change':
            map_plot1(u_avg, f'{v1_label}', yr_pick, sn, lonmin, lonmax, latmin, latmax, grid_change=uslope*(yr_pick-1980))
            map_plot1(v_avg, f'{v2_label}', yr_pick, sn, lonmin, lonmax, latmin, latmax, grid_change=vslope*(yr_pick-1980))
            map_plot1(uslope*(yr_pick-1980), f'{v1_label} Change', yr_pick, sn, lonmin, lonmax, latmin, latmax, grid_change=u_avg)
            map_plot1(vslope*(yr_pick-1980), f'{v2_label} Change', yr_pick, sn, lonmin, lonmax, latmin, latmax, grid_change=v_avg)
        # if you just want to plot the slope of each coordinate's regression:
        elif trend == 'slope':
            map_plot1(u_avg, f'{v1_label}', yr_pick, sn, lonmin, lonmax, latmin, latmax, grid_change=uslope)
            map_plot1(v_avg, f'{v2_label}', yr_pick, sn, lonmin, lonmax, latmin, latmax, grid_change=vslope)
            map_plot1(uslope, f'{v1_label} Change', yr_pick, sn, lonmin, lonmax, latmin, latmax, grid_change=u_avg)
            map_plot1(vslope, f'{v2_label} Change', yr_pick, sn, lonmin, lonmax, latmin, latmax, grid_change=v_avg)
    