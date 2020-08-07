"""
Created on Wed Jun 24 16:17:46 2020

@author: Anna Murphree (annammurphree@gmail.com)
For Summer 2020 Internship with Dr. Yaping Zhou

This code opens statistics files for EPEs and makes a quick histogram of a chosen statistic.
(Refer to Event_parameter_list.docx for parameter indices.) 
"""
def stat_histo(statfile, statname, statind):
    from netCDF4 import Dataset
    import numpy as np
    import matplotlib.pyplot as plt
    
    # open the event stats from IMERG:
    stats = Dataset(f'{statfile}', mode='r')
    #variables = stats.variables.keys()
    #print(variables)
    #print(stats.dimensions)
    params = stats.variables['event_parameters'][:]  
    
    data = []
    
    for event in params:
        '''
        dstart = event[0]   # start time (days) of the year?
        dend = event[1]     # end time (days) of the year?
        year = event[22]    # year of event
        #print('Dates: ',dstart,'-', dend, ',',year)
        long = event[2]     # longitude (deg East)
        #print('IMERG long: ',long)
        lati = event[3]     # latitude (deg North)
        #print('IMERG lat: ',lati)
        area = event[7]      # total area (km^2)
        dura = event[4]      # event duration (hrs)
        dura_h = event[5]    # event duration of heavy rain only (hrs)
        vol = event[6]       # event heavy rain volume (m^3)
        '''
        data.append(event[statind])
        
    stats.close()
    #print(len(data))
    #print(data)
    
    plt.hist(x=data, bins=100, range=(0, np.max(data)), color='#0504aa', alpha=0.7)
    plt.grid(axis='y', alpha=0.75)
    plt.xlabel('Value')
    plt.ylabel('Frequency')
    plt.title(f'{statname}')
    plt.show()
 
#stat_histo('extreme_event_stats_2017_MAM.nc', 'Duration', 4)
