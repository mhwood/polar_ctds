
import os
import requests
import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
import gsw
import argparse
import netCDF4 as nc4

def wood_manual_corrections(L2_file_name,max_depth):

    manual_chop_list = ['OMG_Ocean_CTD_L2_20150726022850.nc',
        # 2016
        'OMG_Ocean_AXCTD_L2_20160915163322.nc',
        'OMG_Ocean_AXCTD_L2_20161001142632.nc',
        'OMG_Ocean_AXCTD_L2_20161006123413.nc',
        'OMG_Ocean_AXCTD_L2_20161008123052.nc',
        'OMG_Ocean_AXCTD_L2_20160928153504.nc',
        # 2017
        'OMG_Ocean_AXCTD_L2_20171018122743.nc',
        'OMG_Ocean_AXCTD_L2_20171014125231.nc',
        'OMG_Ocean_AXCTD_L2_20171014140506.nc',
        'OMG_Ocean_AXCTD_L2_20171016142122.nc',
        'OMG_Ocean_AXCTD_L2_20171018122743.nc',
        'OMG_Ocean_AXCTD_L2_20171019125245.nc',
        # 2018
        'OMG_Ocean_AXCTD_L2_20180822145112.nc',
        'OMG_Ocean_AXCTD_L2_20180823150228.nc',
        # 2019
        'OMG_Ocean_AXCTD_L2_20190815132657.nc',
        'OMG_Ocean_AXCTD_L2_20190828174923.nc',
        'OMG_Ocean_AXCTD_L2_20190815165721.nc',
        'OMG_Ocean_AXCTD_L2_20190814173114.nc',
        # 2020
        'OMG_Ocean_AXCTD_L2_20200825171120.nc']

    manual_chop_dict = {'OMG_Ocean_CTD_L2_20150726022850.nc': 200,
         # 2016
        'OMG_Ocean_AXCTD_L2_20160915163322.nc': 921,
        'OMG_Ocean_AXCTD_L2_20161001142632.nc': 654,
        'OMG_Ocean_AXCTD_L2_20161006123413.nc': 415,
        'OMG_Ocean_AXCTD_L2_20161008123052.nc': 310,
        'OMG_Ocean_AXCTD_L2_20160928153504.nc': 581,
        # 2017
        'OMG_Ocean_AXCTD_L2_20171014125231.nc': 260,
        'OMG_Ocean_AXCTD_L2_20171014140506.nc': 362,
        'OMG_Ocean_AXCTD_L2_20171016142122.nc': 926,
        'OMG_Ocean_AXCTD_L2_20171018122743.nc': 309,
        'OMG_Ocean_AXCTD_L2_20171019125245.nc': 500,
        # 2018
        'OMG_Ocean_AXCTD_L2_20180822145112.nc': 650,
        'OMG_Ocean_AXCTD_L2_20180823150228.nc': 597,
        # 2019
        'OMG_Ocean_AXCTD_L2_20190815132657.nc': 322,
        'OMG_Ocean_AXCTD_L2_20190828174923.nc': 44,
        'OMG_Ocean_AXCTD_L2_20190815165721.nc': 553,
        'OMG_Ocean_AXCTD_L2_20190814173114.nc': 252,
        # 2020
        'OMG_Ocean_AXCTD_L2_20200825171120.nc': 752}

    if L2_file_name in manual_chop_dict.keys():
        max_depth = manual_chop_dict[L2_file_name]

    return(max_depth)

def bin_ctd_data(depth,potential_temperature,practical_salinity):

    max_depth = np.max(depth)
    binned_depth = np.arange(1,int(max_depth) + 1) #ignore the 1m surface layer - data is not trustworthy here
    binned_potential_temperature = np.zeros_like(binned_depth).astype(float)
    binned_practical_salinity = np.zeros_like(binned_depth).astype(float)
    binning_count = np.zeros_like(binned_depth).astype(int)
    for di in range(len(binned_depth) - 1):
        depth_indices = np.logical_and(depth >= binned_depth[di], depth < binned_depth[di + 1])
        if np.sum(depth_indices)>0:
            binned_potential_temperature[di] = np.mean(potential_temperature[depth_indices])
            binned_practical_salinity[di] = np.mean(practical_salinity[depth_indices])
            binning_count[di] = np.sum(depth_indices)
        # else:
        #     print('    No data for depth '+str(binned_depth[di]))

    # remove bins that didnt have any data
    has_data = binning_count > 0
    binned_depth = binned_depth[has_data]
    binned_potential_temperature = binned_potential_temperature[has_data]
    binned_practical_salinity = binned_practical_salinity[has_data]

    #filter out some ridiculous values
    salt_indices = np.logical_and(binned_practical_salinity>10,binned_practical_salinity<45)
    binned_depth = binned_depth[salt_indices]
    binned_potential_temperature = binned_potential_temperature[salt_indices]
    binned_practical_salinity = binned_practical_salinity[salt_indices]

    return(binned_depth,binned_potential_temperature,binned_practical_salinity)

def output_file_as_nc(output_file, lon, lat,
                      year,month,day,hour,minute,
                      binned_depth,binned_potential_temperature,binned_practical_salinity,source,project_name):

    ds = nc4.Dataset(output_file,'w')
    ds.createDimension('records', len(binned_depth))
    dep = ds.createVariable('depth', 'f4', 'records')
    dep[:] = binned_depth
    salt = ds.createVariable('practical_salinity', 'f4', 'records')
    salt[:] = binned_practical_salinity
    ptemp = ds.createVariable('potential_temperature', 'f4', 'records')
    ptemp[:] = binned_potential_temperature

    ds.longitude = lon
    ds.latitude = lat
    ds.year = year
    ds.month = month
    ds.day = day
    ds.hour = hour
    ds.minute = minute
    ds.source = source
    ds.project_name = project_name

    ds.close()

def output_file_as_csv(output_file, lon, lat,
                      year,month,day,hour,minute,
                      binned_depth,binned_potential_temperature,binned_practical_salinity,source,project_name):
    raise ValueError('Did not implement csv output yet')

def clean_and_bin_CTDs(data_dir, project_dir, file_path, output_type, sources, printing):

    ########################################################################################################################
    # Read in the list of CTDs

    f = open(file_path)
    lines = f.read()
    f.close()
    lines = lines.split('\n')
    lines.pop(0)
    axctd_set = []
    ctd_set = []
    for line in lines:
        line = line.split(',')
        if len(line) > 2:
            if line[1]=='AXCTDs':
                axctd_set.append(line[0])
            if line[1]=='CTDs':
                ctd_set.append(line[0])

    ########################################################################################################################
    # This section is for the CTDs

    if 'CTDs' in sources and len(ctd_set)>0:

        counter = 0

        suspect_files = ['OMG_Ocean_CTD_L2_20160913091721.nc', 'OMG_Ocean_CTD_L2_20150822141405.nc']

        for file_ID in ctd_set:
            year = file_ID[:4]

            if year not in os.listdir(project_dir):
                os.mkdir(os.path.join(project_dir, year))

            counter += 1
            if printing:
                print('     Processing CTD file ' + file_ID + ' (' + str(counter) + ' of ' + str(len(ctd_set))+ ')')

            ctd_file = 'OMG_Ocean_CTD_L2_'+file_ID.split('_')[0]+file_ID.split('_')[1]+'.nc'

            if ctd_file not in suspect_files:

                # read in the data
                ds = xr.open_dataset(os.path.join(data_dir,'PODAAC','Data','CTDs',year,ctd_file))
                lat = float(ds['lat'])
                lon = float(ds['lon'])
                depth = np.array(ds['depth'][:])
                temperature = np.array(ds['temperature'][:])
                conductivity = np.array(ds['conductivity'][:])
                depth = np.reshape(depth, (np.size(depth),))
                temperature = np.reshape(temperature, (np.size(temperature),))
                conductivity = np.reshape(conductivity, (np.size(conductivity),))

                # remove nans if they exists
                if np.sum(np.isnan(depth))>0:
                    conductivity = conductivity[~np.isnan(depth)]
                    temperature = temperature[~np.isnan(depth)]
                    depth = depth[~np.isnan(depth)]
                raw_depth = np.copy(depth)

                # if a bunch of values are chopped from the down cast, flip the order so the upcast is used
                if depth[0]>100:
                    depth=np.flip(depth)
                    conductivity=np.flip(conductivity)
                    temperature=np.flip(temperature)

                # make some conversions using the gsw toolbox
                pressure = gsw.conversions.p_from_z(-1*depth,lat)
                practical_salinity = gsw.conversions.SP_from_C(conductivity,temperature,pressure)
                absolute_salinity = gsw.conversions.SA_from_SP(practical_salinity,pressure,lon,lat)
                potential_temperature = gsw.conversions.pt_from_t(absolute_salinity,temperature,pressure,0)

                # many of these CTDs hit the bottom and become muddy, causing issues on the up-cast
                # to avoid data degradation, only take the down cast
                max_depth_index = np.argmax(depth)
                max_depth = int(np.max(depth)) - 1  # go one meter from the bottom
                max_depth = wood_manual_corrections(ctd_file, max_depth)
                indices = depth<max_depth
                indices[max_depth_index:]=0
                potential_temperature=potential_temperature[indices]
                practical_salinity=practical_salinity[indices]
                depth = depth[indices]

                # bin the profiles into 1-meter bins
                binned_depth,binned_potential_temperature,binned_practical_salinity = \
                    bin_ctd_data(depth, potential_temperature, practical_salinity)

                # # output the data into a conveniently loopable csv file
                # time = ctd_file.split('_')[-1][:-3]
                # output_file = os.path.join(project_dir,'Processed','Data','CTDs',year,'CTD_'+time[:8]+'_'+time[8:]+'.csv')
                # binned_output_array = output_binned_ctd_csv(output_file, lon, lat, time, binned_practical_salinity,
                #                       binned_potential_temperature, binned_depth)

                year = int(file_ID.split('_')[0][:4])
                month = int(file_ID.split('_')[0][4:6])
                day = int(file_ID.split('_')[0][6:8])
                hour = int(file_ID.split('_')[1][:2])
                minute = int(file_ID.split('_')[1][2:])

                if output_type=='csv':
                    output_file = os.path.join(project_dir, str(year),
                                               'CTD_' + file_ID.split('_')[0] + '_' + file_ID.split('_')[1] + '.csv')
                    output_file_as_csv(output_file, lon, lat,
                                      year, month, day, hour, minute,
                                      binned_depth,binned_potential_temperature,binned_practical_salinity, ctd_file, 'Oceans Melting Greenland')

                if output_type=='nc':
                    output_file = os.path.join(project_dir, str(year),
                                               'CTD_' + file_ID.split('_')[0] + '_' + file_ID.split('_')[1] + '.nc')
                    output_file_as_nc(output_file, lon, lat,
                                      year, month, day, hour, minute,
                                      binned_depth,binned_potential_temperature,binned_practical_salinity, ctd_file, 'Oceans Melting Greenland')

                # # make a plot if desired
                # if plotting:
                #     output_file = os.path.join(project_dir, 'Processed', 'Plots', 'CTDs',year,
                #                                'CTD_' + time[:8] + '_' + time[8:] + '.png')
                #     create_output_plot(output_file, raw_depth, temperature, conductivity, binned_output_array)

            if ctd_file in suspect_files:
                print('    Skipped '+str(ctd_file)+' because its suspect')

    ########################################################################################################################
    # This section is for the AXCTDs



    if 'AXCTDs' in sources and len(axctd_set)>0:

        counter = 0

        suspect_files = ['OMG_Ocean_AXCTD_L2_20171015104006.nc',
                         'OMG_Ocean_AXCTD_L2_20171019145906.nc',
                         'OMG_Ocean_AXCTD_L2_20190823092743.nc',
                         'OMG_Ocean_AXCTD_L2_20160915131214.nc',
                         'OMG_Ocean_AXCTD_L2_20160916142120.nc',
                         'OMG_Ocean_AXCTD_L2_20190829132903.nc',
                         'OMG_Ocean_AXCTD_L2_20190904132624.nc',
                         # 2018
                         'OMG_Ocean_AXCTD_L2_20180824121607.nc',
                         'OMG_Ocean_AXCTD_L2_20180910142654.nc']

        for file_ID in axctd_set:
            year = file_ID[:4]

            if year not in os.listdir(project_dir):
                os.mkdir(os.path.join(project_dir, year))

            counter += 1
            if printing:
                print('     Processing AXCTD file ' + file_ID + ' (' + str(counter) + ' of ' + str(len(axctd_set))+ ')')

            axctd_file = 'OMG_Ocean_AXCTD_L2_'+file_ID.split('_')[0]+file_ID.split('_')[1]+'.nc'

            if axctd_file not in suspect_files:

                # read in the data
                ds = xr.open_dataset(
                    os.path.join(data_dir, 'PODAAC', 'Data', 'AXCTDs', year, axctd_file))
                lat = float(ds['lat'])
                lon = float(ds['lon'])
                depth = np.array(ds['depth'][:])
                temperature = np.array(ds['temperature'][:])
                conductivity = np.array(ds['conductivity'][:])
                depth = np.reshape(depth, (np.size(depth),))
                temperature = np.reshape(temperature, (np.size(temperature),))
                conductivity = np.reshape(conductivity, (np.size(conductivity),))

                # remove nans if they exists
                if np.sum(np.isnan(depth)) > 0:
                    conductivity = conductivity[~np.isnan(depth)]
                    temperature = temperature[~np.isnan(depth)]
                    depth = depth[~np.isnan(depth)]
                raw_depth = np.copy(depth)

                # make some conversions using the gsw toolbox
                pressure = gsw.conversions.p_from_z(-1 * depth, lat)
                practical_salinity = gsw.conversions.SP_from_C(conductivity, temperature, pressure)
                absolute_salinity = gsw.conversions.SA_from_SP(practical_salinity, pressure, lon, lat)
                potential_temperature = gsw.conversions.pt_from_t(absolute_salinity, temperature,
                                                                  pressure, 0)

                # many of these CTDs hit the bottom and become muddy,
                # to avoid data degradation, only take the down cast
                max_depth_index = np.argmax(depth)
                max_depth = int(np.max(depth)) - 1  # go one meter from the bottom
                #use_file, wood_manual_corrections(L2_file_name, profile)
                max_depth = wood_manual_corrections(axctd_file,max_depth)
                indices = depth < max_depth
                indices[max_depth_index:] = 0
                potential_temperature = potential_temperature[indices]
                practical_salinity = practical_salinity[indices]
                depth = depth[indices]

                # bin the profiles into 1-meter bins
                binned_depth, binned_potential_temperature, binned_practical_salinity = \
                    bin_ctd_data(depth, potential_temperature, practical_salinity)

                year = int(file_ID.split('_')[0][:4])
                month = int(file_ID.split('_')[0][4:6])
                day = int(file_ID.split('_')[0][6:8])
                hour = int(file_ID.split('_')[1][:2])
                minute = int(file_ID.split('_')[1][2:])

                if output_type=='csv':
                    output_file = os.path.join(project_dir, str(year),
                                               'CTD_' + file_ID.split('_')[0] + '_' + file_ID.split('_')[1] + '.csv')
                    output_file_as_csv(output_file, lon, lat,
                                      year, month, day, hour, minute,
                                      binned_depth,binned_potential_temperature,binned_practical_salinity, axctd_file, 'Oceans Melting Greenland')

                if output_type=='nc':
                    output_file = os.path.join(project_dir, str(year),
                                               'CTD_' + file_ID.split('_')[0] + '_' + file_ID.split('_')[1] + '.nc')
                    output_file_as_nc(output_file, lon, lat,
                                      year, month, day, hour, minute,
                                      binned_depth,binned_potential_temperature,binned_practical_salinity, axctd_file, 'Oceans Melting Greenland')


                # # make a plot if desired
                # if plotting:
                #     output_file = os.path.join(project_dir, 'Processed', 'Plots', 'AXCTDs', year,
                #                                'CTD_' + time[:8] + '_' + time[8:] + '.png')
                #     create_output_plot(output_file, raw_depth, temperature, conductivity,
                #                        binned_output_array)
            if axctd_file in suspect_files:
                print('    Skipped '+str(ctd_file)+' because its suspect')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--data_directory", action="store",
                        help="Directory where CTD data will be stored",
                        dest="data_dir", type=str, required=True)

    parser.add_argument("-o", "--output_directory", action="store",
                        help="Directory where metadata and processed data will be stored",
                        dest="project_dir", type=str, required=True)

    parser.add_argument("-f", "--file_path", action="store",
                        help="File with list of CTDs to process",
                        dest="file_path", type=str, required=True)

    parser.add_argument("-t", "--output_type", action="store",
                        help="The type of data to output (options: csv, nc).", dest="output_type",
                        type=str, required=True)

    parser.add_argument("-s", "--sources", action="store",
                        help="(Optional) List of sources. A value of all will choose both CTDs and AXCTDs. "
                             "Default value is all.", default='all', dest="sources", type=str, nargs='+',
                        required=False)

    parser.add_argument("-r", "--printing", action="store",
                        help="Choose whether to print output status messages. Default is 1 (true).",
                        default=1,
                        dest="printing_int", type=int, required=False)


    args = parser.parse_args()

    data_dir = args.data_dir
    project_dir = args.project_dir
    file_path = args.file_path
    output_type = args.output_type
    sources = args.sources

    printing_int = args.printing_int
    if printing_int == 1:
        printing = True
    else:
        printing = False

    if sources=='all':
        sources = ['CTDs','AXCTDs']

    print('Processing CTDs from the OMG Campaign')
    print('Processing parameters: ')
    print('    Data directory:   ' + data_dir)
    print('    Output directory: '+project_dir)
    print('    CTD Metadata File:    ' + str(file_path))
    print('    Print Messages:   ' + str(printing))
    print(' ')

    if printing:
        print(' Processing the CTD data (removing up-casts, bottom hits, obviously bad data, and binning)')
    clean_and_bin_CTDs(data_dir, project_dir, file_path, output_type, sources, printing)
