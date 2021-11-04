
import os
import numpy as np
import argparse
import netCDF4 as nc4

def output_file_as_nc(output_file, lon, lat,
                      year,month,day,hour,minute,
                      profile,source,project_name):

    ds = nc4.Dataset(output_file,'w')
    ds.createDimension('records', np.shape(profile)[0])
    dep = ds.createVariable('depth', 'f4', 'records')
    dep[:] = profile[:,2]
    salt = ds.createVariable('practical_salinity', 'f4', 'records')
    salt[:] = profile[:,1]
    ptemp = ds.createVariable('potential_temperature', 'f4', 'records')
    ptemp[:] = profile[:,0]

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

def format_individual_ctds(output_dir,data_dir,file_path,output_type,version,correction,output_years):

    print('Formatting individual data from the Hadley EN Dataset')
    print('    Version: '+version)
    print('    Correction: '+correction)

    # read in the file list
    f = open(file_path)
    lines = f.read()
    f.close()
    lines = lines.split('\n')
    lines.pop(0)

    indices = np.zeros((len(lines),)).astype(int)
    years = np.zeros((len(lines),)).astype(int)
    months = np.zeros((len(lines),)).astype(int)
    days = np.zeros((len(lines),)).astype(int)
    hours = np.zeros((len(lines),)).astype(int)
    minutes = np.zeros((len(lines),)).astype(int)
    longitudes = np.zeros((len(lines),))
    latitudes = np.zeros((len(lines),))
    year_months = []

    for ll in range(len(lines)):
        index = lines[ll].split(',')[1]
        year = lines[ll].split(',')[2]
        month = lines[ll].split(',')[3]
        day = lines[ll].split(',')[4]
        hour = lines[ll].split(',')[5]
        minute = lines[ll].split(',')[6]
        longitude = lines[ll].split(',')[7]
        latitude = lines[ll].split(',')[8]

        indices[ll] = index
        years[ll] = year
        months[ll] = month
        days[ll] = day
        hours[ll] = hour
        minutes[ll] = minute
        longitudes[ll] = longitude
        latitudes[ll] = latitude

        if str(year)+'-'+str(month) not in year_months:
            year_months.append(str(year)+'-'+str(month))

    for year_month in year_months:
        yr = int(year_month.split('-')[0])
        if yr>=np.min(np.array(output_years)) and yr<=np.max(np.array(output_years)):
            if str(yr) not in os.listdir(output_dir):
                os.mkdir(os.path.join(output_dir,str(yr)))
            mnth = int(year_month.split('-')[1])
            print('Working on files in '+str(yr)+'/'+str(mnth))
            ym_indices = np.where(np.logical_and(years==yr,months==mnth))[0]

            subfolder = os.path.join(data_dir, 'EN.' + version + '.profiles.' + correction + '.' + str(yr))
            month_file = 'EN.'+version + '.f.profiles.' + correction + '.' + str(yr) + '{:02d}'.format(mnth) + '.nc'
            monthData = nc4.Dataset(os.path.join(subfolder, month_file))


            for ym in range(np.size(ym_indices)):
                ym_index = ym_indices[ym]

                index = int(indices[ym_index])

                parms = monthData.variables['PROJECT_NAME'][:]
                tem = monthData.variables['POTM_CORRECTED'][:]
                sal = monthData.variables['PSAL_CORRECTED'][:]
                dep = monthData.variables['DEPH_CORRECTED'][:]

                profile = np.hstack([np.reshape(tem[index, :], (np.shape(dep)[1], 1)),
                                     np.reshape(sal[index, :], (np.shape(dep)[1], 1)),
                                     np.reshape(dep[index, :], (np.shape(dep)[1], 1))])
                goodIndices = np.logical_and(np.logical_and(profile[:, 0] < 99999, profile[:, 1] < 99999),
                                             profile[:, 2] < 99999)
                profile = profile[goodIndices, :]

                year = years[ym_index]
                month = months[ym_index]
                day = days[ym_index]
                hour = hours[ym_index]
                minute = minutes[ym_index]
                longitude = longitudes[ym_index]
                latitude = latitudes[ym_index]

                project_name = ''
                for letter in parms[index]:
                    project_name += letter.decode('utf-8')
                project_name = '_'.join(project_name.strip().split())

                file_name = 'CTD_'+str(year)+'{:02d}'.format(month)+'{:02d}'.format(day)+'_'+'{:02d}'.format(hour)+'{:02d}'.format(minute)+'_'+project_name

                if output_type=='nc':
                    output_file = os.path.join(output_dir, str(yr),file_name+'.nc')
                    source = month_file

                    output_file_as_nc(output_file, longitude, latitude,
                                      year, month, day, hour, minute,
                                      profile, source, project_name)
                else:
                    raise ValueError('output_type must be nc - no other types have yet been implemented')

            monthData.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-o", "--output_directory", action="store",
                        help="The directory where the data will be stored.", dest="output_dir",
                        type=str, required=True)

    parser.add_argument("-d", "--data_directory", action="store",
                        help="The directory where the data is located.", dest="data_dir",
                        type=str, required=True)

    parser.add_argument("-f", "--file_path", action="store",
                        help="The path to the file with the list of ctds to output.", dest="file_path",
                        type=str, required=True)

    parser.add_argument("-t", "--output_type", action="store",
                        help="The type of data to output (options: csv, nc).", dest="output_type",
                        type=str, required=True)

    parser.add_argument("-v", "--version", action="store",
                        help="The version of the dataset (default is 4.2.2).", dest="version",
                        type=str, required=False, default='4.2.2')

    parser.add_argument("-c", "--correction", action="store",
                        help="Corrections applied to the data (default is l09).", dest="correction",
                        type=str, required=False, default='l09')

    parser.add_argument("-y", "--years", action="store",
                        help="The years to download (default is 1960 - 2021).", dest="years",
                        type=int, nargs='+', required=False, default=-1)



    args = parser.parse_args()
    output_dir = args.output_dir
    data_dir = args.data_dir
    file_path = args.file_path
    output_type = args.output_type
    version = args.version
    correction = args.correction
    years = args.years

    if years == -1 or -1 in years:
        years = np.arange(1960, 2020).tolist()

    format_individual_ctds(output_dir,data_dir,file_path,output_type,version,correction,years)