
import os
import netCDF4 as nc4
import numpy as np
import shapefile
import matplotlib.pyplot as plt
import matplotlib.path as mplPath
from pyproj import Proj, Transformer
import argparse

def JDtoYMD(JD,refYear=1950):

    printDebug=False

    yearFound=False
    for yr in range(refYear,2050):
        if not yearFound:
            dayAddition=365
            if yr%4==0:
                dayAddition+=1
            if JD-dayAddition<=0:
                year=yr
                yearFound = True
            else:
                JD-=dayAddition

    if printDebug:
        print('Year = '+str(year)+', days remaining: '+str(JD))

    if year % 4 != 0:
        daysBeforeMonthsDictionary = {1: 0, 2: 31, 3: 59, 4: 90, 5: 120, 6: 151, 7: 181, 8: 212, 9: 243, 10: 273,
                                      11: 304, 12: 334}
    else:
        daysBeforeMonthsDictionary = {1: 0, 2: 31, 3: 60, 4: 91, 5: 121, 6: 152, 7: 182, 8: 213, 9: 244, 10: 274,
                                      11: 305,
                                      12: 335}
    monthFound=False
    for mo in range(1,13):
        if not monthFound:
            monthDays=daysBeforeMonthsDictionary[mo]
            if JD - monthDays <= 0:
                month = mo-1
                monthFound = True
                JD-=daysBeforeMonthsDictionary[month]

    if not monthFound:
        month=mo
        JD -= daysBeforeMonthsDictionary[month]

    if printDebug:
        print('Month = ' + str(month) + ', days remaining: ' + str(JD))

    day=int(np.floor(JD))
    hrs=24*(JD-day)
    hour=int(np.floor(hrs))
    minute=int(np.floor(60*(hrs-hour)))
    day+=1

    return(year,month,day,hour,minute)

def YMDtoDecYr(year,month,day):
    DOY=0
    if year%4!=0:
        daysBeforeMonthsDictionary={1:0,2:31,3:59,4:90,5:120,6:151,7:181,8:212,9:243,10:273,11:304,12:334}
    else:
        daysBeforeMonthsDictionary = {1: 0, 2: 31, 3: 60, 4: 91, 5: 121, 6: 152, 7: 182, 8: 213, 9: 244, 10: 274, 11: 305,
                                  12: 335}
    DOY+=daysBeforeMonthsDictionary[month]
    DOY+=day
    if year%4!=0:
        decYr = year+DOY/365.0
    else:
        decYr = year+DOY/366.0
    return(decYr)

def reproject_polygon(polygon_array,inputCRS,outputCRS,x_column=0,y_column=1):
    transformer = Transformer.from_crs('EPSG:' + str(inputCRS), 'EPSG:' + str(outputCRS))
    # inProj = Proj(init='epsg:'+str(inputCRS))
    # outProj = Proj(init='epsg:'+str(outputCRS))
    x2, y2 = transformer.transform(polygon_array[:,y_column], polygon_array[:,x_column])
    output_polygon=np.copy(polygon_array)
    output_polygon[:,x_column] = x2
    output_polygon[:,y_column] = y2
    return output_polygon

def search_hadley_data(output_dir,data_dir,shapefile_path,pole,version,correction,years,min_depth,location_name):
    min_year = np.min(np.array(years))
    max_year = np.max(np.array(years))
    total_output = 'File_ID,File_Index,Year,Month,Day,Hour,Minute,Longitude,Latitude'
    total_records = []
    total_points = []

    if pole=='north':
        epsg = 3413
    elif pole=='south':
        epsg = 3976
    else:
        raise ValueError('Pole must be north or south')


    sf = shapefile.Reader(shapefile_path)
    shape = sf.shapes()[0]
    outline = np.array(shape.points)
    outline_reprojected = reproject_polygon(outline, 4326, epsg)
    path = mplPath.Path(outline_reprojected)

    for year in range(min_year, max_year + 1):
        print('    Searching for profiles in ' + str(year))
        subfolder = data_dir + '/EN.' + version + '.profiles.' + correction + '.' + str(year)
        # year_output = 'File_ID,Year,Month,Day,Hour,Minute,Longitude,Latitude'
        # year_records = []
        # year_points = []
        for month in range(1, 13):
            # print('    Searching for profiles in ' + str(year) + '/' + str(month))

            month_file = 'EN.'+version + '.f.profiles.' + correction + '.' + str(year) + '{:02d}'.format(month) + '.nc'

            if month_file in os.listdir(subfolder):
                monthData = nc4.Dataset(os.path.join(subfolder,month_file))
                lon = monthData.variables['LONGITUDE'][:]
                lat = monthData.variables['LATITUDE'][:]
                indices = np.arange(len(lon))
                locations = np.hstack([np.reshape(lon, (len(lon), 1)),
                                       np.reshape(lat, (len(lon), 1))])

                locations = reproject_polygon(locations, 4326, epsg)
                inside = path.contains_points(locations)
                indices = indices[inside]
                locations = locations[inside, :]

                # indices = indices[locations[:, 1] < np.max(outline[:, 1])]
                # locations = locations[locations[:, 1] < np.max(outline[:, 1]), :]

                if np.shape(locations)[0]>0:

                    # plt.plot(locations[:,0],locations[:,1],'b.')
                    # plt.plot(outline_reprojected[:, 0], outline_reprojected[:, 1], 'k-')
                    # plt.show()
                    # locations = locations[locations[:,1]<-60,:]

                    if len(indices) > 0:
                        parms = monthData.variables['PROJECT_NAME'][:]
                        jds = monthData.variables['JULD'][:]
                        tem = monthData.variables['POTM_CORRECTED'][:]
                        sal = monthData.variables['PSAL_CORRECTED'][:]
                        dep = monthData.variables['DEPH_CORRECTED'][:]

                        for ii in range(len(indices)):
                            longitude = lon[indices[ii]]
                            latitude = lat[indices[ii]]

                            parm = parms[indices[ii]]
                            project = ''
                            for ll in range(len(parm)):
                                project += parm[ll].decode('utf-8')
                            project = project.strip()
                            project = '_'.join(project.split())

                            Y, Mo, D, H, Mi = JDtoYMD(jds[indices[ii]])
                            # dec_yr = YMDtoDecYr(Y, Mo, D)
                            output_id = str(Y) + '{:02d}'.format(Mo) + '{:02d}'.format(D) + '_' + '{:02d}'.format(
                                H) + '{:02d}'.format(Mi) + '_' + project

                            profile = np.hstack([np.reshape(tem[indices[ii], :], (np.shape(dep)[1], 1)),
                                                 np.reshape(sal[indices[ii], :], (np.shape(dep)[1], 1)),
                                                 np.reshape(dep[indices[ii], :], (np.shape(dep)[1], 1))])
                            goodIndices = np.logical_and(np.logical_and(profile[:, 0] < 99999, profile[:, 1] < 99999),
                                                         profile[:, 2] < 99999)
                            profile = profile[goodIndices, :]

                            if len(profile) > 0:
                                if np.max(profile[:,2])>min_depth:
                                #     ctd_output = 'Longitude,Latitude,Decimal_Year,Pot.Temp,Salinity,Depth'
                                #     for ll in range(len(profile)):
                                #         ctd_output += '\n' + str(longitude) + ',' + str(latitude)
                                #         ctd_output += ',' + str(dec_yr)
                                #         ctd_output += ',' + str(profile[ll, 0]) + ',' + str(profile[ll, 1]) + ',' + str(
                                #             profile[ll, 2])
                                #     if str(Y) not in os.listdir(output_folder + '/Data'):
                                #         os.mkdir(output_folder + '/Data/' + str(Y))
                                #     output_file = output_folder + '/Data/' + str(Y) + '/' + output_id + '.csv'
                                #     f = open(output_file, 'w')
                                #     f.write(ctd_output)
                                #     f.close()

                                    output_line = '\n' + output_id +','+str(indices[ii])+ ',' + str(Y) + ',' + str(Mo) + ',' + str(D) + ',' + str(
                                        H) + ',' + str(Mi) + ',' + str(longitude) + ',' + str(latitude)
                                    # year_output += output_line
                                    total_output += output_line
                                    # year_records.append([output_id, Y, Mo, D, H, Mi])
                                    # year_points.append([longitude, latitude])
                                    total_records.append([output_id, Y, Mo, D, H, Mi])
                                    total_points.append([longitude, latitude])

        # output_file = output_folder + '/Metadata/File Lists/Hadley_CTD_Locations_' + str(year) + '.csv'
        # f = open(output_file, 'w')
        # f.write(year_output)
        # f.close()

        # output_file = output_folder + '/Metadata/Shapefiles/Hadley_CTD_Locations_' + str(year)
        # fields = ['ID', 'Year', 'Month', 'Day', 'Hour', 'Minute']
        # fieldTypes = ['C', 'N', 'N', 'N', 'N', 'N']
        # espg = 4326
        # createPointShapefile(output_file, fields, fieldTypes, year_records, year_points, espg)

    if location_name!='':
        output_file = output_dir + '/Hadley_CTD_Locations_'+'_'.join(location_name.split())+'.csv'
    else:
        output_file = output_dir + '/Hadley_CTD_Locations.csv'
    f = open(output_file, 'w')
    f.write(total_output)
    f.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-o", "--output_directory", action="store",
                        help="The directory where the file list will be stored.", dest="output_dir",
                        type=str, required=True)

    parser.add_argument("-d", "--data_directory", action="store",
                        help="The directory where the data is located.", dest="data_dir",
                        type=str, required=True)

    parser.add_argument("-s", "--shapefile_path", action="store",
                        help="The shapefile used to outline the search area.", dest="shapefile_path",
                        type=str, required=True)

    parser.add_argument("-p", "--pole", action="store",
                        help="The porjection of the shapefile.", dest="pole",
                        type=str, required=True)

    parser.add_argument("-v", "--version", action="store",
                        help="The version of the dataset (default is 4.2.2).", dest="version",
                        type=str, required=False, default='4.2.2')

    parser.add_argument("-c", "--correction", action="store",
                        help="Corrections applied to the data (default is l09).", dest="correction",
                        type=str, required=False, default='l09')

    parser.add_argument("-y", "--years", action="store",
                        help="The years to search through (default is 1960 - 2021).", dest="years",
                        type=int, nargs='+', required=False, default=-1)

    parser.add_argument("-m", "--min_depth", action="store",
                        help="The minimum depth the CTD must reach to be counted (default = 50 m).", dest="min_depth",
                        type=int, required=False, default=50)

    parser.add_argument("-l", "--location_name", action="store",
                        help="(Optional) Location name to add to the metadata file. Default is none.",
                        default='',
                        dest="location_name", type=str, required=False)


    args = parser.parse_args()
    output_dir = args.output_dir
    data_dir = args.data_dir

    shapefile_path = args.shapefile_path
    pole = args.pole
    version = args.version
    correction = args.correction
    location_name = args.location_name

    years = args.years
    min_depth = args.min_depth

    if years == -1 or -1 in years:
        years = np.arange(1960,2022).tolist()

    search_hadley_data(output_dir,data_dir,shapefile_path,pole,version,correction,years,min_depth,location_name)