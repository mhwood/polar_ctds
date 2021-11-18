
import os
import numpy as np
import matplotlib.path as mplPath
import netCDF4 as nc4
import argparse
# import shapefile
from pyproj import Proj, Transformer

def reproject_polygon(polygon_array,inputCRS,outputCRS,x_column=0,y_column=1):
    transformer = Transformer.from_crs('EPSG:' + str(inputCRS), 'EPSG:' + str(outputCRS))
    # inProj = Proj(init='epsg:'+str(inputCRS))
    # outProj = Proj(init='epsg:'+str(outputCRS))
    x2, y2 = transformer.transform(polygon_array[:,y_column], polygon_array[:,x_column])
    output_polygon=np.copy(polygon_array)
    output_polygon[:,x_column] = x2
    output_polygon[:,y_column] = y2
    return output_polygon


def search_omg_data(output_dir,data_dir,search_with_bbox,bbox,search_with_shapefile,shapefile_path,years,sources,printing):

    min_year = np.min(np.array(years))
    max_year = np.max(np.array(years))
    total_output = 'File_ID,Type,Year,Month,Day,Hour,Minute,Second,Longitude,Latitude'

    epsg = 3413

    if search_with_bbox:
        outline = np.array([[bbox[0],bbox[1]],
                            [bbox[2],bbox[1]],
                            [bbox[2],bbox[3]],
                            [bbox[0],bbox[3]],
                            [bbox[0],bbox[1]]])
        outline_reprojected = reproject_polygon(outline, 4326, epsg)
        path = mplPath.Path(outline_reprojected)

    if search_with_shapefile:
        sf = shapefile.Reader(shapefile_path)
        shape = sf.shapes()[0]
        outline = np.array(shape.points)
        outline_reprojected = reproject_polygon(outline, 4326, epsg)
        path = mplPath.Path(outline_reprojected)

    for source in sources:
        for year in range(min_year, max_year + 1):
            subfolder = os.path.join(data_dir,'Data',source,str(year))

            print('    Searching for profiles in ' + str(year) + ' for source ' + source)

            for file_name in os.listdir(subfolder):
                if file_name[-2:]=='nc':
                    ds = nc4.Dataset(os.path.join(subfolder,file_name))
                    lon = ds.variables['lon'][:][0]
                    lat = ds.variables['lat'][:][0]
                    locations = np.array([[lon,lat]])
                    locations = reproject_polygon(locations, 4326, epsg)
                    if path.contains_point(locations[0,:]):
                        Y = int(file_name.split('_')[-1][:4])
                        Mo = int(file_name.split('_')[-1][4:6])
                        D = int(file_name.split('_')[-1][6:8])
                        H = int(file_name.split('_')[-1][8:10])
                        Mi = int(file_name.split('_')[-1][10:12])
                        S = int(file_name.split('_')[-1][12:14])
                        output_id = str(Y) + '{:02d}'.format(Mo) + '{:02d}'.format(D) + '_' + '{:02d}'.format(
                            H) + '{:02d}'.format(Mi)  + '{:02d}'.format(S)
                        output_line = '\n' + output_id  +','+source+ ',' + str(Y) + ',' + str(
                            Mo) + ',' + str(D) + ',' + str(
                            H) + ',' + str(Mi) + ','+ str(S)+',' + str(lon) + ',' + str(lat)
                        total_output += output_line

        output_file = output_dir + '/OMG_CTD_Locations.csv'
        f = open(output_file, 'w')
        f.write(total_output)
        f.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--data_directory", action="store",
                        help="Directory where CTD data will be stored",
                        dest="data_dir", type=str, required=True)

    parser.add_argument("-o", "--output_directory", action="store",
                        help="Directory where metadata and processed data will be stored",
                        dest="project_dir", type=str, required=True)

    parser.add_argument("-b", "--bounding_box", action="store",
                        help="(Optional) A bounding box used to search the database. "
                             "Default value is no box (in this case, a shapefile must be provided).",
                        default=-1, dest="bounding_box", type=float, nargs='+',
                        required=False)

    parser.add_argument("-f", "--shapefile_path", action="store",
                        help="(Optional) A path to a shapefile used to search the database. "
                             "Default is no path (in this case, a bounding box must be provided).", default='', dest="shapefile_path", type=str,
                        required=False)

    parser.add_argument("-y", "--years", action="store",
                        help="(Optional) List of years to process. A value of -1 will choose all years (2015-2020 as available). "
                             "Default value is -1.", default=-1, dest="years", type=int, nargs='+',
                        required=False)

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
    bounding_box = args.bounding_box
    shapefile_path = args.shapefile_path

    if bounding_box!=-1:
        search_with_bbox = True
        search_with_shapefile = False
        if len(bounding_box)!=4:
            raise ValueError ('Bounding box must be [min_x, min_y, max_x, max_y]')
    elif shapefile_path!='':
        search_with_bbox = False
        search_with_shapefile = True
    else:
        raise ValueError('Must provide a bounding box or a shapefile for searching')

    years = args.years
    if years == -1 or -1 in years:
        years = [2015, 2016, 2017, 2018, 2019, 2020]
    sources = args.sources
    if sources=='all':
        sources = ['CTDs','AXCTDs']
    printing_int = args.printing_int
    if printing_int == 1:
        printing = True
    else:
        printing = False

    print('Searching for CTDs in the domain provided')
    print('Processing parameters: ')
    print('    Data directory: ' + data_dir)
    print('    Output directory: '+project_dir)
    if search_with_bbox:
        print('    Bounding box:     ['+str(bounding_box[0])+', '+str(bounding_box[1])+', '+str(bounding_box[2])+', '+str(bounding_box[3])+']')
    if search_with_shapefile:
        print('    Shapefile path:   '+shapefile_path)
    print('    Years:            ' + str(years))
    print('    Sources:          ' + str(sources))
    print('    Print Messages:   ' + str(printing))
    print(' ')

    search_omg_data(project_dir,data_dir,search_with_bbox,bounding_box,search_with_shapefile,shapefile_path,years,sources,printing)


