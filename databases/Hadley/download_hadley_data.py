
import os
import numpy as np
import requests
import argparse
import zipfile

def download_hadley_data(output_dir,version,correction,years):

    print('Downloading data from the Hadley EN Dataset')
    print('    Version: '+version)
    print('    Correction: '+correction)

    for y in range(len(years)):
        year = years[y]
        file_name = 'EN.'+version+'.profiles.'+correction+'.'+str(year)+'.zip'
        print(' Downloading ' + file_name +' ('+str(y+1)+' of '+str(len(years))+')')

        url = 'https://www.metoffice.gov.uk/hadobs/en4/data/en4-2-1/EN.'+version+'/'+file_name
        output_file = os.path.join(output_dir,file_name)

        session = requests.Session()
        adapter = requests.adapters.HTTPAdapter(max_retries=20)
        session.mount('https://', adapter)
        session.mount('http://', adapter)
        resp = session.get(url)
        open(output_file, 'wb').write(resp.content)

        print('   Unzipping... ')
        with zipfile.ZipFile(output_file, "r") as zip_ref:
            zip_ref.extractall(os.path.join(output_dir,file_name[:-4]))

        os.remove(output_file)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-o", "--output_directory", action="store",
                        help="The directory where the data will be stored.", dest="output_dir",
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
    version = args.version
    correction = args.correction
    years = args.years

    if years == -1 or -1 in years:
        years = np.arange(1960,2020).tolist()

    download_hadley_data(output_dir,version,correction,years)