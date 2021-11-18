
import os
import requests
import argparse


def download_ctds_from_podaac(project_dir,auth,years,sources,printing):

    podaac_url_base = 'https://podaac-tools.jpl.nasa.gov/drive/files/allData/omg/L2/CTD'

    if 'PODAAC' not in os.listdir(project_dir):
        os.mkdir(os.path.join(project_dir,'PODAAC'))
    if 'Data' not in os.listdir(os.path.join(project_dir,'PODAAC')):
        os.mkdir(os.path.join(project_dir,'PODAAC','Data'))
    if 'CTDs' in sources:
        if 'CTDs' not in os.listdir(os.path.join(project_dir,'PODAAC','Data')):
            os.mkdir(os.path.join(project_dir,'PODAAC','Data','CTDs'))
    if 'AXCTDs' in sources:
        if 'AXCTDs' not in os.listdir(os.path.join(project_dir, 'PODAAC','Data')):
            os.mkdir(os.path.join(project_dir, 'PODAAC','Data','AXCTDs'))

    for year in years:
        if 'AXCTDs' in sources:
            if str(year) not in os.listdir(os.path.join(project_dir, 'PODAAC','Data','AXCTDs')):
                os.mkdir(os.path.join(project_dir, 'PODAAC','Data','AXCTDs',str(year)))
        if 'CTDs' in sources:
            if str(year) not in os.listdir(os.path.join(project_dir, 'PODAAC','Data','CTDs')):
                os.mkdir(os.path.join(project_dir, 'PODAAC','Data','CTDs',str(year)))

    downloaded_files = 0

    ########################################################################################################################
    # This section is for the CTDs

    if 'CTDs' in sources:

        url = 'https://podaac-opendap.jpl.nasa.gov/opendap/allData/omg/L2/CTD/CTD/'

        ctd_urls = []
        for year in years:
            for month in range(1, 13):
                resp = requests.get(url + '/' + str(year) + '/' + '{:02d}'.format(month))
                resp_text = resp.text
                lines = resp_text.split('\n')
                for line in lines:
                    line = line.strip()
                    line = line.split(': ')
                    if line[0] == '"sameAs"':
                        ctd_url = line[1][1:-6]
                        if ctd_url.split('/')[-1][:3] == 'OMG' and int(ctd_url.split('/')[-1].split('_')[-1][:4]) in years:
                            podaac_url = podaac_url_base + '/CTD/' + str(year) + '/' + '{:02d}'.format(month) + '/' + \
                                         ctd_url.split('/')[-1]
                            ctd_urls.append(podaac_url)

        counter = 1
        for ctd_url in ctd_urls:
            counter += 1
            year = ctd_url.split('/')[-1].split('_')[4][:4]
            if ctd_url.split('/')[-1] not in os.listdir(project_dir + '/PODAAC/Data/CTDs/' + str(year)):
                downloaded_files +=1
                if printing:
                    print('Downloading ' + ctd_url.split('/')[-1] + ' (' + str(counter) + ' of ' + str(len(ctd_urls)) + ')')
                outputFile = os.path.join(project_dir,'PODAAC','Data','CTDs', str(year),ctd_url.split('/')[-1])
                with requests.get(ctd_url, auth=auth, stream=True) as r:
                    r.raise_for_status()
                    open(outputFile, 'wb').write(r.content)

    ########################################################################################################################
    # This section is for the AXCTDs

    if 'AXCTDs' in sources:

        url = 'https://podaac-opendap.jpl.nasa.gov/opendap/allData/omg/L2/CTD/AXCTD/'

        resp = requests.get(url)
        resp_text = resp.text
        lines = resp_text.split('\n')

        ctd_urls = []
        for line in lines:
            line = line.strip()
            line = line.split(': ')
            if line[0] == '"sameAs"':
                ctd_url = line[1][1:-6]
                if ctd_url.split('/')[-1][:3] == 'OMG' and int(ctd_url.split('/')[-1].split('_')[-1][:4]) in years:
                    podaac_url = podaac_url_base + '/AXCTD/' + ctd_url.split('/')[-1]
                    ctd_urls.append(podaac_url)

        counter = 1
        for ctd_url in ctd_urls:
            counter += 1
            year = ctd_url.split('/')[-1].split('_')[4][:4]
            if ctd_url.split('/')[-1] not in os.listdir(project_dir + '/PODAAC/Data/AXCTDs/' + str(year)):
                downloaded_files += 1
                if printing:
                    print('Downloading ' + ctd_url.split('/')[-1] + ' (' + str(counter) + ' of ' + str(len(ctd_urls)) + ')')
                outputFile = os.path.join(project_dir,'PODAAC','Data','AXCTDs',str(year),ctd_url.split('/')[-1])
                with requests.get(ctd_url, auth=auth, stream=True) as r:
                    r.raise_for_status()
                    open(outputFile, 'wb').write(r.content)

    if printing:
        if downloaded_files>0:
            print('    Downloaded '+str(downloaded_files)+' total files')
        else:
            print('    Didn\'t download any files - all files already obtained')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--data_directory", action="store",
                        help="Directory where CTD data will be stored",
                        dest="project_dir", type=str, required=True)

    parser.add_argument("-u", "--uname", action="store",
                        help="Username for PO.DAAC account",
                        dest="podaac_uname", type=str, required=True)

    parser.add_argument("-p", "--pword", action="store",
                        help="Password for PO.DAAC account",
                        dest="podaac_pword", type=str, required=True)

    parser.add_argument("-y", "--years", action="store",
                        help="(Optional) List of years to process. A value of -1 will choose all years (2015-2021 as available). "
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

    project_dir = args.project_dir
    podaac_uname = args.podaac_uname
    podaac_pword = args.podaac_pword
    years = args.years
    if years == -1 or -1 in years:
        years = [2015, 2016, 2017, 2018, 2019, 2020, 2021]
    sources = args.sources
    if sources=='all':
        sources = ['CTDs','AXCTDs']
    printing_int = args.printing_int
    if printing_int == 1:
        printing = True
    else:
        printing = False

    print('Downloading CTDs from the OMG Campaign')
    print('Processing parameters: ')
    print('    Data directory: '+project_dir)
    print('    PO.DAAC Username: '+podaac_uname)
    print('    PO.DAAC Password: ' + podaac_pword)
    print('    Years:            ' + str(years))
    print('    Sources:          ' + str(sources))
    print('    Print Messages:   ' + str(printing))
    print(' ')

    if printing:
        print(' Step 1: Downloading data from PO.DAAC')
    auth = (podaac_uname,podaac_pword)
    download_ctds_from_podaac(project_dir,auth,years,sources,printing)
