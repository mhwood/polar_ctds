# Oceans Melting Greenland
These CTDs were collected during NASA's Oceans Melting Greenland (OMG) Mission 2015-2021.

## Data Storage
The OMG data are stored are JPL's Physical Oceanograph Distributed Active Archive Center. To use these codes, you will need your credentials for this system, which are different than those used to sign into Earthdata. First, navigate to the [PO.DAAC page](https://podaac.jpl.nasa.gov/dataaccess) and log in using your Earthdata credentials. Then, follow the link to "Access API Credentials" to view your PO.DAAC username and password.

## Using these scripts
To download and access data from the OMG Mission using the scripts provided here, first identify the following keywords to be used in the command line:

| Term | Short Key | Description|
| ---- | --------- | -----------|
| data_directory | -d | The directory where the OMG data will be stored, such as an external drive. |
| output_directory | -o | The directory where the metadata and local CTD files will be stored. |
| uname | -u | Username for PO.DAAC (see above) |
| pword | -p | Password for PO.DAAC (see above) |
| sources | -s | (Optional) CTDs, AXCTDs, or both (default is both). |
| correction | -c | (Optional) The set of corrections that are applied to the data (default is l09). |
| years | -y | (Optional) The years of interest (default is 2015 to 2021). |
| shapefile_path | -s | A path to a shapefile which delineates the region of interest. |
| pole | -p | The pole for the region of interest (north or south). |

Start by downloading the OMG data from PO.DAAC, over your timespan of interest
```
python3 download_omg_data.py -d data_directory -u uname -p pword -s [sources] -y [years]
```
Note that the options in brackets are optional.

Next, search through the downloaded files to identify measurements in your region of interest. The search area is defined by a shapefile in longitude-latitude coordinates (EPSG: 4326).
```
python3 search_omg_data.py -o output_directory -d data_directory -s shapefile_path -y [years]
```
This script will create a metdata file which stores the times and locations of all measurements available in your region of interest.

Finally, create individual CTDs files using the data stored in the metadata file:
```
python3 format_individual_ctds_from_file_list.py -o output_directory -d data_directory -f file_path -t output_type -v [version] -c [correction] -y [years]
```

## Citation
To cite this data (described [HERE](https://omg.jpl.nasa.gov/portal/browse/OMGEV-AXCTD/), use the following:

AXCTDs:

OMG. 2020. OMG Airborne eXpendable Conductivity Temperature Depth (AXCTD) Profiles. Ver. 1. PO.DAAC, CA, USA. Dataset accessed [YYYY-MM-DD] at https://dx.doi.org/10.5067/OMGEV-AXCT1.

CTDs:

 acquired from PO.DAAC	OMG. 2020. OMG Conductivity, Temperature and Depth (CTD) profiles. Ver. 1. PO.DAAC, CA, USA. Dataset accessed [YYYY-MM-DD] at https://dx.doi.org/10.5067/OMGEV-CTDS1.
