# Hadley Centre
This oceanograhic data is collected and compiled by the UK Met Office Hadley Center. 


## Using these scripts
To download and access data from Hadley Centre using the scripts provided here, first identify the following keywords to be used in the command line:

| Term | Short Key | Description|
| ---- | --------- | -----------|
| data_directory | -d | The directory where the Hadley data will be stored, such as an external drive. |
| output_directory | -o | The directory where the metadata and local CTD files will be stored. |
| version | -v | (Optional) The version of the Hadley set to acquire (default is 4.2.2). |
| correction | -c | (Optional) The set of corrections that are applied to the data (default is l09). |
| years | -y | (Optional) The years of interest (default is 1960 to 2022). |
| shapefile_path | -s | A path to a shapefile which delineates the region of interest. |
| pole | -p | The pole for the region of interest (north or south). |

Start by downloading data from the Hadley Centre, over your timespan of interest
```
python3 download_hadley_data.py -d data_directory -v [version] -c [correction] -y [years]
```
Note that the options in brackets are optional.

Next, search through the downloaded files to identify measurements in your region of interest. The search area is defined by a shapefile in longitude-latitude coordinates (EPSG: 4326).
```
python3 search_hadley_data_using_shapefile.py -o output_directory -d data_directory -s shapefile_path -p pole -v [version] -c [correction] -y [years]
```
This script will create a metdata file which stores the times and locations of all measurements available in your region of interest.

Finally, create individual CTDs files using the data stored in the metadata file:
```
python3 format_individual_ctds_from_file_list.py -o output_directory -d data_directory -f file_path -t output_type -v [version] -c [correction] -y [years]
```

## Citation
To cite this data (described [HERE](https://www.metoffice.gov.uk/hadobs/en4/terms_and_conditions.html)), use the following:

Good, S. A., M. J. Martin and N. A. Rayner, 2013. EN4: quality controlled ocean temperature and salinity profiles and monthly objective analyses with uncertainty estimates, Journal of Geophysical Research: Oceans, doi:10.1002/2013JC009067.

and mention

"EN.4.2.2 data were obtained from https://www.metoffice.gov.uk/hadobs/en4/ and are British Crown Copyright, Met Office, \[year of first publication\], provided under a Non-Commercial Government Licence http://www.nationalarchives.gov.uk/doc/non-commercial-government-licence/version/2/."
