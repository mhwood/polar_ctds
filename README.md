# Temperature and Salinity measurements in polar regions
A collection of scripts to find, organize, and process temperature and salinity data from the Arctic and Antarctic regions

## Data sources:
The datasets accessed in this collection are from the following sources:
|Short Name|Description|Temporal Coverage|Spatial Coverage|
|----------|-----------|-----------------|----------------|
|[Hadley Centre](https://www.metoffice.gov.uk/hadobs/en4/)| Collection of oceanographic measurements | 1900-present | Global |
| [Oceans Melting Greenland](https://omg.jpl.nasa.gov/portal/)|NASA's Oceans Melting Greenland Mission | 2015 - 2021 | Greenland Shelf |

## Requirements
It is recommended that a fresh anaconda environement be created when using these scripts. The following packages are required:
```
conda install -c anaconda numpy
conda install -c anaconda matplotlib
conda install -c anaconda requests
conda install -c anaconda xarray
conda install -c conda-forge gsw
conda install -c conda-forge pyproj
```

