
# AviSpace

Avian niche specificity and overlap in the Sierra Nevada

AviSpace is a Shiny application which uses principal components analysis
and remotely sensed geospatial habitat variables in combination with
passive acoustic monitoring and machine learning species classifiers to
quantify the prefered habitat niche of 90 bird species in the Sierra
Nevada, USA.

![AviSpace cover image](images/AviSpaceCoverImage.png)

# Functionalities

1.  View prefered habitat niche polygons and overlap of any two selected
    species
2.  View range map and range overlap of any two selected species
3.  View boxplots of any single habitat variable vs any number of
    species

# Accessibility

The AviSpace app can be accessed online at
(<https://kalmialatifolia.shinyapps.io/Avispace4/>)

Alternatively, you can download the shiny app to run locally: [View the
Shiny app folder](AviSpace4.0/)

# File descriptions and variable definitions

**AviSpace4.0** - Self contained Shiny app for download

**fullProtocol.R** - R script containing the complete data processing
pipeline.

**speciesDetails.csv** - CSV file containing total number of detections,
total number of sites where the species was detected, and niche area for
each of the 90 species.

| Variable | Description |
|----------|-------------| 
| species | Species common name | 
| nicheArea | Species habitat niche area in unitless PC space. Habitat specialists have smaller habitat niche areas, generalists have larger ones | 
| totalDetections | The total number of times this species was detected with greater than 90% certainty across 577 monitoring sites for the duration of the study | 
| totalSites | The total number of sites where this species was detected at least once with greater than 90% certainty across the duration of the study |


**GeospatialMetadata.csv** - CSV file containing details on each of the
geospatial variables used as eigenvectors in the PCA, including variable
descriptions and links to the original publications and datasets.

| Variable | Description |
|----------|-------------| 
| Variable Name | Shorthand variable name used in R scripts and column headers | 
| Category | One of six possible BioCube categories: Climate, Disturbance, Function, Phenology, Physiology, Terrain | 
| Description | Longhand variable definition |
| Source | Original academic article where this dataset was published, First author and publication year | | DOI | DOI link to access the original publication | 
| Dataset | Link to access the repository where the dataset is stored | 
| Used in PCA | Was this variable used in the final principal component analysis? Yes or no | 
| Notes | Notes. If this variable was not used in the final analysis, why? |


**siteDetections_foliarTraits_Biocube.rds** - Full dataset with site locations, species detection rates, and spatial variables. One row per site.

| Variable | Description |
|----------|-------------| 
| cell_unit | ID of geographic cell in which the ARU was placed
| long | Longitude | 
| lat | Latitude | 
| burnYear | Year of the most recent burn. NA if no fires on record | 
| Species Columns [5:96] | Column names are bird species names. Column values are detection rates at that site |
| Variable Columns [97:182] | Column names are spatial variables as described in *GeospatialMetadata.csv*. Column values are the weighted average of that variable within a 150 m radius buffer around the recording location |
| geometry | Point geometry of lat/lon coordinates |
| landcover_type | Land cover type according to USGS CONUS GAP 2011 (https://doi.org/10.5066/F7ZS2TM0) |


**RangeOverlaps.rds** - The overlapping geographic ranges of pairs of bird species based on breeding range maps from eBird status and trends (https://science.ebird.org/en/status-and-trends)

| Variable | Description |
|----------|-------------| 
| species1 | a bird species (common name) |
| species2 | a bird species (common name) |
| area1 | the range area of species 1 in kilometers squared |
| area2 | the range area of species 2 in kilometers squared |
| overlap_area | the overlapping range area of species1 and species2 in kilometers squared |
| overlap_percent | the percent range overlap of species1 and species2 |

**speciesNicheAreas.rds** - Species habitat niche areas based on PCA of preferred habitat.

| Variable | Description |
|----------|-------------|
| species | bird species (common name) |
| area | niche area |

**NicheRangeOverlap.rds** - Percent overlaps in geographic range and habitat preference of all 4,002 pairwise species combinations.

| Variable | Description |
|----------|-------------| 
| species1 | a bird species (common name) |
| species2 | a bird species (common name) |
| niche_overlap_percent | the percent overlap in habitat niche space of species1 and species2 |
| range_overlap_percent |the percent geographic range overlap of species1 and species2 |
| genus_species1 | the genus name of species 1 |
| genus_species2 | the genus name of species 2 |
| congeneric | do species 1 and 2 belong to the same genus? |
| sisterSpecies | are species 1 and 2 sister species? |


**speciesPolygons.rds** - Multipolygon geometry describing the preferred niche of each species in principal component space

| Variable | Description |
|----------|-------------| 
| species | Species common name | 
| geometry | Multipolygon geometry of the species' niche area in PC space |




# Citation

If you use AviSpace in your published work, please cite the original
manuscript:

Berman, LM; Schneider, FD; Pavlick, RP; Peery, MZ; Wood, CM; Zheng, T;
Shafron, E; Ye, Z; Queally, N; Dean, M; Tagliabue, G; Winiarski, J;
Kramer, A; Townsend, P. (in prep) Remote sensing and bioacoustics reveal
avian niche partitioning and habitat filtering
