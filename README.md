Description:

The LandscapeMetrics version package was developed in order to calculate certain landscape ecology metrics, assisting in research on forest fragmentation analysis.They were based on the methodology of McGarigal & Marks (1995).


Installation the LandscapeMetrics:

The package is available in the GitHub repository (<github.com/LandscapeMetrics>) and you can download and install it from the source site, you need the devtools package and use the following commands:

library("devtools")
devtools::install_github("LandscapeMetrics/LandscapeMetrics")

Using functions:
library("LandscapeMetrics")

Import a raster:
data(class)

Area calculation:
area <- CA(class)

References
McGarigal, K., Marks, B. J. (1995). Fragstats: Spatial pattern analysis program for quantifying landscape structure. Portland: USDA Forest Service. Retrieved from https://fs.usda.gov/treesearch/pubs/3064
