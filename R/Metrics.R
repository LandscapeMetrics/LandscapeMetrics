#' Area in hectare of patches
#'
#'Get the area in hectares of all raster classes.
#'@param x is a matrix of data with patches identified as classes.
#'@return A table with the values in hectares of each class and their proportion in percentage in the area.
#'@import raster
#'@import rgdal
#'@import rgeos
#'@export
CA <- function(x){
  class_forest=x==1
  fragment_forest=clump(class_forest)
  cell_fragment=freq(fragment_forest)
  cell_fragment_table=data.frame(cell_fragment)
  cell_fragment_table=cell_fragment_table[!is.na(cell_fragment_table[,1]),]
  area_pixel=res(raster)[1] * res(raster)[2]
  sum_area=sum(cell_fragment_table$count)
  total_area_class_ha=sum_area * area_pixel/10000
  return(total_area_class_ha)
}


