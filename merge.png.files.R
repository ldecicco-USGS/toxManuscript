merge.png.pdf <- function(pdfFile, pngFiles, path, deletePngFiles=FALSE) {
  
  library(grid)
  #### Package Install ####
  pngPackageExists <- require ("png")
  
  if ( !pngPackageExists ) {
    install.packages ("png")
    library ("png")
    
  }
  #########################
  
  pdf(pdfFile,width=7,height=7)
  
  n <- length(pngFiles)
  
  for( i in 1:n) {
    
    pngFile <- paste0(path,"/",pngFiles[i])
    
    pngRaster <- readPNG(pngFile)
    
    grid.raster(pngRaster, width=unit(0.8, "npc"), height= unit(0.8, "npc"))
    
    if (i < n) plot.new()
    
  }
  
  dev.off()
  
  if (deletePngFiles) {
    
    unlink(pngFiles)
  }
  
}