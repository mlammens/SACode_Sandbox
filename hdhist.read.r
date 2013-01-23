hdhist.read <- function( hdhistFile ) {
  ########################################################################################
  ## Read data from habdyn history TXT file (hd-hist), produced by the RAMAS GIS - 
  ## Habitat Dynamics module.  The history (HIST) file provides a key
  ## for determining which populaitons in individual *.ptc files match up
  ## with populations in the *.mp file that is produced when several *.ptc
  ## files are combined by the Habitat Dyanimcs (HabDyn) module.
  ##
  ## Args:
  ##  hdhistFile: The path and name of a habdyn history file to be read.
  ##
  ########################################################################################
  
  # Print the name of the hd-hist file to the console
  print( paste( "Begin hdhist.read function with file: ", hdhistFile ) )
  # Read the hd-hist file in as a long unsorted 'list' structure
  Hdhist <- readLines( hdhistFile )
  
  ## Create an empty list structure to store the hd-hist file information 
  hdh.list <- vector("list",length=0)
  
  ## Linked Metapop (*.mp) file - Line 1 of the hd-hist file - the 
  ## name and file path for the *.mp file created by Habitat Dynaimcs, which
  ## combines multiple *.ptc files to incorporate habitat change into meta-
  ## populaiton simulations.
  hdh.list$mp.link <- Hdhist[1]
  
  ## Number of *.ptc files in HabDyn run - Line 2
  hdh.list$PtcFileN <- as.numeric( Hdhist[2] )
  
  ## *.ptc files used in HabDyn run - Line 3 to Line 2 + PtcFileN - each
  ## file path is preceeded by the HabDyn 'iteration' number.
  hdh.list$PtcFiles <- read.table( hdhistFile, header=FALSE, skip=2, nrows=hdh.list$PtcFileN, col.names=c('iter','ptc.file') )
  
  ## Habitat Dynamics Patch Information Matrix - Line 3 + PtcFileN to End-of-file
  ## This information matrix contains the following columns:
  ## iter - the HabDyn iteration number, each iteration corresponds to a 
  ##  particular *.ptc file
  ## Fyear - First year of simulation corresponding to the iteration
  ## Lyear - Last year of simulation corresponding to the iteration
  ## new2old - new2old[ new ] = old - the corrsponding metapopulation number 
  ##  for the present patch-iteration combination
  ## old2new - old2new[ old ] = new 
  hdh.list$patch.info.mat <- read.table( hdhistFile, header=TRUE, skip=(2+hdh.list$PtcFileN))
  
  return(hdh.list)
}