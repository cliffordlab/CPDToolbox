# 
# Overview
#   Generating assessment results for BCP in order to compare with other CP detection methods 
#   in MATLAB
#

# Get working directory ~ this should be where CPD is located ex.: setwd('~/repos/CPDToolbox');
pathToRepo = getwd()

library(bcp)
library(stringr)
source(str_c(pathToRepo, '/CPDMethods/bcp/assess_CP.R'))

# Test on RR data
# Get names of RR and TCP files
rr_files <-list.files(str_c(pathToRepo, '/realDataAnalysis/RRIntervals'))

for (f in 1:length(rr_files)) 
{
  # Load RR interval data
  rr = read.table(str_c(pathToRepo, '/realDataAnalysis/RRIntervals/',rr_files[f]))

  # Convert list to vector
  rr <- as.vector(rr[[1]])

  # Do changepoint detection
  bcp.out <- bcp(rr, w0 = 0.2, p0 = 0.3)
  plot(bcp.out, main="NN")
  plot(bcp.out)
  
  # Find where posterior probablity = 1, these are change point
  posteriorProb <- bcp.out$posterior.prob
  ecp <- which(posteriorProb >= 0.6)
  
  write.table(ecp, file=str_c(pathToRepo,'/realDataAnalysis/bcpResultsRDA/',f), row.names=FALSE, col.names=FALSE)
}
