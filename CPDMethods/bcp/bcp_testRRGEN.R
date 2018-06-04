# 
# Overview
#   Generating assessment results for BCP in order to compare with other CP detection methods 
#   in MATLAB.
#   Set working directory to where you keep toolbox
#   ex.: setwd('~/repos/CPDToolbox')
#

# Get working directory ~ this should be where CPD is located ex.: setwd('~/repos/CPDToolbox');
pathToRepo = getwd()

library(bcp)
library(stringr)
source(str_c(pathToRepo, '/CPDMethods/bcp/assess_CP.R'))

# Test on RR data
# Get names of RR and TCP files
rr_files <-list.files(str_c(pathToRepo,"/artificialDataAnalysis/generatedData/rr"))
rrTcp_files <-list.files(str_c(pathToRepo,"/artificialDataAnalysis/generatedData/rr_tcp"))

out <- matrix(,nrow=length(rr_files),ncol=4)

for (f in 1:length(rr_files)) 
{
  # Load RR interval data and TCPs
  rr = read.table( str_c(pathToRepo, "/artificialDataAnalysis/generatedData/rr/", rr_files[f]))
  tcp = read.table(str_c(pathToRepo, "/artificialDataAnalysis/generatedData/rr_tcp/", rrTcp_files[f]))
  
  # Convert lists to vectors
  rr <- as.vector(rr[[1]])
  tcp <- as.vector(tcp[[1]])
  
  # Do changepoint detection
  bcp.out <- bcp(rr, w0 = 0.2, p0 = 0.3)
  plot(bcp.out, main="RRGEN")
  plot(bcp.out)
  
  # Find where posterior probablity = 1, these are change points
  posteriorProb <- bcp.out$posterior.prob
  ecp <- which(posteriorProb >= 0.6)
  
  # Convert ECP and TCP to time 
  rrTime <- cumsum(rr)
  ecpTime <- rrTime[ecp]
  tcpTime <- rrTime[tcp]
  
  # Find accuracy of CP detection
  tolerance <- 5
  result <- accuracy(ecpTime, tcpTime, tolerance)
  
  out[f,1] <- result$TPR
  out[f,2] <- result$f1Score
  out[f,3] <- result$PPV
  out[f,4] <- result$falsePos
}

write.table(out, file=str_c(pathToRepo,"/artificialDataAnalysis/bcpResultsADA/rrgenTestResults"), row.names=FALSE, col.names=FALSE)
