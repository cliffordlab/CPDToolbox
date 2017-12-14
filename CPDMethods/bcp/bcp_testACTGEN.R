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

# Test on ACTGEN data
# Get names of RR and TCP files
act_files <-list.files(str_c(pathToRepo,"/artificialDataAnalysis/generatedData/act"))
actTcp_files <-list.files(str_c(pathToRepo,"/artificialDataAnalysis/generatedData/act_tcp"))

out <- matrix(,nrow=length(act_files),ncol=4)

for (f in 1:length(act_files)) 
{
  # Load RR interval data and TCPs
  act = read.table(str_c(pathToRepo,"/artificialDataAnalysis/generatedData/act/",act_files[f]))
  tcp = read.table(str_c(pathToRepo, "/artificialDataAnalysis/generatedData/act_tcp/",actTcp_files[f]))
  
  # Convert lists to vectors
  act <- as.vector(act[[1]])
  tcp <- as.vector(tcp[[1]])
  
  # Do changepoint detection
  bcp.out <- bcp(act, w0 = 0.3, p0 = 0.3)
  plot(bcp.out, main="ACTGEN")
  plot(bcp.out)
  
  # Find where posterior probablity = 1, these are change points
  posteriorProb <- bcp.out$posterior.prob
  ecp <- which(posteriorProb >= 0.5)
  
  # Convert ECP and TCP to time 
  actTime <- cumsum(act)
  ecpTime <- actTime[ecp]
  tcpTime <- actTime[tcp]
  
  # Find accuracy of CP detection
  tolerance <- 5
  result <- accuracy(ecpTime, tcpTime, tolerance)
  
  out[f,1] <- result$TPR
  out[f,2] <- result$FNR
  out[f,3] <- result$PPV
  out[f,4] <- result$falsePos
}

write.table(out, file=str_c(pathToRepo,"/artificialDataAnalysis/bcpResultsADA/actgenTestResults"), row.names=FALSE, col.names=FALSE)
