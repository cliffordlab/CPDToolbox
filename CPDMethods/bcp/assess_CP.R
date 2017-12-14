# 
# Overview
#   Calculates accuracy of changepoint detection method in terms of TPR,
#   FNR, PPV and F1 score
#  
# Input
#   ecpTime, tcpTime: Time of estimated and true changepoints
#   tolerance: Tolerance in seconds
#

accuracy <- function(ecpTime, tcpTime, tolerance) 
{
  fp = 0;
  tp = 0; 
  fn = 0;
  tcpFound = rep(0, length(ecpTime));
  
  for (iCP in 1:length(ecpTime))
  {
    time_diff <- abs(tcpTime - ecpTime[iCP])
    idx_detected <- min(which(time_diff <= tolerance))
    
    # If idx_detected is empty, ECP is not whithin tolerance of TCPs so we have a FP
    if (is.infinite(idx_detected))
    {
      fp <- fp + 1
      tcpFound[iCP] <- 0
    }
    else
    {
      tcpFound[iCP] <- idx_detected
      if (iCP > 1)
      {
        if (tcpFound[iCP] != tcpFound[iCP-1])
        {
          tp <- tp + 1;
        }
      }
      else
      {
        tp <- tp + 1;
      }
    }
  }
  
  for (iCP in 1:length(tcpTime))
  {
    if (is.infinite(min(which(tcpFound == iCP))))
    {
      fn <- fn + 1;
    }
  }  
  
  # Calculate and assign metrics
  tpr <- tp / (tp + fn)
  ppv <- tp / (tp + fp)
  f1 <-  (2 * tp) / (2 * tp + fp + fn)
  fnr <- fn / (fn + tp)
  
  result <- list(truePos = tp, falseNeg = fn, falsePos =fp, f1Score = f1, TPR = tpr, PPV = ppv, FNR = fnr)
  
  return(result)
}