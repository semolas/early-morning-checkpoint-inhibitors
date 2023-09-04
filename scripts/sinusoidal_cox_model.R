library(survival)
library(survminer)
library(readxl)
library(openxlsx)
library(msm) # contains deltamethod

################
Useful functions
################

#### Calculate CI of amplitude A and phase Phi with the Delta Method

calculate_APhi <- function(coxfit){
  # Get coefficients of the Cox model
  b1 <- coxfit$coeff[1]
  b2 <- coxfit$coeff[2]
  
  # Calculate amplitude A and phase Phi
  A <- sqrt(b1^2 + b2^2)
  phi_rad <- if ((b1 > 0 & b2 > 0)|(b1 < 0 & b2 > 0)) atan(b2/b1) else pi + atan(b2/b1) # in radians
  #if (b1 > 0 & b2 >0) print('1st quadrant') else if (b1 < 0 & b2 < 0) print('3rd quadrant') else if (b1 < 0 & b2 > 0) print('2nd quadrant') else if (b1 > 0 & b2 < 0) print('4th quadrant')
  
  # Calculate standard error for A and Phi with deltamethod, by mapping (b1,b2) -> (A=f(b1,b2), phi=g(b1,b2))
  mean <- coefficients(coxfit) # estimated means
  cov <- vcov(coxfit)
  se_APhi <- 
    # if Phi in 1st or 3rd quadrant, Phi=atan is positive
    if ((b1 > 0 & b2 > 0)|(b1 < 0 & b2 > 0))
    {deltamethod(list(~ sqrt(x1^2+x2^2), ~ atan(x2/x1)), mean, cov)} 
    # if Phi in 2nd or 4th quadrant, add pi to atan to get a positive phase
    else
      {deltamethod(list(~ sqrt(x1^2+x2^2), ~ pi+atan(x2/x1)), mean, cov)}
  
  # Calculate confidence intervals for A and Phi
  ci95_A <- list(A - 1.96*se_APhi[1], A + 1.96*se_APhi[1])
  ci95_Phi <- list(phi_rad - 1.96*se_APhi[2], phi_rad + 1.96*se_APhi[2])
  
  # Calculate pvalues for A and Phi
  pval_A <- 2 * (1 - pnorm(A/se_APhi[1]))
  pval_Phi <- 2 * (1 - pnorm(phi_rad/se_APhi[2]))
  
  # Put all in a table
  tab <- matrix(c(A, phi_rad,
                  se_APhi[1], se_APhi[2],
                  ci95_A[1], ci95_Phi[1],
                  ci95_A[2], ci95_Phi[2],
                  pval_A, pval_Phi), 
                ncol=5)
  colnames(tab) <- c('estimated','SE','95%CI_lower','95%CI_upper','pvalue')
  rownames(tab) <- c('A','Phi (rad)')
  names(A) <- 'A'
  names(phi_rad) <- 'phi_rad'
  print(tab)
  # Print worst infusion time
  phi_hour_float <- phi_rad * 24/(2*pi) # in hours
  phi_hour_clock <- tfloat_to_tclock(phi_hour_float)
  cat('\nPhi (floating hours) = ', phi_hour_float,'\nPhi (hour clock) = ', phi_hour_clock)
  return(c(A, phi_rad))
}


#### Convert floating hours to hour clock

tfloat_to_tclock <- function(tfloat){
  hh <- as.integer(tfloat)
  mm <- as.integer(round((tfloat - as.integer(tfloat))*60))
  # Deal with exceptions for intersection point > 24:00, setting bound to 23:59
  if (hh == 24){
    hh <- 23
    mm <- 59
  }
  tclock <- paste(hh,':',mm)
  return(tclock)
}


################ 
Main
################

#### Cox model for the infusion timing as a periodic variable

# Read data
df <- read_excel("/Users/simo/AtLeast20percInfus/new_db2.xlsx")

# Change units for infusion timing variable, from hour_float to radians, to input in our sinusoidal Cox model
df$t_median_rad <- df$t_median_float * 2*pi/24

# Sinusoidal Cox model for PS0-1 subpopulation
coxmod_ps01 <- coxph(Surv(surv_days, status) ~ cos(t_median_rad) + sin(t_median_rad), data=df, subset=(ps_immuno1_cat=='0-1'))
summary(coxmod_ps01)

# Get estimates, CIs, and pvalues for amplitude A and phase Phi
vals_ps01 <- calculate_APhi(coxmod_ps01)

# Calculate hazard ratio at worst, best and average infusion time
phi_rad <- vals_ps01[["phi_rad"]]
b1 <- coxmod_ps01$coeff[["cos(t_median_rad)"]]
b2 <- coxmod_ps01$coeff[["sin(t_median_rad)"]]

HR_worst <- exp(b1*cos(phi_rad) + b2*sin(phi_rad))
HR_best  <- exp(b1*cos(phi_rad+pi) + b2*sin(phi_rad+pi))
HR_avg   <- exp(b1*cos(phi_rad+pi/2) + b2*sin(phi_rad+pi/2))

cat(HR_worst, HR_best, HR_avg)

HR_worstOVERavg  <- HR_worst/HR_avg  # == 5.514201
HR_worstOVERbest <- HR_worst/HR_best # == 30.40641

