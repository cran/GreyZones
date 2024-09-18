#' @export
#' @importFrom irrCAC kappa2.table
#' @importFrom irrCAC gwet.ac1.table
#' @importFrom irrCAC bp2.table
#' @importFrom irrCAC linear.weights
#' @importFrom irrCAC quadratic.weights
detectGreyZones.main <- function(table ){
  R <- nrow(table)
  n <- sum(table)
  
  kappa <- kappa2.table(table)
  kappa.linear <- kappa2.table(table, weights = linear.weights(1:R))
  kappa.quadratic <- kappa2.table(table, weights = quadratic.weights(1:R))
  AC2.linear <- gwet.ac1.table(table, weights = linear.weights(1:R))
  AC2.quadratic <- gwet.ac1.table(table, weights = quadratic.weights(1:R))
  BP.linear <- bp2.table(table, weights = linear.weights(1:R))
  BP.quadratic <- bp2.table(table, weights = quadratic.weights(1:R))
  
  expf <- array(0,dim=c(R,R))
  st.res <- array(0,dim=c(R,R))
  
  for (j in 1:R){
    for (k in 1:R){
      expf[j,k]= (table[j,k]+table[k,j])/2
      st.res[j,k]= (table[j,k]-expf[j,k])/sqrt(expf[j,k])
    }
  }
  st.res[st.res=="NaN"] <- 0
  
  delta <- st.res/kappa$coeff.val 
  delta2 <- max(delta)
  
  #-------------- Threshold ------------
  tau_Delta <- (-0.0080+0.4090*kappa$coeff.val^{2}+3.331*10^{-5}*n-2.467*10^{-8}*n^{2})^{-0.6266}
  
  if (delta2 > tau_Delta){# Grey zone detected, report AC2 and BP with quadratic weights.
    agreement <- data.frame(Weights = c("Quadratic"),
                            AC2 = AC2.quadratic$coeff.val,
                            BP = BP.quadratic$coeff.val)    
  } else { # No grey zone detected, report kappa, AC2 and BP with linear and quadratic weights.
    agreement <- data.frame(Weights=c("Linear","Quadratic"),
                            Kappa=c(kappa.linear$coeff.val,kappa.quadratic$coeff.val),
                            AC2=c(AC2.linear$coeff.val,AC2.quadratic$coeff.val),
                            BP=c(BP.linear$coeff.val,BP.quadratic$coeff.val)) 
  }
  
  return(list(delta = delta, Delta = delta2, tau_Delta = tau_Delta, 
              result = ifelse(delta2 > tau_Delta, "There is a gray zone.", 
                              "There is no gray zone."),
              agreement = agreement))
  }
