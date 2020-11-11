
load.lib=c('lavaan', 'dplyr')
install.lib<-load.lib[!load.lib %in% installed.packages()] # Install missing libraries
sapply(load.lib,require,character=TRUE) # Load libraries


####One polygenic score####
gsensY = function(rxy, rgx,rgy,n,h2,constrain=NULL,print=FALSE) { 
  mat=matrix(c(1,rgx,rgy,rgx,1,rxy,rgy,rxy,1),ncol=3,nrow=3)
  colnames(mat)=c("G","X","Y"); rownames(mat)=c("G","X","Y")
  
  model1 <- paste('
                  Y ~ bxy*X+bgy*GG # Y depends on X and true polygenic score
                  X ~ bgx*GG # X depends on true polygenic score
                  GG =~ l*G # true polygenic score is estimated by G
                  GG ~~ 1*GG # true polygenic score will be standardised
                  Y ~~ Y # residual error of Y
                  X ~~ X # residual error of X
                  G ~~ vg*G # measurement error in G due to SNP selection & sampling error
                  bgy+bgx*bxy == sqrt(',h2,') # heritability constraint
                  conf:=bgx*bgy                  #Confounding effect in standardized model
                  total:=conf+bxy
                                ',constrain,' #optional constraints
                  ')
 
   
fit1=lavaan(model1,sample.cov=mat,sample.nobs=n,estimator="GLS")
if (print) {  summary(fit1)}
pe=parameterestimates(fit1)
pe$pvalue=formatC(2*pnorm(-abs(pe$z)), digits=5)
results= data.frame(rbind(
                    pe[pe$label=="bxy",],
                    pe[pe$label=="conf",],
                    pe[pe$label=="total",]
                                        ))[,5:10]
results = results %>% 
  mutate_if(is.numeric, round, 3) # round all numeric variables



rownames(results)=c("Adjusted Bxy","Genetic confounding","Total effect")
results
}


gsensX = function(rxy, rgx,rgy,n,h2,constrain=NULL,print=FALSE) { 
  mat=matrix(c(1,rgx,rgy,rgx,1,rxy,rgy,rxy,1),ncol=3,nrow=3)
  colnames(mat)=c("G","X","Y"); rownames(mat)=c("G","X","Y")
  
  model1 <- paste('
                  Y ~ bxy*X+bgy*GG # Y depends on X and true polygenic score
                  X ~ bgx*GG # X depends on true polygenic score
                  GG =~ l*G # true polygenic score is estimated by G
                  GG ~~ 1*GG # true polygenic score will be standardised
                  Y ~~ Y # residual error of Y
                  X ~~ X # residual error of X
                  G ~~ vg*G # measurement error in G due to SNP selection & sampling error
                bgx==sqrt(',h2,')
                  conf:=bgx*bgy                  #Confounding effect in standardized model
                  total:=conf+bxy
                  ',constrain,' #optional constraints
                  ')
 
   
fit1=lavaan(model1,sample.cov=mat,sample.nobs=n,estimator="GLS")
if (print) {  summary(fit1)}
pe=parameterestimates(fit1)
pe$pvalue=formatC(2*pnorm(-abs(pe$z)), digits=5)
results= data.frame(rbind(
                    pe[pe$label=="bxy",],
                    pe[pe$label=="conf",],
                    pe[pe$label=="total",]
                                        ))[,5:10]
results = results %>% 
  mutate_if(is.numeric, round, 3) # round all numeric variables
  
rownames(results)=c("Adjusted Bxy","Genetic confounding","Total effect")
results
}


gsensYl = function(rxy, rgx,rgy,n,h2,constrain=NULL,print=FALSE) { 
  mat=matrix(c(1,rgx,rgy,rgx,1,rxy,rgy,rxy,1),ncol=3,nrow=3)
  colnames(mat)=c("G","X","Y"); rownames(mat)=c("G","X","Y")
  
  model1 <- paste('
                  Y ~ bxy*X+bgy*GG # Y depends on X and true polygenic score
                  X ~ bgx*GG # X depends on true polygenic score
                  GG =~ l*G # true polygenic score is estimated by G
                  GG ~~ 1*GG # true polygenic score will be standardised
                  Y ~~ Y # residual error of Y
                  X ~~ X # residual error of X
                  G ~~ vg*G # measurement error in G due to SNP selection & sampling error
                  bgy+bgx*bxy == sqrt(',h2,') # heritability constraint
                  conf:=bgx*bgy                  #Confounding effect in standardized model
                  total:=conf+bxy
                  rgxm:=l*bgx             
                   rgym:=l*(bgy+bgx*bxy)    
                    rxym:=bxy+bgx*bgy
                  ',constrain,' #optional constraints
                   ')

  fit1=lavaan(model1,sample.cov=mat,sample.nobs=n,estimator="ML",optim.method = "NLMINB.CONSTR",optim.dx.tol=0.01)
  if (print) {  summary(fit1)}
  pe=parameterestimates(fit1)
  pe$pvalue=formatC(2*pnorm(-abs(pe$z)), digits=5)
  results= data.frame(rbind(
    pe[pe$label=="bxy",],
    pe[pe$label=="conf",],
    pe[pe$label=="total",]
  ))[,5:10]
  results = results %>% 
    mutate_if(is.numeric, round, 3) # round all numeric variables
  
  rownames(results)=c("Adjusted Bxy","Genetic confounding","Total effect")
  results
}


####Two polygenic scores####
gsensXY = function(rxy,rg1x,rg2y,rg1y,rg2x,rg1g2,n,h2.x,h2.y,constrain=NULL,print=FALSE) {
  lower = c(
    1,
    rg1g2,1,
    rg1x,rg2x,1,
    rg1y,rg2y,rxy,1)
  mat=getCov(lower,names=c("G1","G2","X","Y"))
  
  model1 <- paste('
                  Y ~ bxy*X +  bg1y*GG1 + bg2y*GG2 # Y depends on X and true polygenic scores
                  X ~ bg1x*GG1 + bg2x*GG2  # X depends on true polygenic scores
                  GG1 =~ lg1*G1 # true polygenic score is estimated by G1
                  GG2 =~ lg2*G2 # true polygenic score is estimated by G2
                  GG1 ~~ 1*GG1 # true polygenic score will be standardised
                  GG2 ~~ 1*GG2 # true polygenic score will be standardised
                  Y ~~ vy*Y # residual error of Y
                  X ~~ vx*X # residual error of X
                  GG1 ~~ bg1g2*GG2 # correlation of true polygenic scores
                  G1 ~~ vg1*G1 # measurement error in G1 due to SNP selection & sampling error
                  G2 ~~ vg2*G2 # measurement error in G2 due to SNP selection & sampling error
                  #heritability constraints
                  bg1x+bg1g2*bg2x == sqrt(',h2.x,')
                  bg2y+(bg2x+bg1g2*bg1x)*bxy + bg1g2*bg1y== sqrt(',h2.y,') # SNP heritability constraints
                  #Confounding and total effects
                  conf:=bg1x*bg1y+bg2x*bg2y+bg1x*bg2y*bg1g2+bg2x*bg1y*bg1g2
                  total:=conf +bxy
        
                  ',constrain,' #optional constraints
                  
                  ')
  fit1=lavaan(model1,sample.cov=mat,sample.nobs=n,estimator="GLS")
  if (print) {  summary(fit1)}
  pe=parameterestimates(fit1)
  pe$pvalue=formatC(2*pnorm(-abs(pe$z)), digits=5)
  results= data.frame(rbind(
    pe[pe$label=="bxy",],
    pe[pe$label=="conf",],
    pe[pe$label=="total",]
  ))[,5:10]
  
  results = results %>% 
    mutate_if(is.numeric, round, 3) # round all numeric variables
  
  rownames(results)=c("Adjusted Bxy","Genetic confounding","Total effect")
  results
  
}



