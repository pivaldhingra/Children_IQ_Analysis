library(glmnet)
set.seed(1000)
exp <- load("/Users/pivaldhingra/Downloads/exposome_NA.RData")
exposome <- exposomeNA[order(exposomeNA$ID),]
covariates <- covariatesNA[order(covariatesNA$ID),]
phenotype <- phenotypeNA[order(phenotypeNA$ID),]

a <- phenotype[5]
b <- covariates[12]
c <- exposome[,86:105]
expm <- data.frame(a,b,c)
expm <- expm[sample(nrow(expm)),]

M0 <- lm(hs_correct_raven ~ 1,data=expm)
Mfull <- lm(hs_correct_raven ~.,data=expm)

Mfull_summary <- summary(Mfull)

Mfull_R2 <- Mfull_summary$r.squared # R2
Mfull_R2
Mfull_adj_R2 <- Mfull_summary$adj.r.squared #adj R2
Mfull_adj_R2

length(coef(Mfull))
anyNA(coef(Mfull))
df.penalty <- 2
# forward 
system.time({
  Mfwd <- step(object = M0, # base model
               scope = list(lower = M0, upper = Mfull), # smallest and largest model 
               trace = 1, # trace prints out information
               direction = "forward" )
})
# backward 
system.time({
  Mback <- step(object = Mfull, # base model 
                scope = list(lower = M0, upper = Mfull), 
                direction = "backward", trace = 1)
})

system.time({
  Mstep <- step(object = M0,
                scope = list(lower = M0, upper = Mfull), 
                direction = "both", trace = 1)
})
Mfull_summary <- summary(Mfull)
Mfull_R2 <- Mfull_summary$r.squared # R2
Mfull_R2
Mfull_adj_R2 <- Mfull_summary$adj.r.squared #adj R2
Mfull_adj_R2

Mfwd_summary <- summary(Mfwd)
Mfwd_R2 <- Mfwd_summary$r.squared # R2
Mfwd_R2
Mfwd_adj_R2 <- Mfwd_summary$adj.r.squared #adj R2
Mfwd_adj_R2

Mstep_summary <- summary(Mstep)
Mstep_R2 <- Mstep_summary$r.squared # R2
Mstep_R2
Mstep_adj_R2 <- Mstep_summary$adj.r.squared #adj R2
Mstep_adj_R2

Mback_summary <- summary(Mback)
Mback_R2 <- Mback_summary$r.squared # R2
Mback_R2
Mback_adj_R2 <- Mback_summary$adj.r.squared #adj R2
Mback_adj_R2

Mfull_pval <- pf(Mfull_summary$fstatistic[1],
                 Mfull_summary$fstatistic[2],
                 Mfull_summary$fstatistic[3], lower.tail = FALSE) # p-val of F-statistics
Mfull_pval

Mfwd_pval <-  pf(Mfwd_summary$fstatistic[1],
                 Mfwd_summary$fstatistic[2],
                 Mfwd_summary$fstatistic[3], lower.tail = FALSE)
Mfwd_pval

Mback_pval <- pf(Mfwd_summary$fstatistic[1],
                 Mfwd_summary$fstatistic[2],
                 Mfwd_summary$fstatistic[3], lower.tail = FALSE)

Mback_pval

Mstep_pval <- pf(Mstep_summary$fstatistic[1],
                 Mstep_summary$fstatistic[2],
                 Mstep_summary$fstatistic[3], lower.tail = FALSE)

Mstep_pval
# compare the three different models
beta.fwd <- coef(Mfwd)
beta.back <- coef(Mback)
beta.step <- coef(Mstep)
c(fwd = length(beta.fwd), back = length(beta.back),
  step = length(beta.step)) # number of coefficients in each
# check if models are nested
names(beta.fwd)[!names(beta.fwd) %in% names(beta.back)]
names(beta.back)[!names(beta.back) %in% names(beta.fwd)]

## AIC
n <- nrow(expm)
ll_fwd <- -n/2 * (1 + log(sum(resid(Mfwd)^2)/n) + log(2*pi))
aic_fwd <- -2*ll_fwd + 2*(n - Mfwd$df + 1) # total number of parameters includes sigma
aic_fwd - AIC(Mfwd)
aic_step <- AIC(Mstep)
aic_back <- AIC(Mback)

aic_all <- round(c(aic_fwd,aic_step,aic_back),2)
names(aic_all) <- c("FWD","Step","Back")
aic_all

## BIC
bic_fwd <- -2*ll_fwd + log(n)*(n - Mfwd$df + 1) # total number of parameters includes sigma
bic_fwd - BIC(Mfwd)
bic_step <- BIC(Mstep)
bic_back <- BIC(Mback)

bic_all <- round(c(bic_fwd,bic_step,bic_back),2)
names(bic_all) <- c("FWD","Step","Back")
bic_all

## Random subsets cross-validation
# compare Mfull to Mstep
M1 <- Mfull
M2 <- Mstep
Mnames <- expression(M[FULL], M[STEP])

# number of cross-validation replications
nreps <- 1e3

ntot <- nrow(expm) # total number of observations
ntrain <- 800 # for fitting MLE's
ntest <- ntot-ntrain # for out-of-sample prediction

# storage space
mspe1 <- rep(NA, nreps) # mspe for M1
mspe2 <- rep(NA, nreps) # mspe for M2

system.time({
  for(ii in 1:nreps) {
    if(ii%%100 == 0) message("ii = ", ii)
    
    train.ind <- sample(ntot, ntrain) # training observations
    # long-form cross-validation
    M1.cv <- update(M1, subset = train.ind)
    M2.cv <- update(M2, subset = train.ind)
    # cross-validation residuals
    M1.res <- expm$hs_correct_raven[-train.ind] - # test observations
      predict(M1.cv, newdata = expm[-train.ind,]) # prediction with training data
    M2.res <- expm$hs_correct_rave[-train.ind] -predict(M2.cv, newdata = expm[-train.ind,])
    # mspe for each model
    mspe1[ii] <- mean(M1.res^2)
    mspe2[ii] <- mean(M2.res^2)
    
  }
})

# compare
par(mfrow = c(1,2))
cex <- 1
boxplot(x = list(mspe1, mspe2), names = Mnames,
        main = "MSPE using Random subsets cross-validation",
        #ylab = expression(sqrt(bar(SSE)[CV])),
        ylab = expression(MSPE),
        col = c("yellow", "orange"),
        cex = cex, cex.lab = cex, cex.axis = cex, cex.main = cex)
boxplot(x = list(sqrt(mspe1), sqrt(mspe2)), names = Mnames,
        main = "Root MSPE using Random subsets cross-validation",
        ylab = expression(sqrt(MSPE)),
        ## ylab = expression(SSE[CV]),
        col = c("yellow", "orange"),
        cex = cex, cex.lab = cex, cex.axis = cex, cex.main = cex)

# compare predictions by training set
par(mfrow=c(1,1))
plot(mspe1, mspe2, pch = 16,
     xlab = Mnames[1], ylab = Mnames[2],
     main = "Compare MSPE of MFull and Mstep")
abline(a = 0, b = 1, col= "red", lwd = 2)

## K-fold cross-validation 

# compare Mfwd to Mstep
set.seed(1000)
M1 <- Mfwd
M2 <- Mstep
Mnames <- expression(M[FWD], M[STEP])

expm <- expm[-1300,]

ntot <- nrow(expm) # total number of observations

# number of cross-validation replications
Kfolds <- 5


expm <- expm[sample(ntot),] # permute rows
expm$index <- rep(1:Kfolds,each=ntot/Kfolds)


# storage space
mspe1 <- rep(NA, Kfolds) # mspe for M1
mspe2 <- rep(NA, Kfolds) # mspe for M2

system.time({
  for(ii in 1:Kfolds) {
    
    train.ind <- which(expm$index!=ii) # training observations
    
    
    # using R functions
    M1.cv <- update(M1, subset = train.ind)
    M2.cv <- update(M2, subset = train.ind)
    # cross-validation residuals
    M1.res <- expm$hs_correct_raven[-train.ind] - # test observations
      predict(M1.cv, newdata = expm[-train.ind,]) # prediction with training data
    M2.res <- expm$hs_correct_raven[-train.ind] -predict(M2.cv, newdata = expm[-train.ind,])
    # mspe for each model
    mspe1[ii] <- mean(M1.res^2)
    mspe2[ii] <- mean(M2.res^2)
    
  }
})

# compare
par(mfrow = c(1,2))
cex <- 1
boxplot(x = list(mspe1, mspe2), names = Mnames,
        main = "MSPE using K-fold cross-validation",
        #ylab = expression(sqrt(bar(SSE)[CV])),
        ylab = expression(MSPE),
        col = c("yellow", "orange"),
        cex = cex, cex.lab = cex, cex.axis = cex, cex.main = cex)
boxplot(x = list(sqrt(mspe1), sqrt(mspe2)), names = Mnames,
        main = "Root MSPE using K-fold cross-validation",
        ylab = expression(sqrt(MSPE)),
        ## ylab = expression(SSE[CV]),
        col = c("yellow", "orange"),
        cex = cex, cex.lab = cex, cex.axis = cex, cex.main = cex)
mean(mspe1)
mean(mspe2)
PRESS1 <- resid(M1)/(1-hatvalues(M1))
PRESS2 <- resid(M2)/(1-hatvalues(M2))
# should match if doing LOO-CV
mean(PRESS1^2)
mean(PRESS2^2)

## LASSO method
set.seed(1000) 

## add 20 useless predictors
expm <- data.frame(cbind(expm,
                         matrix(rnorm(20*nrow(expm),0,1),ncol=20)))

X <- model.matrix(Mfull)[,-1] ## get covariates
y <- expm$hs_correct_raven  ## get outcome
Y <- na.omit(y)

## split into test and train
expm <- expm[sample(nrow(expm)),]
ntrain <- 500 
train_id <- 1:ntrain
X_train <- X[train_id,] 
X_test <- X[-train_id,]
Y_train <- Y[train_id]
Y_test <- Y[-train_id]
### LASSO
## fit models
M_lasso <- glmnet(x=X_train,y=Y_train,alpha = 1)

## plot paths
plot(M_lasso,xvar = "lambda",label=TRUE, main= "plotting paths according to LASSO method")
## fit with crossval
cvfit_lasso <-  cv.glmnet(x=X_train,y=Y_train,alpha = 1)

## plot MSPEs by lambda
plot(cvfit_lasso, main="plotting Mean Sqaured Prediction Error by lambda for LASSO method")

## estimated betas for minimum lambda 
coef(cvfit_lasso, s = "lambda.min")## alternatively could use "lambda.1se"

## predictions
pred_lasso <- predict(cvfit_lasso,newx=X_test,  s="lambda.min")
MSPE_lasso <- mean((pred_lasso-Y_test)^2)

## RIDGE
## fit models
M_ridge <- glmnet(x=X_train,y=Y_train,alpha = 0)

## plot paths
plot(M_ridge,xvar = "lambda",label=TRUE, main= "plotting paths accordinng to RIDGE method")

## fit with crossval
cvfit_ridge <-  cv.glmnet(x=X_train,y=Y_train,alpha = 0)

## plot MSPEs by lambda
plot(cvfit_ridge, main ="plotting Mean Sqaured Prediction Error by lambda for RIDGE method")

## estimated betas for minimum lambda 
coef(cvfit_ridge, s = "lambda.min")## alternatively could use "lambda.1se"

## predictions
pred_ridge <- predict(cvfit_ridge,newx=X_test,  s="lambda.min")

## MSPE in test set
MSPE_ridge <- mean((pred_ridge-Y_test)^2)



## compare prediction error for lasso and ridge
MSPE_lasso
MSPE_ridge




## Create a table to display the numbers
library(kableExtra)
tab <- matrix(c(0.5206107, 0.5200204, 0.5200204, 0.5200204,0.5126775, 0.5143736, 0.5143736, 
                0.5143736, 68.24673, 7626.89, 7558.64, 7558.64, 120.0868, 7766.50, 7646.41, 
                7646.41, 41.57002, NA, NA, NA, 41.7584, NA, NA, NA, 2.052856e-185, 8.152393e-191, 
                8.152393e-191, 8.152393e-191  ), ncol=4, byrow=TRUE)
colnames(tab) <- c('Full_mod', 'FWD', 'STEP', 'BACK')
rownames(tab) <- c('R2', 'adj_R2', 'AIC', "BIC", "MSPE for LASSO method ", "MSPE for RIDGE method", "F pval")
tab <- as.table(tab)
kbl(tab,
    caption = "The comparison between differnet methods to analyse the data",
    booktabs = T, digits = 5) %>%
  kable_styling(latex_options = c("striped", "hold_position"))



