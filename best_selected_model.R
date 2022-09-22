# Load package
library(coda)
library(tidyverse)
library(runjags)
library(rjags)

# Set working directory
setwd("...")

# Load data
myDataALL = read.csv("example_data.csv")

# Create dataframe only with the values needed in the models
myDataSOME = cbind.data.frame(Speed=myDataALL$Speed, Hours=myDataALL$Hours, ID=myDataALL$ID,
                              Sex=myDataALL$Sex, Season=myDataALL$Season, DistRoad1=myDataALL$DistRoad1,
                              DistRoad2=myDataALL$DistRoad2, DistSettl=myDataALL$DistSettl)

# Create matrix of data
DataArray = data.matrix(myDataSOME)

# Transform the hour 0 to 24
w=which(DataArray[,"Hours"]==0)
DataArray[w,"Hours"] <- 24

# Transform the response variable through the square root function
DataArray[,"Speed"] = DataArray[,"Speed"]^(1/2)

# If necessary, transform distances to anthropic features that are zero to values near zero
# to avoid NA in the models. E.g., in our case the problem was only in the distances to
# settlements, although with the example data this step is not necessary:
# w=which(DataArray[,"DistSettl"]==0)
# DataArray[w,"DistSettl"] <- 10^(-100)

# Create a vector of zeros with 24 entries which will be used in the model
zerovec = array(0, dim=c(24))

# Save the number of rows of the data
nr = nrow(DataArray)

# Create a vector for each season and sex that represents the positions in the data
# for the corrisponding season and sex. These will be used for the interaction terms in
# the model.
# SP = spring = 1, ES = early summer = 2, LS = late summer = 3, A = autumn = 4
# F = female = 1, M = male = 2
vec_SP_F = array(0, dim=c(nr))
vec = which(DataArray[,"Season"] == 1 & DataArray[,"Sex"] == 1)
vec_SP_F[vec] = 1

vec_ES_F = array(0, dim=c(nr))
vec = which(DataArray[,"Season"] == 2 & DataArray[,"Sex"] == 1)
vec_ES_F[vec] = 1

vec_LS_F = array(0, dim=c(nr))
vec = which(DataArray[,"Season"] == 3 & DataArray[,"Sex"] == 1)
vec_LS_F[vec] = 1

vec_A_F = array(0, dim=c(nr))
vec = which(DataArray[,"Season"] == 4 & DataArray[,"Sex"] == 1)
vec_A_F[vec] = 1

vec_SP_M = array(0, dim=c(nr))
vec = which(DataArray[,"Season"] == 1 & DataArray[,"Sex"] == 2)
vec_SP_M[vec] = 1

vec_ES_M = array(0, dim=c(nr))
vec = which(DataArray[,"Season"] == 2 & DataArray[,"Sex"] == 2)
vec_ES_M[vec] = 1

vec_LS_M = array(0, dim=c(nr))
vec = which(DataArray[,"Season"] == 3 & DataArray[,"Sex"] == 2)
vec_LS_M[vec] = 1

vec_A_M = array(0, dim=c(nr))
vec = which(DataArray[,"Season"] == 4 & DataArray[,"Sex"] == 2)
vec_A_M[vec] = 1

# Create a vector that represents the positions in the data for the male sex.
# This will be used for the interaction terms in the model.
# M = male = 2
VEC_M = array(0, dim=c(nr))
sexM = which(DataArray[,"Sex"] == 2)
VEC_M[sexM] = 1

# Create a list with the information needed in the JAGS model
dataList = list(
  DataArray = DataArray,
  zerovec = zerovec,
  nr = nr,
  vec_SP_F = vec_SP_F,
  vec_ES_F = vec_ES_F,
  vec_LS_F = vec_LS_F,
  vec_A_F = vec_A_F,
  vec_SP_M = vec_SP_M,
  vec_ES_M = vec_ES_M,
  vec_LS_M = vec_LS_M,
  vec_A_M = vec_A_M,
  VEC_M = VEC_M
)

# Define the model:
# - h_x_x represent the circadian effects in interaction with season and sex
# - i represents the random effects of the individual
# - s, t and t_m represent the sex, season and interaction term between sex and season,
# respectively
# - b, c, d and the exponentials represent the effects of primary and secondary roads,
# and human settlements, respectively
modelString = "
model {

for(z in 1:nr){
DataArray[z,1] ~ dnorm(
vec_SP_F[z]*h_SP_F[DataArray[z,2]] + vec_ES_F[z]*h_ES_F[DataArray[z,2]] +
vec_LS_F[z]*h_LS_F[DataArray[z,2]] + vec_A_F[z]*h_A_F[DataArray[z,2]] +
vec_SP_M[z]*h_SP_M[DataArray[z,2]] + vec_ES_M[z]*h_ES_M[DataArray[z,2]] +
vec_LS_M[z]*h_LS_M[DataArray[z,2]] + vec_A_M[z]*h_A_M[DataArray[z,2]] +
i[DataArray[z,3]] + s[DataArray[z,4]] + t[DataArray[z,5]] +
VEC_M[z]*t_m[DataArray[z,5]] +
b[DataArray[z,5]]*exp(-B[DataArray[z,5]]*DataArray[z,6]) +
c[DataArray[z,5]]*exp(-C[DataArray[z,5]]*DataArray[z,7]) +
d[DataArray[z,5]]*exp(-D[DataArray[z,5]]*DataArray[z,8]),
rigma)
}

rigma ~ dgamma(1,1)

for(n in 1:4){
b[n] ~ dnorm(0, 1/1000)
c[n] ~ dnorm(0, 1/1000)
d[n] ~ dnorm(0, 1/1000)
}

for(n in 1:4){
B[n] ~ dunif(3/4,3/0.2)
C[n] ~ dunif(3/4,3/0.2)
D[n] ~ dunif(3/4,3/0.2)
}

for(k1 in 1:23){
for(k2 in (k1+1):24){
Sigma_SP_F[k1,k2] = exp(-nu_SP_F*min(k2-k1, 24+k1-k2))/psi_SP_F
Sigma_SP_F[k2,k1] = Sigma_SP_F[k1,k2]
}
}
for(k1 in 1:24){
Sigma_SP_F[k1,k1] = 1/psi_SP_F
}
SigmaMat_SP_F = inverse(Sigma_SP_F)
h_SP_F[1:24] ~ dmnorm(zerovec,SigmaMat_SP_F)
nu_SP_F ~ dunif(3/24,3/1)
psi_SP_F ~ dgamma(1,1)


for(k1 in 1:23){
for(k2 in (k1+1):24){
Sigma_ES_F[k1,k2] = exp(-nu_ES_F*min(k2-k1, 24+k1-k2))/psi_ES_F
Sigma_ES_F[k2,k1] = Sigma_ES_F[k1,k2]
}
}
for(k1 in 1:24){
Sigma_ES_F[k1,k1] = 1/psi_ES_F
}
SigmaMat_ES_F = inverse(Sigma_ES_F)
h_ES_F[1:24] ~ dmnorm(zerovec,SigmaMat_ES_F)
nu_ES_F ~ dunif(3/24,3/1)
psi_ES_F ~ dgamma(1,1)


for(k1 in 1:23){
for(k2 in (k1+1):24){
Sigma_LS_F[k1,k2] = exp(-nu_LS_F*min(k2-k1, 24+k1-k2))/psi_LS_F
Sigma_LS_F[k2,k1] = Sigma_LS_F[k1,k2]
}
}
for(k1 in 1:24){
Sigma_LS_F[k1,k1] = 1/psi_LS_F
}
SigmaMat_LS_F = inverse(Sigma_LS_F)
h_LS_F[1:24] ~ dmnorm(zerovec,SigmaMat_LS_F)
nu_LS_F ~ dunif(3/24,3/1)
psi_LS_F ~ dgamma(1,1)


for(k1 in 1:23){
for(k2 in (k1+1):24){
Sigma_A_F[k1,k2] = exp(-nu_A_F*min(k2-k1, 24+k1-k2))/psi_A_F
Sigma_A_F[k2,k1] = Sigma_A_F[k1,k2]
}
}
for(k1 in 1:24){
Sigma_A_F[k1,k1] = 1/psi_A_F
}
SigmaMat_A_F = inverse(Sigma_A_F)
h_A_F[1:24] ~ dmnorm(zerovec,SigmaMat_A_F)
nu_A_F ~ dunif(3/24,3/1)
psi_A_F ~ dgamma(1,1)

for(k1 in 1:23){
for(k2 in (k1+1):24){
Sigma_SP_M[k1,k2] = exp(-nu_SP_M*min(k2-k1, 24+k1-k2))/psi_SP_M
Sigma_SP_M[k2,k1] = Sigma_SP_M[k1,k2]
}
}
for(k1 in 1:24){
Sigma_SP_M[k1,k1] = 1/psi_SP_M
}
SigmaMat_SP_M = inverse(Sigma_SP_M)
h_SP_M[1:24] ~ dmnorm(zerovec,SigmaMat_SP_M)
nu_SP_M ~ dunif(3/24,3/1)
psi_SP_M ~ dgamma(1,1)


for(k1 in 1:23){
for(k2 in (k1+1):24){
Sigma_ES_M[k1,k2] = exp(-nu_ES_M*min(k2-k1, 24+k1-k2))/psi_ES_M
Sigma_ES_M[k2,k1] = Sigma_ES_M[k1,k2]
}
}
for(k1 in 1:24){
Sigma_ES_M[k1,k1] = 1/psi_ES_M
}
SigmaMat_ES_M = inverse(Sigma_ES_M)
h_ES_M[1:24] ~ dmnorm(zerovec,SigmaMat_ES_M)
nu_ES_M ~ dunif(3/24,3/1)
psi_ES_M ~ dgamma(1,1)


for(k1 in 1:23){
for(k2 in (k1+1):24){
Sigma_LS_M[k1,k2] = exp(-nu_LS_M*min(k2-k1, 24+k1-k2))/psi_LS_M
Sigma_LS_M[k2,k1] = Sigma_LS_M[k1,k2]
}
}
for(k1 in 1:24){
Sigma_LS_M[k1,k1] = 1/psi_LS_M
}
SigmaMat_LS_M = inverse(Sigma_LS_M)
h_LS_M[1:24] ~ dmnorm(zerovec,SigmaMat_LS_M)
nu_LS_M ~ dunif(3/24,3/1)
psi_LS_M ~ dgamma(1,1)


for(k1 in 1:23){
for(k2 in (k1+1):24){
Sigma_A_M[k1,k2] = exp(-nu_A_M*min(k2-k1, 24+k1-k2))/psi_A_M
Sigma_A_M[k2,k1] = Sigma_A_M[k1,k2]
}
}
for(k1 in 1:24){
Sigma_A_M[k1,k1] = 1/psi_A_M
}
SigmaMat_A_M = inverse(Sigma_A_M)
h_A_M[1:24] ~ dmnorm(zerovec,SigmaMat_A_M)
nu_A_M ~ dunif(3/24,3/1)
psi_A_M ~ dgamma(1,1)

alpha ~ dnorm(0, 1/1000)
sigma ~ dgamma(1,1)
for(n in 1:18){ 
i[n] ~ dnorm(alpha, sigma)
}

s[1] <- 0
s[2] ~ dnorm(0, 1/1000)

t[1] <- 0
for(i in 2:4){
t[i] ~ dnorm(0, 1/1000)
}

t_m[1] <- 0
for(v in 2:4){ 
t_m[v] ~ dnorm(0, 1/1000)
}

}
" # close quote for modelString

# Save model
writeLines( modelString , con="TEMPmodel.txt" )

# Run two parallel chains and save variable values
runJagsOut <- run.jags( method="parallel" ,
                        model="TEMPmodel.txt" ,
                        monitor=c("alpha","h_SP_F", "h_ES_F", "h_LS_F", "h_A_F",
                                  "h_SP_M", "h_ES_M", "h_LS_M", "h_A_M","i",
                                  "sigma","s","t","t_m","rigma","b",
                                  "c","d","B","C","D","nu_SP_F", "psi_SP_F",
                                  "nu_ES_F", "psi_ES_F", "nu_LS_F", "psi_LS_F",
                                  "nu_A_F", "psi_A_F","nu_SP_M", "psi_SP_M",
                                  "nu_ES_M", "psi_ES_M", "nu_LS_M", "psi_LS_M",
                                  "nu_A_M", "psi_A_M") ,
                        data=dataList ,
                        n.chains= 2 ,
                        burnin= 50000 , 
                        sample= 2000,
                        thin= 10,
                        summarise=FALSE ,
                        plots=FALSE )

codaSamples = as.mcmc.list( runJagsOut ) # converts the runjags object to coda format, useful to
# operate on posterior distributions

# Save the DIC value for model selection
dic <- runjags::extract(runJagsOut, what = "dic")

# Save workspace
save.image("workspace.Rdata")

