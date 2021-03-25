#########################################################################################
#                                                                                       #
# The following code performs fixed-effects aster analyses on data examining the effects#
# of density and effective genetic population size (Ne) on female (seeds seet) fitness. #
# LM analysis of biomass (above and below ground) and number of aborted ovules follows  #
# the aster analyses.                                                                   #
#                                                                                       #    
#                                                                                       #
#                                                                                       #
# For any questions, please contact Mason Kulbaba: mason.kulbaba@gmail.com              #
#                                                                                       #
#                                                                                       #
#########################################################################################

#set your working directory
setwd("C:/Users/Mason Kulbaba/Dropbox/git/density-gen-Ne")

#Begin with analysis of female fitenss (seeds set)


#Load data
fin<- read.csv("data/aster.dat.csv")

#ensure class of factor variables
fin$Den<- as.factor(fin$Den)
fin$Gen<- as.factor(fin$Gen)
fin$plotID<- as.factor(fin$plotID)
fin$plantID<- as.factor(fin$plantID)
fin$familyID<- as.factor(fin$familyID)

#load aster
library(aster)

#set graphical model: Flower Number -> Fruit Number -> Subsampled Fruit -> Seed Number

#set vars
vars<- c( "flw", "frt", "frt.2","seeds")

#reshape data so that all response variables are located in a single vector in a new data
#set called "redata"
redata <- reshape(fin, varying = list(vars), direction = "long",timevar = "varb", times = as.factor(vars), v.names = "resp")


#Designation of terminal fitness variable (aka "artifice" variable...I don't like that name)
fit <- grepl("seeds", as.character(redata$varb))
fit<- as.numeric(fit)

#add fit to redata file
redata$fit <- fit

#check
with(redata, sort(unique(as.character(varb)[fit == 0])))
with(redata, sort(unique(as.character(varb)[fit == 1])))

#add a variable "root" to redata, where value is 1. This is the "starting point"
# of the aster graphical model (i.e. a seed planted)
redata<- data.frame(redata, root=1)


#set graph. model and family for each node
pred<- c(0,1,2,3) # 4 nodes
fam<- c(1,2,1,2)

#check family 
sapply(fam.default(), as.character)[fam]

#run aster with only fitness
aout<- aster(resp~varb , pred, fam, varb, id, root, data=redata)

summary(aout, show.graph=T)

#add density
aout.d<- aster(resp~varb + fit:(Den), pred, fam, varb, id, root, data=redata)

summary(aout.d, show.graph=T)

#liklihood ratio test
anova(aout, aout.d) #density is significant

#add density and Ne
aout.dg<- aster(resp~varb + fit:(Den + Gen), pred, fam, varb, id, root, data=redata)

summary(aout.dg, show.graph=T)

anova(aout.d, aout.dg)#Ne is significant

#add interaction between Den & Gen

aout.dg2<- aster(resp~varb + fit:(Den + Gen + Den*Gen), pred, fam, varb, id, root, data=redata)

summary(aout.dg2)

anova(aout.dg, aout.dg2)# interaction is significant


#add plot ID 


aoutc<- aster(resp~varb + fit:(plotID + Den + Gen + Den*Gen), pred, fam, varb, id, root, data=redata)

summary(aoutc, show.graph=T)


anova(aout.dg2, aoutc)# plotID significant


###################################################################################
#Generate Fitness Estimates for High and Low Ne  with respect to Density Treatment#
###################################################################################

#The aster analyses for high and low Ne are performed in parallel below (i.e. each step
# is performed twice, once for high Ne and once for low Ne analses).

#First, isloate high (HG) and low (LG) data from the main redata file. Therefore, don't 
# have to do "reshape" data step. HG/hg = High Ne, LG/lg = Low Ne

redataHG<- subset(redata, Gen=="HG")

redataLG<- subset(redata, Gen=="LG")

#drop unused levels in each file
redataHG<- droplevels(redataHG)

redataLG<- droplevels(redataLG)

#aster analysis with just fitness
aoutHG<- aster(resp~varb, pred, fam, varb, id, root, data=redataHG)


aoutHG2<- aster(resp~varb + fit:plotID, pred, fam, varb, id, root, data=redataHG)

aoutHG3<- aster(resp~varb + fit:(Den+ plotID), pred, fam, varb, id, root, data=redataHG)

aoutHG4<- aster(resp~varb + fit:(Den+ plotID +familyID), pred, fam, varb, id, root, data=redataHG)


summary(aoutHG, show.graph = T)
summary(aoutHG2, show.graph=T)
summary(aoutHG3, show.graph=T)
summary(aoutHG4, show.graph=T)


anova(aoutHG, aoutHG2, aoutHG3, aoutHG4)# density is significant in high Ne


#LG with just fitness
aoutLG<- aster(resp~varb, pred, fam, varb, id, root, data=redataLG)

aoutLG1<- aster(resp~varb + fit:plotID, pred, fam, varb, id, root, data=redataLG)

aoutLG2<- aster(resp~varb + fit:(plotID + Den), pred, fam, varb, id, root, data=redataLG)

aoutLG3<- aster(resp~varb + fit:(plotID + Den + familyID), pred, fam, varb, id, root, data=redataLG)


summary(aoutLG, show.graph = T)
summary(aoutLG1, show.graph=T)
summary(aoutLG2, show.graph=T)
summary(aoutLG3, show.graph=T)

anova(aoutLG, aoutLG1, aoutLG2, aoutLG3)# density significant in low Ne

######################################################
#####Fitness estimates for high and low Ne plants#####
#####across three density treatmetns (H, M, L)########
######################################################

#generate MLE of saturated model mean value parameter vector: mu
pout.HG<- predict(aoutHG3, se.fit =TRUE)

pout.LG<- predict(aoutHG, se.fit=TRUE)

##############



#make design matrix data.frame of indivudals for each density level (low, med., high)
fred.hg <- data.frame(Den=levels(redataHG$Den), plotID=levels(redataHG$plotID), flw=1, frt=1,frt.2=1, seeds=1,root = 1)

fred.lg <- data.frame( Den=levels(redataLG$Den), flw=1, frt=1,frt.2=1, seeds=1,root = 1)


#reshape the design matrix just as the actual data
renewdata.hg <- reshape(fred.hg, varying = list(vars),
                     direction = "long", timevar = "varb",
                     times = as.factor(vars), v.names = "resp")

renewdata.lg <- reshape(fred.lg, varying = list(vars),
                        direction = "long", timevar = "varb",
                        times = as.factor(vars), v.names = "resp")

#make character string from "varb" of renewdata,
#without actual values (i.e., the layers of varb in renewdata)
layer<- gsub("[0-9]", "", as.character(renewdata.hg$varb))


#add layer to renewdata
renewdata.hg<- data.frame(renewdata.hg, layer= layer)

renewdata.lg<- data.frame(renewdata.lg, layer= layer)

#seed seed.ct in new layer col of renewdata as numeric, called fit
# note, only need one "fit" object as it is the same for both High and Low Ne
fit<- as.numeric(layer=="seeds")


########################################################
# This part is new

#Generate fintess estimates and standard errors for each treatment
nDen<- nrow(fred.hg)
nnode<- length(vars)
amat<- array(0, c(nDen, nnode, nDen))
dim(amat)# makes an 3 x 8 x 3 matrix (3 region types and 8 nodes of graphicla model)

#only want means for k'th individual that contribute to expected
#fitness, and want to add only sm.2 entries

foo<- grepl("seeds", vars)
for(k in 1:nDen)
  amat[k, foo, k]<- 1

#check
foo #yes, only last node is "true"; corresponds to seeds (sm.2)

#generate predicted valuses using aout object, with renewdata, and amat format
pout.amat<- predict(aoutHG3, newdata= renewdata.hg, varvar= varb,
                    idvar= id, root = root, se.fit=TRUE, amat = amat)

foo<- cbind(pout.amat$fit, pout.amat$se.fit)

rownames(foo)<- as.character(fred.hg$Den)
colnames(foo)<- c("Estimate", "Standard Error")
foo

#add fit to renewdata
renewdata.hg<- data.frame(renewdata.hg, fit = fit)

renewdata.lg<- data.frame(renewdata.lg, fit = fit)

#rerun prediction of aout, with "made up" renewdata
pout.hg<- predict(aoutHG2, newdata= renewdata.hg, varvar= varb,
               idvar = id, root = root, se.fit = TRUE)

pout.lg<- predict(aoutLG2, newdata= renewdata.lg, varvar= varb,
               idvar = id, root = root, se.fit = TRUE)


sapply(pout.hg, class)
sapply(pout.lg, class)
# lengths of fit and se.fit (12) match row number of renewdata 
#(as should be with predict.aster)
sapply(pout.hg, length)
sapply(pout.lg, length)

#therefore, can make 12 CIs, one for each of 4 nodes of graphical model, 
#and 3 density treatmetns (4 x 3 =12).


#put the parameter estimates into a matrix with individuals in rows
#and nodes in columns


#extract HG resutls
nnode<- length(vars)
sally.hg<- matrix(pout.hg$fit, ncol = nnode)
dim(sally.hg)# makes 3 x 4 matrix: 3 densities by 4 nodes

#name the rows (by Den Treat) and columns (as nodes)
rownames(sally.hg)<- as.character(fred.hg$Den) #is this correct?
colnames(sally.hg)<- unique(as.character(renewdata.hg$varb))


#view matrix 
round(sally.hg, 3)

#now generate matrix of standard errors
nnode2<- length(vars)
sally2<- matrix(pout.hg$se.fit, ncol = nnode)
dim(sally2)# makes 3 x 4 matrix: 3 densities by 4 nodes

#name the rows (by Den Treat) and columns (as nodes)
rownames(sally2)<- as.character(fred.hg$Den) #is this correct?
colnames(sally2)<- unique(as.character(renewdata.hg$varb))

#view matrix 
round(sally2, 3)

#combine estimates with standard errors for only final node: seeds
ests<- sally.hg[,grepl("seeds", colnames(sally.hg))]
se<-  sally2[,grepl("seeds", colnames(sally2))]

HG<- cbind(ests, se)# these are the fitness ests (and SE) for high Ne, across 3 densities

HG
#now extract LG results

#extrac HG resutls
nnode<- length(vars)
sally.lg<- matrix(pout.lg$fit, ncol = nnode)
dim(sally.lg)# makes 3 x 4 matrix: 3 densities by 4 nodes

#name the rows (by Den Treat) and columns (as nodes)
rownames(sally.lg)<- unique(as.character(renewdata.lg$Den))
colnames(sally.lg)<- unique(as.character(renewdata.lg$varb))

#view matrix 
round(sally.lg, 3)

#now generate matrix of standard errors

nnode2<- length(vars)
sally.lg2<- matrix(pout.lg$se.fit, ncol = nnode)
dim(sally.lg2)# makes 3 x 4 matrix: 3 densities by 4 nodes

#name the rows (by Den Treat) and columns (as nodes)
rownames(sally.lg2)<- unique(as.character(renewdata.lg$Den))
colnames(sally.lg2)<- unique(as.character(renewdata.lg$varb))

#view matrix 
round(sally.lg2, 3)

#combine estimates with standard errors for only final node: seeds
ests<- sally.lg[,grepl("seeds", colnames(sally.lg))]
se<-  sally.lg2[,grepl("seeds", colnames(sally.lg2))]

LG<- cbind(ests, se)# these are the fitness ests (and SE) for high Ne, across 3 densities

#fitness and standard errors for HG and LG treatments (across densities)
HG

LG

#Estiamte of male fitness (seeds sired)

#See markdown file for estimates of density-specific male fitness

fin<- read.csv("data/aster.sire.dat.csv")

fin$Den<- as.factor(fin$Den)
fin$Gen<- as.factor(fin$Gen)
fin$plotID<- as.factor(fin$plotID)
fin$plantID<- as.factor(fin$plantID)
fin$familyID<- as.factor(fin$familyID)


vars<- c( "flw", "frt", "sires")

redata <- reshape(fin, varying = list(vars), direction = "long",timevar = "varb", times = as.factor(vars), v.names = "resp")

fit <- grepl("sires", as.character(redata$varb))
fit<- as.numeric(fit)

redata$fit <- fit

with(redata, sort(unique(as.character(varb)[fit == 0])))
with(redata, sort(unique(as.character(varb)[fit == 1])))

redata<- data.frame(redata, root=2)

pred<- c(0,1,2)
fam<- c(1,2,2)

sapply(fam.default(), as.character)[fam]

aout<- aster(resp~varb, pred, fam, varb, id, root, data=redata)

summary(aout)

aout2<- aster(resp~varb +fit:(plotID), pred, fam, varb, id, root, data=redata)

summary(aout2)

aout3<- aster(resp~varb +fit:(plotID + Den), pred, fam, varb, id, root, data=redata)

summary(aout3)

aout4<- aster(resp~varb +fit:(plotID + Den + familyID), pred, fam, varb, id, root, data=redata)

summary(aout4)

anova(aout, aout2, aout3, aout4)

##################################################################################
#Calculate mean seed set for each treat x Ne treat to "relativize" female fitness#
##################################################################################

#calculate mean seed set for each individual plot

aggregate(fin$seeds, by=list(fin$plotID), mean)


#identify families in each plot

hg<- subset(fin, Gen=="HG")
lg<- subset(fin, Gen=="LG")


aggregate(hg$familyID, by=list(hg$Den), unique)

aggregate(lg$familyID, by=list(lg$Den), unique)


###############################################################################################
###############################################################################################
#To compliment familyID estimates of fitness via siring success, produce familyID estimates of 
#fitness via seed set, on L, M, H density for High Ne
###############################################################################################
###############################################################################################

#Load data
fin<- read.csv("C:/Users/Mason Kulbaba/Dropbox/git/density-gen-Ne/data/aster.sire.dat.csv")

#sum of number of seeds that we successfully assigned paternity to, per density treatment
aggregate(fin$sires, by=list(fin$Den), sum)

#change class of factor variables
fin$Den<- as.factor(fin$Den)
fin$Gen<- as.factor(fin$Gen)
fin$plotID<- as.factor(fin$plotID)
fin$plantID<- as.factor(fin$plantID)
fin$familyID<- as.factor(fin$familyID)

#load aster
library(aster)

#set vars
vars<- c( "flw", "frt", "frt.2","seeds")


datH<- subset(fin, Den=="H")
datH$familyID<- droplevels(datH$familyID)
datM<- subset(fin, Den=="M")
datM$familyID<- droplevels(datM$familyID)
datL<- subset(fin, Den=="L")
datL$familyID<- droplevels(datL$familyID)

#Perform the same data reshaping steps as before

#reshape data so that all response variables are located in a single vector in a new data
#set called "redata"
redataH <- reshape(datH, varying = list(vars), direction = "long",timevar = "varb", times = as.factor(vars), v.names = "resp")
redataM <- reshape(datM, varying = list(vars), direction = "long",timevar = "varb", times = as.factor(vars), v.names = "resp")
redataL <- reshape(datL, varying = list(vars), direction = "long",timevar = "varb", times = as.factor(vars), v.names = "resp")


#Designation of fitness variable
fit <- grepl("seeds", as.character(redataH$varb))
fit<- as.numeric(fit)

#add fit to each of three redata files (one for each ensity treat)
redataH$fit <- fit
redataM$fit <- fit
redataL$fit <- fit

#check
with(redataH, sort(unique(as.character(varb)[fit == 0])))
with(redataH, sort(unique(as.character(varb)[fit == 1])))

with(redataM, sort(unique(as.character(varb)[fit == 0])))
with(redataM, sort(unique(as.character(varb)[fit == 1])))

with(redataL, sort(unique(as.character(varb)[fit == 0])))
with(redataL, sort(unique(as.character(varb)[fit == 1])))

#add a variable "root" to redata, where value is 2
# value is 2 here to compliment the male fitness estimates (see density_male.R)
# Only working with High Ne (2 full-sib individuals from 6 families), and can't always 
# asign paternity between two full-sibs. So, we assigned paternity to families

redataH<- data.frame(redataH, root=2)
redataM<- data.frame(redataM, root=2)
redataL<- data.frame(redataL, root=2)

#set graph. model and family for each node
pred<- c(0,1,2,3)
fam<- c(1,2,2,2)


#aster analyses for three density treatments
aoutH<- aster(resp~varb + fit:(familyID), pred, fam, varb, id, root, data=redataH)
summary(aoutH, show.graph = T)

aoutM<- aster(resp~varb + fit:(familyID), pred, fam, varb, id, root, data=redataM)
summary(aoutM, show.graph = T)

aoutL<- aster(resp~varb + fit:(familyID), pred, fam, varb, id, root, data=redataL)
summary(aoutL, show.graph = T)

#generate estimates: High -> Med -> Low

#####################
####HIGH DENSITY#####
#####################


#generate MLE of saturated model mean value parameter vector: mu
pout<- predict.aster(aoutH, se.fit=TRUE)


#make design matrix
fred <- data.frame(familyID=levels(redataH$familyID), flw=1, frt=1, frt.2=1, seeds=1,root = 2)

#reshape the "made up data" just as the actual data
renewdata <- reshape(fred, varying = list(vars),
                     direction = "long", timevar = "varb",
                     times = as.factor(vars), v.names = "resp")

#make character string from "varb" of renewdata,
#without actual values (e.g., the layers of varb in renewdata)
layer<- gsub("[0-9]", "", as.character(renewdata$varb))

#add layer to renewdata
renewdata<- data.frame(renewdata, layer= layer)

#add fit (seeds) in new layer col of renewdata as numeric, called fit
fit<- as.numeric(layer=="seeds")

#add fit to renewdata
renewdata<- data.frame(renewdata, fit = fit)

#rerun prediction of aout, with "made up" renewdata
pout<- predict(aoutH, newdata= renewdata, varvar= varb,
               idvar = id, root = root, se.fit = TRUE)

sapply(pout, class)

sapply(pout, length)

#put the parameter estimates into a matrix with individuals in rows
#and nodes along columns
nnode<- length(vars)
sally<- matrix(pout$fit, ncol = nnode)
dim(sally)# makes 10 x 4 matrix: 10 families by 4 nodes

#name the rows (by familyID) and columns (as nodes)
rownames(sally)<- unique(as.character(renewdata$familyID))
colnames(sally)<- unique(as.character(renewdata$varb))

#view matrix 
round(sally, 3)

#use just totalseeds as predicted (expected) fitneses
herman<- sally[,grepl("seeds", colnames(sally))]

#Generate Standard Errors for these estimates

nFam<- nrow(fred)
nnode<- length(vars)
amat<- array(0, c(nFam, nnode, nFam))
dim(amat)# makes an 10 x 4 x 10 matrix


foo<- grepl("seeds", vars)
for(k in 1:nFam)
  amat[k, foo, k]<- 1

#use aout object, with renewdata, and amat format
pout.amat<- predict(aoutH, newdata= renewdata, varvar= varb,
                    idvar= id, root = root, se.fit=TRUE, amat = amat)

#pout.amat$fit should be the same as file "herman"
herman
pout.amat$fit #they are the same.  Good.

#combine std.err with estimates, and then round
#to three decimal places
foo<- cbind(pout.amat$fit, pout.amat$se.fit)
rownames(foo)<- as.character(fred$familyID)
colnames(foo)<- c("High Den Fitness", "SE")
round(foo, 3) 
H_estimates<- round(foo,3)

H_estimates

########################
#####MEDIUM DENSITY#####
########################

#generate MLE of saturated model mean value parameter vector: mu
pout<- predict.aster(aoutM, se.fit=TRUE)


#make design matrix
fred <- data.frame(familyID=levels(redataM$familyID), flw=1, frt=1, frt.2=1, seeds=1,root = 2)

#reshape the "made up data" just as the actual data
renewdata <- reshape(fred, varying = list(vars),
                     direction = "long", timevar = "varb",
                     times = as.factor(vars), v.names = "resp")

#make character string from "varb" of renewdata,
#without actual values (e.g., the layers of varb in renewdata)
layer<- gsub("[0-9]", "", as.character(renewdata$varb))

#add layer to renewdata
renewdata<- data.frame(renewdata, layer= layer)

#seed seed.ct in new layer col of renewdata as numeric, called fit
fit<- as.numeric(layer=="seeds")

#add fit to renewdata
renewdata<- data.frame(renewdata, fit = fit)

#rerun prediction of aout, with "made up" renewdata
pout<- predict(aoutM, newdata= renewdata, varvar= varb,
               idvar = id, root = root, se.fit = TRUE)

sapply(pout, class)
sapply(pout, length)

#put the parameter estimates into a matrix with familyID in rows
#and nodes along columns
nnode<- length(vars)
sally<- matrix(pout$fit, ncol = nnode)
dim(sally)# makes 9 x 4 matrix: 9 families by 4 nodes

#name the rows (by Den Treat) and columns (as nodes)
rownames(sally)<- unique(as.character(renewdata$familyID))
colnames(sally)<- unique(as.character(renewdata$varb))

#view matrix 
round(sally, 3)

#use just totalseeds as predicted (expected) fitneses
herman<- sally[,grepl("seeds", colnames(sally))]

#Generate Standard Errors for these estimates

nFam<- nrow(fred)
nnode<- length(vars)
amat<- array(0, c(nFam, nnode, nFam))
dim(amat)# makes an 9 x 4 x 9 matrix


foo<- grepl("seeds", vars)
for(k in 1:nFam)
  amat[k, foo, k]<- 1

#use aout object, with renewdata, and amat format
pout.amat<- predict(aoutM, newdata= renewdata, varvar= varb,
                    idvar= id, root = root, se.fit=TRUE, amat = amat)

#pout.amat$fit should be the same as file "herman"
herman
pout.amat$fit #they are the same.  Good.

#combine std.err with estimates, and then round
#to three decimal places
foo<- cbind(pout.amat$fit, pout.amat$se.fit)
rownames(foo)<- as.character(fred$familyID)
colnames(foo)<- c("Med Den Fitness", "SE")
round(foo, 3) 
M_estimates<- round(foo, 3)

M_estimates

#####################
#####LOW DENSITY#####
#####################

#generate MLE of saturated model mean value parameter vector: mu
pout<- predict.aster(aoutL, se.fit=TRUE)


#make data.frame of indivudals for each block (1-8)
fred <- data.frame(familyID=levels(redataL$familyID), flw=1, frt=1, frt.2=1, seeds=1,root = 2)

#reshape the "made up data" just as the actual data
renewdata <- reshape(fred, varying = list(vars),
                     direction = "long", timevar = "varb",
                     times = as.factor(vars), v.names = "resp")

#make character string from "varb" of renewdata,
#without actual values (e.g., the layers of varb in renewdata)
layer<- gsub("[0-9]", "", as.character(renewdata$varb))

#add layer to renewdata
renewdata<- data.frame(renewdata, layer= layer)

#seed seed.ct in new layer col of renewdata as numeric, called fit
fit<- as.numeric(layer=="seeds")

#add fit to renewdata
renewdata<- data.frame(renewdata, fit = fit)

#rerun prediction of aout, with "made up" renewdata
pout<- predict(aoutM, newdata= renewdata, varvar= varb,
               idvar = id, root = root, se.fit = TRUE)

sapply(pout, class)
sapply(pout, length)

#put the parameter estimates into a matrix with familyID in rows
#and nodes along columns
nnode<- length(vars)
sally<- matrix(pout$fit, ncol = nnode)
dim(sally)# makes 8 x 4 matrix: 8 families by 4 nodes

#name the rows (by Den Treat) and columns (as nodes)
rownames(sally)<- unique(as.character(renewdata$familyID))
colnames(sally)<- unique(as.character(renewdata$varb))

#view matrix 
round(sally, 3)

#use just totalseeds as predicted (expected) fitneses
herman<- sally[,grepl("seeds", colnames(sally))]

#Generate Standard Errors for these estimates

nFam<- nrow(fred)
nnode<- length(vars)
amat<- array(0, c(nFam, nnode, nFam))
dim(amat)# makes an 8 x 4 x 8 matrix


foo<- grepl("seeds", vars)
for(k in 1:nFam)
  amat[k, foo, k]<- 1

#use aout object, with renewdata, and amat format
pout.amat<- predict(aoutL, newdata= renewdata, varvar= varb,
                    idvar= id, root = root, se.fit=TRUE, amat = amat)

#pout.amat$fit should be the same as file "herman"
herman
pout.amat$fit #they are the same.  Good.

#combine std.err with estimates, and then round
#to three decimal places
foo<- cbind(pout.amat$fit, pout.amat$se.fit)
rownames(foo)<- as.character(fred$familyID)
colnames(foo)<- c("Low Den Fitness", "SE")
round(foo, 3) 
L_estimates<- round(foo, 3)

L_estimates


#make rownames first column and combine all estimates
library(data.table)

L_estimates<- as.data.frame(L_estimates)
M_estimates<- as.data.frame(M_estimates)
H_estimates<- as.data.frame(H_estimates)

setDT(L_estimates, keep.rownames = TRUE)[]
setDT(M_estimates, keep.rownames = TRUE)[]
setDT(H_estimates, keep.rownames = TRUE)[]

#write.table(L_estimates, file="C:/Users/Mason Kulbaba/Dropbox/Rscripts/density/Results/L_fem_est.csv", sep=",", row.names = F, quote = F)
#write.table(M_estimates, file="C:/Users/Mason Kulbaba/Dropbox/Rscripts/density/Results/M_fem_est.csv", sep=",", row.names = F, quote = F)
#write.table(H_estimates, file="C:/Users/Mason Kulbaba/Dropbox/Rscripts/density/Results/H_fem_ests.csv", sep=",", row.names = F, quote = F)


###################################################################################
##### The following code conducs LM analsysi of above and below ground biomass#####
###################################################################################

#reload data
fin<- read.csv("data/aster.dat.csv")




library(emmeans)

# subset high and low Ne data

hi<- subset(fin, Gen=="HG")
lo<- subset(fin, Gen=="LG")

#test for equal variances

library(car)


leveneTest(mass.b ~ Treat, data = lo)


# make a function for standard error

stderr<- function(x) sd(x)/sqrt(length(x))

#Generate some summary statistics first

#High Ne
#Above ground biomass
aggregate(hi$mass.a, by=list(hi$Den), mean)
aggregate(hi$mass.a, by=list(hi$Den), stderr)

#Below ground biomass
aggregate(hi$mass.b, by=list(hi$Den), mean)
aggregate(hi$mass.b, by=list(hi$Den), stderr)

#Low Ne
#Above ground biomass
aggregate(lo$mass.a, by=list(lo$Den), mean)
aggregate(lo$mass.a, by=list(lo$Den), stderr)

#Below ground biomass
aggregate(lo$mass.b, by=list(lo$Den), mean)
aggregate(lo$mass.b, by=list(lo$Den), stderr)



###############################################
#total model - above ground biomass

#log transform biomass

##############################################################
#above ground biomass

f.lm<- lm(log(mass.a) ~ (Den) , data=fin)

f.lm2<- lm(log(mass.a) ~ (Den + Gen) , data=fin)

f.lm3<- lm(log(mass.a) ~ (Den + Gen + Den*Gen) , data=fin)

anova(f.lm, f.lm2, f.lm3)

summary(f.lm)
summary(f.lm2)
summary(f.lm3)

vec<- c("Gen")

emmeans(f.lm2, "Den", type='response', by=vec)

pairs(emmeans(f.lm2, "Den", type='response', by=vec))
test(emmeans(f.lm2, "Den", type='response', by=vec))

################################################################
#Now check for position effects based on three classes 

#First create classes based on plantID value
  #Two middle position = 1 &12

  #Abutting hexagons = 6 & 7

  #remaining outer plants = 2, 3, 4, 5, 8, 9, 10, 11

#subset data based on above criteria

#reload data
fin<- read.csv("C:/Users/Mason Kulbaba/Dropbox/git/density-gen-Ne/data/aster.dat.csv")


fin$pos<- 0

fin$pos[fin$plantID == 1 | fin$plantID==12]="mid"

fin$pos[fin$plantID == 6 | fin$plantID==7]="abs"

fin$pos[fin$plantID == 2 | fin$plantID==3| fin$plantID == 4 | fin$plantID==5|fin$plantID == 8| fin$plantID==9| fin$plantID == 10 | fin$plantID==11]="outer"

fin$pos<- as.factor(fin$pos)

#analysis with position
l<- lm(log(mass.a) ~ plotID, data=fin)

l2<- lm(log(mass.a) ~ plotID + Den, data=fin)

l3<- lm(log(mass.a) ~ plotID + Den +Gen, data=fin)

l4<- lm(log(mass.a) ~ plotID + Den + Gen +Den*Gen, data=fin)


lb<- lm(log(mass.b) ~ plotID, data=fin)

lb2<- lm(log(mass.b) ~ plotID + Den, data=fin)

lb3<- lm(log(mass.b) ~ plotID + Den + Gen, data=fin)

lb4<- lm(log(mass.b) ~ plotID + Den + Gen +Den*Gen, data=fin)

anova(lb, lb2, lb3, lb4)

#check for positional effect

lm1<- lm(log(mass.a) ~ (Den) , data=fin)

lm2<- lm(log(mass.a) ~ (Den + pos) , data=fin)

lm3<- lm(log(mass.a) ~ Den + pos + Den*pos , data=fin)

anova(lm1, lm2, lm3)

summary(lm1)
summary(lm2)
summary(lm3)

library(lmtest)

coeftest(lm3)


lm2<- lm(log(mass.a) ~ (Den + Gen) , data=fin)

lm3<- lm(log(mass.a) ~ (Den + Gen + Den*Gen) , data=fin)

lm4<- lm(log(mass.a) ~ (Den + Gen + Den*Gen + pos) , data=fin)


anova(lm, lm2, lm3, lm4)

summary(l)
summary(lm)
summary(lm2)
summary(lm3)
summary(lm4)

vec<- c("Gen")

emmeans(lm4, "Den", type='response', by=vec)

pairs(emmeans(lm4, "Den", type='response', by=vec))
test(emmeans(lm4, "Den", type='response', by=vec))


#redo with analysis of below ground biomass

lm<- lm(log(mass.b) ~ (Den) , data=fin)

lm2<- lm(log(mass.b) ~ (Den + Gen) , data=fin)

lm3<- lm(log(mass.b) ~ (Den + Gen + Den*Gen) , data=fin)

lm4<- lm(log(mass.b) ~ (Den + Gen + Den*Gen + pos) , data=fin)


anova(lm, lm2, lm3, lm4)

summary(lm)
summary(lm2)
summary(lm3)
summary(lm4)

vec<- c("Den")

emmeans(lm4, "Den", type='response', by=vec)

pairs(emmeans(lm4, "Den", type='response', by=vec))
test(emmeans(lm4, "Den", type='response', by=vec))



#below ground biomass

b.lm<- lm(log(mass.b) ~ (Den) , data=fin)

b.lm2<- lm(log(mass.b) ~ (Den + Gen) , data=fin)

b.lm3<- lm(log(mass.b) ~ (Den + Gen + Den*Gen) , data=fin)

anova(b.lm, b.lm2, b.lm3)

summary(b.lm)
summary(b.lm2)
summary(b.lm3)

emmeans(b.lm2, "Den", type='response', by=vec)

pairs(emmeans(b.lm2, "Den", type='response', by=vec))
test(emmeans(b.lm2, "Den", type='response', by=vec))


##########################################################################
##### The following code describes the LM analysis of aborted embryos#####
##########################################################################

#Load data
fin<- read.csv("data/aster.dat.csv")

#change class of factor variables
fin$Den<- as.factor(fin$Den)
fin$Gen<- as.factor(fin$Gen)
fin$plotID<- as.factor(fin$plotID)
fin$plantID<- as.factor(fin$plantID)
fin$familyID<- as.factor(fin$familyID)

#################################################################
#aboted embryo analysis with linear models
lm<- lm(log(aborted +1) ~ plotID , data=fin)

lm1<- lm(log(aborted +1) ~ plotID+Den , data=fin)

lm2<- lm(log(aborted +1) ~ plotID+Den + Gen, data=fin)

lm3<- lm(log(aborted +1) ~ plotID+Gen + Den + Den*Gen, data=fin)#interaction not significant       

anova(lm, lm1, lm2, lm3)

summary(lm)
summary(lm1)
summary(lm2)
summary(lm3)



emmeans(lm3, "Den", "Gen", type="response")
#plot(emmeans(lm3, "Den", "Gen", type="response"))
pairs(emmeans(lm3, "Den","Gen", type="response"))
test(emmeans(lm3, "Den","Gen", type="response"))

###########################
#Upon advice of Reviewer 2: any association between fitness (seeds set/seeds sired) and flower number?

# note, frt represents all flowers
lm1<- lm(log(seeds +1) ~ frt, data=fin)

lm2<- lm(log(seeds +1) ~ Den +frt, data=fin)

lm3<- lm(log(seeds +1) ~ Den +Gen+frt, data=fin)


anova(lm1, lm2, lm3)

summary(lm3)


library(MASS)

flw<- fin$frt

sd<- fin$seeds

flw1<- fitdistr(flw, densfun = "normal")
flw2<- fitdistr(flw, densfun = "negative binomial")
flw3<- fitdistr(flw, densfun = "Poisson")

AIC(flw1, flw2, flw3)


seed1<- fitdistr(sd, densfun = "normal")
seed2<- fitdistr(sd, densfun = "negative binomial")
seed3<- fitdistr(sd, densfun = "Poisson")

AIC(seed1, seed2, seed3)


gs<- glm(seeds ~ Den + Gen + flw, family=poisson, data=fin)

summary(gs)

#quick plots by den

high<- subset(fin, Den=="H")
med<- subset(fin, Den=="M")
low<- subset(fin, Den=="L")


plot(high$frt, high$seeds)

plot(med$frt, med$seeds)

plot(low$frt, low$seeds)


###########################################################################
# Analysis of proportion of aborted embryos


#Load data
fin<- read.csv("data/aster.dat.csv")

#first, make variable of sum of aborted ovules and seeds: fertilized ovules

fin$fert.ovules<- (fin$aborted + fin$seeds)

#name make column of proportion of aborted embryos

fin$prop.abort<- fin$aborted / fin$fert.ovules

#Change NA to zeros

#fin[is.na(fin)] <- 0


#try removing na
#fin<- fin[complete.cases(fin), ]



#plantID and plotID as factor

fin$plantID<- as.factor(fin$plantID)

fin$plotID<- as.factor(fin$plotID)

#make two column resoonse object for binomial

#These models yield probability of successfull seed set

resp<- cbind(fin$seeds, fin$fert.ovules)

resp2<- cbind(fin$seeds, fin$aborted)

t<- glm(resp2 ~ Den,family='binomial', data=fin)

t1<- glm(resp2 ~ Den + Gen,family='binomial', data=fin)

t2<- glm(resp2 ~ Den + Gen + (Den*Gen),family='binomial', data=fin)

t3<- glm(resp2 ~ plotID + Den + Gen + (Den*Gen),family='binomial', data=fin)


summary(t)
summary(t1)
summary(t2)
summary(t3)

AIC(t, t1, t2, t3)
anova(t, t1, t2, t3)

library(emmeans)


emmeans(t3,"Den","Gen", type="response")

pairs(emmeans(t3, "Gen", "Den", type='response'))
test(emmeans(t3,"Den", "Gen", type='response'))

# following models yield probability of embryo abortion

resp.abort<- cbind(fin$aborted, fin$seeds)

t<- glm(resp.abort ~ Gen,family='binomial', data=fin)

t1<- glm(resp.abort ~  Den + Gen,family='binomial', data=fin)

t2<- glm(resp.abort ~  Gen + Den + (Gen*Den) ,family='binomial', data=fin)

t3<- glm(resp.abort ~  Gen + Den + (Gen*Den) + plotID ,family='binomial', data=fin)


summary(t)
summary(t1)
summary(t2)
summary(t3)

AIC(t, t1, t2, t3)
anova(t, t1, t2, t3)

library(emmeans)


emmeans(t3,"Den","Gen", type="response")

pairs(emmeans(t3, "Gen", "Den", type='response'))
test(emmeans(t3,"Den", "Gen", type='response'))



