#########################################################################################
#                                                                                       #
# The following code performs fixed-effects aster analyses on data examining the effects#
# of density and effective genetic population size (Ne) on male (seeds sired) fitness.  #
#                                                                                       #
#                                                                                       #    
#                                                                                       #
#                                                                                       #
# For any questions, please contact Mason Kulbaba: mason.kulbaba@gmail.com              #
#                                                                                       #
#                                                                                       #
#########################################################################################



#Load data
fin<- read.csv("C:/Users/Mason Kulbaba/Dropbox/git/density-Ne/data/aster.sire.dat.csv")

#change class of factor variables
fin$Den<- as.factor(fin$Den)
fin$Gen<- as.factor(fin$Gen)
fin$plotID<- as.factor(fin$plotID)
fin$plantID<- as.factor(fin$plantID)
fin$familyID<- as.factor(fin$familyID)

#load aster
library(aster)

#set graphical model for male fitness: survive to flower ->  total flowers -> number of seed sired

#total flower number was determined by persistent pedicels and all fruit fruit prodiced,
# that is represented by variable "frt"


#set vars
vars<- c( "flw", "frt", "sires")

#reshape data so that all response variables are located in a single vector in a new data
#set called "redata"
redata <- reshape(fin, varying = list(vars), direction = "long",timevar = "varb", times = as.factor(vars), v.names = "resp")


#Designation of fitness variable
fit <- grepl("sires", as.character(redata$varb))
fit<- as.numeric(fit)

#add fit to redata
redata$fit <- fit

#check
with(redata, sort(unique(as.character(varb)[fit == 0])))
with(redata, sort(unique(as.character(varb)[fit == 1])))

#add a variable "root" to redata, where value is 2
redata<- data.frame(redata, root=2)

#make familyID a factor

redata$familyID<- as.factor(redata$familyID)

#set graph. model and family for each node
pred<- c(0,1,2)
fam<- c(1,2,2)

#check family 
sapply(fam.default(), as.character)[fam]

#analysis with only fintess data
aout<- aster(resp~varb, pred, fam, varb, id, root, data=redata)

#add family ID
aout2<- aster(resp~varb +fit:familyID, pred, fam, varb, id, root, data=redata)

#add density
aout3<- aster(resp~varb +fit:Den, pred, fam, varb, id, root, data=redata)

summary(aout)
summary(aout2)
summary(aout3)

anova(aout, aout2)# family not significant

anova(aout, aout3)# density is VERY not significant


###########################################################################
#####divide into Density treatments and drop unused levels of familyID#####
###########################################################################

#Note: to compare with relativized female fitness, still want male fitness estimates
# for each density treatment of High Ne plots

redataL<-subset(redata, Den=="L")

redataL$familyID<- droplevels(redataL$familyID)

redataM<- subset(redata, Den=="M")

redataM$familyID<- droplevels(redataM$familyID)

redataH<- subset(redata, Den=="H")

redataH$familyID<- droplevels(redataH$familyID)



##############################
#####Estimate Low Density#####
##############################

#run aster
aoutL<- aster(resp~varb + fit:(familyID), pred, fam, varb, id, root, data=redataL)

summary(aout, show.graph = T)

#generate MLE of saturated model mean value parameter vector: mu
pout<- predict.aster(aoutL, se.fit=TRUE)


# make design matrix
fred <- data.frame(familyID=levels(redataL$familyID), flw=1, frt=1, sires=1,root = 2)

#reshape the "made up data" just as the actual data
renewdata <- reshape(fred, varying = list(vars),
                     direction = "long", timevar = "varb",
                     times = as.factor(vars), v.names = "resp")

#make character string from "varb" of renewdata,
#without actual values (e.g., the layers of varb in renewdata)
layer<- gsub("[0-9]", "", as.character(renewdata$varb))

#add layer to renewdata
renewdata<- data.frame(renewdata, layer= layer)

#add seeds dired in new layer col of renewdata as numeric, called fit
fit<- as.numeric(layer=="sires")

#add fit to renewdata
renewdata<- data.frame(renewdata, fit = fit)

#rerun prediction of aout, with "made up" renewdata
pout<- predict(aoutL, newdata= renewdata, varvar= varb,
               idvar = id, root = root, se.fit = TRUE)

sapply(pout, class)
sapply(pout, length)

#put the parameter estimates into a matrix with individuals in rows
#and nodes along columns

nnode<- length(vars)
sally<- matrix(pout$fit, ncol = nnode)
dim(sally)# makes 8 x 3 matrix: 8 families by 3 nodes

#name the rows (by familyID Treat) and columns (as nodes)
rownames(sally)<- unique(as.character(renewdata$familyID))
colnames(sally)<- unique(as.character(renewdata$varb))

#view matrix 
round(sally, 3)


#use just totalseeds as predicted (expected) fitneses
herman<- sally[,grepl("sires", colnames(sally))]

#generate sandard error for predictions

nFam<- nrow(fred)
nnode<- length(vars)
amat<- array(0, c(nFam, nnode, nFam))
dim(amat)# makes an 8 x 3 x 8 matrix

#only want means for k'th individual that contributes to expected
#fitness, and want to add only 'sired' entries


foo<- grepl("sires", vars)
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

#################################
#####Redo for Medium Density#####
#################################

aoutM<- aster(resp~varb + fit:(familyID), pred, fam, varb, id, root, data=redataM)

summary(aoutM, show.graph=T)

#generate MLE of saturated model mean value parameter vector: mu
pout<- predict.aster(aoutM, se.fit=TRUE)

# make design matrix
fred <- data.frame(familyID=levels(redataM$familyID), flw=1, frt=1, sires=1,root = 2)

#reshape the "made up data" just as the actual data
renewdata <- reshape(fred, varying = list(vars),
                     direction = "long", timevar = "varb",
                     times = as.factor(vars), v.names = "resp")

#make character string from "varb" of renewdata,
#without actual values (e.g., the layers of varb in renewdata)
layer<- gsub("[0-9]", "", as.character(renewdata$varb))

#add layer to renewdata
renewdata<- data.frame(renewdata, layer= layer)

#add seeds sired in new layer col of renewdata as numeric, called fit
fit<- as.numeric(layer=="sires")

#add fit to renewdata
renewdata<- data.frame(renewdata, fit = fit)

#rerun prediction of aout, with "made up" renewdata
pout<- predict(aoutM, newdata= renewdata, varvar= varb,
               idvar = id, root = root, se.fit = TRUE)

sapply(pout, class)
sapply(pout, length)


#put the parameter estimates into a matrix with individuals in rows
#and nodes along columns

nnode<- length(vars)
sally<- matrix(pout$fit, ncol = nnode)
dim(sally)# makes 9 x 3 matrix: 9 families by 3 nodes

#name the rows (by family ID) and columns (as nodes)
rownames(sally)<- unique(as.character(renewdata$familyID))
colnames(sally)<- unique(as.character(renewdata$varb))

#view matrix 
round(sally, 3)

#use just totalseeds as predicted (expected) fitneses
herman<- sally[,grepl("sires", colnames(sally))]

#generate standard errors for the predicted values
nFam<- nrow(fred)
nnode<- length(vars)
amat<- array(0, c(nFam, nnode, nFam))
dim(amat)# makes an 9 x 3 x 9 matrix

#only want means for k'th individual that contributes to expected
#fitness, and want to add only total.pods entries


foo<- grepl("sires", vars)
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
colnames(foo)<- c("Medium Mean Fitness", "SE")
round(foo, 3) 
M_estimates<- round(foo, 3)
M_estimates


###############################
#####Redo for High Density#####
###############################

#aster analysis

aoutH<- aster(resp~varb + fit:(familyID), pred, fam, varb, id, root, data=redataH)

summary(aoutH, show.graph=T)


#generate MLE of saturated model mean value parameter vector: mu
pout<- predict.aster(aoutH, se.fit=TRUE)


#make design matrix
fred <- data.frame(familyID=levels(redataH$familyID), flw=1, frt=1, sires=1,root = 2)

#reshape the "made up data" just as the actual data
renewdata <- reshape(fred, varying = list(vars),
                     direction = "long", timevar = "varb",
                     times = as.factor(vars), v.names = "resp")

#make character string from "varb" of renewdata,
#without actual values (e.g., the layers of varb in renewdata)
layer<- gsub("[0-9]", "", as.character(renewdata$varb))

#add layer to renewdata
renewdata<- data.frame(renewdata, layer= layer)

#add seeds sired in new layer col of renewdata as numeric, called fit
fit<- as.numeric(layer=="sires")

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
dim(sally)# makes 10 x 3 matrix: 10 families by 3 nodes

#name the rows (by family ID) and columns (as nodes)
rownames(sally)<- unique(as.character(renewdata$familyID))
colnames(sally)<- unique(as.character(renewdata$varb))

#view matrix 
round(sally, 3)


#use just totalseeds as predicted (expected) fitneses
herman<- sally[,grepl("sires", colnames(sally))]

#generate standard errors for these predicted values
nFam<- nrow(fred)
nnode<- length(vars)
amat<- array(0, c(nFam, nnode, nFam))
dim(amat)# makes an 10 x 3 x 10 matrix

#only want means for k'th individual that contributes to expected
#fitness, and want to add only total.pods entries


foo<- grepl("sires", vars)
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
colnames(foo)<- c("High Mean Fitness", "SE")
round(foo, 3) 
H_estimates<- round(foo, 3)
H_estimates


#make rownames first column and combine all estimates
library(data.table)

L_estimates<- as.data.frame(L_estimates)
M_estimates<- as.data.frame(M_estimates)
H_estimates<- as.data.frame(H_estimates)

setDT(L_estimates, keep.rownames = TRUE)[]
setDT(M_estimates, keep.rownames = TRUE)[]
setDT(H_estimates, keep.rownames = TRUE)[]