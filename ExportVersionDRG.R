#This code is for the analysis of ED patients wirth seizure / not seizure status
#Email: davidrobert.grimes@tcd.ie for any queries, code presented as is... 

library(epiR)
library(readxl)

#CHANGE your location to whereever you store the analysis file: 
location <- "Z:/Dropbox/RStudioStuff/Epilepsy/AnalysisDataCleaned.xlsx"

epildata <- read_excel(location,sheet = "DxData")
eegdata <- read_excel(location,sheet = "EEGData")
#Note: EEGData appears bigger due to spaces and duplicates, but these are excluded in regression!

#First, we create DFs of existing AND preexisting epilepsy diagnoses 
epiexisting <- epildata[epildata[, 3] == 1, ]
newpresents <- epildata[epildata[, 3] == 0, ]

#For new presentations, how many were diagnosed as seizure disorder and as mimics? 
#Note that in Lines 17 - 33 we present "unclear" is not seizure disorder! We vary this next.. 
seiznew <- newpresents[grep("Epileptic Sz|Alcohol/Drug related Sz|First Sz|Symptomatic Sz|Seizure", newpresents$`Final Dx`, ignore.case = TRUE), ]
seizmimic <- newpresents[!grepl("Epileptic Sz|Alcohol/Drug related Sz|First Sz|Symptomatic Sz|Seizure", newpresents$`Final Dx`, ignore.case = TRUE), ]


#The true positive and true negative rate for seizure disorders in ED: True positives and negatives
truepos <- seiznew[grep("Epileptic Sz|Alcohol/Drug related Sz|First Sz|Symptomatic Sz|Seizure", seiznew$`Initial Dx`, ignore.case = TRUE), ]
trueneg <- seizmimic[!grepl("Epileptic Sz|Alcohol/Drug related Sz|First Sz|Symptomatic Sz|Seizure", seizmimic$`Initial Dx`, ignore.case = TRUE), ]

#How many WITH seizures were misdiagnosed at ED? how many WITHOUT seizure misdiagnosed?
falseneg <- seiznew[!grepl("Epileptic Sz|Alcohol/Drug related Sz|First Sz|Symptomatic Sz|Seizure", seiznew$`Initial Dx`, ignore.case = TRUE), ]
falsepos <- seizmimic[grep("Epileptic Sz|Alcohol/Drug related Sz|First Sz|Symptomatic Sz|Seizure", seizmimic$`Initial Dx`, ignore.case = TRUE), ]

#Run a sensitivity / specificity analysis for new presentations! 
dat.v01 <- c(nrow(truepos),nrow(seizmimic) - nrow(trueneg),nrow(seiznew) - nrow(truepos),nrow(trueneg))
rval.tes01 <- epi.tests(dat.v01, method = "exact", digits = 2, conf.level = 0.95)
print(rval.tes01)

#Now we run a similar analysis, with new presumption (worse case) - all unclears are true seizures. 
#now, what is we recode all unclears as seizures? 
seiznew2 <- newpresents[grep("Epileptic Sz|Alcohol/Drug related Sz|First Sz|Symptomatic Sz|Seizure|Unclear", newpresents$`Final Dx`, ignore.case = TRUE), ]
seizmimic2 <- newpresents[!grepl("Epileptic Sz|Alcohol/Drug related Sz|First Sz|Symptomatic Sz|Seizure|Unclear", newpresents$`Final Dx`, ignore.case = TRUE), ]


#true positives and true negatives! 
truepos2 <- seiznew2[grep("Epileptic Sz|Alcohol/Drug related Sz|First Sz|Symptomatic Sz|Seizure", seiznew2$`Initial Dx`, ignore.case = TRUE), ]
trueneg2 <- seizmimic2[!grepl("Epileptic Sz|Alcohol/Drug related Sz|First Sz|Symptomatic Sz|Seizure", seizmimic2$`Initial Dx`, ignore.case = TRUE), ]

#How many WITH seizures were misdiagnosed at ED? how many WITHOUT seizure misdiagnosed?
falseneg <- seiznew2[!grepl("Epileptic Sz|Alcohol/Drug related Sz|First Sz|Symptomatic Sz|Seizure", seiznew2$`Initial Dx`, ignore.case = TRUE), ]
falsepos <- seizmimic2[grep("Epileptic Sz|Alcohol/Drug related Sz|First Sz|Symptomatic Sz|Seizure", seizmimic2$`Initial Dx`, ignore.case = TRUE), ]

#New sensitivity / specificity analysis 
dat.v02 <- c(nrow(truepos2),nrow(seizmimic2) - nrow(trueneg2),nrow(seiznew2) - nrow(truepos2),nrow(trueneg2))
rval.tes02 <- epi.tests(dat.v02, method = "exact", digits = 2, conf.level = 0.95)
print(rval.tes02)


#Next Question: does Prior epilepsy diagnosis lead to more accurate outcomes? 
#now we add a chi2 table for ALL data! 
# Finding rows where col1 matches col2 and col3 is 1
matchingdata_hx <- subset(epildata, `Initial Dx` == `Final Dx` & `PriorEpilepsy` == 1)
nomatchingdata_hx <- subset(epildata, `Initial Dx` != `Final Dx` & `PriorEpilepsy` == 1)

matchingdata_nohx <- subset(epildata, `Initial Dx` == `Final Dx` & `PriorEpilepsy` == 0)
nomatchingdata_nohx <- subset(epildata, `Initial Dx` != `Final Dx` & `PriorEpilepsy` == 0)

# Create a 2x2 contingency table of the above:
datachi <- matrix(c(nrow(matchingdata_hx),nrow(matchingdata_nohx), nrow(nomatchingdata_hx), nrow(nomatchingdata_nohx)), nrow = 2, byrow = TRUE)
colnames(datachi) <- c("Epilepsy Hx", "No Hx")
rownames(datachi) <- c("Correct", "Incorrect")
datachi <- as.table(datachi)
print(datachi)

#Chi 2 test on above table: 
chi_squared_test <- chisq.test(datachi)
print(chi_squared_test)


#Exclude unclear cases..
#a version where unclear endings are removed! 
matchingdata_nohx_noclear <- matchingdata_nohx[matchingdata_nohx$`Final Dx` != "Unclear", ]
matchingdata_hx_noclear <- matchingdata_hx[matchingdata_hx$`Final Dx` != "Unclear", ]
nomatchingdata_hx_noclear <- nomatchingdata_hx[nomatchingdata_hx$`Final Dx` != "Unclear", ]
nomatchingdata_nohx_noclear <- nomatchingdata_nohx[nomatchingdata_nohx$`Final Dx` != "Unclear", ]

# Create a 2x2 contingency table
datachinoclear <- matrix(c(nrow(matchingdata_hx_noclear),nrow(matchingdata_nohx_noclear), nrow(nomatchingdata_hx_noclear), nrow(nomatchingdata_nohx_noclear)), nrow = 2, byrow = TRUE)
colnames(datachinoclear) <- c("Epilepsy Hx", "No Hx")
rownames(datachinoclear) <- c("Correct", "Incorrect")
datachinoclear <- as.table(datachinoclear)
print(datachinoclear)

#test again!
chi_squared_test_noclear <- chisq.test(datachinoclear)
print(chi_squared_test_noclear)


#EEG Analysis! 

#Model disagreement as a logistic regression with just "had EEG"
disagreement_model <- glm(Disagree ~ EEG, data = eegdata, family = binomial)
summary(disagreement_model)

#Model disagreement as a logistic regression with just "had EEG" and "Abnormal" predictors! 
disagreement_model2 <- glm(Disagree ~ EEG + Abnormal, data = eegdata, family = binomial)
summary(disagreement_model2)

#In patients who had an EEG only, was abnormal result associated with discordance? 
eegonly <- eegdata[eegdata$EEG == 1 & !is.na(eegdata$EEG), ]
disagree_eeg<- glm(Disagree ~ Abnormal, data = eegonly, family = binomial)
summary(disagree_eeg)
