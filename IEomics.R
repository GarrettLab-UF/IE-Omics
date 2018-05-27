##remove all R objects from memory, this program can be memory intensive as it is dealing with huge datasets
rm(list = ls())

##User inputs
RTWindow = .3 # For a selected ion in current injection, exclude ion within this retention time window of that ion in following injections
noiseCount = 15 # Number of times precursor is selected before it is considered noise and removed across all retention times
MZWindow = .02 #leave the same, combines m/z values within this Da window

if("gWidgets" %in% rownames(installed.packages()) == FALSE) {install.packages("gWidgets")}
if("gWidgetstcltk" %in% rownames(installed.packages()) == FALSE) {install.packages("gWidgetstcltk")}
require(gWidgets)
require(gWidgetstcltk)
options(guiToolkit="tcltk") 

runNumber = ginput(message="What injections is this? starting at 1 (numeric value) \n(add 1 everytime you import a new MS/MS \nfor the same experiment)", title="Injection Number",icon="question") # MUST INCREMENT THIS EACH TIME YOU DO A RUN
runNumber = as.numeric(runNumber)
polarity = ginput(message="What is the polarity? \ninput: 'Positive' or 'Negative' without parenthesis\n(Must be spelled exactly as above \nwith first letter capitalized)", title="Polarity",icon="question") # "Positive"
RTMargin<-as.numeric(RTWindow)/2
userMaxRT = ginput(message="maximum retention time in chromatography", title="maximum retention time",icon="question") # maximum RT time
userMaxRT<-as.numeric(userMaxRT)
#import ms2
inputMS2 <- choose.files(caption="import .ms2 file with MS/MS data \nto be used to generate an exclusion list",multi=FALSE)
#export the exclusion list
outputCSVDirectory <- file.path(dirname(inputMS2), paste("ExclusionList_injection", runNumber, ".csv", sep = ""))

if(runNumber > 1){
  #import already made exclusion csv to combine with new ms2 data
  inputDirectoryCSV <- choose.files(caption="import previously generated exclusion list \n(.csv file) The new selected precursors will be appended \nto this exclusion list",multi=FALSE)
}

#header for matrix
exclusionMatrix = matrix(c("Mass [m/z]", "Formula [M]", "Formula type", "Species", "CS [z]", "Polarity", "Start [min]", "End [min]", "Comment"), nrow=1)
##Converting MS2 file into a bunch of characters
MSMS = scan(file = inputMS2, what = 'character')
MSMS = MSMS[19:length(MSMS)]#cuts out header

##Append rows with the respective column information
i = 1
while(i<length(MSMS)){ #could be changed to for loop later
  if(MSMS[i] == "S"){ #start of each block of MS2 file
    #ScurrentIndex = i
    MZindex = i + 3
    RTindex = i + 6
    minRT = (as.numeric(MSMS[RTindex]) - RTMargin)
    if(i > 1 & i < 700){print(paste("minRT: ", minRT, " at i: ", i))}
    if(minRT < 0){ minRT = 0 }
    maxRT = (as.numeric(MSMS[RTindex]) + RTMargin)
    if(maxRT > userMaxRT){ maxRT = userMaxRT}
    comment = paste("run",runNumber, sep="")
    b = matrix(c(MSMS[MZindex], '', '', '', '', polarity, minRT, maxRT, comment), nrow=1)
    exclusionMatrix = rbind(exclusionMatrix, b)
  }
  i = i + 1
}

## Moved this function down

# exclusionDF = as.data.frame(exclusionMatrix)
# for (h in length(MZvalues)){
#    tempUniqueMZwithinRange = subset(exclusionDF, as.numeric(exclusionDF[h, 1]) <= (MZvalues[h] + MZMargin) & as.numeric(exclusionDF[h, 1]) >= (MZvalues[h] - MZMargin))
# }

##Combine old exclusion with new exclusion list
if(runNumber > 1){
  exclusionCSV = read.csv(inputDirectoryCSV, sep=",", na.strings="NA", dec=".", strip.white=TRUE,header=FALSE)
  exclusionCSV = as.matrix(exclusionCSV)
  exclusionMatrix = exclusionMatrix[2:nrow(exclusionMatrix), ]
  exclusionMatrix = rbind(exclusionCSV, exclusionMatrix)
}

##CHANGED jpk Pulled out of IF statement above, the data cleaning needs to happen even the first time
header = exclusionMatrix[1,]   
col1 = as.numeric(exclusionMatrix[, 1])
col7 = as.numeric(exclusionMatrix[, 7])
exclusionMatrix = exclusionMatrix[order(col1, col7), ] #sort by MZ, then minRT
exclusionMatrix = exclusionMatrix[-nrow(exclusionMatrix),]#delete last row (adds the header to bottom when it sorts)
exclusionMatrix = rbind(header, exclusionMatrix)

##Code to replace mass to charges within the same mass error with the same mass value
z <- 2
while (z < nrow(exclusionMatrix)) {
  if(as.numeric(exclusionMatrix[(z+1),1])<(as.numeric(exclusionMatrix[z,1])+MZWindow)){
    exclusionMatrix[(z+1),1]<-exclusionMatrix[z,1]
  }
  z<-z+1
}

##CHANGED jpk (moved down here, you want this on the appended dataset, plus you can get errors later if the background of to similar m/z is replaced by a max RT (e. g. 0-18 become 0-3.4)) Data cleaning - filters background artifacts
filterComment = "background"
MZvalues = exclusionMatrix[2:nrow(exclusionMatrix),1]
uniqueVals = unique(MZvalues)
for (i in 1:length(uniqueVals)) {
  MZvalues = exclusionMatrix[2:nrow(exclusionMatrix),1]
  if(length(MZvalues[MZvalues == uniqueVals[i]]) > noiseCount){
    exclusionMatrix = exclusionMatrix[c(TRUE, (MZvalues != uniqueVals[i])), ]
    exclusionMatrix = rbind(exclusionMatrix, matrix(c(uniqueVals[i], '', '', '', '', polarity, '0', userMaxRT, filterComment), nrow=1))
  }
}

col1 = as.numeric(exclusionMatrix[, 1])
col7 = as.numeric(exclusionMatrix[, 7])
exclusionMatrix = exclusionMatrix[order(col1, col7), ] #sort by MZ, then minRT
exclusionMatrix = exclusionMatrix[-nrow(exclusionMatrix),]
exclusionMatrix = rbind(header, exclusionMatrix)
# 
##Data cleaning - compresses MZs based on overlapping RT for a more compact and readable exclusion List
MZvalues = as.numeric(exclusionMatrix[2:nrow(exclusionMatrix),1])
row = 3
for (j in unique(MZvalues)){
  loopSize = length(which(MZvalues == j))
  #print(loopSize)
  k = 1
  while(k < loopSize){
    if((as.numeric(exclusionMatrix[row, 7]) <= as.numeric(exclusionMatrix[(row-1), 8]))){ #RT Min
      exclusionMatrix[(row-1), 8] = as.numeric(exclusionMatrix [row, 8])
      exclusionMatrix = exclusionMatrix[-row, ]
      row = row - 1
    }
    k = k + 1
    row = row + 1
  }
  row = row + 1
}

#Necessary for fixing comments (e.g. sometimes the comment says background but the RT is not userMaxRT.
#                                    I spot-checked and confirmed that it was okay to do this, I wasn't 
#                                    ruining original data)
for(l in 1:nrow(exclusionMatrix)){ 
  if(exclusionMatrix[l,9] == "background"){
    exclusionMatrix[l,8] = userMaxRT
  }
}

write.table(exclusionMatrix, outputCSVDirectory, sep=",", col.names=FALSE, row.names=FALSE, quote=TRUE, na="NA")

