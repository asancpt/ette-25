#### Import, process, and export simulated data ##############
# The simulated data need to be processed to remove observations following the . . .
# . . . first simulated remedication time, if any (i.e. all observartions following . . .
# . . . an observation in which QUIT=1 are removed from the data set)
importDir <- "./models-clinical-trial/"
fileName <- paste(importDir, "analgesicSim.tab", sep="")
simtab <- read.table(fileName, header=T, skip=1)
#
# Determine number of simulated trials, and sample size
nSim <- sum(!duplicated(simtab$ISIM))
nSubj <- sum(!duplicated(simtab$ID))
# Add var to hold logical vector . . . indicating obs to be removed
simtab$REMOVE <- rep(F, nrow(simtab))

##### Mark observations to be removed
for(iSim in 1:nSim){
for(iSubj in 1:nSubj){
TF.subj <- simtab$ISIM==iSim & simtab$ID==iSubj
TF.quit <- simtab$QUIT[TF.subj]==1
if(any(TF.quit)){
Time.v <- simtab$TIME[TF.subj]
qTime <- min(Time.v[TF.quit])
simtab$REMOVE[TF.subj] <- simtab$TIME[TF.subj] > qTime
} # end-if
} # end-for nSubj
} # end-for iSim
# Remove marked observations
simtab <- simtab[!simtab$REMOVE,]
#Removerecordsforwhichtheobservationismissing(MDV==1forTIME==0)
TF.v <- simtab$MDV==1
simtab$PRLS[TF.v] <- rep(NA,sum(TF.v))
keepCols <- c( "ISIM", "ID", "TIME", "MDV", "PRLS", "QUIT", "PTIM", "DOSE",
"CL", "VC", "Q", "VP", "KA", "ALAG", "FBIO")
exportDir <- "./data/"
fileName <- paste(exportDir, "analgesicSim.csv", sep="")
z.export.csv(simtab[, keepCols], fileName, MissingVal=".",AddHash=T)

