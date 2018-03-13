######################### Generate NONMEM data set ########################
##### Specify parameters for clinical trial
nArm <- 3
nSubjArm <- 200
dose.v <- c(0, 500, 1000) # ug
pdObsTime.v <- c(0, 15, 30, 45, 60, 90, 120)/60 # [hr]


##### Calculate clin trial design parameters
nSubj <- nArm*nSubjArm
subj.v <- 1:nSubj
nDose <- length(dose.v)
nPDObs <- length(pdObsTime.v)
# Vars needed for NONMEM
# ID TIME MDV PRLF=DV QUIT PTIM DOSE
# ID . . . NONMEM ID
# TIME . . . Dose or Observation Time [hr]
# PRLF . . . Pain relief score
# QUIT . . . Remedication indicator
# (0=subject stays in trial, 1=subject remedicates and quits trial)
# PTIM . . . Time of previous observation
# DOSE . . . Dose [ug]
##### Create PD observation records
pdObs.tmp <- data.frame(ID=subj.v,
TIME=rep(NA, nSubj),
MDV=rep(NA, nSubj),
PRLF=rep(NA, nSubj),
QUIT=rep(NA, nSubj),
PTIM=rep(NA, nSubj),
DOSE=rep(dose.v, each=nSubjArm))
##### Create PD observation records
pdObs <- pdObs.tmp
for (iObs in 2:nPDObs){
pdObs <- rbind(pdObs, pdObs.tmp)
}
#
pdObs$TIME <- rep(pdObsTime.v, each=nSubj)
pdObs$MDV <- rep(0, nrow(pdObs))
pdObs$PRLF <- rep(NA, nrow(pdObs))
pdObs$QUIT <- rep(NA, nrow(pdObs))
TF.v <- pdObs$TIME==0
pdObs$MDV[TF.v] <- rep(1,sum(TF.v))
pdObs <- pdObs[order(pdObs$ID, pdObs$TIME), ]
##### specify "previous time"
pTime <- pdObs$TIME[pdObs$ID==1]
pTime <- c(0, pTime[-length(pTime)])
pdObs$PTIM <- rep(pTime, nSubj)
pdObs <- pdObs[order(pdObs$ID, pdObs$TIME, -pdObs$MDV), ]

exportDir <- "./data/"
fileName <- paste(exportDir, "analgesicTemplate.csv", sep="")
##### Define function to export NONMEM datasets as CSV file
z.export.csv <- function(data.df, fileName, MissingVal=".", AddHash=F){
if(AddHash){
names(data.df)[1] <- paste("#", names(data.df)[1])
}
cat(names(data.df), sep=", ", file=fileName)
cat("\n",file=fileName, append=T)
write.table(data.df, file=fileName, sep = ",", na = MissingVal, col.names = FALSE, quote=FALSE, row.names = FALSE,
            append=T)
}
z.export.csv(pdObs, fileName, MissingVal=".", AddHash=T)
pdObs