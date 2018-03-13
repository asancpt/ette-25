#APPENDIX 25.4 S-PLUS CODE CREATED BARPLOTS OF
#PAIN RELIEF SCORES VERSUS TIME BY DOSE
##################################################################
##### Plot Number of Subjects by Time and Pain Relief Score (conditioned on Dose)
##### summarize number ofsubjects with a givenlevel ofpain relief (by dose and time)
n0 <- tapply(simtab$PRLS==0, list(simtab$DOSE, simtab$TIME), sum, na.rm=T)
n1 <- tapply(simtab$PRLS==1, list(simtab$DOSE, simtab$TIME), sum, na.rm=T)
n2 <- tapply(simtab$PRLS==2, list(simtab$DOSE, simtab$TIME), sum, na.rm=T)
n3 <- tapply(simtab$PRLS==3, list(simtab$DOSE, simtab$TIME), sum, na.rm=T)
n4 <- tapply(simtab$PRLS==4, list(simtab$DOSE, simtab$TIME), sum, na.rm=T)
##### Combine summaries of pain relief scores into a single dataframe
n.all <- rbind(n0, n1, n2, n3, n4)
#
nPRLS <- as.data.frame(n.all)
names(nPRLS) <- dimnames(n1)[[2]]
#
nPRLS$DOSE <- rep(as.numeric(dimnames(n1)[[1]]), 5)
nPRLS$DOSE <- nPRLS$DOSE/1000 # Convert Dose to mg
#### Create stacked barplot of pain relief scores
#graphsheet(color.style="black and white")
par(mfrow=c(2,2))
par(cex=1)
#
# Get number of doses, and col that specifies DOSE
Dose.v <- unique(nPRLS$DOSE)
rmCol <- match("DOSE", names(nPRLS))
#
for (iDose in Dose.v){
TF.v <- nPRLS$DOSE==iDose
xvals <- barplot(as.matrix(nPRLS[TF.v,-c(1,rmCol)]), ylim=c(0,200),
col=seq(5,1), density=500)
text(xvals, -10, names(nPRLS)[-c(1,rmCol)])
mtext("Number of Subjects", side=2, line=2.5, cex=1)
mtext("Time [hr]", side=1, line=1)
mtext(paste("Dose =", iDose, "mg"))
}
#
#key(10,200, rectangles=list(size=5, col=seq(5,1), density=500, angle=20),
#text=c("0: none", "1: a little", "2: medium", "3: a lot", "4: complete"),
#title="Pain Relief", cex.title=1.1,
#border=F
#)
#export.graph("./plots/PRLS.barplot.wmf", ExportType="WMF")
