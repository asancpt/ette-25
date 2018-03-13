
#APPENDIX 25.6 S-PLUS CODE CREATED PLOTS OF
#PROBABILITY/PROPORTION PAIN RELIEF VERSUS
#TIME BY PAIN RELIEF AND DOSE
##################################################################
#### Import PK simulated PK data ##############
importDir <- "./models-pain/"
fileName <- paste(importDir, "analgesicEst.tab", sep="")
esttab <- read.table(fileName, header=T, skip=1)
#
#### Calculate Proportion of subjects with a given level (or
# . . . greater) Pain Relief . . . by Dose and Time
prop1<- tapply(esttab$PRLF>=1, list(esttab$TIME,esttab$DOSE),sum,na.rm=T)/
tapply(!is.na(esttab$PRLF), list(esttab$TIME,esttab$DOSE),sum,na.rm=T)
prop2<- tapply(esttab$PRLF>=2, list(esttab$TIME,esttab$DOSE),sum,na.rm=T)/
tapply(!is.na(esttab$PRLF), list(esttab$TIME,esttab$DOSE),sum,na.rm=T)
prop3<- tapply(esttab$PRLF>=3, list(esttab$TIME,esttab$DOSE),sum,na.rm=T)/
tapply(!is.na(esttab$PRLF), list(esttab$TIME,esttab$DOSE),sum,na.rm=T)
prop4<- tapply(esttab$PRLF>=4, list(esttab$TIME,esttab$DOSE),sum,na.rm=T)/
tapply(!is.na(esttab$PRLF), list(esttab$TIME,esttab$DOSE),sum,na.rm=T )
#
tmp <- data.frame(rbind(prop1, prop2, prop3, prop4))
dose.c <- sort(unique(esttab$DOSE))
names(tmp) <- as.character(dose.c)
#
prop.c <- c("prop1", "prop2", "prop3", "prop4")
#
tmp$PRLF <- rep(prop.c, each=nrow(prop1))
tmp$TIME <- rep(as.numeric(dimnames(prop1)[[1]]), length(prlf.c))
#
est.prop <- data.frame(make.groups("0"=tmp[,1], "0.5"=tmp[,2], "1"=tmp[,3]))
names(est.prop) <- c("CPROB", "DOSE")
est.prop$DOSE.F <- paste(as.character(est.prop$DOSE), "mg")
est.prop$TIME <- rep(tmp$TIME, length(dose.c))
est.prop$PRLF <- rep(tmp$PRLF, length(dose.c))
####CalculateProbabilityofagivenlevel(orgreater)PainReliefbyDoseandTime
P1 <- tapply(esttab$P1, list(esttab$TIME, esttab$DOSE), mean, na.rm=T)
P2 <- tapply(esttab$P2, list(esttab$TIME, esttab$DOSE), mean, na.rm=T)
P3 <- tapply(esttab$P3, list(esttab$TIME, esttab$DOSE), mean, na.rm=T)
P4 <- tapply(esttab$P4, list(esttab$TIME, esttab$DOSE), mean, na.rm=T)
tmp <- data.frame(rbind(P1, P2, P3, P4))
dose.c <- sort(unique(esttab$DOSE))
names(tmp) <- as.character(dose.c)
prlf.c <- c("P1", "P2", "P3", "P4")
tmp$PRLF <- rep(prlf.c, each=nrow(P1))
tmp$TIME <- rep(as.numeric(dimnames(P1)[[1]]), length(prlf.c))
est.cProb <- data.frame(make.groups("0"=tmp[,1], "0.5"=tmp[,2], "1"=tmp[,3]))
names(est.cProb) <- c("CPROB", "DOSE")
est.cProb$DOSE.F <- paste(as.character(est.cProb$DOSE), "mg")
est.cProb$TIME <- rep(tmp$TIME, length(dose.c))
est.cProb$PRLF <- rep(tmp$PRLF, length(dose.c))
est.prop.cProb <- rbind(est.prop, est.cProb)
################ Plot Probability/Proportion PRLF vs Time ###################
trellis.device(color=F)
superpose.line.lst <- trellis.par.get("superpose.line")
superpose.line.lst$col <- rep(3:6, 2)
superpose.line.lst$col <- rep(1, 8)
superpose.line.lst$lty <- rep(1:4, 2)
superpose.line.lst$lwd <- rep(3, 8)
trellis.par.set("superpose.line", superpose.line.lst)
superpose.symbol.lst <- trellis.par.get("superpose.symbol")
superpose.symbol.lst$cex <- rep(1, 8)
superpose.symbol.lst$col <- rep(3:6, 2)
superpose.symbol.lst$col <- rep(1, 8)
superpose.symbol.lst$font <- rep(1, 8)
superpose.symbol.lst$pch <- rep(c("1", "2", "3", "4"), 2)
trellis.par.set("superpose.symbol", superpose.symbol.lst)
tmp.plt <- xyplot(CPROB ~ TIME|DOSE.F, data=est.prop.cProb, groups=PRLF,
type=rep(c("l","p"), each=4), as.table=T, layout=c(3,1),
scales=list(cex=1),
ylab=list("Probability/Proportion Pain Relief >= m", cex=1.2),
xlab=list("Time [hr]", cex=1.2),
par.strip.text=list(cex=0.8),
panel=panel.superpose)
key.lst <- list(text=list(c("m=1", "m=2", "m=3", "m=4"), cex=1, adj=1),
lines=list(type="o",
col=superpose.symbol.lst$col[1:4],
lty=superpose.line.lst$lty[1:4],
lwd=superpose.line.lst$lwd[1:4],
pch=superpose.symbol.lst$pch[1:4],
cex=superpose.symbol.lst$cex[1:4]),
# x=0.5, y=0.5, corner=c(0,0.5),
space="top", columns=4, between.columns=2,
border=1, transparent=F,
title="Pain Relief", cex.title=1,
border=1,
between=1)
# update(tmp.plt, key=key.lst)
update(tmp.plt)
export.graph("./plots/ObsPred.PRLF.TIME.wmf", ExportType="WMF")

