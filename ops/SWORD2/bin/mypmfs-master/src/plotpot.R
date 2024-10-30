# Copyright (C) 2020 Guillaume Postic (guillaume.postic@u-paris.fr)

options(warn=-1) # because "lwd.ticks" works but returns a warning

args = commandArgs(trailingOnly=TRUE)

if (length(args) < 5) {
  stop("Missing argument(s)
        Arguments:
        - args[1]    The *.nrg file
        - args[2]    Bin width (Å)
        - args[3]    Upper limit of the distances (Å)
        - args[4]    Name of the residue pair
        - args[5]    Output dir
        \n", call.=FALSE)
}

mydata = read.table(args[1], h=F)

halfbin = as.numeric(args[2])/2
mylimit = as.numeric(args[3])-halfbin
outfile = paste(args[5], '/', args[4], '.svg', sep="")
svg(outfile)
plot(seq(from=halfbin, to=mylimit, by=as.numeric(args[2])), mydata$V1, type="l", xlab="Interatomic distance (Å)", ylab="Interaction score", col="darkblue", las=1, cex.lab=1.75, cex.axis=1.5, cex.main=1.5, cex.sub=1.5, tck=0.02, lwd = 4, lwd.ticks=3, xaxp = c(0, 100, 50))
abline(h=0, lwd=3, lty=2)
box(lwd=3)
dev.off()
