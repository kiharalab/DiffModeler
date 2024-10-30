# Copyright (C) 2017 Guillaume Postic (guillaume.postic@upmc.fr)

args = commandArgs(trailingOnly=TRUE)

if (length(args) < 3) {
  stop("Missing argument(s)
        Arguments:
        - args[1]    The *.spline file
        - args[2]    Name of the residue pair
        - args[3]    Output dir
        \n", call.=FALSE)
}

mydata = read.table(args[1], h=F)

outfile = paste(args[3], '/', args[2], '.svg', sep="")
svg(outfile)
plot(mydata, type="l", xlab="Distance (Å)", ylab="Pairwise score")
abline(h=0, col="red", lwd=2, lty=2)
dev.off()
