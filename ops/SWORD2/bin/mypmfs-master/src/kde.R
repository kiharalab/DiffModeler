# Copyright (C) 2017 Guillaume Postic (guillaume.postic@upmc.fr)

args = commandArgs(trailingOnly=TRUE)
cmd_args <- commandArgs(TRUE)

if (length(args) < 7) {
  stop("Missing argument(s)
        Arguments:
        - args[1]    Distances file (*.dat): the interatomic distances observed for a given residue pair
        - args[2]    Bin width (Å)
        - args[3]    Upper limit of the distances (Å)
        - args[4]    Kernel
        - args[5]    Bandwidth
        - args[6]    Adjust
        - args[7]    Cut
        \n", call.=FALSE)
}

totalargs = length(args)

infile = args[1]
mydata = read.table(infile, h=F)
binwidth = as.numeric(args[2])
mymax = as.numeric(args[3])

mybreaks = seq(from=0, to=mymax, by=binwidth)
datahist = hist(mydata$V1, breaks=mybreaks, plot=F)

halfbinwidth = binwidth/2
mymax2 = mymax-halfbinwidth
mybins = seq(from=halfbinwidth, to=mymax2, by=binwidth)
mydensity = density(mydata$V1, kernel=args[4], bw=args[5], adjust=as.numeric(args[6]), cut=as.numeric(args[7]))
densfunction = approxfun(mydensity)
datakde = densfunction(mybins)

# Correction
totaldata = mymax/binwidth
for(i in seq(from=1, to=totaldata, by=1)){
    if (datahist$density[i] == 0){ # if hist() says 0 for this bin...
        datakde[i] = 0 # ...then it's 0 for this bin!
    }
}

if (totalargs < 8){ # i.e. Plot only if there is no extra argument
    # Plotting the estimated density
    outfile2 = gsub(".dat", "_dens.svg", infile)
    svg(outfile2)
    plot(mydensity, xlab="Distance (Å)", main=" ")
    invisible(dev.off())
}

cat(datakde, sep=" ")

if (totalargs < 8){
    # Plotting the histogram
    outfile1 = gsub(".dat", "_hist.svg", infile)
    svg(outfile1)
    plot(datahist, xlab="Distance (Å)", main=" ")
}
