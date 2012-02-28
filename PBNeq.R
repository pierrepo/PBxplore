#! /usr/bin/env Rscript

# 2012 - P. Poulain, A. G. de Brevern

#====================================================================
# load data
#====================================================================
args = commandArgs(TRUE)
filename = args[1]
lower = as.integer(args[2])
upper = as.integer(args[3])

neq = read.table(filename, header = TRUE)

if (length(colnames(neq)) != 2) {
    cat("Bad format. Two columns expected.\nBye bye\n")
    quit()
}

# lower / upper bound
if (is.na(lower)) lower = min(neq[,1])
if (is.na(upper)) upper = max(neq[,1])

# select data range
neq = neq[lower:upper,]

#====================================================================
# graphical parameters
#====================================================================
name="PBNeq.png"
png(filename=name, width = 1600, height = 1200)
par(
    # default margins are: 5.1 4.1 4.1 2.1
    # extend bottom margin for text (+5 line)
    mar = c(5.1, 5.1, 4.1, 2.1),
    oma = c(2,0,0,0), # 2 lines for comments: 0 to 4
    lwd=3,            # line width
    bty = 'o',        # type of box around graphic
    font.lab = 2,     # axis label font (bold)
    font.axis = 2,    # axis font (bold)
    cex.lab=2.5,      # axis label width
    cex.axis=2.0      # axis width
)

#====================================================================
# plot Neq 
#====================================================================
plot(neq, type = 'l', 
    xlab = 'Residue number', ylab = 'Neq', 
    xlim=c(lower, upper), ylim=c(1,max(round(neq[,2]))+2))

#====================================================================
# close output file
#====================================================================
cat(paste("wrote", name, "\n"))
void = dev.off()
# avoid the ouput of "null device"

