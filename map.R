#! /usr/bin/env Rscript

# 2012 - P. Poulain, A. G. de Brevern

#====================================================================
# graphical parameters
#====================================================================
# name of output file:
name = "show.pdf"

# define graphical parameters
definePar = function() {
par(
    # default margins are: 5.1 4.1 4.1 2.1
    # extend bottom margin for text (+5 line)
    mar = c(5.1, 5.1, 4.1, 2.1),
    oma = c(2,0,0,0), # 2 lines for comments: 0 to 4
    lwd=3,          # line width
    bty = 'o',      # type of box around graphic
    font.lab = 2,   # axis label font (bold)
    font.axis = 2,  # axis font (bold)
    cex.lab=1.7,    # axis label width
    cex.axis=1.5    # axis width
)}

# print filename and date on top
printNameDate = function() {
mtext( paste(name, Sys.Date()), side = 3, line = 3)
}

#====================================================================
# load data
#====================================================================
args = commandArgs(TRUE)
filename = args[1]
data = read.table(filename, header = TRUE, row.names = 1)

# select data range
#map = data[27:37,]
map = data

# define row and col names
xnames = rownames(map)
ynames = colnames(map)

total = sum(map[3,])
# convert 0 to NA
map[map==0] = NA
# normalize by the number of points
map = map / total

#====================================================================
# start output file
#====================================================================
png(filename="map.png", width = log(length(xnames))*250, height = 800)
definePar()

colorpal = rev(heat.colors(10))
#====================================================================
# plot
#====================================================================

image(as.matrix(map), xaxt="n", yaxt="n", xlab="Residue number", ylab="PB", col=colorpal)
box()
axis(1, seq(0, 1, 1/(length(xnames)-1)), xnames)
axis(2, seq(0, 1, 1/(length(ynames)-1)), ynames, font = 4)



#====================================================================
# close output file
#====================================================================
toto = dev.off()
# avoid the ouput of "null device"


