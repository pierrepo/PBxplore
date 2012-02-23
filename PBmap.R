#! /usr/bin/env Rscript

# 2012 - P. Poulain, A. G. de Brevern

#====================================================================
# load data
#====================================================================
args = commandArgs(TRUE)
filename = args[1]
lower = as.integer(args[2])
upper = as.integer(args[3])

data = read.table(filename, header = TRUE, row.names = 1)

if (length(colnames(data)) != 16) {
    cat("Bad format. 16 + 1 columns expected.\nBye bye\n")
    quit()
}

# lower / upper bound
if (is.na(lower)) lower = min(as.numeric(rownames(data)))
if (is.na(upper)) upper = max(as.numeric(rownames(data)))

# select data range
map = data[lower:upper,]

# define row and col names
xnames = rownames(map)
ynames = colnames(map)

# get total number of counts
total_count = sum(map[3,])
# convert 0 to NA (to have a white background)
#map[map==0] = NA
# normalize by the total number of counts
map = map / total_count


#====================================================================
# graphical parameters
#====================================================================
name="PBmap.png"
png(filename=name, width = log(length(xnames))*250, height = 800)
par(
    # default margins are: 5.1 4.1 4.1 2.1
    # extend bottom margin for text (+5 line)
    mar = c(5.1, 5.1, 4.1, 2.1),
    oma = c(2,0,0,0), # 2 lines for comments: 0 to 4
    lwd=3,            # line width
    bty = 'o',        # type of box around graphic
    font.lab = 2,     # axis label font (bold)
    font.axis = 2,    # axis font (bold)
    cex.lab=1.7,      # axis label width
    cex.axis=1.5      # axis width
)

# color gradient
# color goes from light yeallow to red
colorpal = rev(heat.colors(10))

# new color gradient
# color goes from dark blue to green/yellow to red
grad = matrix(nrow=848, ncol=3)
grad[1,1] = 20
grad[1,2] = 20
grad[1,3] = 232
for(i in 2:212){
grad[i,1] = grad[i-1,1]
grad[i,2] = grad[i-1,2]+1
grad[i,3] = grad[i-1,3]
}
for(i in 213:424){
grad[i,1] = grad[i-1,1]
grad[i,2] = grad[i-1,2]
grad[i,3] = grad[i-1,3]-1
}
for(i in 425:636){
grad[i,1] = grad[i-1,1]+1
grad[i,2] = grad[i-1,2]
grad[i,3] = grad[i-1,3]
}
for(i in 637:848){
grad[i,1] = grad[i-1,1]
grad[i,2] = grad[i-1,2]-1
grad[i,3] = grad[i-1,3]
}
colorpal = rgb(grad[,1]/255,grad[,2]/255,grad[,3]/255)

#====================================================================
# plot
#====================================================================

image(as.matrix(map), axes=FALSE, xlab="Residue number", ylab="PB", col=colorpal)
box()
axis(1, seq(0, 1, 1/(length(xnames)-1)), xnames)
axis(2, seq(0, 1, 1/(length(ynames)-1)), ynames, font = 4)

mtext(paste(name, Sys.Date()), side = 3, line = 3)

cat(paste("wrote", name, "\n"))
#====================================================================
# close output file
#====================================================================
void = dev.off()
# avoid the ouput of "null device"


