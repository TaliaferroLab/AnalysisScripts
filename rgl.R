require(rgl)
require(RColorBrewer)
#####Prepare data#####
trips_pulldown <- read.table('~/Documents/MIT/AncientExonShape/MutatedOligos/Broadrun/Oligo1/Oligo1_triplemutations_InputvsPulldown.txt', header=T)
trips_control <- read.table('~/Documents/MIT/AncientExonShape/MutatedOligos/Broadrun/Oligo1/Oligo1_triplemutations_InputvsNoProt.txt', header=T)
trips <- merge(trips_pulldown, trips_control, by = c('pos1', 'pos2', 'pos3'), suffixes = c('.pulldown', '.control'))
#Remove duplicates and convert to log2
#trips <- trips[trips$pos1 > trips$pos2 & trips$pos2 > trips$pos3,]
trips <- cbind(trips, log2(trips$relativerate.pulldown), log2(trips$relativerate.control))
colnames(trips)[14:15] <- c('log2relrate.pulldown', 'log2relrate.control')
#Only trips with at least 2000 counts in each file
trips <- trips[trips$fileAcount.pulldown >= 2000 & trips$fileBcount.pulldown>= 2000 & trips$fileAcount.control >= 2000 & trips$fileBcount.control >= 2000,]
#Take delta log2relrate
trips <- cbind(trips, trips$log2relrate.pulldown - trips$log2relrate.control)
colnames(trips)[16] <- 'deltalog2relrate'
#Take top n most enriched and depleted trips
trips_enriched <- head(trips[with(trips, order(-deltalog2relrate)),], 300)
trips_depleted <- tail(trips[with(trips, order(-deltalog2relrate)),], 300)
trips <- rbind(trips_enriched, trips_depleted)
rm(trips_enriched)
rm(trips_depleted)
minenrich <- min(trips$deltalog2relrate)
maxenrich <- max(trips$deltalog2relrate)
normfactor <- max(maxenrich, abs(minenrich))
normalizedenrich <- function(x) {
  enrich <- x[16]
  if(enrich > 0){
    norm <- enrich / normfactor
  } else if ( enrich < 0) {
    norm <- (enrich / normfactor)
  }
}

trips <- cbind(trips, apply(trips, 1, normalizedenrich))
colnames(trips)[ncol(trips)] <- 'deltalog2relrate.normalized'

#Sort by normalized enrichment so that color matches enrichment
trips <- trips[with(trips, order(-deltalog2relrate.normalized)),]



####SPHERES####
x <- rep(1:3, each=9)
y <- rep(rep(1:3, each=3), 3)
z <- rep(1:3, 9)
colours <- terrain.colors(27)
plot3d(x,y,z,col=colours, type="s", size=1, alpha = rate*0.1, box = FALSE, axes = TRUE)
plot3d(trips$pos1, trips$pos2, trips$pos3,col=colours, type="s", size=abs(trips$deltalog2relrate.normalized), alpha = 1, box = FALSE, axes = TRUE)


###CUBES#####
lim <- function(x){c(0, max(abs(x))) * 1}

rgl_init()
x <- trips$pos1
y <- trips$pos2
z <- trips$pos3
deltalog2relrate.normalized <- trips$deltalog2relrate.normalized

###COLORS####
#colors <- colorRampPalette(brewer.pal(11, 'RdBu'))
#colors <- colorRampPalette(rev(c('#0000FF','#D8D8FF','white','#FED6AE','#FF4000')))
myColorRamp <- function(colors, values) {
  v <- (values - min(values)) / diff(range(values))
  x <- colorRamp(colors, space = 'rgb')(v)
  rgb(x[,1], x[,2], x[,3], maxColorValue = 255)
}
colors <- myColorRamp(c('#0099ff', '#ff9933'), deltalog2relrate.normalized) #Blue to Orange


###PLOT###
shapelist3d(cube3d(alpha = abs(1)), x,y,z, size=abs(deltalog2relrate.normalized**2), col=colors, alpha = abs(1))  #FIX THESE COLORS
rgl.lines(c(0, 150), c(0, 0), c(0, 0), color = "darkgrey")
rgl.lines(c(0, 0), c(0, 150), c(0, 0), color = "red")
rgl.lines(c(0, 0), c(0, 0), c(0, 150), color = "darkgreen")
axis3d('x', pos=c( NA, 0, 0 ), col = "darkgrey")
axis3d('y', pos=c( 0, NA, 0 ), col = "red")
axis3d('z', pos=c( 0, 0, NA ), col = "darkgreen")
title3d(xlab='Nucleotide 1', col = 'darkgrey')
title3d(ylab='Nucleotide 2', col = 'red', angle = 90)
title3d(zlab='Nucleotide 3', col = 'darkgreen')

#Make lines that are 4 units long in x, y, and z from the 3 most enriched and depleted trips
#Put labels on the ends of those lines
#Because we are now allowing all trips, not just x > y > z, the top 3 unique trips are lines 1, 7, and 13 (and last, last - 6, last - 12)
rgl.lines(x = c(c(trips[1,1], trips[1,1] + 4), c(trips[7,1], trips[7,1] + 4), c(trips[13,1], trips[13,1] + 4)), 
          y = c(c(trips[1,2], trips[1,2] + 4), c(trips[7,2], trips[7,2] + 4), c(trips[13,2], trips[13,2] + 4)), 
          z = c(c(trips[1,3], trips[1,3] + 4), c(trips[7,3], trips[7,3] + 4), c(trips[13,3], trips[13,3] + 4)), col = '#FF4000')
rgl.lines(x = c(c(trips[nrow(trips),1], trips[nrow(trips),1] + 4), c(trips[nrow(trips) - 6,1], trips[nrow(trips) - 6,1] + 4), c(trips[nrow(trips) - 12,1], trips[nrow(trips) - 12,1] + 4)), 
          y = c(c(trips[nrow(trips),2], trips[nrow(trips),2] + 4), c(trips[nrow(trips) - 6,2], trips[nrow(trips) - 6,2] + 4), c(trips[nrow(trips) - 12,2], trips[nrow(trips) - 12,2] + 4)), 
          z = c(c(trips[nrow(trips),3], trips[nrow(trips),3] + 4), c(trips[nrow(trips) - 6,3], trips[nrow(trips) - 6,3] + 4), c(trips[nrow(trips) - 12,3], trips[nrow(trips) - 12,3] + 4)), col = '#0000FF')
rgl.texts(x = trips[1,1] + 4, y = trips[1,2] + 4, z = trips[1,3] + 4, paste('(', trips[1,1], ',', trips[1,2], ',', trips[1,3], ')', sep = ''), col = '#FF4000', adj = c(0,0), cex = 0.8)
rgl.texts(x = trips[7,1] + 4, y = trips[7,2] + 4, z = trips[7,3] + 4, paste('(', trips[7,1], ',', trips[7,2], ',', trips[7,3], ')', sep = ''), col = '#FF4000', adj = c(0,0), cex = 0.8)
rgl.texts(x = trips[13,1] + 4, y = trips[13,2] + 4, z = trips[13,3] + 4, paste('(', trips[13,1], ',', trips[13,2], ',', trips[13,3], ')', sep = ''), col = '#FF4000', adj = c(0,0), cex = 0.8)
rgl.texts(x = trips[nrow(trips),1] + 4, y = trips[nrow(trips),2] + 4, z = trips[nrow(trips),3] + 4, 
          paste('(', trips[nrow(trips),1], ',', trips[nrow(trips),2], ',', trips[nrow(trips),3], ')', sep = ''), col = '#0000FF', adj = c(0,0), cex = 0.8)
rgl.texts(x = trips[nrow(trips) - 6,1] + 4, y = trips[nrow(trips) - 6,2] + 4, z = trips[nrow(trips) - 6,3] + 4, 
          paste('(', trips[nrow(trips) - 6,1], ',', trips[nrow(trips) - 6,2], ',', trips[nrow(trips) - 6,3], ')', sep = ''), col = '#0000FF', adj = c(0,0), cex = 0.8)
rgl.texts(x = trips[nrow(trips) - 12,1] + 4, y = trips[nrow(trips) - 12,2] + 4, z = trips[nrow(trips) - 12,3] + 4, 
          paste('(', trips[nrow(trips) - 12,1], ',', trips[nrow(trips) - 12,2], ',', trips[nrow(trips) - 12,3], ')', sep = ''), col = '#0000FF', adj = c(0,0), cex = 0.8)

#rgl.bbox(color = c('#E6E6E6', '#585858'), emission = 'gray', specular = 'gray', shininess = 5, alpha = 0.8, xlen = 0, ylen = 0, zlen = 0)
aspect3d(1,1,1)



#Make two movies individually, then combine them using ImageMagick
movie3d(spin3d(axis = c(0,3,0), rpm = 2), duration = 23.6, dir = getwd(), clean = FALSE, convert = FALSE, movie = 'A', type = 'gif')

#Draw lines through common coordinates (this will have to be manually changed)
#4 common pairs: (39, 37, z), (x, 22, 19), (39, y, 19), (39, 30, z), if letter, draw from 0 to 143
rgl.lines(x = c(c(62, 62), c(0, 143), c(62, 62), c(60, 60), c(0, 143), c(60, 60)),
          y = c(c(60, 60), c(62, 62), c(0, 143), c(62, 62), c(60, 60), c(0, 143)),
          z = c(c(0, 143), c(60, 60), c(60, 60), c(0, 143), c(62, 62), c(62, 62)), col = '#ff9933', lwd = 3)
#Label lines. If line is missing a coord (x, y or z), the coord of the label will be 120
#Add 1 to the other two coords so that the label rests on the line
rgl.texts(x = 63, y = 61, z = 143, paste('(', '62', ',', '60', ',', 'z)', sep = ''), col = '#FF4000', cex = 0.8, adj = c(0,0))
rgl.texts(x = 143, y = 63, z = 61, paste('(', 'x', ',', '62', ',', '60', ')', sep = ''), col = '#FF4000', cex = 0.8, adj = c(0,0))
rgl.texts(x = 63, y = 143, z = 61, paste('(', '62', ',', 'y', ',', '60', ')', sep = ''), col = '#FF4000', cex = 0.8, adj = c(0,0))
rgl.texts(x = 61, y = 63, z = 133, paste('(', '60', ',', '62', ',', 'z)', sep = ''), col = '#FF4000', cex = 0.8, adj = c(0,0))
rgl.texts(x = 133, y = 61, z = 63, paste('(', 'x', ',', '60', ',', '62)', sep = ''), col = '#FF4000', cex = 0.8, adj = c(0,0))
rgl.texts(x = 61, y = 133, z = 63, paste('(', '60', ',', 'y', ',', '62)', sep = ''), col = '#FF4000', cex = 0.8, adj = c(0,0))

movie3d(spin3d(axis = c(0,0,3), rpm = 2), duration = 23.6, dir = getwd(), clean = FALSE, convert = FALSE, movie = 'B', type = 'gif')

system('convert -delay 1x10 *.png movie.gif')
system('rm A*.png')
system('rm B*.png')


#' @param new.device a logical value. If TRUE, creates a new device
#' @param bg the background color of the device
#' @param width the width of the device
rgl_init <- function(new.device = FALSE, bg = "white", width = 900) { 
  if( new.device | rgl.cur() == 0 ) {
    rgl.open()
    par3d(windowRect = 50 + c( 0, 0, width, width ) )
    rgl.bg(color = bg )
  }
  rgl.clear(type = c("shapes", "bboxdeco"))
  rgl.viewpoint(theta = 15, phi = 20, zoom = 1)
}