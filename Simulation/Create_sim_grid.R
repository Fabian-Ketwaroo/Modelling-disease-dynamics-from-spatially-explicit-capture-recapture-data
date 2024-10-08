rm(list=ls())
library(nimble)
library(nimbleSCR)
# load("~/SCR/disease_dpp_data_2014_2018.RData")
# load("~/SCR/koscrdata_T2014_2018.Rdata")
# load("~/SCR/ktraplocs_2014_2018.Rdata")
# load("~/SCR/off_den_2014_2018.Rdata")

# load("C:/Users/Fabian Ketwaroo/Documents/SCR/R Models y1/Rmodsy1/koscrdata_T2014_2018.Rdata")
# load("C:/Users/Fabian Ketwaroo/Documents/SCR/R Models y1/Rmodsy1/ktraplocs_2014_2018.Rdata")
# load("C:/Users/Fabian Ketwaroo/Documents/SCR/R Models y1/Rmodsy1/disease_dpp_data_2014_2018.RData")
# load("C:/Users/Fabian Ketwaroo/Documents/SCR/R Models y1/Rmodsy1/off_den_2014_2018.Rdata")

load("C:/Users/fke/OneDrive - Vogelwarte/SCR/Revisions/ktraplocs_2014_2018.Rdata")


# Present set up
buffer = 0 #  0.5km which is 0.75 km buffer size 

ntraps <- nrow(traplocs)
Xl <- min(traplocs[, 1] - buffer)
Xu <- max(traplocs[, 1] + buffer )
Yl <- min(traplocs[, 2] - buffer)
Yu <- max(traplocs[, 2] + buffer ) # these values here are already scaled.
area <- (Xu - Xl) * (Yu - Yl)/4; area 

plot(traplocs) # trapping grid
plot(traplocs, xlim = c(Xl,Xu), ylim = c(Yl, Yu)) # trapping grid
#abline(h = c(min(traplocs[, 2]),  max(traplocs[, 2]) ) )
#abline(v = c(min(traplocs[, 1]),  max(traplocs[, 1]) ) )

Xu;Xl
Yl;Yu

t1 = cbind( traplocs[,1]+1, traplocs[,2])
buffer = 1 #  0.5km which is 0.75 km buffer size 

ntraps <- nrow(t1)
Xl <- min(t1[, 1] - buffer)
Xu <- max(t1[, 1] + buffer )
Yl <- min(t1[, 2] )
Yu <- max(t1[, 2]  ) # these values here are already scaled.
area <- (Xu - Xl) * (Yu - Yl)/4; area 
plot(t1, xlim = c(Xl,Xu), ylim = c(Yl, Yu) )

Xu;Xl
Yl;Yu



# CREATE HABITAT GRID 
coordsHabitatGridCenter <- cbind(rep(seq(Xu, Xl, by= -1), 10),  sort(rep(seq(Xu, Xl, by= -1), 10)))
colnames(coordsHabitatGridCenter) <- c("x","y")
head(coordsHabitatGridCenter)

# Create trap grids I wont use 
trapCoords <- cbind(rep(seq(-3, 3,by=1),7),
                    sort(rep(seq(-3, 3,by=1),7)))
colnames(trapCoords) <- c("x","y")

# PLOT CHECK
plot(coordsHabitatGridCenter[,"y"] ~ coordsHabitatGridCenter[,"x"], pch = 1, cex = 1.5) #pch=16) 
points(t1, col="red", pch = 15)
par(xpd=TRUE)
legend(x = 2, y = 8.3,
       legend=c("Habitat grid centers", "Traps"),
       pt.cex = c(1.5,1),
       horiz = T,
       pch=c(1,16),
       col=c("black", "red"),
       bty = 'n')


habitatMask <- matrix(1, nrow = 11, ncol= 11, byrow = TRUE)


## Rescale coordinates
scaledObjects <- scaleCoordsToHabitatGrid(
  coordsData = trapCoords,
  coordsHabitatGridCenter = coordsHabitatGridCenter)


## Get lower and upper cell coordinates
lowerAndUpperCoords <- getWindowCoords(
  scaledHabGridCenter = scaledObjects$coordsHabitatGridCenterScaled,
  plot.check = T)

#points(scaledObjects$coordsDataScaled, pch= 15) # trap location
#points(lowerAndUpperCoords$lowerObsCoords, pch = 16, col = "green", cex =1.5)
#points(lowerAndUpperCoords$upperObsCoords, pch = 17, , col = "blue", cex = 1)




HabitatGid_dis <- getWindowCoords(
  scaledHabGridCenter = scaledObjects$coordsHabitatGridCenterScaled,
  plot.check = T)

points(t1, pch=16)

(habitatGrid = lowerAndUpperCoords$habitatGrid)
# What about I just use 0-11 on xaxis and 0-8 on yaxis. In this way I wont need to compute the cells above 8 on the x-axis
#(habitatGrid = habitatGrid[4:11,1:11])

#save(habitatGrid, file = "BadgerHabitatGrid.Rdata")
#traplocs = t1
#save(traplocs, file = "rescale_traplocs_2014_2018.Rdata")




#(habitatGrid = habitatGrid[4:11,1:11])


(habitatGrid = matrix(45:1, 5, 9, byrow = T))

save(habitatGrid, file = "simBadgerHabitatGrid.Rdata")

traplocs = t1
save(traplocs, file = "simrescale_traplocs_2014_2018.Rdata")
