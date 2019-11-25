library(sp)
library(raster)
library(rgdal)
library(dismo)
library(plyr)

# Make a list of directories that can be looped over that can be looped over
dirs = list.dirs("D:/FinalModelRunsFixed/sod", full.names = TRUE)
#open Original raster (OR) to get extent and projection
OR = raster("D:/FinalModelRunsFixed/sod/1028/initialcommunitiesmap.img")

# All data using visual symptoms and unburned (comparing burned plots doesn't work as well )
s2006PS <- readOGR("D:/sod_data/2006DiseasestatusPS.shp")
s2007PS <- readOGR("D:/sod_data/2007DiseasestatusPS.shp")
s2009PS <- readOGR("D:/sod_data/2009DiseasestatusPSUB.shp")
s2010PS <- readOGR("D:/sod_data/2010DiseasestatusPSUB.shp")
s2011PS <- readOGR("D:/sod_data/2011DiseasestatusPSUB.shp")
s2013PS <- readOGR("D:/sod_data/2013DiseasestatusPSUB.shp")

# Create a holding data frame for the output to build confusion matrix and calculate odds ratio
oddR2006 <- data.frame("Number" = 0,"OddsRatio1" = 0, "pospos" = 0,"posneg" = 0, "negpos" = 0, "negneg" = 0, "total" = 0) 
oddR2007 <- data.frame("Number" = 0,"OddsRatio1" = 0, "pospos" = 0,"posneg" = 0, "negpos" = 0, "negneg" = 0, "total" = 0) 
oddR2009 <- data.frame("Number" = 0,"OddsRatio1" = 0, "pospos" = 0,"posneg" = 0, "negpos" = 0, "negneg" = 0, "total" = 0) 
oddR2010 <- data.frame("Number" = 0,"OddsRatio1" = 0, "pospos" = 0,"posneg" = 0, "negpos" = 0, "negneg" = 0, "total" = 0) 
oddR2011 <- data.frame("Number" = 0,"OddsRatio1" = 0, "pospos" = 0,"posneg" = 0, "negpos" = 0, "negneg" = 0, "total" = 0) 
oddR2013 <- data.frame("Number" = 0,"OddsRatio1" = 0, "pospos" = 0,"posneg" = 0, "negpos" = 0, "negneg" = 0, "total" = 0) 

## Get simulated years that match plot data years since we start in 1991 then year 16 of the simulation is 2006 and year 23 is 2013
file6 <- list.files(path = dirs, pattern = "ramorum-16.img", full.names = TRUE)
file7 <- list.files(path = dirs, pattern = "ramorum-17.img", full.names = TRUE)
file9 <- list.files(path = dirs, pattern = "ramorum-19.img", full.names = TRUE)
file10 <- list.files(path = dirs, pattern = "ramorum-20.img", full.names = TRUE)
file11 <- list.files(path = dirs, pattern = "ramorum-21.img", full.names = TRUE)
file13 <- list.files(path = dirs, pattern = "ramorum-23.img", full.names = TRUE)

nums <- c(seq(1,12,1), seq(14, length(file6), 1))
rcl <- c(2, Inf, 2, 1, 1.99, 1, 0, 0.99, NA)
rclmat <- matrix(rcl, ncol=3, byrow=TRUE)
#Loop over dirs to create confusion tables and 
for (i in 1:length(file6)){
  #open raster that needs shifted, projected and rescaled
  rast6 <- raster(file6[i])
  rast7 <- raster(file7[i])
  rast9 <- raster(file9[i])
  rast10 <- raster(file10[i])
  rast11 <- raster(file11[i])
  rast13 <- raster(file13[i])
  
  #Create raster same size and projection as the Original then set the values to that of the raster of interest and scale up 30 m resolution to 60 m resolution to account for plot size and gps uncertainity
  rast26 <- raster(ext = OR@extent, resolution = res(OR),crs = projection(OR))
  rast26[] <- getValues(raster(file6[i]))
  rast26 <- reclassify(rast26, rclmat)
  rast26 <- aggregate(rast26, fact = 2, fun = max, na.rm = TRUE)
  
  rast27 <- raster(ext = OR@extent, resolution = res(OR),crs = projection(OR))
  rast27[] <- getValues(rast7)
  rast27 <- reclassify(rast27, rclmat)
  rast27 <- aggregate(rast27, fact = 2, fun = max, na.rm = TRUE)
  
  rast29 = raster(ext=OR@extent, resolution=res(OR),crs=projection(OR))
  rast29[] = getValues(rast9)
  rast29 <- reclassify(rast29, rclmat)
  rast29 <- aggregate(rast29, fact = 2, fun = max, na.rm = TRUE)
  
  rast210 = raster(ext=OR@extent, resolution=res(OR),crs=projection(OR))
  rast210[] = getValues(rast10)
  rast210 <- reclassify(rast210, rclmat)
  rast210 <- aggregate(rast210, fact = 2, fun = max, na.rm = TRUE)
  
  rast211 = raster(ext=OR@extent, resolution=res(OR),crs=projection(OR))
  rast211[] = getValues(rast11)
  rast211 <- reclassify(rast211, rclmat)
  rast211 <- aggregate(rast211, fact = 2, fun = max, na.rm = TRUE)
  
  rast213 = raster(ext=OR@extent, resolution=res(OR),crs=projection(OR))
  rast213[] = getValues(rast13)
  rast213 <- reclassify(rast213, rclmat)
  rast213 <- aggregate(rast213, fact = 2, fun = max, na.rm = TRUE)
  
  #Extract the data from the raster to the points
  s2006PS$ModelSOD <- extract(rast26, s2006PS)
  s2007PS$ModelSOD <- extract(rast27, s2007PS)
  s2009PS$ModelSOD <- extract(rast29, s2009PS)
  s2010PS$ModelSOD <- extract(rast210, s2010PS)
  s2011PS$ModelSOD <- extract(rast211, s2011PS)
  s2013PS$ModelSOD <- extract(rast213, s2013PS)
  
  # Use with other datasets
  actual <- as.vector(as.numeric(s2006PS$pram2006))
  model<- as.vector(s2006PS$ModelSOD)
  eda <- data.frame("model" = model, "actual" = actual)
  eda$model[eda$model==3] <- 2
  Pos <- eda[eda$actual==2,]
  Neg <- eda[eda$actual==1,]
  comp <- data.frame(pos = NA, neg = NA)
  Pos$comp <- Pos$actual - Pos$model
  Neg$comp <- Neg$actual - Neg$model
  comp$neg <- sum(Pos$comp)
  comp$pos <- nrow(Pos)-sum(Pos$comp)
  comp[2,1] <- abs(sum(Neg$comp))
  comp[2,2] <- nrow(Neg)-abs(sum(Neg$comp))
  row.names(comp) <- c("actual positive", "actual negative")
  names(comp) <- c("model positive", "model negative")
  comp$rate <- c(comp[1,1]/sum(comp[1,])*100,comp[2,2]/sum(comp[2,])*100)
  oddsratio <-(comp[1,1]*comp[2,2])/(comp[1,2]*comp[2,1])
  comps <- data.frame(comp[1,1], comp[1,2], comp[2,1], comp[2,2], nrow(eda))
  
  oddR2006[i,] <- c(file6[i], oddsratio, comps[1,1], comps[1,2], comps[1,3], comps[1,4], comps[1,5])
  
  actual <- as.vector(as.numeric(s2007PS$pram2006))
  model <- as.vector(s2007PS$ModelSOD)
  eda <- data.frame("model" = model, "actual" = actual)
  eda$model[eda$model==3] <- 2
  Pos <- eda[eda$actual==2,]
  Neg <- eda[eda$actual==1,]
  comp <- data.frame(pos = NA, neg = NA)
  Pos$comp <- Pos$actual-Pos$model
  Neg$comp <- Neg$actual-Neg$model
  comp$neg <- sum(Pos$comp)
  comp$pos <- nrow(Pos)-sum(Pos$comp)
  comp[2,1] <- abs(sum(Neg$comp))
  comp[2,2] <- nrow(Neg)-abs(sum(Neg$comp))
  row.names(comp) <- c("actual positive", "actual negative")
  names(comp) <- c("model positive", "model negative")
  comp$rate <- c(comp[1,1]/sum(comp[1,])*100,comp[2,2]/sum(comp[2,])*100)
  oddsratio <-(comp[1,1]*comp[2,2])/(comp[1,2]*comp[2,1])
  comps <- data.frame(comp[1,1],comp[1,2],comp[2,1],comp[2,2],nrow(eda))
  
  oddR2007[i,] <- c(file7[i], oddsratio, comps[1,1], comps[1,2], comps[1,3], comps[1,4], comps[1,5])
  
  actual <- as.vector(as.numeric(s2009PS$pram2006))
  model <- as.vector(s2009PS$ModelSOD)
  eda <- data.frame("model"= model, "actual"=actual)
  eda$model[eda$model==3]=2
  Pos <- eda[eda$actual==2,]
  Neg <- eda[eda$actual==1,]
  comp <- data.frame(pos=NA, neg=NA)
  Pos$comp = Pos$actual-Pos$model
  Neg$comp = Neg$actual-Neg$model
  comp$neg = sum(Pos$comp)
  comp$pos = nrow(Pos)-sum(Pos$comp)
  comp[2,1] = abs(sum(Neg$comp))
  comp[2,2] = nrow(Neg)-abs(sum(Neg$comp))
  row.names(comp) <- c("actual positive", "actual negative")
  names(comp) <- c("model positive", "model negative")
  comp$rate = c(comp[1,1]/sum(comp[1,])*100,comp[2,2]/sum(comp[2,])*100)
  oddsratio =(comp[1,1]*comp[2,2])/(comp[1,2]*comp[2,1])
  comps = data.frame(comp[1,1],comp[1,2],comp[2,1],comp[2,2],nrow(eda))
  
  oddR2009[i,]= c(file9[i], oddsratio, comps[1,1], comps[1,2], comps[1,3], comps[1,4], comps[1,5])
  
  actual = as.vector(as.numeric(s2010PS$pram2006))
  model= as.vector(s2010PS$ModelSOD)
  eda <- data.frame("model"= model, "actual"=actual)
  eda$model[eda$model==3]=2
  Pos <- eda[eda$actual==2,]
  Neg <- eda[eda$actual==1,]
  comp <- data.frame(pos=NA, neg=NA)
  Pos$comp = Pos$actual-Pos$model
  Neg$comp = Neg$actual-Neg$model
  comp$neg = sum(Pos$comp)
  comp$pos = nrow(Pos)-sum(Pos$comp)
  comp[2,1] = abs(sum(Neg$comp))
  comp[2,2] = nrow(Neg)-abs(sum(Neg$comp))
  row.names(comp) <- c("actual positive", "actual negative")
  names(comp) <- c("model positive", "model negative")
  comp$rate = c(comp[1,1]/sum(comp[1,])*100,comp[2,2]/sum(comp[2,])*100)
  oddsratio =(comp[1,1]*comp[2,2])/(comp[1,2]*comp[2,1])
  comps = data.frame(comp[1,1],comp[1,2],comp[2,1],comp[2,2],nrow(eda))
  
  oddR2010[i,]= c(file10[i], oddsratio, comps[1,1], comps[1,2], comps[1,3], comps[1,4], comps[1,5])
  
  actual = as.vector(as.numeric(s2011PS$pram2006))
  model= as.vector(s2011PS$ModelSOD)
  eda <- data.frame("model"= model, "actual"=actual)
  eda$model[eda$model==3]=2
  Pos <- eda[eda$actual==2,]
  Neg <- eda[eda$actual==1,]
  comp <- data.frame(pos=NA, neg=NA)
  Pos$comp = Pos$actual-Pos$model
  Neg$comp = Neg$actual-Neg$model
  comp$neg = sum(Pos$comp)
  comp$pos = nrow(Pos)-sum(Pos$comp)
  comp[2,1] = abs(sum(Neg$comp))
  comp[2,2] = nrow(Neg)-abs(sum(Neg$comp))
  row.names(comp) <- c("actual positive", "actual negative")
  names(comp) <- c("model positive", "model negative")
  comp$rate = c(comp[1,1]/sum(comp[1,])*100,comp[2,2]/sum(comp[2,])*100)
  oddsratio =(comp[1,1]*comp[2,2])/(comp[1,2]*comp[2,1])
  comps = data.frame(comp[1,1],comp[1,2],comp[2,1],comp[2,2],nrow(eda))
  
  oddR2011[i,]= c(file11[i], oddsratio, comps[1,1], comps[1,2], comps[1,3], comps[1,4], comps[1,5])
  
  actual = as.vector(as.numeric(s2013PS$pram2006))
  model= as.vector(s2013PS$ModelSOD)
  eda <- data.frame("model"= model, "actual"=actual)
  eda$model[eda$model==3]=2
  Pos <- eda[eda$actual==2,]
  Neg <- eda[eda$actual==1,]
  comp <- data.frame(pos=NA, neg=NA)
  Pos$comp = Pos$actual-Pos$model
  Neg$comp = Neg$actual-Neg$model
  comp$neg = sum(Pos$comp)
  comp$pos = nrow(Pos)-sum(Pos$comp)
  comp[2,1] = abs(sum(Neg$comp))
  comp[2,2] = nrow(Neg)-abs(sum(Neg$comp))
  row.names(comp) <- c("actual positive", "actual negative")
  names(comp) <- c("model positive", "model negative")
  comp$rate = c(comp[1,1]/sum(comp[1,])*100,comp[2,2]/sum(comp[2,])*100)
  oddsratio =(comp[1,1]*comp[2,2])/(comp[1,2]*comp[2,1])
  comps = data.frame(comp[1,1],comp[1,2],comp[2,1],comp[2,2],nrow(eda))
  
  oddR2013[i,] <- c(file13[i], oddsratio, comps[1,1], comps[1,2], comps[1,3], comps[1,4], comps[1,5])
}
## set up odds ratio data frame
oddR <- oddR2006

for (i in 2:7) {
  oddR2006[,i] <- as.numeric(oddR2006[,i])
  oddR2007[,i] <- as.numeric(oddR2007[,i])
  oddR2009[,i] <- as.numeric(oddR2009[,i])
  oddR2010[,i] <- as.numeric(oddR2010[,i])
  oddR2011[,i] <- as.numeric(oddR2011[,i])
  oddR2013[,i] <- as.numeric(oddR2013[,i])
}
## use all years of data to manage 
oddR[,3:7] <- oddR2006[,3:7] + oddR2007[,3:7] + oddR2009[,3:7] + oddR2010[,3:7] + oddR2011[,3:7] + oddR2013[,3:7]
oddR[,2] <- (oddR[,3]*oddR[,6])/(oddR[,4]*oddR[,5])
oddR$posaccurracy <- oddR$pospos / (oddR$pospos + oddR$posneg)
oddR$negaccurracy <- oddR$negneg / (oddR$negneg + oddR$negpos)
oddR$accuraccy <- (oddR$pospos + oddR$negneg) / (oddR$total)

# write.csv(oddR, './data/OddsRatioandAccuracy.csv')

## Start here if the odds ratio csv is already in the data folder
Odds_ratio <- read.csv('./data/OddsRatioandAccuracy.csv', header = T)
odds_ratio_mean <- colMeans(Odds_ratio[,3:ncol(Odds_ratio)])

odd <- Odds_ratio[,3:ncol(Odds_ratio)]
odds_ratio_sd <- colwise(sd)(odd)
odds_ratio_mean
odds_ratio_sd

