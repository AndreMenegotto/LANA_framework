
# ------------------------------------------------------------
# Script:   2_NameStats.R
# Purpose:  Check similarity among valid names of marine species
# Updated:  09-12-2025
# ------------------------------------------------------------


#. Load packages ----
library(stringdist)
library(vioplot)
library(dplyr)



#. Loading and preparing the data ----

## Read the data set
# To access WoRMS data it is necessary to apply a request at 'https://www.marinespecies.org/usersrequest.php'
patch <- "./Data/raw/WoRMS_download_2025-04-01/"
taxon <- data.table::fread(paste(patch,"taxon.txt",sep=""), na.strings=c("","NA"))
taxon <- as.data.frame(taxon)


## Select only valid species names
species <- taxon %>% filter(taxonRank %in% "Species" & taxonomicStatus %in% "accepted")
species <- species %>% select(scientificName, phylum) %>% distinct(scientificName, .keep_all = T)
any(duplicated(species$scientificName)) #FALSE

# Remove hybrid or unidentified species
rmSP <- grep("\\[|\\]|[0-9]|×", species$scientificName)
species <- species[-rmSP,]
species$scientificName <- gsub("['\"]", "", species$scientificName)

speciesList <- strsplit(species$scientificName, " ")
rmSP <- which(sapply(speciesList, function(x) any(nchar(x) == 1)))
speciesList <- speciesList[-rmSP]
species <- species[-rmSP,]

table(lengths(speciesList))
# 3 = subgenus
# 4 = incertae sedis
# >4 = virus description

pos <- which(lengths(speciesList) < 4)
speciesList <- speciesList[pos]
length(speciesList) #229,578

species <- species[pos,]
length(table(species$phylum)) #85

# Select only genus and epithet (ignore subgenus)
genus <- sapply(speciesList,"[[",1)
epithet <- sapply(speciesList,tail,1L)

tmpSpecies <- paste(genus, epithet, sep=" ")
length(tmpSpecies) #229,578
length(unique(tmpSpecies)) #229,508 (removing subgenus created duplicates)

# Select unique names
tmpSpecies <- unique(tmpSpecies)
speciesList <- strsplit(tmpSpecies, " ")
genus <- sapply(speciesList,"[[",1)
epithet <- sapply(speciesList,tail,1L)



#. Comparing similarity (general names) ----

## Get the lowest number of edits among all valid names

# Original vector
vec <- tmpSpecies
n <- length(vec)

# Output to save the lowest number of edits and the closest name
min_dist <- rep(Inf, n)
closest  <- rep(NA_character_, n)

# Loop using blocks
chunk_size <- 500
starts <- seq(1, n, by = chunk_size)

for(a in seq_along(starts))
{
  # Select the first names
  i_start <- starts[a]
  i_end   <- min(i_start + chunk_size - 1, n)
  block_i <- vec[i_start:i_end]
  
  # Compare these names with all the next ones
  j_start <- i_end + 1
  if (j_start > n)
  {
    break
  }
  block_j <- vec[j_start:n]
  
  cat(sprintf("Comparando blocos %d/%d (%d-%d)\n",
              a, length(starts), i_start, i_end))
  flush.console()
  
  # Calculate the distance matrix for block_i × block_j
  m <- stringdistmatrix(block_i, block_j, method = "dl")
  
  # Update the minimal distance for i (rows)
  new_min_i <- apply(m, 1, min)
  pos_i <- apply(m, 1, which.min)
  better_i <- new_min_i < min_dist[i_start:i_end]
  min_dist[i_start:i_end][better_i] <- new_min_i[better_i]
  closest[i_start:i_end][better_i]  <- block_j[pos_i[better_i]]
  
  # Update the minimal distance for j (columns)
  new_min_j <- apply(m, 2, min)
  pos_j <- apply(m, 2, which.min)
  j_idx <- j_start:n
  better_j <- new_min_j < min_dist[j_idx]
  min_dist[j_idx][better_j] <- new_min_j[better_j]
  closest[j_idx][better_j]  <- block_i[pos_j[better_j]]
}
nearNames <- data.frame(vec, closest, min_dist)
#save(nearNames, file = "./Data/tmpRData/nearNames.RData")
#load("./Data/tmpRData/nearNames.RData")


## Check results
sum(nearNames$min_dist==1) #768
sum(nearNames$min_dist==2) #5,528
sum(nearNames$min_dist==3) #21,637
round(c(768,5528,21637)/nrow(nearNames)*100, 2)


## Bar plot
barplot(table(nearNames$min_dist)/nrow(nearNames), space=0, col="grey80", ylim=c(0,.25), xlim=c(0,21), las=2, axisnames=F, cex.axis=1.25)
axis(side=1, at=seq(1,21,2)-.5, labels=NA, tck=-.035)
axis(side=1, at=seq(2,20,2)-.5, labels=NA, tck=-.02, lwd=0, lwd.ticks = 1)
axis(side=1, at=seq(1,21,2)-.5, labels=seq(1,21,2), tick=F, cex.axis=1.25)


## Get names with three adits or less
edt3 <- nearNames[nearNames$min_dist<4,]

speciesList1 <- strsplit(edt3$vec, " ")
speciesList2 <- strsplit(edt3$closest, " ")

# Split genus and epithet
genus1 <- sapply(speciesList1,"[[",1)
epithet1 <- sapply(speciesList1,tail,1L)
genus2 <- sapply(speciesList2,"[[",1)
epithet2 <- sapply(speciesList2,tail,1L)

# Get the number of edits for each component
dlGE <- stringdist(genus1, genus2, method = "dl")
dlEP <- stringdist(epithet1, epithet2, method = "dl")

# Check where the edits are concentrated
View(cbind(edt3, dlGE, dlEP))
round(sum(dlGE==0 | dlEP==0)/nrow(edt3)*100, 2) #98.94%
round(sum(dlGE==0)/nrow(edt3)*100, 2) #75.7%
round(sum(dlEP==0)/nrow(edt3)*100, 2) #23.24%



#. Comparing similarity (when shared elements) ----

## same genera (how distant are the epithets when names share the genus?)
dG <- table(genus)
posDupG <- which(dG>1) # which genera have more than one species...

length(dG) #33,926 genera
length(posDupG) #20,846 non-exclusive
sum(dG==1) #13,080

# Get all the names with non unique genera
tmpSpecies <- speciesList[which(genus %in% names(dG[posDupG]))]
tmpSpecies <- as.data.frame(do.call(rbind, tmpSpecies), stringsAsFactors = FALSE)
dim(tmpSpecies) #216,428

# For each genus
dupGen <- names(dG[posDupG])
minDLep <- numeric()
for(i in 1:length(dupGen)) #20,846
{
  # Get the minimum distance among its epithets
  pos <- which(tmpSpecies$V1 == dupGen[i])
  epi <- tmpSpecies$V2[pos]
  
  dlEp <- stringdistmatrix(epi, epi, method = "dl")
  dlEp[upper.tri(dlEp, diag = TRUE)] <- NA
  dlEp <- na.omit(as.vector(dlEp))
  minDLep[i] <- min(dlEp)
}


## same epithet (how distant are the genera when names share epithet?)
dE <- table(epithet)
posDupE <- which(dE>1) # which epithet is used in more than one species...

length(dE) #85,865 epithets
length(posDupE) #24,655 non-exclusive
sum(dE==1) #61,210

# Most common epithets in WoRMS
sort(dE, decreasing = T)[1:10]
#gracilis=634, australis=600, elegans=489, japonica=488, pacifica=450, antarctica=368, simplex=361, orientalis=349, elongata=344, minuta=333

# Get all the names with non unique epithet
tmpSpecies <- speciesList[which(epithet %in% names(dE[posDupE]))]
tmpSpecies <- as.data.frame(do.call(rbind, tmpSpecies), stringsAsFactors = FALSE)
dim(tmpSpecies) #168,298

# For each epithet
dupEps <- names(dE[posDupE])
minDLge <- numeric()
for(i in 1:length(dupEps)) #24,655
{
  # Get the minimum distance among its genera
  pos <- which(tmpSpecies$V2 == dupEps[i])
  gen <- tmpSpecies$V1[pos]
  
  dlGe <- stringdistmatrix(gen, gen, method = "dl")
  dlGe[upper.tri(dlGe, diag = TRUE)] <- NA
  dlGe <- na.omit(as.vector(dlGe))
  minDLge[i] <- min(dlGe)
}


## Summary
sum(minDLge<2)/length(minDLge)*100 #57   | 0.23
sum(minDLge<3)/length(minDLge)*100 #385  | 1.56
sum(minDLge<4)/length(minDLge)*100 #1372 | 5.56

sum(minDLep<2)/length(minDLep)*100 #303  | 1.45
sum(minDLep<3)/length(minDLep)*100 #1581 | 7.58
sum(minDLep<4)/length(minDLep)*100 #4031 | 19.34


## Bar plot
barplot(table(minDLge)/length(minDLge), space=0, col="grey80", ylim=c(0,.2), xlim=c(0,21), las=2, axisnames=F, cex.axis=1.25)
barplot(table(minDLep)/length(minDLep), space=0, col=adjustcolor("darkred",alpha.f=.1), axes=F, axisnames=F, add=T)
axis(side=1, at=seq(1,21,2)-.5, labels=NA, tck=-.035)
axis(side=1, at=seq(2,20,2)-.5, labels=NA, tck=-.02, lwd=0, lwd.ticks = 1)
axis(side=1, at=seq(1,21,2)-.5, labels=seq(1,21,2), tick=F, cex.axis=1.25)


#### END ####