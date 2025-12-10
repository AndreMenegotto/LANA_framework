
# ------------------------------------------------------------
# Script:   1_ApplicationCrete.R
# Purpose:  Test the LANA framework in practice
# Updated:  09-12-2025
# ------------------------------------------------------------

# Let's retrieve some marine records around Crete from GBIF...
#...and harmonize it with WoRMS


#. Load packages ----
library(terra)
library(rgbif)
library(sf)
library(dplyr)
library(worrms)



#. Get GBIF records ----

## Download land polygons
land10 <- rnaturalearth::ne_download(scale = 10, type = "land", category = "physical")
land10 <- vect(land10)
land10 <- disagg(land10)


## Find Crete island and create a 100km buffer
x <- expanse(land10, "km")
land10 <- land10[order(x, decreasing = T)]
crete <- land10[90]
creteBuffer <- buffer(crete, 100*1000) #100 km


## Check it
plot(creteBuffer)
plot(land10, add=T)
plot(crete, col="grey", add=T)


## Convert the polygons to wkt
crete100 <- st_as_sf(creteBuffer)
geome <- st_as_text(st_make_valid(st_sfc(crete100$geometry)))
nchar(geome)


## Get number of records in GBIF
occ_data(hasCoordinate=T, occurrenceStatus="PRESENT", taxonKey=c(1,6,4),
         geometry=geome, limit=0) #343,877


## Download the records from GBIF
# (gbif_download <- occ_download(
#   pred("hasCoordinate", TRUE),
#   pred("occurrenceStatus", "PRESENT"),
#   pred_within(geome),
#   pred_in("taxonKey",c(1,6,4)),
# 
#   user = "", pwd = "", email = "",
#   format = "SIMPLE_CSV"
# ))


## Read the data
if (!dir.exists("./Data/Raw/")) {dir.create("./Data/Raw/", recursive = TRUE)}
getGBIF <- occ_download_get(key='0029179-251120083545085', path = "./Data/Raw/")
creteGBIF <- occ_download_import(getGBIF)
dim(creteGBIF) #343,276


## Clean dataset
creteSP <- creteGBIF %>% filter(taxonRank %in% c("SPECIES","SUBSPECIES","VARIETY","FORM") & !basisOfRecord %in% "FOSSIL_SPECIMEN")
dim(creteSP) #316,780
length(unique(creteSP$scientificName)) #11,141

creteSP_points <- vect(creteSP, geom=c("decimalLongitude","decimalLatitude"), crs="epsg:4326")
creteInvBuffer <- disagg(buffer(crete, -5*1000))[1,]

creteSP_marine <- erase(creteSP_points, creteInvBuffer)
length(creteSP_marine) #234,850


## plot points
# plot(creteBuffer, axes=F)
# plot(creteSP_marine, add=T, col=rgb(0,0,1,.1), cex=.25)
# plot(crete, add=T)


## Get spatially filtered data set
creteMarine <- values(creteSP_marine)
length(unique(creteMarine$scientificName)) #9,335
nrow(creteMarine) #234,850



#. Applying the LANA framework ----

## Load functions
source("./Scripts/LANA_Functions.R")


## Get species names
zz <- creteMarine %>% distinct(scientificName, .keep_all = T)
dim(zz) #9,335


## Step 1: Preparation
targetNames <- preFilter(taxInput = zz)
nrow(targetNames) #9,310

# Note that hybrid species were removed
setdiff(zz$scientificName, targetNames$scientificname)


## Step 2: Search
x <- wTaxaMatch(taxTarget = targetNames, fuzzyMatch = F)
length(unique(x$SN_gbif)) #4,190
nrow(x) #4,509

# Intermediate step to ensure that habitat and rank information are filled
# and standardize the use of 'Linnaeus' abbreviation
x2 <- alignData(wormsTab = x, taxTarget = targetNames, kRank = "kingdom", verbose = F)
length(unique(x2$SN_gbif)) #4,190
nrow(x2) #4,509

# Intermediate step to pre-filter suggestions above species level
# or not belonging to the indicated ranks
x3 <- filterTaxa(wormsTab = x2, kRank = "kingdom", kTaxa = c("Animalia","Plantae","Chromista"))
length(unique(x3$SN_gbif)) #4,189
nrow(x3) #4,508


## Step 3: Score
x4 <- taxaScoring(wormsTab = x3, taxTarget = targetNames)



# Define new search (fuzzy match)
# Look for names not found or with low score
f1 <- setdiff(targetNames$scientificname, x3$SN_gbif)
f2 <- x4$SN_gbif[which(x4$orthScore<1 | x4$authScore<1 | x4$systScore<1)]
buscar <- unique(c(f1,f2))

# but avoid cases where an alternative perfect suggestion is already present
goodMatches <- unique(x4$SN_gbif[which(x4$orthScore==1 & x4$authScore==1 & x4$systScore==1)])
buscar <- setdiff(buscar, goodMatches)

# Names to search
buscar <- targetNames[targetNames$scientificname %in% buscar,]
nrow(buscar) #6,052 (notice that some names where found but don't have the highest score, so they'll be queried again)

if(nrow(buscar)>0)
{
  # Search B
  xb <- wTaxaMatch(taxTarget = buscar, fuzzyMatch = T)
  length(unique(xb$SN_gbif)) #1,056
  nrow(xb) #1,258
  
  # Alignment B
  x2b <- alignData(wormsTab = xb, taxTarget = targetNames, kRank = "kingdom", verbose = F)
  length(unique(x2b$SN_gbif)) #1,056
  nrow(x2b) #1,258
  
  # Filter 1 B
  x3b <- filterTaxa(wormsTab = x2b, kRank = "kingdom", kTaxa = c("Animalia","Plantae","Chromista"))
  length(unique(x3b$SN_gbif)) #1,055
  nrow(x3b) #1,257
  
  # Scoring B
  x4b <- taxaScoring(wormsTab = x3b, taxTarget = targetNames)
  length(unique(x4b$SN_gbif)) #1,055
  nrow(x4b) #1,257
  
  x4c <- rbind(x4, x4b)
}

x4U <- x4c %>% distinct(AphiaID, SN_gbif, .keep_all = T)
length(unique(x4U$SN_gbif)) #4,313
nrow(x4U) #4,687


## Step 4: Filter
x5 <- filterNames(wormsTab = x4U, taxTarget = targetNames)
length(unique(x5$SN_gbif)) #2,681
nrow(x5) #2,681


## Step 5: Classification
x6 <- confMatch(wormsTab = x5)
nrow(x6) #2,681


## Step 6: Correction
x7 <- taxaRefresh(wormsTab = x6, taxTarget = targetNames, verbose = T)
nrow(x7) #2,672 #Some names were removed after correction they were found non-marine or invalid names
#save.image(file = "./Data/tmpRData/Crete.RData")



#. Summary of the match ----
confLevels <- c("Keep [H-CL]","Keep [M-CL]","Keep [L-CL]","Drop [H-CL]","Drop [M-CL]","Drop [L-CL]")

summMatch <- table(x7$confMatch)
pos <- match(confLevels, names(summMatch))
summMatch <- summMatch[pos]
summMatch[is.na(names(summMatch))] <- 0

tabSumm <- matrix(summMatch, nrow = 3, ncol = 2, byrow = F)
rownames(tabSumm) <- c("H-CL","M-CL","L-CL")
colnames(tabSumm) <- c("Keep","Drop")

{
  cat("\n\n"); msg<-"Matching confidence summary:\n"; cat(msg); cat(strrep("-", nchar(msg)), "\n")
  
  cat(sprintf("%-5s | %-6s | %-6s\n", "", "Keep", "Drop"))
  cat(strrep("-", 29), "\n")
  
  for (i in 1:nrow(tabSumm))
  {
    cat(sprintf("%-5s | %-6d | %-6d\n",
                rownames(tabSumm)[i],
                tabSumm[i, 1],
                tabSumm[i, 2]))
  }
  cat(strrep("-", nchar(msg)), "\n")
}



#. Summary of the search and additional outputs ----

# Searched
NM_searched <- unique(targetNames$scientificname)
length(NM_searched) #9,310

# Found
NM_found <- unique(x4U$SN_gbif)
length(NM_found) #4,313

# Found and kept
NM_kept <- unique(x7$SN_gbif)
length(NM_kept) #2,672


# Additional outputs
NM_notFound <- setdiff(NM_searched, NM_found)
length(NM_notFound) #4,997

NM_Excluded <- setdiff(NM_found, NM_kept)
WoRMS_Excluded <- x4U[which(x4U$SN_gbif %in% NM_Excluded),]
length(NM_Excluded) #1,641

{
  cat("\n\n"); msg<-"Searching summary:\n"; cat(msg); cat(strrep("-", nchar(msg)+10), "\n")
  
  cat(paste("Names classified: ", length(NM_kept), sep=""), "\n")
  cat(paste("Names excluded:   ", length(NM_Excluded), sep=""), "\n")
  cat(paste("Names not found:  ", length(NM_notFound), sep=""), "\n")
  cat(strrep("-", nchar(msg)+10), "\n")
}



#. Output exploration ----


## Mismatch summary: Moderate confidence
temp <- x7 %>% filter(confMatch == "Keep [M-CL]")

(1-round(sum(temp$orthScore >= .9)/nrow(temp), 2))*100 #0
(1-round(sum(temp$authScore >= 0.9, na.rm = T)/nrow(temp), 2))*100 #84%
(1-round(sum(temp$systScore == 1)/nrow(temp), 2))*100 #15%

# None=0, Kingdom=0.2, Phylo=0.4, Class=0.6, Order=0.8, Family=1.0
sys <- table(temp$systScore[temp$systScore < 1])
round(sys/sum(sys), 2)

# -1=NA, 0=mismatch, 0.5=partial match
aut <- temp$authScore[temp$authScore < .9 | is.na(temp$authScore)]
aut <- ifelse(is.na(aut), -1, aut)
aut <- ifelse(aut > 0 & aut < .9, .5, aut)
aut <- table(aut)
round(aut/sum(aut), 2)


## Mismatch summary: Low confidence
temp <- x7 %>% filter(confMatch == "Keep [L-CL]")# %>% select(orthScore,authScore,systScore)

(1-round(sum(temp$orthScore >= .9)/nrow(temp), 2))*100 #9%
(1-round(sum(temp$authScore > 0.9, na.rm = T)/nrow(temp), 2))*100 #91%
(1-round(sum(temp$systScore == 1)/nrow(temp), 2))*100 #100%

# None=0, Kingdom=0.2, Phylo=0.4, Class=0.6, Order=0.8, Family=1.0
sys <- table(temp$systScore[temp$systScore < 1])
round(sys/sum(sys), 2)

# -1=NA, 0=mismatch, 0.5=partial match
aut <- temp$authScore[temp$authScore < .9 | is.na(temp$authScore)]
aut <- ifelse(is.na(aut), -1, aut)
aut <- ifelse(aut > 0 & aut < .9, .5, aut)
aut <- table(aut)
round(aut/sum(aut), 2)


## View recommendations to be discarded
temp <- x7[grep("Drop", x7$confMatch),]
View(temp)


## count valid species
length(unique(x7$scientificname[grepl("Keep",x7$confMatch)])) #2,519
length(unique(x7$scientificname[x7$confMatch %in% "Keep [H-CL]"])) #2,189



#. Mapping points ----

## Create the raster
r <- rast(ext(creteBuffer), resolution=.1, crs="epsg:4326")
p <- as.polygons(r)
p <- terra::intersect(p, creteBuffer)
centR <- centroids(p)


## Get geographical information
cc <- geom(creteSP_marine)
creteFinalSP <- data.frame(creteMarine, long=cc[,"x"], lat=cc[,"y"])


## Filter data set by marine records with valid name
# creteFinalSP <- creteFinalSP %>%
#   select(scientificName, species, lat, long) %>%
#   right_join(x7[grepl("Keep",x7$confMatch),] %>% select(SN_gbif,AphiaID,scientificname,kingdom,phylum,class,order,family,confMatch), by = c("scientificName" = "SN_gbif")) %>%
#   rename(valid_name = scientificname)

creteFinalSP <- creteFinalSP %>%
  select(scientificName, species, lat, long) %>%
  left_join(x7 %>% select(SN_gbif,AphiaID,scientificname,kingdom,phylum,class,order,family,confMatch), by = c("scientificName" = "SN_gbif")) %>%
  rename(valid_name = scientificname)


# Spatialize the occurrence records and relate to the grid cells
occPoints <- vect(creteFinalSP, geom=c("long","lat"), crs="epsg:4326")
xRelate <- relate(p, occPoints, relation="intersects", pairs=T)

creteFinalSP$ID <- 1:nrow(creteFinalSP)
xRelate <- as.data.frame(xRelate)

creteFinalSP2 <- creteFinalSP %>%
  left_join(xRelate, by = c("ID" = "id.y")) %>%
  rename(cell=id.x) %>% select(!ID)


## Plot original records OR...
res <- creteFinalSP2 %>%
  group_by(cell) %>%
  summarise(n_recs = n(), .groups = "drop")
sum(res$n_recs) #234,850

# Plot records recommended to be kept
res <- creteFinalSP2 %>%
  filter(grepl("Keep", creteFinalSP2$confMatch)) %>%
  group_by(cell) %>%
  summarise(n_recs = n(), .groups = "drop")
sum(res$n_recs) #41,982

# Do the plot
resRecs <- rep(0, length(p))
resRecs[res$cell] <- res$n_recs

colPall <- c("#f7fbff","#deebf7","#c6dbef","#9ecae1","#6baed6","#4292c6","#2171b5","#08519c","#08306b","black")
myCol <- colPall[cut(log10(resRecs), seq(0, 5, length.out=11))]
max(log10(resRecs))

plot(creteBuffer, border="transparent", axes=F, box=F, lty=2)
plot(centR, add=T, pch=21, cex=log10(resRecs)/2, bg=myCol, col=NA)
plot(land10, col="grey", add=T)
plot(creteBuffer, box=T, lty=2, add=T)

  
## Plot richness map

# Quantify observed richness by grid cell
res <- creteFinalSP2 %>%
  filter(grepl("Keep", creteFinalSP2$confMatch)) %>%
  group_by(cell) %>%
  summarise(n_species = n_distinct(valid_name))

resKeep <- rep(0, length(p))
resKeep[res$cell] <- res$n_species

res <- creteFinalSP2 %>%
  filter(creteFinalSP2$confMatch %in% "Keep [H-CL]") %>%
  group_by(cell) %>%
  summarise(n_species = n_distinct(valid_name))

resKHigh <- rep(0, length(p))
resKHigh[res$cell] <- res$n_species

# Spatial correlation
cor(resKeep, resKHigh) #.99

# Plot
mm <- plyr::round_any(max(resKeep), 10, f = ceiling)
myCol <- rev(viridis::magma(length(seq(0,mm,5))-1))[cut(resKeep, seq(0,mm,5))]

plot(creteBuffer, border="transparent", axes=F, box=F, lty=2)
plot(p, col=myCol, border="grey90", add=T)
plot(land10, col="grey", add=T)
plot(creteBuffer, box=T, lty=2, add=T)



#. Comparing species names ----

## Filter data set by marine records with valid name
TabFinalSP <- creteMarine %>%
  select(scientificName, species) %>%
  right_join(x7[grepl("Keep",x7$confMatch),] %>% select(SN_gbif,AphiaID,scientificname,kingdom,phylum,class,order,family,confMatch), by = c("scientificName" = "SN_gbif")) %>%
  rename(valid_name = scientificname) %>%
  distinct(scientificName, .keep_all = T)


## Compare total richness
length(unique(TabFinalSP$valid_name)) #2519
length(unique(TabFinalSP$species)) #2525

length(unique(TabFinalSP$valid_name[TabFinalSP$confMatch == "Keep [H-CL]"])) #2189
length(unique(TabFinalSP$species[TabFinalSP$confMatch == "Keep [H-CL]"])) #2193


## Compare valid names (for each original scientific name)
sum(TabFinalSP$species==TabFinalSP$valid_name)/nrow(TabFinalSP) #2478 | 93%


## Compare valid names (between databases)
library(VennDiagram)

# Make diagram
lista <- list(GBIF = unique(TabFinalSP$species), WoRMS = unique(TabFinalSP$valid_name))
venn.plot <- venn.diagram(
  x = lista,
  category.names = c("GBIF", "WoRMS"),
  disable.logging = T,
  filename = NULL,
  output = F,
  scaled = F,
  
  # Circles
  lwd = 1,
#  lty = 'blank',
#  col = c("#33a02c","#1f78b4"),
  fill = c("#33a02c","#1f78b4"),
  alpha = 0.5,
  
  # Numbers
  cex = 1.5,
  fontface = "bold",
  fontfamily = "sans",
  
  # Set names
  cat.cex = 1.5,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.fontfamily = "sans",
  cat.col = c("#006d2c","#045a8d")
)

# Plot diagramm
grid::grid.newpage()
grid::grid.draw(venn.plot)


#### END ####