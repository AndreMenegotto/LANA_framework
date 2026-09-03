
# ------------------------------------------------------------
# Script:   LANA_Functions.R
# Purpose:  Execute the LANA framework
# Updated:  13-07-2026
# ------------------------------------------------------------


## Required packages
require(stringdist)
require(stringr)
require(stringi)
require(dplyr)



#. Prepraration ----

# Prepare input data set
preFilter <- function(taxInput)
{
  # Ensure all columns are lowercase
  colnames(taxInput) <- tolower(colnames(taxInput))
  
  # Remove ranks above species level
  taxInput <- taxInput[which(taxInput$taxonrank %in% c("SPECIES","SUBSPECIES","VARIETY","FORM")),]

  # Remove hybrid species (general cases)
  posG <- grep("(^|\\s)×", taxInput$scientificname)
  posE <- grep(" x | × ", taxInput$scientificname)
  delPos <- unique(c(posG, posE))
  
  if(length(delPos)>0)
  {
    taxInput <- taxInput[-delPos,]
  }

  # Remove hybrid species (when "X" indicate hybrid and is not an author's name first letter)
  #Note: Identify strings where " X " is followed by any name starting with lower case (indicative of epithet)
  delPos <- grep(" X .*\\b\\p{Ll}", taxInput$scientificname, perl = TRUE)
  if(length(delPos)>0)
  {
    taxInput <- taxInput[-delPos,]
  }

  # Remove unidentified species
  delPos <- grep(" sp\\. ", taxInput$scientificname)
  if(length(delPos)>0)
  {
    taxInput <- taxInput[-delPos,]
  }

  # Split species from authors names
  spAut <- t(apply(taxInput, 1, get_SciNames_Author))
  
  # Standardize encoding / remove other characters
  spAut[,1] <- stringi::stri_trans_general(spAut[,1], "Latin-ASCII")
  spAut[,1] <- stringr::str_squish(gsub("[^a-zA-Z ().-]", "", spAut[,1]))
  
  # Return data set
  taxInput <- data.frame(taxInput, spAut)
  return(taxInput)
}


# Isolate taxon name from author
get_SciNames_Author <- function(dataset)
{
  ## Split the name
  nome <- unlist(str_split(string = dataset["scientificname"], pattern = " "))
  
  # Check if second name is epithet or subgenus
  primeiro_caractere <- substr(nome[2], 1, 1)
  subgenus <- grepl("^[A-Z\\(]", primeiro_caractere)
  
  ## Check the taxon rank
  if(dataset["taxonrank"] == "SPECIES")
  {
    # String name length for species
    n=2
    
    # Longer if subgenus
    if(subgenus)
    {
      n=3
    }
  }
  else
  {
    # Extend name length if additional particle is present
    part <- ifelse(any(nome=="subsp." | nome=="var." | nome=="f."), 1, 0)
    
    # String name length for subspecies
    n=3+part
    
    # Longer if subgenus
    if(subgenus)
    {
      n=4+part
    }
  }
  
  resu <- c(
    sciName=paste(nome[1:n], collapse = " "),
    Authors=ifelse(length(nome) < (n+1), NA, paste(nome[(n+1):length(nome)], collapse = " "))
  )
  return(resu)
}



#. Search ----

# Search for names in WoRMS
wTaxaMatch <- function(taxTarget, fuzzyMatch=F, wormsTabLocal=NULL)
{
  # If a WoRMS offline data set is available
  if(!is.null(wormsTabLocal) & fuzzyMatch==F)
  {
    # Match the suggested names with the WoRMS offline data set
    pos <- match(spTrim(taxTarget$sciName, delAbv=T, delSubGen=T), spTrim(wormsTabLocal$scientificname, delAbv=T, delSubGen=T))
    if(any(!is.na(pos)))
    {
      wormsTab <- data.frame(wormsTabLocal[na.omit(pos),], SN_gbif=taxTarget$scientificname[!is.na(pos)])
    }
    else
    {
      wormsTab <- numeric()
    }
    
    # If any name was not found yet
    pos <- which(!taxTarget$scientificname %in% wormsTabLocal$SN_gbif)
    if(length(pos)>0)
    {
      # Prepare input objects for the online search
      inputNames <- taxTarget[pos,"sciName"]
      originalNames <- taxTarget[pos,"scientificname"]
    }
    else
    {
      # Everything was found offline
      return(wormsTab)
    }
  }
  else
  {
    # Prepare input and output objects for the online search
    inputNames <- taxTarget[,"sciName"]
    originalNames <- taxTarget[,"scientificname"]
    wormsTab <- numeric()
  }
  
  # Prepare loop controls
  tot <- length(inputNames)

  # Number of names to be searched at each loop
  chunkSize <- 50
  
  # Start the search...
  pb <- txtProgressBar(min = 0, max = tot, style = 3)
  while(length(inputNames)>0)
  {
    nmFind <- inputNames[1:chunkSize]
    snRef <- originalNames[1:chunkSize]
    
    # Safety control (only search names when there is communication with the API)
    if(worms_connection())
    {
      # Procedure to get error and enter the while loop
      #Note: It's hard to not find any match within 50 names, so an error might be an internal problem in the server
      #In this case, it's safer to wait a few seconds and try again
      #Note: This 'dummy_function()' doesn't exist. It's here just to simulate the first 'try-error' for tmp1 enter the loop
      tmp1 <- try(dummy_function(), silent = TRUE)
      tryMatch <- 0
      
      while(any(class(tmp1)=="try-error") & tryMatch < 5)
      {
        if(tryMatch>0){Sys.sleep(3)}
        
        # Search names
        if(fuzzyMatch)
        {
          # Notice that we use a modification of the original 'wm_records_taxamatch' to include 'extant_only'
          tmp1 <- try(wm_records_taxamatch2(name = nmFind, marine_only = F, extant_only = F), silent = T)
        }
        else
        {
          # Notice that we use a modification of the original 'wm_records_names' to include 'extant_only'
          tmp1 <- try(wm_records_names2(name = nmFind, marine_only = F, extant_only = F), silent = T)
        }
        
        tryMatch <- tryMatch+1
      }
      
      # If any result was found
      if(any(class(tmp1)!="try-error"))
      {
        if(fuzzyMatch)
        {
          posX1 <- which(lengths(tmp1)>0)
        }
        else
        {
          # Prioritize true match (replaced names can be incorporated in the fuzzy match option)
          posX1 <- which(lengths(tmp1)>0 & grepl("replaced", tmp1)==F)
        }
        
        # Get the valid results
        if(length(posX1)>0)
        {
          # Add the original scientific name
          tmp1 <- add_sciNames(tmpList = tmp1[posX1], names50 = snRef[posX1])
          
          # Ignore suggestions above species level
          #Note: Species ID is 220, but the searched name may be Aggr. and Coll. sp. (which are unassessed and excluded further)
          #Keep them is necessary to avoid searching again in a second round
          tmp1 <- tmp1[tmp1$taxonRankID > 210 | is.na(tmp1$taxonRankID), ]
          
          # Bind the result with the output data frame
          if(nrow(tmp1)>0)
          {
            wormsTab <- rbind(wormsTab, tmp1)
          }
        }
      }
      
      # Loop control
      inputNames <- inputNames[-c(1:chunkSize)]
      originalNames <- originalNames[-c(1:chunkSize)]
      
      setTxtProgressBar(pb, tot-length(inputNames))
    }
  }
  close(pb)
  
  # Return output
  return(wormsTab)
}


# Include original name as a additional column in the WoRMS output
add_sciNames <- function(tmpList, names50)
{
  # Include the names
  lista_modificada <- Map(function(df, nome)
  {
    df$SN_gbif <- nome
    return(df)
  }, tmpList, names50)
  
  # Ensure all data sets have the same columns
  if(length(unique(lengths(lista_modificada)))>1)
  {
    all_cols <- colnames(lista_modificada[[which.max(lengths(lista_modificada))]])
    
    padroniza <- function(df)
    {
      missing <- setdiff(all_cols, names(df))
      df[missing] <- NA
      
      df <- df[all_cols]
      return(df)
    }
    lista_modificada <- lapply(lista_modificada, padroniza)
  }
  
  # Convert the list to a single data frame
  lista_modificada <- do.call(rbind, lista_modificada)
  return(lista_modificada)
}


# Check user connection and WoRMS API status
worms_connection <- function(api_url = "https://www.marinespecies.org/rest", timeout = 5)
{
  # Check if the user has internet access
  if(curl::has_internet())
  {
    # Check of the server is responding
    response <- tryCatch({httr::HEAD(api_url, httr::timeout(timeout))}, error = function(e) {return(NULL)})
    
    if (!is.null(response) && httr::status_code(response) == 200)
    {
      return(T)
    }
    else
    {
      return(F)
    }
  }
}



#. Alignment ----

# Integrates all function of alignment required to filter and compare the data sets
alignData <- function(wormsTab, taxTarget, kRank, verbose = T)
{
  if(verbose){message("Filling marine habitat information")}
  wormsTab <- fill_MarHabitat(wormsTab = wormsTab)
  
  if(verbose){message("Filling taxonomic rank information")}
  wormsTab <- fill_taxRank(wormsTab = wormsTab, kRank = kRank)
  
  # if(verbose){message("Standardizing Linnaeus' name abbreviation")}
  # stTabs <- stLinnaeus(wormsTab = wormsTab, taxTarget = taxTarget)
  
  # Return output
  return(wormsTab)
}


# Update habitat
fill_MarHabitat <- function(wormsTab)
{
  # Fill in information based on other habitats
  # Logic here: if marine is NA but any other is 1, so marine is 0
  habClass <- apply(wormsTab[,c("isMarine","isBrackish","isFreshwater","isTerrestrial")], 1, sum, na.rm=T)
  habMar <- wormsTab$isMarine
  wormsTab$isMarine <- ifelse(is.na(habMar) & habClass>0, 0, habMar)
  
  habBrack <- wormsTab$isBrackish
  wormsTab$isBrackish <- ifelse(is.na(habBrack) & habClass>0, 0, habBrack)
  
  # Check the presence of persistent NAs in valid names
  posNA <- which(is.na(wormsTab$isMarine) & !is.na(wormsTab$scientificname) & !wormsTab$status %in% taxStatus("invNames") & !wormsTab$isExtinct %in% 1)
  rankLevel <- 0
  
  # Number of names to be searched at each loop
  chunkSize <- 50
  
  # Repeat this procedure at different levels while there are NAs in the data set
  while(length(posNA)>0)
  {
    # Define the rank to be searched in the inner While() loop below
    #Note: At each outer While() loop the search is done at a specific taxonomic level
    rankLevel <- rankLevel+1
    
    # Repeat this procedure while there are NAs in the current pNA vector
    while(length(posNA)>0)
    {
      # Use a temporary data frame to allow for a search in loop
      tmp1 <- wormsTab[na.omit(posNA[1:chunkSize]),]
      
      # Control de taxonomic level to search
      cont <- 0
      while(cont < rankLevel)
      {
        if(worms_connection())
        {
          tmp1 <- try(worrms::wm_record(id = tmp1$parentNameUsageID), silent = T)
          cont <- cont+1
        }
      }
      
      # Fill in information based on other habitats
      habClass <- apply(tmp1[,c("isMarine","isBrackish","isFreshwater","isTerrestrial")], 1, sum, na.rm=T)
      habMar <- tmp1$isMarine
      x <- ifelse(is.na(habMar) & habClass>0, 0, habMar)
      
      # Save the information and update the vector
      wormsTab$isMarine[na.omit(posNA[1:chunkSize])] <- x
      posNA <- posNA[-c(1:chunkSize)]
    }
    
    # After finish the past vector, check again the presence of persistent NAs
    posNA <- which(is.na(wormsTab$isMarine) & !is.na(wormsTab$scientificname) & !wormsTab$status %in% taxStatus("invNames") & !wormsTab$isExtinct %in% 1)
  }
  
  # Return output
  return(wormsTab)
}


# Status definition
taxStatus <- function(st)
{
  if(st=="accNames")
  {
    return(c("accepted","unreplaced junior homonym"))
  }
  else if(st=="invNames")
  {
    return(c("unassessed","taxon inquirendum","nomen dubium","quarantined","uncertain","deleted","nomen nudum","unavailable name","interim unpublished","temporary name"))
  }
}


# Update rank
# This is necessary if any higher rank is a criteria of inclusion
fill_taxRank <- function(wormsTab, kRank=NULL)
{
  # If no rank is provided
  if(is.null(kRank))
  {
    return(wormsTab)
  }
  
  # Required because 'wm_records_names' uses lower case
  kRank <- tolower(kRank)
  
  # Target column names in 'wm_records_names'
  recRank <- c("kingdom","phylum","class","order","family","genus")
  
  # Find for empty information in the target rank
  #Note: Records with no scientific name will be excluded in filter 1 and do not need to be filled in.
  #Other unwanted names (e.g. extinct, invalid, etc.) are still filled in because this information may be necessary before excluding them in filter 2.
  posNA <- which(is.na(wormsTab[,kRank]) & !is.na(wormsTab$scientificname))
  
  # Number of names to be searched at each loop
  chunkSize <- 100
  
  # Start the search
  while(length(posNA)>0)
  { 
    # IDs to search
    IDpos <- na.omit(posNA[1:chunkSize])
    idFind <- unique(wormsTab$AphiaID[IDpos])
    
    if(worms_connection())
    {
      # Get the classification
      tmp <- worrms::wm_classification_(id = idFind)
      tmp <- tmp[,c("id","rank","scientificname")]
      
      # Convert output from design matrix to a community matrix type
      tmp <- tmp %>% tidyr::pivot_wider(names_from = rank, values_from = scientificname)
      colnames(tmp) <- tolower(colnames(tmp))
      
      # Match rows from output with rows from searched IDs
      #Note: A match is required because search is done over unique IDs
      #This is necessary because the matrix conversion is sensible to two entries with the same ID
      rowsM <- match(wormsTab$AphiaID[IDpos], tmp$id)
      colsM <- which(colnames(tmp) %in% recRank)
      tmp <- tmp[rowsM, colsM]
      
      # Ignore ranks above kRank
      colsM <- which(colnames(tmp) %in% kRank)
      if(length(colsM)>0)
      {
        # Select the columns (from kRank to genus), match position and fill the information (if the column is available)
        tmp <- tmp[, colsM:ncol(tmp)]
        colsM <- match(colnames(tmp), colnames(wormsTab))
        
        if(length(colsM)>0)
        {
          wormsTab[IDpos,colsM] <- tmp
        }
      }
      
      # Loop control
      posNA <- posNA[-c(1:chunkSize)]
    }
  }
  
  # Return output
  return(wormsTab)
}


# Standardize Linnaeus' name in the WoRMS output
# This is necessary because many records in GBIF use the abbreviated version
# stLinnaeus <- function(wormsTab, taxTarget)
# {
#   # Match rows of the preFiltered input data set with the WoRMS output
#   pos <- match(wormsTab$SN_gbif, taxTarget$scientificname)
#   taxTarget <- taxTarget[pos,]
#   
#   # Find species with authority described as 'Linnaeus' in the WoRMS output
#   # Find species with authority described as 'Linnaeus' in both data sets
#   posLn1 <- which(grepl("Linnaeus", wormsTab$authority) & grepl("L\\.", taxTarget$Authors))
#   if(length(posLn1)>0)
#   {
#     taxTarget$Authors[posLn1] <- gsub("L.","Linnaeus" , taxTarget$Authors[posLn1])
#   }
#   
#   posLn2 <- which(grepl("L\\.", wormsTab$authority) & grepl("Linnaeus", taxTarget$Authors))
#   if(length(posLn2)>0)
#   {
#     wormsTab$authority[posLn2] <- gsub("L.","Linnaeus" , wormsTab$authority[posLn2])
#   }
#   
#   posLn3 <- which(grepl("\\(L\\.\\)", wormsTab$authority) & grepl("\\(L\\.\\)", taxTarget$Authors))
#   if(length(posLn3)>0)
#   {
#     wormsTab$authority[posLn3] <- gsub("\\(L\\.\\)","Linnaeus" , wormsTab$authority[posLn3])
#     taxTarget$Authors[posLn3] <- gsub("\\(L\\.\\)","Linnaeus" , taxTarget$Authors[posLn3])
#   }
# 
#   # posLn <- grep("Linnaeus", wormsTab$authority)
#   # 
#   # if(length(posLn)>0)
#   # {
#   #   # Start the standardization
#   #   for(i in 1:length(posLn))
#   #   {
#   #     # Change Linnaeus for its abbreviation versions
#   #     vList <- list()
#   #     vList[[1]] <- gsub("Linnaeus(, \\d{4})?","L." , wormsTab$authority[posLn[i]])
#   #     vList[[2]] <- gsub("Linnaeus","L." , wormsTab$authority[posLn[i]])
#   #     vList[[3]] <- gsub("\\([^)]*Linnaeus[^)]*\\)", "(L.)", wormsTab$authority[posLn[i]])
#   #     
#   #     # If the authority is provided in the input data set
#   #     if(!is.na(taxTarget$Authors[posLn[i]]))
#   #     {
#   #       # For each version of the name...
#   #       for(j in 1:length(vList))
#   #       {
#   #         # If the name's version matches
#   #         if(taxTarget$Authors[posLn[i]]==vList[[j]])
#   #         {
#   #           # Update Linnaeus name version in the authority column
#   #           wormsTab$authority[posLn[i]] <- vList[[j]]
#   #           break
#   #         }
#   #       }
#   #     }
#   #   }
#   # }
#   
#   # Return output
#   return(list(wormsTab=wormsTab, taxTarget=taxTarget))
# }



#. Filter 1----

# First filter (taxonomic based)
filterTaxa <- function(wormsTab, kRank, kTaxa)
{
  inDF <- wormsTab
  kTaxa <- c(kTaxa, NA)
  
  # Exclude matches above species level
  #Note: We keep quarantined names to ensure the best match when checking duplicates.
  #All quarantined and other invalid namesare excluded in Filter 2
  pos <- which(wormsTab$taxonRankID >= 220 | is.na(wormsTab$taxonRankID))
  wormsTab <- wormsTab[pos,]
  
  # Exclude suggestions from different taxa
  pos <- which(wormsTab[[kRank]] %in% kTaxa)
  wormsTab <- wormsTab[pos,]
  
  #Note: Name as in 'wm_records_names'
  #c("kingdom","phylum","class","order","family","genus")
  
  # Return output
  return(wormsTab)
}



#. Score ----

# Get scores
taxaScoring <- function(wormsTab, taxTarget)
{
  # Get systematic, authorial and orthographic scores for each suggestion
  sScore <- systScore(wormsTab = wormsTab, taxTarget = taxTarget)
  aScore <- authScore(wormsTab = wormsTab, taxTarget = taxTarget)
  oScore <- orthScore(wormsTab = wormsTab, taxTarget = taxTarget)
  
  # Merge scores to each suggestion
  wormsTab <- bind_cols(wormsTab, orthScore=oScore, authScore=aScore, systScore=sScore)
  
  # Return output
  return(wormsTab)
}


# Confidence score of the systematic axis
systScore <- function(wormsTab, taxTarget)
{
  # Match rows of the preFiltered input data set with the current WoRMS output
  pos <- match(wormsTab$SN_gbif, taxTarget$scientificname)
  taxTarget <- taxTarget[pos,]
  
  # Compare the taxonomic classification between data sets
  EqualKingdom <- wormsTab$kingdom==taxTarget$kingdom
  EqualPhylum <- wormsTab$phylum==taxTarget$phylum
  EqualClass <- wormsTab$class==taxTarget$class
  EqualOrder <- wormsTab$order==taxTarget$order
  EqualFamily <- wormsTab$family==taxTarget$family
  
  # Combine results and get the minimum rank with the same classification
  matchTaxon <- cbind(EqualFamily, EqualOrder, EqualClass, EqualPhylum, EqualKingdom, None=1)*1
  resScore <- apply(matchTaxon, 1, which.max)
  
  # Calculate the index and return results
  #Note: Score ranges from 0 (different kingdom) to 1 (same family)
  resScore <- 1-(resScore-1)/(ncol(matchTaxon)-1)
  return(resScore)
}


# Confidence score of the authorial axis
authScore <- function(wormsTab, taxTarget, YearDiff_tol=15)
{
  # Match rows of the preFiltered input data set with the current WoRMS output
  pos <- match(wormsTab$SN_gbif, taxTarget$scientificname)
  taxTarget <- taxTarget[pos,]
  
  # Get a ("highly") trimmed version of the authors' names
  #Note: This version remove isolated characters - usually initials of the names
  wAuthDate1 <- authTrim(wormsTab$authority, delNC1 = T)
  gAuthDate1 <- authTrim(taxTarget$Authors, delNC1 = T)
  
  # Compare the names and get the score (binary: equal or different)
  resScore <- (wAuthDate1==gAuthDate1)*1
  
  # # Get a ("lesser") trimmed version of the authors' names and update score
  # #Note: This is important because removing single letters as above inhibit comparisons of abbreviated authorities
  # posNA <- which(is.na(resScore)) ## NA? Yes, because L. abbreviation will return as NA and only these ones are compared here
  # if(length(posNA)>0)
  # {
  #   wAuthDate <- authTrim(wormsTab$authority[posNA])
  #   gAuthDate <- authTrim(taxTarget$Authors[posNA])
  #   resScore[posNA] <- (wAuthDate==gAuthDate)*1
  # }
  
  ## Split authors' name from publication year
  pos0 <- which(resScore %in% 0 & !is.na(wAuthDate1) & !is.na(gAuthDate1))
  if(length(pos0)>0)
  {
    ## Split
    wAuthDateS <- SplitNameYear(wAuthDate1[pos0])
    gAuthDateS <- SplitNameYear(gAuthDate1[pos0])

    ## For names
    # First check if names (isolated from year) match
    Names <- mapply(identical, wAuthDateS$names, gAuthDateS$names)*1

    # Correct typos and compare names
    name0 <- which(Names %in% 0 & !is.na(wAuthDateS$names) & !is.na(gAuthDateS$names))
    if(length(name0)>0)
    {
      tmp <- mapply(correct_typo, wAuthDateS$names[name0], gAuthDateS$names[name0], USE.NAMES = F)
      wAuthDateS$names[name0] <- strsplit(tmp[1,], " ")
      gAuthDateS$names[name0] <- strsplit(tmp[2,], " ")
      
      Names[name0] <- mapply(identical, wAuthDateS$names[name0], gAuthDateS$names[name0])*1
    }
    
    # Check if names match when ignoring order or repetition
    #Note: Sometimes an authorship is simplified, e.g. (L.) L. can be shown as L. and represent the same species 
    #However, these cases are penalized (*.95) to give more weight to cases of identical representation
    name0 <- which(Names %in% 0)
    if(length(name0)>0)
    {
      tmp <- mapply(function(a, b) setequal(unique(a), unique(b)), wAuthDateS$names[name0], gAuthDateS$names[name0])
      Names[name0] <- ifelse(tmp == 1, .95, tmp)
    }

    # Get "similarity" of the names compounding the authorities
    #Note1: If one authority is totally within the other, like an abbreviation, this is considered similar
    #However, these similarities are penalized (*.9) to give more weight to cases of truly similar authors
    #Note2: Only cases where no author is available are considered NA
    name0 <- which(Names %in% 0)
    if(length(name0)>0)
    {
      tmp <- mapply(stringAbbr, wAuthDateS$names[name0], gAuthDateS$names[name0])
      Names[name0] <- ifelse(tmp == 1, .9, tmp)
    }
    
    ## For Dates
    # First check if year (isolated from names) match
    Year <- mapply(function(a, b){any(a == b)}, wAuthDateS$year, gAuthDateS$year)*1

    # Get temporal distance using an exponential decay curve
    year0 <- which(Year %in% 0 & Names > 0 & !is.na(Names))
    if(length(year0)>0)
    {
      anos_dif <- mapply(function(a, b){min(abs(a - b))}, wAuthDateS$year[year0], gAuthDateS$year[year0])
      Year[year0] <- exp(log(.01)/YearDiff_tol * anos_dif)
    }
    
    # When authors' name mismatch or the years are not available, return 0
    #This is because 'year' is ignored when authors mismatch but not when they match
    Year[Names %in% 0 | is.na(Year)] <- 0
    
    ## Get total score
    resScore[pos0] <- ifelse(wormsTab$kingdom[pos0] %in% c("Animalia","Protozoa"),
                             rowSums(cbind(Names*0.8, Year*0.2)), #Year is only considered for those following the ICZN
                             Names)                               #For those under the Botanical code (Plants, fungi, algae) and bacterias, the year is ignored
  }
  
  # Return results
  return(round(resScore, 2))
}


# Trim authority names to improve comparison
authTrim <- function(vec, delNC1=F)
{
  # Standardize Linnaeus
  posLn <- which(grepl("\\(L\\.\\)", vec))
  if(length(posLn)>0)
  {
    vec[posLn] <- gsub("\\(L\\.\\)","Linnaeus" , vec[posLn])
  }
  
  posLn <- which(grepl("L\\.$", vec))
  if(length(posLn)>0)
  {
    vec[posLn] <- gsub("L\\.$","Linnaeus" , vec[posLn])
  }
  
  # Remove 'not authors'
  xNOT <- grep("\\b(not|non)\\b", vec)
  if(length(xNOT)>0)
  {
    vec[xNOT] <- gsub("\\b(not|non)\\b.*", "", vec[xNOT])
  }
  
  # Remove isolated characters - name initials (if required)
  vec <- mapply(strsplit,vec, '[[:punct:] ]+', USE.NAMES = F)
  if(delNC1)
  {
    vec <- lapply(vec, function(x) x[nchar(x)>1])
  }
  
  # Format coding, convert to lower case, and remove extra spaces
  vec <- lapply(vec, paste, collapse = " ")
  vec <- do.call(c, vec)
  
  vec <- stringr::str_squish(gsub("(?i)\\b(et al|in|de|sensu|ex|jr)\\b", "", vec))
  vec <- stringr::str_squish(gsub("(?i)\\b(XXX|XX|X)?(IX|IV|V?I{1,3})\\b\\.?", "", vec, perl = T))
  vec <- stringr::str_squish(tolower(stringi::stri_trans_general(vec, "Latin-ASCII")))
  vec <- ifelse(vec=="na" | vec=="", NA, vec)
  
  # Return output
  return(vec)
}


# Split the name of the authors and year of publication in the authorship
SplitNameYear <- function(x)
{
  # Extract the years (4 digits) into a single string (if more than one)
  ano <- regmatches(x, gregexpr("\\b\\d{4}\\b", x))
  #ano <- lapply(ano, function(x){as.numeric(x)})
  ano <- sapply(ano, function(a) if(length(a)==0) NA else paste(a, collapse = " "))
  ano <- lapply(strsplit(ano, " "), function(x){as.numeric(x)})
  
  # Remove years from the string
  nome <- gsub("\\b\\d{4}\\b", "", x)
  nome <- trimws(gsub("\\s{2,}", " ", nome))
  nome <- sapply(nome, function(a) if(nchar(a)==0) NA else a, USE.NAMES = F)
  nome <- strsplit(nome, " ")
  
  # Return results
  out <- list(names=nome, year=ano)
  return(out)
}


# Correct typos in the name of the authors
correct_typo <- function(string1, string2)
{
  # Split strings
  string1 <- unlist(string1)
  string2 <- unlist(string2)
  
  # Convert to list and find the shortest string
  ll <- list(string1,string2)
  if(length(unique(lengths(ll)))>1)
  {
    posMin <- which.min(sapply(ll, length))
  }
  else
  {
    posMin <- 1
  }
  
  # Find the best match and the shortest distance between the names in each list
  best_match <- sapply(ll[[posMin]], function(x) which.min(stringdist(x, ll[[-posMin]], method = "dl")))
  best_dist <- sapply(ll[[posMin]], function(x) min(stringdist(x, ll[[-posMin]], method = "dl")))
  
  if(any(best_dist>0 & best_dist < 3))
  {
    # Replace the names for the version mirrored in the other string
    nch <- nchar(names(best_match))
    posD <- which(best_dist < 3 & (best_dist/nch) < .21)
    
    if(length(posD)>0)
    {
      ll[[posMin]][posD] <- ll[[-posMin]][best_match[posD]]
    }
  }
  
  # Return list
  out <- sapply(ll, function(v) paste(unlist(v), collapse = " "))
  return(out)
}


# Get Damerau-Levenshtein distance between strings
calcDists <- function(vecW, vecG)
{
  vecW <- unlist(strsplit(vecW, " "))
  vecG <- unlist(strsplit(vecG, " "))
  
  # Remove subgenus to compare only genus and epithet
  if(any(grepl("\\(", vecW))){vecW <- vecW[-c(grep("\\(", vecW))]}
  if(any(grepl("\\(", vecG))){vecG <- vecG[-c(grep("\\(", vecG))]}
  
  # Ignore cases where names have different lengths
  if(length(vecW) != length(vecG))
  {
    resu <- c(10, 10)
    names(resu) <- c("dlGe", "dlEp")
    return(resu)
  }
  
  # Levenshtein distance
  # Mede quantas operacoes (insercao, delecao, substituicao) sao necessarias para transformar uma palavra na outra.
  # distGe <- adist(vecW[1], vecG[1])[1,1]
  # distEp <- adist(vecW[2], vecG[2])[1,1]
  
  # Damerau-Levenshtein distance
  # Conta o numero minimo de operacoes para transformar uma palavra em outra
  dlGe <- stringdist::stringdist(vecW[1], vecG[1], method = "dl")
  dlEp <- stringdist::stringdist(paste(vecW[2:length(vecW)], collapse=""), paste(vecG[2:length(vecG)], collapse=""), method = "dl")
  
  # Return output
  resu <- c(dlGe, dlEp)
  names(resu) <- c("dlGe", "dlEp")
  
  return(round(resu, 2))
}


# Get the proportion of shared elements between strings (publication info)
stringAbbr <- function(string1, string2, range=F)
{
  # Split strings
  vW <- unlist(string1)
  vG <- unlist(string2)
  
  # Get the number of matches, considering abbreviations
  matches1 <- sapply(vG, function(m) any(sapply(vW, function(n) is_abbreviation(m, n) | is_abbreviation(n, m))))
  matches2 <- sapply(vW, function(m) any(sapply(vG, function(n) is_abbreviation(m, n) | is_abbreviation(n, m))))
  
  # Calculate the proportion of shared elements
  r1 <- sum(matches1) / length(vG)
  r2 <- sum(matches2) / length(vW)
  
  # Provide the maximum number of shared elements (considering both directions)
  if(range==F)
  {
    r <- max(c(r1, r2))
  }
  else
  {
    # Or provide the number of shared elements for each direction
    r <- c(r1, r2)
  }
  
  # Return output
  return(r)
}


is_abbreviation <- function(abbr, words)
{
  # Verify if any item (word or part of the word) of one strong appears in the other
  any(grepl(abbr, words))
}


# Confidence score of the orthographic axis
orthScore <- function(wormsTab, taxTarget)
{
  # Match rows of the preFiltered input data set with the current WoRMS output
  pos <- match(wormsTab$SN_gbif, taxTarget$scientificname)
  taxTarget <- taxTarget[pos,]
  
  # First comparison (direct comparison)
  resScore <- (spTrim(wormsTab$scientificname) == spTrim(taxTarget$sciName))*1
  
  # Second comparison (excluding infraspecific rank indicators or subgenus)
  #Note: This comparison is penalized (*.95 or *.90) to give more weight to truly equal names
  pos0 <- which(resScore == 0)
  if(length(pos0)>0)
  {
    x1 <- (spTrim(wormsTab$scientificname[pos0] , delSubGen = T) == spTrim(taxTarget$sciName[pos0], delSubGen = T))*.95
    x2 <- (spTrim(wormsTab$scientificname[pos0], delAbv = T) == spTrim(taxTarget$sciName[pos0], delAbv = T))*.95
    x3 <- (spTrim(wormsTab$scientificname[pos0], delAbv = T, delSubGen = T) == spTrim(taxTarget$sciName[pos0], delAbv = T, delSubGen = T))*.9 #Note: More penalized because has two alterations
    tmp <- cbind(x1, x2, x3)
    resScore[pos0] <- apply(tmp, 1, max)
  }
  
  # Third comparison (considering canonical roots, i.e. excluding suffix)
  #Note: This comparison is penalized (*.95) to give more weight to truly equal names
  pos0 <- which(resScore == 0)
  if(length(pos0)>0)
  {
    wNameStem <- mapply(canonicalName, spTrim(wormsTab$scientificname[pos0], delAbv = T))
    gNameStem <- mapply(canonicalName, spTrim(taxTarget$sciName[pos0], delAbv = T))
    resScore[pos0] <- (wNameStem==gNameStem)*.95
  }
  
  # Fourth comparison (considering the accepted names at species level)
  #Note: This comparison is penalized (*.90) to give more weight to comparisons between original names
  pos0 <- which(resScore == 0 & !is.na(wormsTab$valid_name))
  if(length(pos0)>0)
  {
    # Get the rank level of the valid names
    spW_val <- spTrim(wormsTab$valid_name[pos0], delAbv = T, delSubGen = T)
    posSub <- which(str_count(spW_val, "\\b\\w+\\b") > 2)
    
    # Get the valid species name for those below species level
    spW_val <- wormsTab$valid_name[pos0]
    if(length(posSub)>0)
    {
      tmp <- getValid_SP(wormsTab = wormsTab[pos0[posSub],])
      spW_val[posSub] <- tmp$scientificname
    }
    
    # Compare the names and update the score
    x1 <- (spTrim(spW_val) == spTrim(taxTarget$species[pos0]))*.90
    x2 <- (spTrim(spW_val, delSubGen = T) == spTrim(taxTarget$species[pos0], delSubGen = T))*.85
    
    resScore[pos0] <- apply(cbind(x1, x2), 1, max)
  }
  
  # # Fifth comparison (calculate similarity in strings composition)
  # #Note: This is necessary because a simple var. in GBIF may be a var. of a subsp. in WoRMS (i.e. it includes one more name)
  # pos0 <- which(resScore %in% 0 & wormsTab$taxonRankID > 220 & !taxTarget$taxonrank %in% "SPECIES")
  # if(length(pos0)>0)
  # {
  #   resScore[pos0] <- mapply(stringAbbr, spTrim(wormsTab$scientificname[pos0], delAbv=T), spTrim(taxTarget$sciName[pos0], delAbv=T))*.9
  # }
  
  # Sixth comparison (calculate orthographic similarity)
  #Note: Get similarity considering the number of typo at genus and epithet separately
  #The weight for differences in the genus is higher than for differences in the epithet
  #pos0 <- which(resScore %in% 0 & wormsTab$taxonRankID %in% 220 & taxTarget$taxonrank %in% "SPECIES")
  pos0 <- which(resScore %in% 0)
  if(length(pos0)>0)
  {
    distis <- t(mapply(calcDists, spTrim(wormsTab$scientificname[pos0], delAbv = T), spTrim(taxTarget$sciName[pos0], delAbv = T)))
    
    # simGen <- 1-(distis[,"dlGe"]/3)
    # simGen <- ifelse(simGen < 0, 0, simGen)
    # 
    # simEpi <- 1-(distis[,"dlEp"]/4)
    # simEpi <- ifelse(simEpi < 0, 0, simEpi)
    # 
    # resScore[pos0] <- simGen*simEpi
    
    # Three edits within the same component are not allowed
    distis[distis[,1]==3,1] <- 10
    distis[distis[,2]==3,2] <- 10
    x1 <- rowSums(distis)
    resScore[pos0] <- c(.75,.5,.25,0)[cut(x1, breaks=c(0,1,2,3,20))]
  }
  
  # Return results
  return(resScore)
}


# Get the canonical version of the name (remove Latin suffix)
canonicalName <- function(sciName)
{
  # Split names into genus and epithet
  fullNam <- unlist(strsplit(sciName, " "))
  posGen <- ifelse(grepl("\\(", fullNam[2]), 2, 1)
  
  namGen <- fullNam[1:posGen]
  namEpi <- fullNam[(posGen+1):length(fullNam)]
  
  # Stem epithet and infraspecific names
  stemEpi <- mapply(latStemming, namEpi)
  
  # Reassemble the name
  canName <- paste(paste(namGen,collapse = " "), paste(stemEpi,collapse = " "), collapse = " ")
  
  # Return output
  return(canName)
}


# Stemming algorithm based on Schinke et al (1996)
latStemming <- function(latName)
{
  # step 1
  nameString <- latName
  
  # step 2
  nameString <- gsub(pattern = "j", replacement = "i", x = nameString)
  nameString <- gsub(pattern = "v", replacement = "u", x = nameString)
  
  # step 3
  que <- c("atque","quoque","neque","itaque","absque","apsque","abusque","adaeque","adusque","denique",
           "deque","susque","oblique","peraeque","plenisque","quandoque","quisque","quaeque",
           "cuiusque","cuique","quemque","quamque","quaque","quique","quorumque","quarumque",
           "quibusque","quosque","quasque","quotusquisque","quousque","ubique","undique","usque",
           "uterque","utique","utroque","utribique","torque","coque","concoque","contorque",
           "detorque","decoque","excoque","extorque","obtorque","optorque","retorque","recoque",
           "attorque","incoque","intorque","praetorque")
  
  if(endsWith(nameString, "que"))
  {
    if(nameString %in% que)
    {
      return(latName)
    }
    else
    {
      nameString <- gsub('.{3}$', '', nameString)
    }
  }
  
  # step 4
  noumSuf <- c("ibus", "ius", "ae", "am", "as", "em", "es", "ia", "is",
               "nt", "os", "ud", "um", "us", "a", "e", "ii", "i", "o", "u") #ii
  
  posSuf <- which(endsWith(nameString, noumSuf))
  
  if(length(posSuf)>0)
  {
    suffs <- noumSuf[posSuf]
    nchars <- nchar(suffs)
    posMax <- which.max(nchars)
    nameString <- gsub(paste('.{',nchars[posMax],'}$',sep=''), '', nameString)
  }
  
  if(nchar(nameString)>2)
  {
    return(nameString)
  }
  else
  {
    return(latName)
  }
  
}


# Trim scinames to improve comparison
spTrim <- function(vec, delAbv=F, delSubGen=F)
{
  # Convert Latin characters, ignore symbols and convert to lower case
  vec <- stringi::stri_trans_general(vec, "Latin-ASCII")
  vec <- stringr::str_squish(gsub("[^a-zA-Z ()]", "", vec))
  vec <- tolower(vec)
  
  # Remove subgenus (if required)
  if(delSubGen)
  {
    vec <- stringr::str_squish(gsub("\\(.*)", "", vec))
  }
  
  # Remove infraspecific rank abbreviation (if required)
  if(delAbv)
  {
    vec <- stringr::str_squish(gsub(pattern = "\\b(var|subvar|subsp|f)\\b", "", vec))
  }
  
  # Return output
  return(vec)
}



#. Filter ----

# Filter suggestions
filterNames <- function(wormsTab, taxTarget)
{
  # Create an temporary index for each row
  wormsTab <- wormsTab %>% mutate(tmp_ID = row_number())
  
  #Note: Duplicated suggestions == more than one WoRMS suggestion for the same GBIF name
  
  # First step: remove duplicates with lower total score
  dups <- unique(wormsTab$SN_gbif[which(duplicated(wormsTab$SN_gbif))])
  if(length(dups)>0)
  {
    delPos <- wormsTab %>%
      filter(SN_gbif %in% dups) %>%
      group_by(SN_gbif) %>%
      mutate(
        scrs = rowSums(cbind(orthScore, systScore, authScore), na.rm = TRUE),
        max_score = max(scrs)
      ) %>%
      ungroup()
    
    # Remove duplicates
    delPos <- delPos[which(delPos$scrs < delPos$max_score),]
    wormsTab <- wormsTab[which(!wormsTab$tmp_ID %in% delPos$tmp_ID),]
  }
  
  
  # Second step: remove duplicates entirely constituted of unwanted species
  dups <- unique(wormsTab$SN_gbif[which(duplicated(wormsTab$SN_gbif))])
  if(length(dups)>0)
  {
    delPos <- wormsTab %>%
      filter(SN_gbif %in% dups) %>%
      group_by(SN_gbif) %>%
      filter(
        (all(isMarine == 0) & all(!isBrackish %in% 1)) | all(isExtinct == 1) | all(status %in% taxStatus("invNames"))
      ) %>%
      ungroup()
    
    wormsTab <- wormsTab[which(!wormsTab$tmp_ID %in% delPos$tmp_ID),]
  }
  
  
  # Third step: for those indicating the same valid species, keep the first accepted or the first non invalid name
  dups <- unique(wormsTab$SN_gbif[which(duplicated(wormsTab$SN_gbif))])
  if(length(dups)>0)
  {
    delPos <- wormsTab %>%
      filter(SN_gbif %in% dups) %>%
      group_by(SN_gbif) %>%
      mutate(
        # Check for single valid species name
        valNam1 = if_else(
          str_detect(word(valid_name, 2, 2), "^\\(.*\\)$"),
          word(valid_name, 1, 3, sep = " "),
          word(valid_name, 1, 2, sep = " ")),
        all_names_equal = n_distinct(valNam1) == 1
      ) %>%
      # Filter the duplicates with single valid species name
      filter(all_names_equal) %>%
      mutate(
        # Get the first accepted or non invalid suggestion
        is_first = row_number() == if (any(status %in% taxStatus("accNames"))) {
          which(status %in% taxStatus("accNames"))[1]  # Linha selecionada como TRUE
        } else {
          which(!status %in% taxStatus("invNames"))[1]
        }
      ) %>%
      ungroup()
    
    delPos <- delPos[which(!delPos$is_first %in% T),]
    wormsTab <- wormsTab[which(!wormsTab$tmp_ID %in% delPos$tmp_ID),]
  }
  
  
  # Fourth step: remove the less plausible duplicate (at this step duplicates indicate different valid species)
  #When scientific name is the same and authorities have score > 0, select the most similar authority
  #When authorities are equally similar, select the most recent suggestion
  #When equally recent, discard all (impossible to distinguish the most plausible suggestion)
  dups <- unique(wormsTab$SN_gbif[which(duplicated(wormsTab$SN_gbif))])
  if(length(dups)>0)
  {
    # Select duplicates with the same scientific name (irrespective of differences in the valid name)
    delPos <- wormsTab %>%
      filter(SN_gbif %in% dups) %>%
      group_by(SN_gbif) %>%
      mutate(
        # Check for single scientific name
        all_names_equal = n_distinct(scientificname) == 1
      ) %>%
      # Filter the duplicates with single scientific name
      filter(all_names_equal) %>%
      # Select duplicates having any similarity in authority
      #Note: Any other duplicates (i.e. with unknown or totally distinct authors) will be entirely excluded
      filter(all(!is.na(authScore)) & all(authScore > 0))
    
    # If the distinction is not possible, remove all
    if(nrow(delPos)==0)
    {
      delPos <- wormsTab %>%
        filter(SN_gbif %in% dups)
      
      wormsTab <- wormsTab[which(!wormsTab$tmp_ID %in% delPos$tmp_ID),]
    }
    else
    {
      # Match rows of the preFiltered input data set with the selected WoRMS output
      pos <- match(delPos$SN_gbif, taxTarget$scientificname)
      taxTarget <- taxTarget[pos,]
      
      # Get the sum of the similarity range (assessed in both directions and not only the highest value as in the score evaluation)
      authSim <- colSums(mapply(stringAbbr, authTrim(delPos$authority, delNC1 = T), authTrim(taxTarget$Authors, delNC1 = T), range=T))
      delPos <- bind_cols(delPos, authSim=authSim)
      
      delPos <- delPos %>%
        mutate(
          # Get the most similar authority as a tiebreaker
          posKeep = authSim == max(authSim),
          posKeep = if(sum(posKeep)>1)
          {
            # If more than 1, get that with the most recent modification (if possible)
            a <- as.Date(substr(modified, 1, 10), '%Y-%m-%d')
            posKeep <- posKeep & (a == max(a))
            if(sum(posKeep) == 1)
            {
              posKeep[a < max(a)] <- FALSE
            }
            else
            {
              posKeep <- rep(FALSE, length(posKeep))
            }
            posKeep
          }
          else
          {
            posKeep
          }
        ) %>%
        ungroup()
      
      delPos <- setdiff(wormsTab$tmp_ID[wormsTab$SN_gbif %in% dups], delPos$tmp_ID[delPos$posKeep])
      wormsTab <- wormsTab[which(!wormsTab$tmp_ID %in% delPos),]
    }
  }
  
  
  # Fifth step: remove undesired names (invalid, extinct or non-marine)
  #Fill in marine habitat information first (for those species that might be missing it)
  wormsTab <- fill_MarHabitat(wormsTab = wormsTab)
  
  delPos <- wormsTab %>%
    filter((isMarine == 0 & !isBrackish %in% 1) | isExtinct == 1 | status %in% taxStatus("invNames"))
  
  if(nrow(delPos)>0)
  {
    wormsTab <- wormsTab[which(!wormsTab$tmp_ID %in% delPos$tmp_ID),]
  }
  
  
  # Return output
  wormsTab <- wormsTab %>% select(-tmp_ID)
  return(wormsTab)
}



#. Classification ----

# Classification of the taxonomic match based on the three scores
confMatch <- function(wormsTab)
{
  # Create the output vector
  flag <- rep(NA, nrow(wormsTab))
  
  # High confidence in (all axes high score)
  HighIn <- which(wormsTab$systScore==1 & wormsTab$orthScore>=.9 & wormsTab$authScore>=.9)
  flag[HighIn] <- "Keep [H-CL]"
  
  # High confidence out (all axes no score)
  HighOut <- which(wormsTab$systScore==0 & wormsTab$orthScore==0 & (wormsTab$authScore==0 | is.na(wormsTab$authScore)))
  flag[HighOut] <- "Drop [H-CL]"
  
  # Medium confidence in (only two axes high score)
  highNamAut <- which(wormsTab$systScore<1 & wormsTab$orthScore>=.9 & wormsTab$authScore>=.9)
  highSysNam <- which(wormsTab$systScore==1 & wormsTab$orthScore>=.9 & (wormsTab$authScore<.9 | is.na(wormsTab$authScore)))
  highSysAut <- which(wormsTab$systScore==1 & (wormsTab$orthScore > 0 & wormsTab$orthScore<.9) & wormsTab$authScore>=.9)
  flag[unique(c(highNamAut, highSysNam, highSysAut))] <- "Keep [M-CL]"
  
  # Medium confidence out (only two axes no score)
  tmp <- cbind(wormsTab$systScore==0, wormsTab$orthScore==0, wormsTab$authScore==0)
  mOut <- which(rowSums(tmp, na.rm = T)>=2 & is.na(flag))
  flag[mOut] <- "Drop [M-CL]"
  
  # Low confidence in (~ upper half of the scores)
  #Note: For this case is required a high orthographic score OR a medium high orthographic score (typo) associated with any authorial score
  tmp <- cbind(wormsTab$systScore, wormsTab$orthScore, wormsTab$authScore)
  lIn <- which(rowSums(tmp, na.rm = T)>1.6 & is.na(flag))
  flag[lIn] <- "Keep [L-CL]"
  
  # Low confidence out (~ lower half of the scores)
  flag[is.na(flag)] <- "Drop [L-CL]"
  
  # Return output
  wormsTab <- bind_cols(wormsTab, confMatch=flag)
  return(wormsTab)
}



#. Correction ----

# Update information and provide final WoRMS table
taxaRefresh <- function(wormsTab, taxTarget, verbose=F)
{
  # Get originally suggested names
  origSugg <- paste(wormsTab$scientificname, wormsTab$authority, sep=" ")
  wormsTab <- data.frame(wormsTab[,1:29],origSugg,wormsTab[,30:33])
  
  # Get valid names at species level
  wormsTab <- getValid_SP(wormsTab = wormsTab, verbose = verbose)
  
  # Exclude unwanted species and revise classification 
  wormsTab <- reviseData(wormsTab = wormsTab, taxTarget = taxTarget)
  
  ## Apply revision separately for In and Out???
  
  # Return output
  return(wormsTab)
}


# Update suggested names to accepted version at species level
getValid_SP <- function(wormsTab, verbose=F)
{
  ## Find names that do not need to be verified
  posOK <- which((wormsTab$status %in% taxStatus("accNames") & wormsTab$taxonRankID == 220) | wormsTab$status %in% taxStatus("invNames"))
  
  if(length(posOK)==nrow(wormsTab))
  {
    return(wormsTab)
  }
  
  # Number of IDs to be searched at each loop
  chunkSize <- 50
  
  ## Start checking for the valid name
  
  # Lets check for the indicated IDs within the own data set before searching online
  cc <- !wormsTab$status %in% taxStatus("accNames") & !wormsTab$status %in% taxStatus("invNames") & !is.na(wormsTab$valid_AphiaID)
  pos <- match(wormsTab$valid_AphiaID[cc], wormsTab$AphiaID)
  if(sum(!is.na(pos))>0)
  {
    wormsTab[which(cc)[!is.na(pos)], 1:27] <- wormsTab[na.omit(pos), 1:27]
  }
  
  # Get remaining unaccepted names
  idUnacc <- unique(wormsTab$valid_AphiaID[which(!wormsTab$status %in% taxStatus("accNames") & !wormsTab$status %in% taxStatus("invNames") & !is.na(wormsTab$valid_AphiaID))])
  
  # Show message and progress bar
  if(verbose & length(idUnacc)>0)
  {
    message("Converting WoRMS suggestions to accepted names")
    tot <- length(idUnacc)
    pb <- txtProgressBar(min = 0, max = length(idUnacc), style = 3)
  }
  
  # Search online those not found in the given data set
  while(length(idUnacc)>0)
  {
    if(worms_connection())
    {
      tmp <- try(worrms::wm_record(id = na.omit(idUnacc[1:chunkSize])), silent = T)
      
      # Update the WoRMS table
      pos <- match(wormsTab$valid_AphiaID, tmp$AphiaID)
      wormsTab[!is.na(pos), 1:ncol(tmp)] <- tmp[na.omit(pos),]
      
      # Check for "valid names" that are quarantined (not returned by wm_record)
      posNA <- which(!na.omit(idUnacc[1:chunkSize]) %in% tmp$AphiaID)
      if(length(posNA)>0)
      {
        pos <- match(wormsTab$valid_AphiaID, idUnacc[posNA])
        wormsTab[!is.na(pos), 1:ncol(tmp)] <- NA
      }
      
      idUnacc <- idUnacc[-c(1:chunkSize)]
    }
    if(verbose){setTxtProgressBar(pb, tot-length(idUnacc))}
  }
  if(verbose & exists("pb")){close(pb)}
  
  
  ## Now update records to species level (and confirm species validity below)
  
  # Lets check for the indicated IDs within the own data set before searching online (only for accepted infraspecific ranks)
  cc <- wormsTab$status %in% taxStatus("accNames")
  pos <- match(wormsTab$parentNameUsageID[cc], wormsTab$AphiaID)
  if(sum(!is.na(pos))>0)
  {
    wormsTab[which(cc)[!is.na(pos)], 1:27] <- wormsTab[na.omit(pos), 1:27]
  }
  
  # Show message and progress bar (controlled only by species level checking)
  if(verbose)
  {
    message("Converting WoRMS suggestions to species level")
    tot <- length(which(wormsTab$taxonRankID > 220  & !is.na(wormsTab$valid_AphiaID)))
    if(tot>0){pb <- txtProgressBar(min = 0, max = tot, style = 3)}
  }
  
  # The procedure below will check in loop until all conditions are met (i.e. accepted species level name)
  # It will search for the valid species name for those valid subspecies or unnacepted species
  
  cont = 0
  while(any((!wormsTab$status %in% taxStatus("accNames") & !wormsTab$status %in% taxStatus("invNames") & !is.na(wormsTab$valid_AphiaID)) | wormsTab$taxonRankID > 220)) #& !wormsTab$status %in% "alternative representation"
  {
    # Which suggestion is not accepted yet (checking again for consistency)
    #Note: "alternative representation" are not excluded from the search because sometimes it indicates a valid accepted name (though it may also just indicate a infraspecies rank)
    #If it indicates a infraspecific rank, the species name is recovered in the last search below
    idUnacc <- unique(wormsTab$valid_AphiaID[which(!wormsTab$status %in% taxStatus("accNames") & !wormsTab$status %in% taxStatus("invNames") & !is.na(wormsTab$valid_AphiaID))])
    
    while(length(idUnacc)>0)
    {
      if(worms_connection())
      {
        tmp <- try(worrms::wm_record(id = na.omit(idUnacc[1:chunkSize])), silent = T)
        
        # Check for "valid names" that are quarantined (not returned by wm_record)
        posNA <- which(!na.omit(idUnacc[1:chunkSize]) %in% tmp$AphiaID)
        if(length(posNA)>0)
        {
          pos <- match(wormsTab$valid_AphiaID, idUnacc[posNA])
          wormsTab[!is.na(pos), 1:ncol(tmp)] <- NA
        }
        
        # Update the WoRMS table
        pos <- match(wormsTab$valid_AphiaID, tmp$AphiaID)
        
        # Ignore cases where the name suggested belongs to a lower taxon rank
        delPos <- which(tmp$taxonRankID[na.omit(pos)] > wormsTab$taxonRankID[!is.na(pos)])
        if(length(delPos)>0)
        {
          tmp <- tmp[-na.omit(pos)[delPos], ]
          pos <- match(wormsTab$valid_AphiaID, tmp$AphiaID)
        }
        
        if(nrow(tmp)>0)
        {
          wormsTab[!is.na(pos), 1:ncol(tmp)] <- tmp[na.omit(pos),]
        }
        
        idUnacc <- idUnacc[-c(1:chunkSize)]
      }
    }
    
    # Which record is below species level (only check for valid names)
    idSubsp <- unique(wormsTab$parentNameUsageID[which(wormsTab$taxonRankID > 220 & (wormsTab$status %in% taxStatus("accNames") | wormsTab$status %in% "alternative representation") )])
    
    while(length(idSubsp)>0)
    {
      if(worms_connection())
      {
        tmp <- try(worrms::wm_record(id = na.omit(idSubsp[1:chunkSize])), silent = T)
        
        # Update the WoRMS table
        pos <- match(wormsTab$parentNameUsageID, tmp$AphiaID)
        wormsTab[!is.na(pos), 1:ncol(tmp)] <- tmp[na.omit(pos),]
        
        # Check for "valid names" that are quarantined (not returned by wm_record)
        posNA <- which(!na.omit(idSubsp[1:chunkSize]) %in% tmp$AphiaID)
        if(length(posNA)>0)
        {
          pos <- match(wormsTab$parentNameUsageID, idSubsp[posNA])
          wormsTab[!is.na(pos), 1:ncol(tmp)] <- NA
        }
        
        idSubsp <- idSubsp[-c(1:chunkSize)]
      }
      if(verbose & exists("pb")){setTxtProgressBar(pb, tot-length(idSubsp))}
    }
    
    
    # Safety control to avoid infinite loop (e.g. unaccepted species suggesting the subspecies as accepted)
    cont <- cont+1
    if(cont==5)
    {
      break
    }
  }
  if(verbose & exists("pb")){close(pb)}
  
  # Return output
  return(wormsTab)
}


# Final check on the data information
reviseData <- function(wormsTab, taxTarget)
{
  # Fill in marine habitat information (for those newly informed valid species that might be missing it)
  wormsTab <- fill_MarHabitat(wormsTab = wormsTab)
  
  # Remove unwanted names (that might have been included during the names update)
  delPos <- wormsTab %>%
    filter( (isMarine == 0 & !isBrackish %in% 1) | (isExtinct == 1) | (status %in% taxStatus("invNames") | (is.na(valid_AphiaID)) | (taxonRankID < 220)) )
  wormsTab <- wormsTab[!wormsTab$AphiaID %in% delPos$AphiaID,]
  
  # Standardize name status
  pos <- which(wormsTab$status %in% taxStatus("accNames"))
  wormsTab$status[pos] <- "accepted"
  
  pos <- which(!wormsTab$status %in% taxStatus("accNames") & !wormsTab$status %in% "alternative representation")
  wormsTab$status[pos] <- "unaccepted"
  
  # Revise taxonomic score
  sScore <- systScore(wormsTab = wormsTab, taxTarget = taxTarget)
  
  posUp <- which(sScore > wormsTab$systScore)
  wormsTab$systScore[posUp] <- sScore[posUp]
  
  # Revise classification
  tmp <- confMatch(wormsTab = wormsTab[posUp, colnames(wormsTab)!="confMatch"])
  wormsTab$confMatch[posUp] <- tmp$confMatch
  
  # Return output
  return(wormsTab)
}



#. Other functions ----

# Wrapping function integrating the entire workflow
zenMatch <- function(taxInput, wormsTabLocal=NULL, fuzzyMatch = T, kRank, kTaxa, verbose = T)
{
  # Preparation
  if(verbose){msg<-"Preparing input data\n"; cat(msg); cat(strrep("-", nchar(msg)), "\n")}
  targetNames <- preFilter(taxInput = taxInput)
  
  # Search
  if(verbose){cat("\n\n"); msg<-"Searching exact names on WoRMS\n"; cat(msg); cat(strrep("-", nchar(msg)), "\n")}
  x <- wTaxaMatch(taxTarget = targetNames, wormsTabLocal = wormsTabLocal)
  
  # Alignment
  if(verbose){cat("\n\n"); msg<-"Aligning information\n"; cat(msg); cat(strrep("-", nchar(msg)), "\n")}
  x2 <- alignData(wormsTab = x, taxTarget = targetNames, kRank = kRank, verbose = F)
  
  # Filter 1
  if(verbose){cat("\n\n"); msg<-"Filtering taxa\n"; cat(msg); cat(strrep("-", nchar(msg)), "\n")}
  x3 <- filterTaxa(wormsTab = x2, kRank = kRank, kTaxa = kTaxa)
  
  # Scoring
  if(verbose){cat("\n\n"); msg<-"Scoring WoRMS suggestions\n"; cat(msg); cat(strrep("-", nchar(msg)), "\n")}
  x4 <- taxaScoring(wormsTab = x3, taxTarget = targetNames)
  
  
  #--- Optional loop adding a fuzzy match search
  if(fuzzyMatch)
  {
    # Define new search
    f1 <- setdiff(targetNames$scientificname, x3$SN_gbif)
    f2 <- x4$SN_gbif[which(x4$orthScore<1 | x4$authScore<1 | x4$systScore<1)]
    buscar <- unique(c(f1,f2))
    
    # But avoid cases where an alternative perfect suggestion is already present
    goodMatches <- unique(x4$SN_gbif[which(x4$orthScore==1 & x4$authScore==1 & x4$systScore==1)])
    buscar <- setdiff(buscar, goodMatches)
    
    # Names to search
    buscar <- targetNames[targetNames$scientificname %in% buscar,]

    if(nrow(buscar)>0)
    {
      # Search B
      if(verbose){cat("\n\n"); msg<-"Applying fuzzy search of the names on WoRMS\n"; cat(msg); cat(strrep("-", nchar(msg)), "\n")}
      xb <- wTaxaMatch(taxTarget = buscar, fuzzyMatch = T)
      
      # Alignment B
      x2b <- alignData(wormsTab = xb, taxTarget = targetNames, kRank = kRank, verbose = F)
      
      # Filter 1 B
      x3b <- filterTaxa(wormsTab = x2b, kRank = kRank, kTaxa = kTaxa)
      
      # Scoring B
      x4b <- taxaScoring(wormsTab = x3b, taxTarget = targetNames)
      
      x4 <- rbind(x4, x4b)
      x4 <- x4 %>% distinct(AphiaID, SN_gbif, .keep_all = T)
    }
  }
  #--- End optional loop
  
  
  # Filter 2
  if(verbose){cat("\n\n"); msg<-"Filtering names\n"; cat(msg); cat(strrep("-", nchar(msg)), "\n")}
  x5 <- filterNames(wormsTab = x4, taxTarget = targetNames)
  
  # Classification
  if(verbose){cat("\n\n"); msg<-"Estimating suggestions confidence\n"; cat(msg); cat(strrep("-", nchar(msg)), "\n")}
  x6 <- confMatch(wormsTab = x5)
  
  # Correction
  if(verbose){cat("\n\n"); msg<-"Updating suggestions to valid species\n"; cat(msg); cat(strrep("-", nchar(msg)), "\n")}
  x7 <- taxaRefresh(wormsTab = x6, taxTarget = targetNames, verbose = verbose)
  
  
  #--- Summary of the match
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
  
  
  #--- Summary of the search and additional outputs
  # Searched
  NM_searched <- unique(targetNames$scientificname)
  
  # Found
  NM_found <- unique(x4$SN_gbif)
  
  # Found and kept
  NM_kept <- unique(x7$SN_gbif)
  
  
  # Additional outputs
  NM_notFound <- setdiff(NM_searched, NM_found)
  
  NM_Excluded <- setdiff(NM_found, NM_kept)
  WoRMS_Excluded <- x4[which(x4$SN_gbif %in% NM_Excluded),]
  
  {
    cat("\n\n"); msg<-"Searching summary:\n"; cat(msg); cat(strrep("-", nchar(msg)+10), "\n")
    
    cat(paste("Names classified: ", length(NM_kept), sep=""), "\n")
    cat(paste("Names excluded:   ", length(NM_Excluded), sep=""), "\n")
    cat(paste("Names not found:  ", length(NM_notFound), sep=""), "\n")
    cat(strrep("-", nchar(msg)+10), "\n")
  }
  
  # Return output
  resu <- list(WoRMS_Matching=x7, WoRMS_Excluded=WoRMS_Excluded, NM_notFound=NM_notFound)
  return(resu)
}


# Standardize local WoRMS data set (offline version) to speed up the search
reshape_DCworms <- function(taxon, speciesProfile)
{
  # Match the two data sets (if necessary)
  zz <- inner_join(taxon, speciesProfile, by = "taxonID")
  
  # Select columns and rows
  zz <- data.frame(zz[,c("scientificNameID","references","scientificName","scientificNameAuthorship","taxonomicStatus")], unacceptreason=NA, taxonRankID=NA, zz[,c("taxonRank","acceptedNameUsageID","acceptedNameUsage")], valid_authority=NA, parentNameUsageID=zz[["parentNameUsageID"]], originalNameUsageID=NA, zz[,c("kingdom","phylum","class","order","family","genus","bibliographicCitation","scientificNameID")], zz[,c("isMarine","isBrackish","isFreshwater","isTerrestrial","isExtinct")], match_type="exact", modified=zz[,c("modified")])
  zz <- zz[!is.na(taxon$specificEpithet),]
  zz <- zz %>% mutate(taxonRankID = recode(taxonRank, "Species"=220, "Subspecies"=230, "Variety"=240, "Subvariety"=250, "Forma"=260, "Subforma"=270, "Mutatio"=280, "Natio"=235))
  
  # Get numeric version of the AphiaIDs
  zz$scientificNameID <- as.numeric(gsub("urn:lsid:marinespecies.org:taxname:","",zz$scientificNameID))
  zz$acceptedNameUsageID <- as.numeric(gsub("urn:lsid:marinespecies.org:taxname:","",zz$acceptedNameUsageID))
  zz$parentNameUsageID <- as.numeric(gsub("urn:lsid:marinespecies.org:taxname:","",zz$parentNameUsageID))
  
  # Remove html expressions
  zz$scientificNameAuthorship <- gsub("<[^>]+>", "", zz$scientificNameAuthorship)
  
  # Standardize column names for the style obtained when using the worrms package
  colnames(zz) <- c("AphiaID","url","scientificname","authority","status","unacceptreason","taxonRankID","rank","valid_AphiaID","valid_name","valid_authority","parentNameUsageID","originalNameUsageID","kingdom","phylum","class","order","family","genus","citation","lsid","isMarine","isBrackish","isFreshwater","isTerrestrial","isExtinct","match_type","modified")
  
  # Fill in NAs in valid_AphiaID (offline search)
  posNA <- which(is.na(zz$valid_AphiaID) & (zz$valid_name==zz$scientificname))
  if(length(posNA)>0)
  {
    zz$valid_AphiaID[posNA] <- zz$AphiaID[posNA]
  }
  
  # Fill in NAs in valid_AphiaID (online search)
  posNA <- which(is.na(zz$valid_AphiaID) & !zz$status %in% taxStatus("invNames"))
  if(length(posNA)>0)
  {
    chunkSize <- 50
    while(length(posNA)>0)
    {
      if(worms_connection())
      {
        # Search the IDs
        tmp <- try(worrms::wm_record(id = zz$AphiaID[na.omit(posNA[1:chunkSize])]), silent = T)
        
        # Update the WoRMS table
        pos <- match(zz$AphiaID, tmp$AphiaID)
        zz[!is.na(pos), 1:ncol(tmp)] <- tmp[na.omit(pos),]
        
        posNA <- posNA[-c(1:chunkSize)]
      }
    }
  }
  
  # Return output
  return(zz)
}


## Modified from worrms (including extant_only)
wm_records_names2 <- function(name, marine_only = TRUE, extant_only = TRUE, ...)
{
  assert(name, "character")
  assert(marine_only, "logical")
  
  args <- cc(list(marine_only = as_log(marine_only), extant_only = as_log(extant_only)))
  args <- c(args,
            stats::setNames(as.list(name),
                            rep('scientificnames[]',
                                length(name))))
  result <- wm_GET(file.path(wm_base(), "AphiaRecordsByNames"),
                   query = args, ...)
  if (identical(result, tibble::tibble()))
  {
    rep(list(tibble::data_frame()), length(name))
  }
  else
  {
    result
  }
}
environment(wm_records_names2) <- asNamespace("worrms")


wm_records_taxamatch2 <- function(name, marine_only = TRUE, extant_only = TRUE, ...)
{
  assert(name, "character")
  assert(marine_only, "logical")
  
  args <- cc(list(marine_only = as_log(marine_only), extant_only = as_log(extant_only)))
  args <- c(args,
            stats::setNames(as.list(name),
                            rep('scientificnames[]',
                                length(name))))
  result <- wm_GET(file.path(wm_base(), "AphiaRecordsByMatchNames"),
                   query = args, ...)
  if (identical(result, tibble::tibble()))
  {
    rep(list(tibble::data_frame()), length(name))
  }
  else
  {
    result
  }
}
environment(wm_records_taxamatch2) <- asNamespace("worrms")


#### END ####