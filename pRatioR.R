(WD <- getwd())
if (!is.null(WD)) setwd(WD)

################################################################################################################

# Install packages if they don't exit

################################################################################################################

list.of.packages <- c("RSQLite", "readr", "stringi", "plyr", "Peptides", "XML")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos="http://cran.us.r-project.org")

# Libraries

library("RSQLite")
library("readr")
library("stringi")
library("plyr")
library("Peptides")
library("XML")

################################################################################################################

# Get Params

################################################################################################################

source(paste0(WD,"/Config_pRatioPD21_12.txt"))

################################################################################################################

# Or put it manually

################################################################################################################

# deltaMassThreshold = 15 # in ppm

# deltaMassAreas = 5 # number of jumps: 1,3 or 5

# tagMass = 229.162932

# tagDecoy = "DECOY_"

################################################################################################################

# Input and output manual

################################################################################################################

# input <- "D:/projects/pRatio/pRatio_R/TMT1/MSF/FR_1/TMT1_PESA8_Fr1.msf"

# output <- input
# 
# output <- substr(output, 1, nchar(output) - 4)
# 
# output <- paste(output,"_results_FDR.txt",sep="")

# output <- "D:/projects/pRatio/pRatio_R/TMT1/MSF/FR_1/TMT1_PESA8_Fr1_results.txt"

################################################################################################################

# Local functions

################################################################################################################

list.dirs <- function(path=".", pattern=NULL, all.dirs=FALSE,
                      full.names=FALSE, ignore.case=FALSE) {
  # use full.names=TRUE to pass to file.info
  all <- list.files(path, pattern, all.dirs,
                    full.names=TRUE, recursive=FALSE, ignore.case)
  dirs <- all[file.info(all)$isdir]
  # determine whether to return full names or just dir names
  if(isTRUE(full.names))
    return(dirs)
  else
    return(basename(dirs))
}

# prepare DATA: mods & prots
ptmAnnotation <- function(x)
{
  pep<-x[1,]$Sequence
  b<-0
  p<-""
  modNum<-1
  modMass<-0
  for (i in x$Position)
  {
    p <- stri_flatten(c(p,substr(pep,b,i+1),"[",x[modNum,]$DeltaMass,"]"),collapse="")
    modMass<-modMass + x[modNum,]$DeltaMass
    #p <- paste(p,substr(pep,b,i+1),"[",x[modNum,]$DeltaMass,"]",sep="")
    b <- i+2
    modNum <- modNum+1
  }
  #p <- paste(p,substr(pep,b,nchar(pep)),sep="")
  p <- stri_flatten(c(p,substr(pep,b,nchar(pep))),collapse="")
  return(c(p,modMass))
}

# filter by deltaMass
filterDeltaMass <- function(x, deltaMassThreshold, deltaMassAreas)
{
  TheoreticalModTag=x[1]
  Mass=x[2]
  ScoreValue=x[3]
  jump1_ppm = abs(TheoreticalModTag - Mass) / TheoreticalModTag * 1e6
  if (jump1_ppm >= deltaMassThreshold)
  {
    if (deltaMassAreas <= 1) { return(0.01) } # jump 1 >= threshold
    else 
    {
      MassCorr <- Mass - 1.0033
      jump23_ppm = abs(TheoreticalModTag - MassCorr) / TheoreticalModTag * 1e6
      if (jump23_ppm >= deltaMassThreshold)
      {  
        if (deltaMassAreas <= 3) { return(0.01) } # jump 23 >= threshold
        else
        {
          MassCorr2 <- Mass - 1.0033
          jump45_ppm = abs(TheoreticalModTag - MassCorr2) / TheoreticalModTag * 1e6
          if (jump45_ppm >= deltaMassThreshold) {return (0.01)} # jump 45 >= threshold
          else {return (ScoreValue)} # jump 45 < threshold
        }
      }
      else
      {
        return (ScoreValue) # jump 23 < threshold
      }
    }
  }
  else
  {
    return(ScoreValue) # jump 1 < threshold
  }
}

# replace modifications to symbols (output to tmp)
parsingMod <- function(weight, symbol) {
  resPratio$Sequence <<- gsub( paste0('\\[',weight,'\\]'), symbol, resPratio$Sequence) # global access
}


################################################################################################################

# Main function

################################################################################################################

# parse the modifications xml file
mod_xml <- xmlParse("modifications.xml")
modifications <- xmlToDataFrame(nodes = xmlChildren(xmlRoot(mod_xml)[["modifSet"]]) )

# list of MSF folders
MSFfolders <- list.dirs(path = paste0(WD,"/",Expto,"/MSF"), pattern=Patern)

for (j in Expto){

  for (k in MSFfolders){

    files <- list.files(path = paste(WD,"/",j,"/MSF/",k,sep=""),pattern="*.msf$", full.names = TRUE)

    for (i in files) {
      
      input <- i
      
      output <- input
      
      output <- substr(output, 1, nchar(output) - 4)
      
      output <- paste(output,"_results_FDR.txt",sep="")
      
      db=dbConnect(SQLite(), dbname=input)
      
      queryMain = "select 
      p.peptideid, 
      fi.filename, 
      sh.firstscan, 
      sh.lastscan, 
      sh.charge, 
      p.sequence,
      sh.mass,  
      ps.scorevalue,  
      sh.retentiontime,  
      p.searchenginerank, 
      p.deltascore  
      from peptides p, 
      peptideScores ps, 
      spectrumHeaders sh, 
      massPeaks mp, 
      workFlowInputFiles fi, 
      processingNodeScores scoreNames 
      where p.peptideid = ps.peptideid 
      and sh.spectrumid = p.spectrumid 
      and (fi.fileid = mp.fileid or mp.fileid = -1) 
      and mp.masspeakid = sh.masspeakid 
      and scoreNames.scoreid = ps.scoreid 
      and scoreNames.ScoreName = 'Xcorr'  
      and p.searchenginerank = 1   
      and ps.scorevalue > 1.5
      order by 
      fi.filename desc,  
      sh.firstscan asc, 
      sh.lastscan asc,  
      sh.charge asc,  
      ps.scorevalue desc" 
      
      data=dbGetQuery(conn = db, queryMain)
      
      queryModifications = "select 
      p.peptideid, 
      paam.aminoacidmodificationid, 
      paam.position,  
      p.sequence,  
      aam.modificationname,  
      aam.deltamass   
      from peptides p, 
      peptideScores ps, 
      spectrumHeaders sh, 
      peptidesaminoacidmodifications paam, 
      aminoacidmodifications aam 
      where p.peptideid = paam.peptideid 
      and sh.spectrumid = p.spectrumid 
      and p.peptideid = ps.peptideid 
      and aam.aminoacidmodificationid = paam.aminoacidmodificationid 
      and p.searchenginerank = 1 
      and ps.scorevalue > 1.5
      order by p.peptideid ASC, paam.position ASC"
      
      dataMod=dbGetQuery(conn = db, queryModifications)
      
      queryProteinInfo = "select 
      pq.peptideid, 
      p.sequence, 
      pq.proteinid,  
      q.description  
      from peptidesProteins pq, 
      spectrumHeaders sh, 
      peptides p, 
      peptideScores ps, 
      proteinAnnotations q 
      where pq.peptideid = p.peptideid 
      and p.peptideid = ps.peptideid 
      and pq.proteinid = q.proteinid 
      and sh.spectrumid = p.spectrumid 
      and p.searchenginerank = 1 
      and ps.scorevalue > 1.5
      order by pq.peptideid asc"
      
      dataProt=dbGetQuery(conn = db, queryProteinInfo)      
      
      
      ## prepare DATA: mods & prots
      dataModAnnotation <- ddply(dataMod,.(PeptideID),ptmAnnotation)
      
      colnames(dataModAnnotation) <- c("PeptideID","Sequence","modMass")
      
      dataModAnnotation$modMass <- as.numeric(dataModAnnotation$modMass)
      
      # mix unmodified and modified
      dataModTmp <- merge(unique(data[,c("PeptideID","Sequence")]),dataModAnnotation,by = "PeptideID",all.x=TRUE)
      dataModTmp[is.na(dataModTmp["Sequence.y"]),"Sequence.y"] <- dataModTmp[is.na(dataModTmp["Sequence.y"]),"Sequence.x"]
      dataModAll <- dataModTmp[,c("PeptideID","Sequence.y","modMass")]
      colnames(dataModAll) <- c("PeptideID","SequenceMod","modMass")
      
      redundances <- aggregate(Description ~ PeptideID, data=dataProt, paste, collapse = " -- ")
      colnames(redundances) <- c("PeptideID","Redundances")
      dataProt.u <- dataProt[!duplicated(dataProt["PeptideID"]),]
      peptideProt <- merge(dataProt.u[,c("PeptideID","Description")], redundances, by="PeptideID", all.x=TRUE)
      
      #*****
      dataAll <- cbind(data, dataModAll$"SequenceMod",  dataModAll$"modMass", peptideProt[ , -which(names(peptideProt) %in% c("PeptideID"))])
      colnames(dataAll) <- c("PeptideID","FileName","FirstScan","LastScan","Charge","Sequence","Mass","ScoreValue","RetentionTime","SearchEngineRank","DeltaScore","SequenceMod","modMass","Description","Redundances")
      
      ## Calculate theoretical mass
      dataAll[is.na(dataAll[,"modMass"]),]$modMass <- 0
      # Be careful!!! The following line of code could print a Warning messages:Sequence 1 has unrecognized amino acid types. Output value might be wrong calculated 
      dataAll <- cbind(dataAll,as.data.frame(unlist(lapply(dataAll[,c("Sequence")], mw, monoisotopic=TRUE))))
      names(dataAll)[length(names(dataAll))]<-"Theoretical" 
      dataAll$Theoretical <- dataAll$Theoretical + 1.00727647
      dataAll <- cbind(dataAll, dataAll$Theoretical + dataAll$modMass + tagMass)
      names(dataAll)[length(names(dataAll))]<-"TheoreticalModTag" 
      dataAll <- cbind(dataAll, abs(dataAll$Mass - dataAll$Theoretical - dataAll$modMass - tagMass) / dataAll$Mass * 1e6)
      names(dataAll)[length(names(dataAll))]<-"deltaMassTargetppm" 
      
      ## Decoy tagging
      isDecoy <- rep(0, dim(dataAll)[1])
      isTarget <- rep(1, dim(dataAll)[1])
      protein <- dataAll[,'Description']
      index <- grep(tagDecoy,protein,fixed=TRUE)
      isDecoy[index] <- 1
      isTarget[index] <- 0
      dataAll <- cbind(dataAll,isDecoy,isTarget)
      
      ## filter by deltaMass
      jump1ScoreValue <-as.data.frame(unlist(apply(dataAll[,c("TheoreticalModTag","Mass","ScoreValue")], 1, filterDeltaMass, deltaMassThreshold=deltaMassThreshold, deltaMassAreas=deltaMassAreas)))
      colnames(jump1ScoreValue) <- "ScoreValueAfterJUMP"
      dataAll$ScoreValue<-jump1ScoreValue$ScoreValueAfterJUMP #Assign the calculated scored after being modified and
                                                  #assign to the column ScoreValue
      
      ## Add xcorr_c
      n = dim(dataAll)[1]
      
      xcorr_c <- function(x) {
        r=1
        if(as.numeric(x[1])>2) {r=1.22}
        xcorr_c = log((as.numeric(x[2]))/r)/log(2*nchar(as.character(x[3])))
        return (xcorr_c)
      }
      
      dataAll <- cbind(dataAll,apply(dataAll[,c("Charge","ScoreValue","Sequence")], 1, xcorr_c))
      
      colnames(dataAll)[ncol(dataAll)] <- "xcorr_c"
      
      # sort by xcorr_c
      #dataAll <- dataAll[order(decreasing = TRUE,dataAll$xcorr_c),]
      ##dataAll <- dataAll[order(decreasing = TRUE,dataAll$ScoreValue),]
      #tmp <- cbind(dataAll[, "xcorr_c"], dataAll[, "isDecoy"])
      ##tmp <- cbind(dataAll[, "ScoreValue"], dataAll[, "isDecoy"])
      #FP <- cumsum(tmp[, 2])
      #tmp <- cbind(tmp, FP)
      #xcorr_cP <- unlist(lapply(1:n, function(x) (tmp[x, 'FP'])/n))
      #dataAll <- cbind(dataAll, xcorr_cP)
      
      ### FDR ScoreValue
      
      dataAll <- dataAll[order(decreasing = TRUE,dataAll$ScoreValue),]
      tmp <- cbind(dataAll[, "ScoreValue"], dataAll[, "isDecoy"], dataAll[, "isTarget"])
      FP <- cumsum(tmp[, 2])
      TP <- cumsum(tmp[, 3])
      tmp <- cbind(tmp, FP, TP)
      xcorr_FDR <- unlist(lapply(1:dim(dataAll)[1], function(x) (tmp[x, 'FP'])/(tmp[x, 'TP'])))
      dataAll <- cbind(dataAll, tmp, xcorr_FDR)
      xcorr_FDRa <- unlist(lapply(1:dim(dataAll)[1], function(x) max(dataAll[1:x,"xcorr_FDR"])))
      dataAll <- cbind(dataAll, xcorr_FDRa)
      
      ### FDR CALC
      dataAll <- dataAll[order(decreasing = TRUE,dataAll$xcorr_c),]
      tmp <- cbind(dataAll[, "xcorr_c"], dataAll[, "isDecoy"], dataAll[, "isTarget"])
      FP <- cumsum(tmp[, 2])
      TP <- cumsum(tmp[, 3])
      tmp <- cbind(tmp, FP, TP)
      xcorr_c_FDR <- unlist(lapply(1:dim(dataAll)[1], function(x) (tmp[x, 'FP'])/(tmp[x, 'TP'])))
      dataAll <- cbind(dataAll, tmp, xcorr_c_FDR)
      xcorr_c_FDRa <- unlist(lapply(1:dim(dataAll)[1], function(x) max(dataAll[1:x,"xcorr_c_FDR"])))
      dataAll <- cbind(dataAll, xcorr_c_FDRa)
      
      res <- dataAll[dataAll$xcorr_c_FDR < 0.01 & dataAll$isTarget == 1,]  
      #res <- dataAll[dataAll$xcorr_c_FDR < 0.01,]            
      
      fileName <- strsplit(data[1,"FileName"], fixed = TRUE, split = "\\")[[1]][length(strsplit(data[1,"FileName"], fixed = TRUE, split = "\\")[[1]])]
      pRatio <- "NA"; pI <- "NA"; Xcorr1Original <- "NA"; Xcorr2Search <- "NA"; Sp <- "NA"; SpRank <- "NA"; ProteinsWithPeptide <- "NA"
      
      resPratio <- cbind(fileName,fileName,res[,c("FirstScan","LastScan","Charge")],pRatio,res[,c("xcorr_c_FDR","Description","SequenceMod")],pI,res[,c("Mass","xcorr_c")],Xcorr1Original,Xcorr2Search,res[,"DeltaScore"],Sp,SpRank,ProteinsWithPeptide,res[,"Redundances"])
      
      colnames(resPratio) <- c("FileName","RAWFile","FirstScan","LastScan","Charge","pRatio","FDR","FASTAProteinDescription","Sequence","pI","PrecursorMass","Xcorr1Search","Xcorr1Original","Xcorr2Search","DeltaCn","Sp","SpRank","ProteinsWithPeptide","Redundances")
      
      #SIMPLYFIED
      resPratio <- resPratio[,c("FileName","RAWFile","FirstScan","LastScan","Charge","Sequence","FASTAProteinDescription","Xcorr1Search","FDR","Redundances")]
      
      # replace modifications to symbols (output to tmp)
      tmp <- apply(modifications[,c('weight', 'symbol')], 1, function(x) parsingMod(x[1],x[2]))

      write.table(resPratio,file = output,col.names = TRUE, row.names = FALSE,sep="\t", quote = FALSE)
      
    } # end files loop
  } # end Expto loop
} # end MSFfolders loop



