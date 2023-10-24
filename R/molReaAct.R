molReaAct <- function(data = data,mol = mol,formula_column,error.term = 0.000010, type = "Dataset",HtoC_ratio_threshold = 1.5, Trans_threshold = c(1,10)) {
  mol$C <- sapply(mol[[formula_column]], function(formula) {
    matches <- gregexpr("C\\d*", formula)
    sums <- sapply(regmatches(formula, matches), function(match) {
      if (length(match) > 0) {
        numeric_part <- sub("C", "", match)
        if (grepl("^\\d+$", numeric_part)) {
          return(as.numeric(numeric_part))
        } else {
          return(1)
        }
      } else {
        return(0)
      }
    })
    sum(sums)
  })
  mol$H <- sapply(mol[[formula_column]], function(formula) {
    matches <- gregexpr("H\\d*", formula)
    sums <- sapply(regmatches(formula, matches), function(match) {
      if (length(match) > 0) {
        numeric_part <- sub("H", "", match)
        if (grepl("^\\d+$", numeric_part)) {
          return(as.numeric(numeric_part))
        } else {
          return(1)
        }
      } else {
        return(0)
      }
    })
    sum(sums)
  })
  mol$O <- sapply(mol[[formula_column]], function(formula) {
    matches <- gregexpr("O\\d*", formula)
    sums <- sapply(regmatches(formula, matches), function(match) {
      if (length(match) > 0) {
        numeric_part <- sub("O", "", match)
        if (grepl("^\\d+$", numeric_part)) {
          return(as.numeric(numeric_part))
        } else {
          return(1)
        }
      } else {
        return(0)
      }
    })
    sum(sums)
  })
  mol$N <- sapply(mol[[formula_column]], function(formula) {
    matches <- gregexpr("N\\d*", formula)
    sums <- sapply(regmatches(formula, matches), function(match) {
      if (length(match) > 0) {
        numeric_part <- sub("N", "", match)
        if (grepl("^\\d+$", numeric_part)) {
          return(as.numeric(numeric_part))
        } else {
          return(1)
        }
      } else {
        return(0)
      }
    })
    sum(sums)
  })
  mol$HtoC_ratio <- mol$H/mol$C
  if(identical(x = colnames(data), y = rownames(mol)) == FALSE){
    stop("Something is incorrect in your row names")
  }
  trans.full = Transformation_Database
  if (type == "Dataset") {
    # Creating a presence matrix for all peaks observed within the dataset
    bulk.peaks = as.data.frame(cbind("Sample_All",row.names(mol)))
    bulk.peaks[,2] = as.numeric(as.character(bulk.peaks[,2]))
    colnames(bulk.peaks) = c("sample", "peak.x")
    
    # Sort molecules from smallest to largest by mass
    bulk.peaks <- bulk.peaks[order(bulk.peaks$peak.x),]
    
    peak.2.peak = NULL
    start_peak = 1
    
    # i = 1
    # Running a loop to compare each peak to each other peak
    for(i in start_peak:(ncol(data)-1)){ # I cannot stress the importance of the "-1" here...
      
      # Creating a data matrix to ensur no repeat or negative differences
      Distance_Results = bulk.peaks[-1:-i,] # Removing all peaks up to, and including the current peak
      Distance_Results$peak.y = bulk.peaks$peak.x[i] # Setting the peak of interest
      Distance_Results$Dist = Distance_Results$peak.x - Distance_Results$peak.y # Finding the difference between all peaks and the peak of interest
      
      # Adding in error terms to the matrix
      Distance_Results$Dist.plus = Distance_Results$Dist + error.term
      Distance_Results$Dist.minus = Distance_Results$Dist - error.term
      Distance_Results$Trans.name = -999
      
      # Reorganizing the data to make it applicable with other scripts
      Distance_Results = Distance_Results[,c("sample", "Dist", "peak.x", "peak.y", "Dist.plus", "Dist.minus", "Trans.name")]
      
      for (current.trans in unique(trans.full$Name)) { # note that for masses with multiple names, only the last name is going to be recorded
        
        mass.diff = trans.full$Mass[which(trans.full$Name == current.trans)]
        if (length(mass.diff) > 1) { break() }
        Distance_Results$Trans.name[ which(Distance_Results$Dist.plus >= mass.diff & Distance_Results$Dist.minus <= mass.diff)] = current.trans
      }
      
      # Removing differences that didn't match any transformation
      Distance_Results = Distance_Results[-which(Distance_Results$Trans.name == -999),]
      
      # Building a larger peak.2.peak file
      peak.2.peak = rbind(peak.2.peak, Distance_Results)
      
      print(paste("Finished running through peak #", i, " on ", date(), sep = ""))
      
    }
    
    # Creating a num.trans file for network generation
    peak.stack = as.data.frame(c(peak.2.peak$peak.x, peak.2.peak$peak.y)); head(peak.stack)
    peak.profile = as.data.frame(tapply(X = peak.stack[,1], INDEX = peak.stack[,1], FUN = 'length' )); dim(peak.profile)
    colnames(peak.profile) = 'num.trans.involved.in'
    peak.profile$peak = row.names(peak.profile)
    peak.profile$type = type
    
    # "num.trans.involved.in" is zero 
    mol.trans0 = mol[!(rownames(mol) %in% rownames(peak.profile)),]
    peak.profile.trans0 = data.frame(num.trans.involved.in = 0, peak = rownames(mol.trans0),type = type)
    rownames(peak.profile.trans0) = peak.profile.trans0$peak
    
    peak.profile.dataset = rbind(peak.profile,peak.profile.trans0)
    peak.profile.dataset = peak.profile.dataset[rownames(mol),]
    mol$peak = rownames(mol)
    peak.profile.dataset.HC = peak.profile.dataset %>% 
      left_join(mol,by = "peak") %>% 
      mutate(Class = case_when(HtoC_ratio >= HtoC_ratio_threshold & num.trans.involved.in > Trans_threshold[2] ~"Labile_Active",
                               HtoC_ratio >= HtoC_ratio_threshold & num.trans.involved.in <= Trans_threshold[1] ~"Labile_Inactive",
                               HtoC_ratio < HtoC_ratio_threshold & num.trans.involved.in > Trans_threshold[2] ~"Recalcitrant_Active",
                               HtoC_ratio < HtoC_ratio_threshold & num.trans.involved.in <= Trans_threshold[1] ~"Recalcitrant_Inactive"))
    
    return(peak.profile.dataset.HC)
    
  }
  else if (type == "Sample"){
    # pull out just the sample names
    samples.to.process = rownames(data)
    
    # current.sample = samples.to.process[1]
    
    for (current.sample in samples.to.process) {
      
      one.sample.matrix = t(data[which(rownames(data) == current.sample),,drop = F])
      one.sample.matrix = data.frame(peak = rownames(one.sample.matrix),one.sample.matrix)
      
      Sample_Peak_Mat <- one.sample.matrix %>% gather("sample", "value", -1) %>% filter(value > 0) %>% select(sample, peak)
      
      Distance_Results <- Sample_Peak_Mat %>% 
        left_join(Sample_Peak_Mat, by = "sample",relationship = "many-to-many") %>% 
        mutate(peak.x = as.numeric(peak.x),peak.y = as.numeric(peak.y)) %>% 
        filter(peak.x > peak.y) %>% mutate(Dist = peak.x - peak.y) %>% 
        select(sample,peak.x,peak.y,Dist)
      
      Distance_Results$Dist.plus = Distance_Results$Dist + error.term
      Distance_Results$Dist.minus = Distance_Results$Dist - error.term
      Distance_Results$Trans.name = -999
      head(Distance_Results)
      
      dist.unique = unique(Distance_Results[,'sample']) #unique samples
      
      for (current.trans in unique(trans.full$Name)) { # note that for masses with multiple names, only the last name is going to be recorded
        
        mass.diff = trans.full$Mass[which(trans.full$Name == current.trans)]
        if (length(mass.diff) > 1) { break() }
        Distance_Results$Trans.name[which(Distance_Results$Dist.plus >= mass.diff & Distance_Results$Dist.minus <= mass.diff)] = current.trans
        
      }
      
      Distance_Results = Distance_Results[-which(Distance_Results$Trans.name == -999),]
      head(Distance_Results)
      
      # find the number of transformations each peak was associated with
      peak.stack = as.data.frame(c(Distance_Results$peak.x,Distance_Results$peak.y)); head(peak.stack)
      peak.profile = as.data.frame(tapply(X = peak.stack[,1],INDEX = peak.stack[,1],FUN = 'length' )); dim(peak.profile)
      colnames(peak.profile) = 'num.trans.involved.in'
      peak.profile$sample = dist.unique
      peak.profile$peak = row.names(peak.profile)
      peak.profile$type = type
      head(peak.profile)
      
      # "num.trans.involved.in" is zero 
      mol.trans0 = Sample_Peak_Mat[!(Sample_Peak_Mat$peak %in% rownames(peak.profile)),]
      peak.profile.trans0 = data.frame(num.trans.involved.in = 0,sample = current.sample, peak = mol.trans0$peak,type = type)
      rownames(peak.profile.trans0) = peak.profile.trans0$peak
      
      peak.profile.all = rbind(peak.profile,peak.profile.trans0)
      
      if (current.sample == samples.to.process[1]) {
        peak.profile.sample = peak.profile.all
      }else{
        peak.profile.sample = rbind(peak.profile.sample,peak.profile.all)
      }
      
      print(paste("Finished running through sample #", current.sample, " on ", date(), sep = ""))
    }
    
    peak.profile.sample.HC = peak.profile.sample %>% 
      left_join(mol.HC,by = "peak") %>% 
      rename(HtoC_ratio = colnames(.)[5]) %>% 
      mutate(Class = case_when(HtoC_ratio >= HtoC_ratio_threshold & num.trans.involved.in > Trans_threshold[2] ~"Labile_Active",
                               HtoC_ratio >= HtoC_ratio_threshold & num.trans.involved.in <= Trans_threshold[1] ~"Labile_Inactive",
                               HtoC_ratio < HtoC_ratio_threshold & num.trans.involved.in > Trans_threshold[2] ~"Recalcitrant_Active",
                               HtoC_ratio < HtoC_ratio_threshold & num.trans.involved.in <= Trans_threshold[1] ~"Recalcitrant_Inactive"))
    
    return(peak.profile.sample.HC)
  }

}
