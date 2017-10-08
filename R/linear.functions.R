.nicifyTerms<-function(term) {
  gsub(":"," ✻ ",term)
}

.getDummiesNames<-function(varname,data) {
  if (is.factor(data[,varname]))
  {
    paste(varname,1:length(levels(data[,varname])[-1]),sep="")
  } else 
    varname
}

.removeFromFormula<-function(aFormula,variables) {
  dep<-all.vars(aFormula)[[1]]
  terms<-attr(terms(aFormula),"factors")
  terms<-terms[variables == rownames(terms),]
  terms<-apply(terms,2,sum)
  list(dep,names(which(terms==0)))
}


.scaleVariables=function(factors,covariates,data) {
  for(factor in factors)
     data[[factor]]<-factor(data[[factor]])
  for(covariate in covariates)
     data[[covariate]] <- jmvcore::toNumeric(data[[covariate]])
  data
}

.createContrasts=function(levels, type) {
  
  nLevels <- length(levels)
  
  if (type == 'simple') {
    
    dummy <- contr.treatment(levels)
    dimnames(dummy) <- NULL
    coding <- matrix(rep(1/nLevels, prod(dim(dummy))), ncol=nLevels-1)
    contrast <- (dummy - coding)
    
  } else if (type == 'deviation') {
    contrast <- matrix(0, nrow=nLevels, ncol=nLevels-1)
    for (i in seq_len(nLevels-1)) {
      contrast[i+1, i] <- 1
      contrast[1, i] <- -1
    }
    
  } else if (type == 'difference') {
    
    contrast <- stats::contr.helmert(levels)
    for (i in 1:ncol(contrast))
      contrast[,i] <- contrast[,i] / (i + 1)
    
    dimnames(contrast) <- NULL
    
  } else if (type == 'helmert') {
    
    contrast <- matrix(0, nrow=nLevels, ncol=nLevels-1)
    
    for (i in seq_len(nLevels-1)) {
      p <- (1 / (nLevels - i + 1))
      contrast[i,i] <- p * (nLevels - i)
      contrast[(i+1):nLevels,i] <- -p
    }
    
  } else if (type == 'polynomial') {
    
    contrast <- stats::contr.poly(levels)
    dimnames(contrast) <- NULL
    
  } else if (type == 'repeated') {
    
    contrast <- matrix(0, nrow=nLevels, ncol=nLevels-1)
    for (i in seq_len(nLevels-1)) {
      contrast[i,  i] <- 1
      contrast[i+1,i] <- -1
    }
    
  } else {
    
    contrast <- NULL
  }
  
  contrast
}

.contrastLabels=function(levels, type) {
  
  nLevels <- length(levels)
  labels <- list()
  
  if (length(levels) <= 1) {
    
    # do nothing
    
  } else if (type == 'simple') {
    
    for (i in seq_len(nLevels-1))
      labels[[i]] <- paste(levels[i+1], '-', levels[1])
    
  } else if (type == 'deviation') {
    
    all <- paste(levels, collapse=', ')
    for (i in seq_len(nLevels-1))
      labels[[i]] <- paste(levels[i+1], '-', all)
    
  } else if (type == 'difference') {
    
    for (i in seq_len(nLevels-1)) {
      rhs <- paste0(levels[1:i], collapse=', ')
      labels[[i]] <- paste(levels[i + 1], '-', rhs)
    }
    
  } else if (type == 'helmert') {
    
    for (i in seq_len(nLevels-1)) {
      rhs <- paste(levels[(i+1):nLevels], collapse=', ')
      labels[[i]] <- paste(levels[i], '-', rhs)
    }
    
  } else if (type == 'repeated') {
    
    for (i in seq_len(nLevels-1))
      labels[[i]] <- paste(levels[i], '-', levels[i+1])
    
  } else if (type == 'polynomial') {
    
    names <- c('linear', 'quadratic', 'cubic', 'quartic', 'quintic', 'sextic', 'septic', 'octic')
    
    for (i in seq_len(nLevels-1)) {
      if (i <= length(names)) {
        labels[[i]] <- names[i]
      } else {
        labels[[i]] <- paste('degree', i, 'polynomial')
      }
    }
  }
  
  labels
}

.scaleContinuous<-function(var,method,by=NULL) {
  
  if (method=="centered") 
          var<-scale(var,scale = F)  
  if (method=="cluster-based centered")    
          var<-unlist(tapply(var,by,scale,scale=F))
  if (method=="standardized") 
          var<-scale(var,scale = T)  
  if (method=="cluster-based standardized")     
          var<-unlist(tapply(var,by,scale,scale=T))
  as.numeric(var)
}


