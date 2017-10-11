.nicifyTerms<-function(term) {
  term<-gsub("`","",term,fixed = T)
  gsub(":"," ✻ ",term)
}

.getDummiesNames<-function(varname,data) {
  if (!(varname %in% names(data)))
    return(varname)
  if (is.factor(data[,varname]))
  {
    paste(varname,1:length(levels(data[,varname])[-1]),sep="")
  } else 
    varname
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
      contrast <- matrix(0, nrow=nLevels, ncol=nLevels-1)
      for (i in seq_len(nLevels-1)) {
        contrast[i+1, i] <- 1
        contrast[1, i] <- -1
      }
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
  } else {
    all <- paste(levels, collapse=', ')
    for (i in seq_len(nLevels-1))
      labels[[i]] <- paste(levels[i+1], '-', all)
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

# This is a crazy piece of software to extract labels of contrasts and
# merge them to the names of the terms of a summary.lm(). Before thinking
# it is too elaborated, please consider that this would attach the right
# labels even in cases two different model terms have the same name (it may happen) 
# or one dummy gets the same name of a conitnuous variable (this may happens to).


.getModelContrastsLabels<-function(contrasts,model) {
  contrasts<-as.data.frame(do.call(rbind,contrasts))
  print(contrasts)
  data=local.getModelData(model)
  facts<-local.getModelFactors(model)  
  labels<-sapply(facts, function(varname) {
    if (varname %in% contrasts$var) {
      contrasttype<-contrasts[contrasts$var==varname,"type"]
      print(paste(contrasts$var,contrasttype))
    } else {
      contrasttype<-"deviation"
    }
    
    var<-data[,varname]
    levels <- base::levels(var)
    .contrastLabels(levels, contrasttype)
  },simplify = F)
  layout<-attr(terms(model),"term.labels")
  vars<-as.character(attr(terms(model),"variables"))
  vars<-vars[c(-1,-2)]
  laylist<-sapply(layout,function(a) strsplit(a,":",fixed = T))
  for (i in seq_along(vars)) {
    if (vars[i] %in% names(labels)) {
      laylist<-sapply(laylist, function(a) {
        print(paste(vars[i],"changed"))
        a[which(a == vars[i])]<-paste(unlist(labels[vars[i]]),collapse = "#")
        a
      })
    }
  }
  laylist<-sapply(laylist,function(a) strsplit(a,"#",fixed = T))
  final<-list()
  nr<-1
  for (j in seq_along(laylist)) {
    lay<-laylist[j]
    records<-expand.grid(lay[[1]],stringsAsFactors = F)
    ready<-apply(records,1,paste,collapse=":")
    for (i in seq_along(ready)) {
      final[nr]<-ready[i]
      nr<-nr+1
    }
  } 
  ### add an empty contrast for the intercept ####
  c("Intercept",unlist(final))
}
