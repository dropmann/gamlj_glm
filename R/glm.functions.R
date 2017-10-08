.getModelDummies<-function(variable,model) {
  nn<-names(model$coefficients)
  nn<-nn[grep(variable,nn)]
  nn[grep(":",nn,invert = T,fixed = T)]     
}

