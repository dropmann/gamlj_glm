
local.getModelTerms<-function(model) names(coefficients(model))

local.getModelData<-function(model) model$model

local.getModelDummies<-function(variable,model) {
  nn<-local.getModelTerms(model)
  nn<-nn[grep(variable,nn)]
  nn[grep(":",nn,invert = T,fixed = T)]     
}

local.getModelFactors<-function(model) names(model$contrasts)



