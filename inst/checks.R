source("R/linear.functions.R")
library(lmerTest)
dat<-read.csv("../gamlj_mixed/data/dat3x2x2_mixed.csv")
dat$x<-as.numeric(scale(dat$x))
dat$wfac<-factor(dat$wfac)
dat$wfac3<-factor(dat$wfac3)
dat$bfac<-factor(dat$bfac)
dat$cluster<-factor(dat$cluster)
contrasts(dat$wfac3)<-.createContrasts(levels(dat$wfac3),"deviation")

#contrasts(dat$wfac3)<-contr.sum(3)


model<-lm(y~wfac3*wfac*x,data=dat)
summ<-summary(model)
rows<-rownames(summ$coefficients)
rows<-rows[rows!="(Intercept)"]

contrasts<-data.frame(1:2)
contrasts$var<-c("wfac","wfac3")
contrasts$type<-c("simple","deviation")
data<-dat
seq_len(3)


.getContrastsList<-function(model) {

labels<-sapply(contrasts$var, function(varname) {
  contrast<-contrasts[contrasts$var==varname,]
  var<-data[,contrast$var]
  levels <- base::levels(var)
  .contrastLabels(levels, contrast$type)
})

layout<-attr(terms(model),"term.labels")
vars<-as.character(attr(terms(model),"variables"))
vars<-vars[c(-1,-2)]
laylist<-sapply(layout,function(a) strsplit(a,":",fixed = T))

for (i in seq_along(vars)) {
  if (vars[i] %in% names(labels)) {
   laylist<-sapply(laylist, function(a) {
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
c(" ",unlist(final))
}

q<-.getContrastsList(model)
q<-unlist(q)

