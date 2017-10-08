library(jmvtools)

instdir<-"/home/marcello/Skinner/Stat/R/Packages/gamlj/gamlj_glm/build/R"
ff<-find.package("coda",lib.loc = instdir,quiet = T)
if (length(ff)==0)
  install.packages("nloptr",lib = instdir)
ff<-find.package("pbkrtest",lib.loc = instdir,quiet = T)
if (length(ff)==0)
  install.packages("pbkrtest",lib = instdir)
#jmvtools::create('galmjglm')
  
install(debug = T)
    
