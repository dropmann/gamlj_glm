dat<-read.csv("data/data3by2.csv")
dat$threegroups<-factor(dat$threegroups)
dat$twogroups<-factor(dat$twogroups)
contrasts(dat$threegroups)<-contr.sum(3)
contrasts(dat$twogroups)<-contr.sum(2)

#contrasts(dat$threegroups)<-contr.treatment(3)
contrasts(dat$threegroups)

nn<-names(mod$coefficients)
form<-y~threegroups
mod<-lm(form,data=dat)

car::Anova(mod,type="III")
(ss<-car::Anova(mod,type="II"))
dim(ss)
ss
names(ss)
ss
anova(mod)
ss[dim(ss)[1],]
(ss<-drop1(mod,.~.,test="F"))
ss[2:dim(ss)[1],c(2,1,5,6)]
as.list(ss[2,])
anova(mod)
modelTerms<-list("a","b","c",c("a","b"))
dep<-"y"
highest<-max(sapply(modelTerms,length))

modelt<-terms(mod)
or<-attr(modelt,"order")
tl<-attr(modelt,"term.labels")
i=1
for( i in 1:highest) {
f<-tl[or<(i+1)]
print(f)
f<-jmvcore::constructFormula(dep, f)
print(f)
print(drop1(mod,scope = .~threegroups+twogroups+x,test="F"))
}
form<-y~x+threegroups+twogroups
car::Anova(lm(form,data=dat),type=2)
drop1(lm(form,data=dat),.~.,test="F")

drop1(mod,c("twogroups","threegroups","x"),test="F")
drop1(mod,form,test="F")
car::Anova(mod,type=2)
anova(mod)
mod$df.residual
summary(mod)
k<-length(coef(mod))
confint(mod,1:k)
ss<-summary(mod)
eresults<-ss$coefficients
cc<-contrasts(dat$threegroups)
variables<-c("threegroups","twogroups")
columns<-rownames(eresults)
for (variable in variables) {
dummies<-paste(variable,colnames(contrasts(dat[,variable])),sep="")
labels<-paste(dummies,"()()",sep="")
for(i in seq_along(dummies))
  columns<-gsub(dummies[i],labels[i],columns,fixed=T)
}
print(columns)
a<-c("a","b","c","d")
b<-c("a","b")
c<-c("t","f")
gsub(b,c,a)
dummies
eresults
termsClasses<-attr(mod$terms,"dataClasses")
termsClasses<-termsClasses[termsClasses=="factor"]

variable<-names(mod$contrasts[1])
.getDummies(variable,mod)
