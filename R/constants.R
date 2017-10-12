WARNS<-list()

WARNS["ph.covariates"]<-"Post-hocs means are estimated keeping constant other independent variable(s) in the model"
WARNS["ph.interactions"]<-"Post-hocs means are estimated averaging across interacting factors (if any) and setting interacting covariates to zero (if any)"

WARNS["se.interactions"]<-"Simple effects are estimated setting higher order moderator (if any) in covariates to zero and averaging across moderating factors levels (if any)"
WARNS["se.covariates"]<-"Simple effects are estimated keeping constant other independent variable(s) in the model"
