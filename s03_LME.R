#library(caret)
library(MuMIn)
setwd("/net/pepper/ABCD/CIFTI/Scripts/pfactor/")

library(lme4)
# getting pvalues from multilevel models is a 
# frequently asked question. code below uses
# lmerTest to generate coefficient estimate pvals
# more discussion here:
# https://rdrr.io/cran/lme4/man/pvalues.html
library(lmerTest)

# this function takes an unstandardized multilevel model
# and returns a dataframe with standardized estimates and 
# standard errors
# source: https://stackoverflow.com/questions/25142901/standardized-coefficients-for-lmer-model
stdCoef.merMod <- function(object) {
    sdy <- sd(getME(object, "y"))
    sdx <- apply(getME(object, "X"), 2, sd)
    sc <- fixef(object)*sdx/sdy
    se.fixef <- coef(summary(object))[,"Std. Error"]
    se <- se.fixef*sdx/sdy
    return(data.frame(stdcoef = sc, stdse = se))
}

# load data - important to specify all inputs that should 
# count as missing values in the 'na.strings' parameter
# read phenotypic data file
pheno = read.csv("Data/ABCD_rest.csv",na.strings=c("NA","NaN"," "),stringsAsFactors=T)

# read PCA expressions
expressions = read.csv("./Results/ABCD_rest_CIFTI_expressions.csv",na.strings=c("NA","NaN"," "))

# merge
dat = merge(pheno,expressions,by="subjectkey")

dat['acs_raked_propensity_score_scaled'] = dat$acs_raked_propensity_score/5000
dat['PF10_lt']=log(dat$PF10+abs(min(dat$PF10))+1) 

controls = c('Gender', 'Age', 'I(Age^2)', 'C(RaceEthnicity)','fd', 'I(fd^2)')

# and multilevel structure
nested_str = '(1| site_id_l/rel_family_id)'

dependent = 'PF10'     # variable being predicted

n_components=250
brain_formula1='A001'
for (i in 2:n_components) {
     icomp <- sprintf("%03d",i) # fix to 3 characters 
     brain_formula1 = paste(brain_formula1, ' + ', 'A', icomp, sep = "")
}

brain_formula=brain_formula1

fmla = paste(dependent, ' ~ ', brain_formula, ' + ', 
             paste(controls, collapse=' + '), ' + ', 
             nested_str, sep='')

fmla_reduced = paste(dependent, ' ~ ', 
             paste(controls, collapse=' + '), ' + ', 
             nested_str, sep='')

model = lmer(formula=fmla, data=dat, na.action='na.exclude')
model_reduced=lmer(formula=fmla_reduced, data=dat, na.action='na.exclude')

r.squaredGLMM(model)
r.squaredGLMM(model_reduced)
delta_rsquared_val=r.squaredGLMM(model)[1]-r.squaredGLMM(model_reduced)[1]
delta_rsquared_val

summary(model)
summary(model_reduced)
anova(model)
anova(model_reduced)
anova(model, model_reduced)

temp<-abs(coef(summary(model))[, 'Pr(>|t|)'])<(0.05/250)
result<-coef(summary(model))[temp, 't value']
length(result)
result

fixedformula <- as.formula(nobars(formula(model))[-2])
mm = model.matrix(fixedformula)
preds1 = mm%*%fixef(model)    
preds2 = predict(model, re.form=NA)
all(preds1 == preds2)
plot(dat$PF10, preds1, xlab='True PF10', ylab='Predicted PF10')


betas = coef(summary(model))[2:251,1]
write.table(betas,"/net/pepper/ABCD/CIFTI/Scripts/pfactor/Results/mm_betas.csv",row.names=F,col.names=F,quote=F,na="NaN")

#with weights
model = lmer(formula=fmla, data=dat, na.action='na.exclude', weights = acs_raked_propensity_score_scaled)
model_reduced=lmer(formula=fmla_reduced, data=dat, na.action='na.exclude', weights = acs_raked_propensity_score_scaled)

r.squaredGLMM(model)
r.squaredGLMM(model_reduced)
delta_rsquared_val=r.squaredGLMM(model)[1]-r.squaredGLMM(model_reduced)[1]
delta_rsquared_val

summary(model)
summary(model_reduced)
anova(model)
anova(model_reduced)
anova(model, model_reduced)
