rm(list=ls())
source("setup.R")

mergecols = c("subjectkey","eventname")

#load family data
fam = read.abcd(file.path(ABCDDataDir,"acspsw03.txt"))
fam = fam[,c(mergecols,"rel_family_id")]

#load site data
site = read.abcd(file.path(ABCDDataDir,"abcd_lt01.txt"))
site = site[,c(mergecols,"site_id_l")]

#load meanFD data and DOF totals
motion = read.csv(file.path(DataDir,"meanfd.csv"))
motion$dofavg = motion$dof/motion$runs
motion$censoredavg = motion$censored/motion$runs
motion$subjectkey = gsub("NDAR","NDAR_",motion$subjectkey)
motion$eventname = "baseline_year_1_arm_1"

#generate nuisance covariates from raw ABCD files
pdem = read.abcd(file.path(ABCDDataDir,"pdem02.txt"))
cols = c(mergecols,"interview_age","sex","demo_sex_v2",
         "demo_comb_income_v2","demo_prnt_marital_v2",
         "demo_Prnt_ed_v2","demo_prtnr_ed_v2","demo_ethn_v2",
         paste0("demo_race_a_p___",c(10:25,77,99)))
pdem = pdem[,cols]
pdem$re_white = as.numeric(pdem$demo_race_a_p___10==1)
pdem$re_black = as.numeric(pdem$demo_race_a_p___11==1)
pdem$re_hisp = as.numeric(pdem$demo_ethn_v2==1)
pdem$re_hisp[is.na(pdem$re_hisp)] = 0

pdem$re_asian = as.numeric(rowSums(pdem[,c(paste0("demo_race_a_p___",c(14:24)))])>0)
pdem$re_other = as.numeric(rowSums(pdem[,c(paste0("demo_race_a_p___",c(12:13,25)))])>0)
pdem$in_lt50 = as.numeric(pdem$demo_comb_income_v2<=6)
pdem$in_50_100 = as.numeric(pdem$demo_comb_income_v2>=7 & pdem$demo_comb_income_v2<=8)
pdem$in_gt100 = as.numeric(pdem$demo_comb_income_v2>=9 & pdem$demo_comb_income_v2<=10)

pdem$ed_prnt = pdem$demo_prnt_ed_v2
pdem$ed_prtnr = pdem$demo_prtnr_ed_v2
pdem$ed_prnt[pdem$ed_prnt>100] = NA
pdem$ed_prtnr[pdem$ed_prtnr>100] = NA

pdem$ed_max = apply(pdem[,c("ed_prnt","ed_prtnr")],1,max,na.rm=T)
pdem$ed_max[is.infinite(pdem$ed_max)] = NA

pdem$ed_lths = as.numeric(pdem$ed_max<13)
pdem$ed_bac = as.numeric(pdem$ed_max==18)
pdem$ed_hs = as.numeric(pdem$ed_max==13 | pdem$ed_max==14)
pdem$ed_grad = as.numeric(pdem$ed_max>=19 & pdem$ed_max<=21)
pdem$ed_sc = as.numeric(pdem$ed_max>=15 & pdem$ed_max<=17)

pdem$ms = pdem$demo_prnt_marital_v2
pdem$ms[pdem$ms>100] = NA
pdem$ms = as.numeric(pdem$ms==1)

temp = with(pdem,re_asian + 2*re_black + 4*re_hisp + 8*re_other + 16*re_white)
temp[temp==1] = "Asian"
temp[temp==2] = "Black"
temp[temp==4] = "Hispanic"
temp[temp==8] = "Other"
temp[temp==16] = "White"
b=1:31
d = c(1:31)[bitwAnd(b,4)>0]
temp[temp %in% d] = "Hispanic"
temp[!(temp %in% c("Asian","Black","Hispanic","Other","White"))] = "Other"
pdem$RaceEthnicity = as.factor(temp)

temp = with(pdem,in_lt50 + 2*in_50_100 + 4*in_gt100)
temp[temp==1] = "[<50K]"
temp[temp==4] = "[>=100K]"
temp[temp==2] = "[>=50K & <100K]"
temp[temp==0] = NA
pdem$HouseholdIncome = as.factor(temp)

temp = with(pdem, ed_lths + 2*ed_bac + 4*ed_hs + 8*ed_grad + 16*ed_sc)
temp[temp==1] = "< HS Diploma"
temp[temp==2] = "Bachelor"
temp[temp==4] = "HS Diploma/GED"
temp[temp==8] = "Post Graduate Degree"
temp[temp==16] = "Some College"
pdem$HighestParentalEducation = temp

pdem$Gender = pdem$sex

pdem$Age = pdem$interview_age/12

pdem$ms = pdem$demo_prnt_marital_v2
pdem$ms[pdem$ms>100] = NA
pdem$ms = as.numeric(pdem$ms==1)

temp = pdem$ms
temp[temp==0] = "no"
temp[temp==1] = "yes"
pdem$HouseholdMaritalStatus = temp

pdem = pdem[,c(mergecols,"Age","Gender","RaceEthnicity","HighestParentalEducation","HouseholdMaritalStatus","HouseholdIncome")]

gp = read.csv(file.path(DataDir,"ABCD_phenotypic_data_BPM.csv"))
gp = gp[,c(mergecols,"PF10","PF10_EXT","PF10_INT")]

#propensity weights
pw = read.abcd(file.path(ABCDDataDir,"acspsw03.txt"))
pw = pw[,c(mergecols,"acs_raked_propensity_score")]

#merge everything together
data = multi.merge(fam,site,motion,pdem,gp,pw,by=mergecols)

data = data[data$eventname=="baseline_year_1_arm_1",]

tr = 0.8
data$GoodTime = (data$TRs - data$censored)*tr/60
data$Include.rest = data$GoodTime >= 8 & data$runs>=2

data$Gender[data$Gender==""] = NA

data$Include.demo = !is.na(data$Age) & !is.na(data$Gender) & data$Gender!="" & !is.na(data$RaceEthnicity) & !is.na(data$HighestParentalEducation) & !is.na(data$HouseholdMaritalStatus) & !is.na(data$HouseholdIncome)

data$Include.p = !is.na(data$PF10)

#exclude for any data problems, like NaN ROIs etc
nansubs = read.csv(file.path(DataDir,"nan_subs.csv"))
data$Include.data = !(data$subjectkey %in% nansubs$subjectkey)

data$Include = data$Include.rest & data$Include.demo & data$Include.p & data$Include.data

sum(data$Include,na.rm=T)

data = data[data$Include==T & !is.na(data$Include),]
data$Subject = gsub("NDAR_","NDAR",data$subjectkey)

sum(data$Include)

data$IncludeUnrelated = 1

t = table(data$rel_family_id)
w = which(t>1)

set.seed(12345)
for (i in 1:length(w)) {
  f = names(w)[i]
  idx = which(data$rel_family_id == f)
  s = sample(length(idx),length(idx)-1)
  data$IncludeUnrelated[idx[s]] = 0
}

#check for families that cross site
t = table(data$rel_family_id,data$site_id_l)
t = t>0
sum(rowSums(t)>1)

data$abcd_site = as.character(data$site_id_l)
table(data$abcd_site)
data$abcd_site_num = as.numeric(gsub('site','',data$abcd_site))

t = table(data$abcd_site)
sites = names(t)[t>=75]
data = data[data$abcd_site %in% sites,]

dim(data)

#merge in loso p factors
pfactors = read.csv("./Data/ABCD_lavaan_pfactor_loso.csv")
pfactors = pfactors[,c("subjectkey",paste0("P",sort(unique(data$abcd_site_num))))]

data = merge(data,pfactors,by="subjectkey",all.x=T)

write.csv(data,file.path(DataDir,"ABCD_rest.csv"),row.names=FALSE,na="NaN")
