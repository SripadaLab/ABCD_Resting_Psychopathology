%%
%utility funcions
addpath /home/slab/users/mangstad/repos/Misc_utils/

Exp = '/net/pepper/ABCD/CIFTI/';

%%
%Setup files
DataFile = [Exp '/Scripts/pfactor/Data/ABCD_rest.csv'];
CorrTemplate = [Exp '/Rest/Gordon_Sub_Cere/[Subject].txt'];
NetsFile = [Exp '/Scripts/pfactor/Data/gordon_sub_cere_parcels.csv'];
Perms = [Exp '/Scripts/pfactor/Data/perms_5880.mat'];

LowMotion = 0;
IncludeTwins = 1;
RemoveSite = 0;

%%
%setup
dat = readtable(DataFile);

if (LowMotion==1)
    dat = dat(dat.fd<0.2,:);
    Perms = [Exp '/Scripts/pfactor/Data/perms_2853.mat'];
end

perms = load(Perms);
perms = perms.pset;

%to drop twins
if (IncludeTwins==0)
    dat = dat(dat.IncludeUnrelated==1,:);
end

N = size(dat,1);
nROI = 418;
P = (nROI*(nROI-1))/2;

netsfile = readtable(NetsFile);
netsfile = netsfile(1:nROI,:);
nets = netsfile.NetworkNumber;

%%
%load connectomes
featuremat = zeros(N,P);

parfor iSubject = 1:N
    Subject = dat.Subject{iSubject};
    fprintf(1,'%d\n',iSubject);
    %file = mc_GenPath(CorrTemplate);
    file = strrep(CorrTemplate,'[Subject]',Subject);
    tmp = load(file);
    tmp = tmp(1:418,1:418);
    tmp = mc_flatten_upper_triangle(tmp);
    featuremat(iSubject,:) = mc_FisherZ(tmp);
end

%%
%setup folds
u = unique(dat.abcd_site_num);
nFold = numel(u);
fold_site_conversion = [[1:nFold]' u];
folds = zeros(size(dat.abcd_site_num));
sitesize = zeros(nFold,1);
for iFold = 1:nFold
    folds(dat.abcd_site_num==u(iFold)) = iFold;
    sitesize(iFold) = sum(folds==iFold);
end

if (RemoveSite==1) 
    dummy_site = mc_Vector2Mask(folds);
    [~,i] = max(sum(dummy_site));
    dummy_site(:,i) = [];
    X = [ones(size(featuremat,1),1) dummy_site];
    b = pinv(X'*X)*X'*featuremat;
    res = featuremat - X*b;
    featuremat_orig = featuremat;
    featuremat = res;
    
    mainpheno_orig = mainpheno;
    b = pinv(X'*X)*X'*mainpheno;
    res = mainpheno - X*b;
    mainpheno = res;
end

%%
%setup leave one site out variables to predict

losopheno = [];
for i = 1:nFold
    name = sprintf('P%d',fold_site_conversion(i,2));
    losopheno = [losopheno dat.(name)];
end
%logmin = abs(min(losopheno))+1;
%logpheno = log(bsxfun(@plus,losopheno,logmin));
mainpheno = [dat.PF10 log(dat.PF10 + abs(min(dat.PF10))+1)];

%%
%get nuisance variables
u = unique(dat.RaceEthnicity);
re = zeros(N,numel(u));
for i = 1:numel(u)
    idx = strcmp(dat.RaceEthnicity,u{i});
    re(idx,i) = 1;
end
s = sum(re);
[~,i] = max(s);
re(:,5) = [];

u = unique(dat.Gender);
gen = zeros(N,numel(u));
for i = 1:numel(u)
    idx = strcmp(dat.Gender,u{i});
    gen(idx,i) = 1;
end
s = sum(gen);
[~,i] = max(s);
gen(:,2) = [];

u = unique(dat.HighestParentalEducation);
hpe = zeros(N,numel(u));
for i = 1:numel(u)
    idx = strcmp(dat.HighestParentalEducation,u{i});
    hpe(idx,i) = 1;
end
s = sum(hpe);
[~,i] = max(s);
hpe(:,i) = [];

u = unique(dat.HouseholdIncome);
hi = zeros(N,numel(u));
for i = 1:numel(u)
    idx = strcmp(dat.HouseholdIncome,u{i});
    hi(idx,i) = 1;
end
s = sum(hi);
[~,i] = max(s);
hi(:,i) = [];

u = unique(dat.HouseholdMaritalStatus);
hms = zeros(N,numel(u));
for i = 1:numel(u)
    idx = strcmp(dat.HouseholdMaritalStatus,u{i});
    hms(idx,i) = 1;
end
s = sum(hms);
[~,i] = max(s);
hms(:,i) = [];

%normal and expanded nuisance
nuisance = [dat.Age dat.Age.^2 gen dat.fd dat.fd.^2 re];
nuisancefull = [dat.Age dat.Age.^2 gen dat.fd dat.fd.^2 re hpe hi hms];


