%%
%utility funcions
addpath /home/slab/users/mangstad/repos/Misc_utils/

NumComp = 250;
Perms = 10000;

%%
%run leave one site out PCA and BBS for all phenotypes with regular
%nuisance
results_p = mc_bbs(featuremat,mainpheno,nuisance,folds,NumComp,'LOSOPheno',0);
results_losop = mc_bbs(featuremat,losopheno,nuisance,folds,NumComp,'LOSOPheno',1,'Scores',results_p.Aa);
%resultsp_log = mc_bbs(featuremat,logpheno,nuisance,folds,NumComp,'LOSOPheno',1,'Scores',resultsp.Aa);

%%
%now run expanded nuisance using PCA from above
results_p_full = mc_bbs(featuremat,mainpheno,nuisancefull,folds,NumComp,'Scores',results_p.Aa,'LOSOPheno',0);
%resultsp_log_full = mc_bbs(featuremat,logpheno,nuisancefull,folds,NumComp,'Scores',resultsp.Aa,'LOSOPheno',1);

%%
%now run 10k perms with regular nuisance
results_p_perms = mc_bbs_perm(featuremat,mainpheno,nuisance,folds,NumComp,Perms,'Scores',results_p.Aa,'LOSOPheno',0,'Perms',perms,'Tolerance',1e-10);

%%
%now run 10k perms with expanded nuisance
results_p_full_perms = mc_bbs_perm(featuremat,mainpheno(:,1),nuisancefull,folds,NumComp,Perms,'Scores',results_p.Aa,'LOSOPheno',0,'Perms',perms,'Tolerance',1e-10);

%%
%now run 10k perms with drop site p
results_p_losop_perms = mc_bbs_perm(featuremat,losopheno,nuisancefull,folds,NumComp,Perms,'Scores',results_p.Aa,'LOSOPheno',1,'Perms',perms,'Tolerance',1e-10);
