%%
%Setup output file
OutputFile = [Exp '/Scripts/pfactor/Results/ABCD_rest_CIFTI_expressions.csv'];
NumComp = 250;

%%
%run PCA over everyone to get 250 components
[coeffall,~,~,~,expall] = pca(featuremat);

mu = mean(featuremat);
x = bsxfun(@minus,featuremat,mu);
Aall = x*coeffall(:,1:NumComp);

%%
%Save out file with first 250 component expressions for lme model in R
subjectkey = dat.subjectkey;
fid = fopen(OutputFile,'w');

fprintf(fid,'subjectkey,');
for i = 1:NumComp
    temp = sprintf('A%03d',i);
    fprintf(fid,'%s,',temp);
end
fprintf(fid,'\n');
for i = 1:numel(subjectkey)
    fprintf(fid,'%s,',subjectkey{i});
    fprintf(fid,'%15.15f,',Aall(i,:));
    fprintf(fid,'\n');
end
fclose(fid);


%%
%run lme in R
%s03_LME.R

%%
%consensus
betas = load([Exp '/Scripts/pfactor/Results/mm_betas.csv']);
cons = coeffall(:,1:NumComp)*betas;

%%
%calculate within network and DMN/TPN connections
zcons = zscore(cons);

inmask = zeros(418,418);
for i = 1:max(nets)
    inmask(nets==i,nets==i) = 1;
end
inmask = inmask==1;

dmntpn = zeros(418,418);
idx = [3,5,8,11,12];
for i = 1:numel(idx)
    for j = 1:numel(idx)
        dmntpn(nets==idx(i),nets==idx(j)) = 1;
        dmntpn(nets==idx(j),nets==idx(i)) = 1;
    end
end
dmntpn = (dmntpn - inmask)>0;

inmask = mc_flatten_upper_triangle(inmask);
dmntpn = mc_flatten_upper_triangle(dmntpn);

plot_jica_component(inmask,1,0,0,nets,'',[1:16]);
plot_jica_component(dmntpn,1,0,0,nets,'',[1:16]);

supra = abs(zcons)>2;
supra = supra';
100*sum(inmask)/numel(inmask)
100*sum(dmntpn)/numel(dmntpn)
100*sum(supra.*inmask)/sum(supra)
100*sum(supra.*dmntpn)/sum(supra)


%need to shift the grid lines in inkscape as matlab for some reason doesn't
%save the svg exactly as displayed in the figure. 1 pixel up, 1 pixel left
close all
plot_jica_component(cons',1,1,2,nets,'Mixed Model Consensus',[1:16]);
print([Exp '/Scripts/pfactor/Results/consensus.svg'],'-dsvg','-r600');
close all
plot_jica_component(sign(betas(5))*coeffall(:,5)',1,1,2,nets,'Component 5',[1:16]);
print([Exp '/Scripts/pfactor/Results/component005.svg'],'-dsvg','-r600');
close all
plot_jica_component(sign(betas(139))*coeffall(:,139)',1,1,2,nets,'Component 139',[1:16]);
print([Exp '/Scripts/pfactor/Results/component139.svg'],'-dsvg','-r600');
