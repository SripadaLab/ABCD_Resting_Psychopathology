%need to add the Misc_utils repo path for visualization code
addpath /home/slab/users/mangstad/repos/Misc_utils

NetsFile = [Exp '/Scripts/pfactor/Data/gordon_sub_cere_parcels.csv'];
nROI = 418;

netsfile = readtable(NetsFile);
netsfile = netsfile(1:nROI,:);
nets = netsfile.NetworkNumber;

comps = load('Results/components1-75.txt');

%visualize a -1/0/1 z score 2 thresholded component
plot_jica_component(comps(:,5)',1,1,2,nets,'component 5',[1:16]);

%visualize an -1/0/1 original value 0.01 thresholded component only in a
%subset of networks
plot_jica_component(comps(:,5)',1,0,0.01,nets,'component 5',[1:12]);

%visualize -1/0/1 top 10% of connections 
plot_jica_component(comps(:,5)',1,2,0.9,nets,'component 5',[1:16]);

%visualize an unthresholded but weighted component
mc_plot_connectome(comps(:,5)',nets);
