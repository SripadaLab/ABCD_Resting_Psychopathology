centroids = readtable([Exp '/Scripts/pfactor/Data/gordon_sub_cere_parcels.csv']);
centroids = centroids(1:418,:);
ROIcoord = [centroids.Centroid_X centroids.Centroid_Y centroids.Centroid_Z];
nROI = size(centroids,1);
distanceMat = zeros(nROI);

for i = 1:(nROI-1)
    for j = (i+1):nROI
        distanceMat(i,j) = sqrt((ROIcoord(i,1)-ROIcoord(j,1))^2+(ROIcoord(i,2)-ROIcoord(j,2))^2+(ROIcoord(i,3)-ROIcoord(j,3))^2);
    end
end

distance=mc_flatten_upper_triangle(distanceMat);
nROIpair=length(distance);

addpath('/home/slab/users/yfang/MethodsCore2/QualityChecks/CheckMotion')
graph.size=ones(1,nROIpair)*4;
graph.size(distanceMask) = 1;
graph.iffit=1;
graph.titleprefix = 'site';

QC_RSFC_all = corr(dat.fd,featuremat);
QC_RSFC_plot(distance,QC_RSFC_all,'ABCD',graph,'all');
