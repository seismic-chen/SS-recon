function [tt, f]=ss_tt_table(dist,depth)
% Calculate the travel-time table for SS precursors at various depths.
% Input: 
% dist: distance vector 
% depth: depth vector
% Output:
% tt: a nx*nd matrix of travel times of SdS phases
% f: an interpolation function
% Example:
% dist = 95:5:170;
% depth = 0:5:1000; 
% [tt, f]=ss_tt_table(dist,depth);
% Jan. 23, 2022, Yunfeng Chen, Global Seismology Group, Zhejiang University
tt=zeros(length(dist),length(depth));
for i = 1:length(dist)
    x=dist(i)/2;
    times=taupTime('prem',0,'S','deg',x);
    indices = find(strcmp({times.phaseName},'S'));
    tt0=times(indices(1)).time;
    for j = 1:length(depth)
        d=depth(j);
        times=taupTime('prem',d,'S','deg',x);
        indices = find(strcmp({times.phaseName},'S'));
        tt(i,j)=2*times(indices(1)).time-2*tt0;
    end
end
[X,Y]=meshgrid(dist,depth);
Z = tt';
f=scatteredInterpolant(X(:),Y(:),Z(:),'linear');