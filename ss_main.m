% Aug. 13, 2021, processing the SS precursor data
% Yunfeng Chen, Global Seismology Group, Zhejiang University
clear; close all; clc;
%% Read in data and calculate bounce point (midpoint)
javaaddpath('./FMI/lib/FMI.jar');
addpath('./FMI/matTaup');
addpath ~/MATLAB/codes/xcorr_example/
%% load the data
datadir = '/home/yunfeng/project/western_pacific/data/*';
event = dir(datadir);
event(1:3)=[];
nevt = length(event);
ss=[];
for i = 1:nevt
    if mod(i,100) == 0
        disp(['Reading ',num2str(i),'/',num2str(nevt),' events']);
    end
    sacfiles=dir(fullfile(event(i).folder,event(i).name,'*.T'));
    for j = 1:length(sacfiles)
        [t,data,SAChdr] = fget_sac(fullfile(sacfiles(j).folder,sacfiles(j).name));
        tmp.d = data;
        tmp.t = t;
        tmp.stla=SAChdr.station.stla;
        tmp.stlo=SAChdr.station.stlo;
        tmp.stel=SAChdr.station.stel;
        tmp.sta =SAChdr.station.kstnm;
        tmp.evla=SAChdr.event.evla;
        tmp.evlo=SAChdr.event.evlo;
        tmp.evdp=SAChdr.event.evdp/1000.; % meter to kilometer
        tmp.mag=SAChdr.event.mag;
        tmp.dist=SAChdr.evsta.dist;
        tmp.az=SAChdr.evsta.az;
        tmp.baz=SAChdr.evsta.baz;
        tmp.gcarc=SAChdr.evsta.gcarc;
        % calculate the bounce point (midpoint)
        [tmp.bplat,tmp.bplon]=gc_midpoint(tmp.evla, tmp.evlo, tmp.stla, tmp.stlo);
        % calculate SS arrival time
        times=taupTime('prem',tmp.evdp,'SS,S^410S,S^660S','sta',[tmp.stla tmp.stlo],...
            'evt',[tmp.evla,tmp.evlo]);
        tmp.t660=times(1).time;
        tmp.t410=times(2).time;
        tmp.tss=times(3).time;
        ss=[ss tmp];
    end
end
% save as mat
save 'ss.mat' 'ss';

figure;
plot([ss.bplon],[ss.bplat],'.')
k=2;
figure;
plot(ss(k).t,ss(k).d/max(ss(100).d)); hold on;
plot([ss(k).t410,ss(k).t410],[-1,1],'--r');
plot([ss(k).t660,ss(k).t660],[-1,1],'--r')
plot([ss(k).tss,ss(k).tss],[-1,1],'--r')
%% calculate SNR
for k=1:length(ss)
    d=ss(k).d;
    t=ss(k).t;
    tss=ss(k).tss;
    ss(k).snr = ss_snr(d,t,tss);
end
%% check polarity reversal
for k=1:length(ss)
    d=ss(k).d;
    t=ss(k).t;
    tss=ss(k).tss;
    [d,is_reversal] = ss_check_polarity(d,t,tss);
    ss(k).is_reversal = is_reversal;
    if is_reversal
        ss(k).d = d;
    end
end
remove = [ss.snr]<=3;
ss(remove) = [];
%% binning
dx=5; dy=5; dh=2;
xmin=110; ymin=20; hmin=100;
xmax=160; ymax=60; hmax=170;
% define grid center
x = xmin+dx/2:dx:xmax;
y = ymin+dy/2:dy:ymax;
h = hmin+dh/2:dh:hmax;
nx=length(x); ny=length(y); nh=length(h); nt = length(t);
disp('Binning')
% 5D case
% d1=zeros(nt, nx, ny, nh, nphi);
% 4D case
d1 = zeros(nt, nx, ny, nh);
fold_map=zeros(nx,ny,nh);
flow=1/75.;
fhigh=1/15.;
for n=1:length(ss)
    j=floor((ss(n).bplat-ymin)/dy)+1;
    i=floor((ss(n).bplon-xmin)/dx)+1;
    k=floor((ss(n).gcarc-hmin)/dh)+1;
%     l=floor(ss(n).phi/dphi)+1;
    fold_map(i,j,k)=fold_map(i,j,k)+1;
    d=ss(n).d;
    % bandpass filter
    d_filt=bandpassSeis(d,1,flow,fhigh,3);
    d1(:,i,j,k)=d1(:,i,j,k)+d_filt(1:nt);
end
% plot fold map
fold_map_xy=sum(fold_map,3);
figure;
set(gcf,'Position',[100 100 1000 800],'color','w')
imagesc(x,y,fold_map_xy'); hold on;
plot([ss.bplon],[ss.bplat],'k.'); 
xlabel('Easting (km)');
ylabel('Northing (km)');
colorbar;
set(gca,'fontsize',14)
axis equal;
%% apply cross-correlation to all stacked traces
din = reshape(d1,nt,nx*ny*nh); % flatten the tensor
N=5; % number of iteration for cross-correlation measurments
t=0:nt-1;
times = repmat(t(:),1,size(din,2));
t0=ones(1,size(din,2))*900; % SS arrival
xwin=[-100 100]; % cross-correlation window
maxlag=50; % maximum time lag
is_plot=1; % flag controls plotting
dout = ss_align(din,times,N,t0,xwin,maxlag,is_plot);
d1 = reshape(dout,nt,nx,ny,nh);
%% plotting
% plot one CMP gather
ix=4;
iy=6;
xaxis=h;
yaxis=0:nt-1;
d2d = squeeze(d1(:,ix,iy,:));
keep = any(d2d);
d2d(:,keep) = d2d(:,keep)*diag(1./max(d2d(:,keep)));
figure;
wigb(d2d(:,keep),1,xaxis(keep),yaxis)