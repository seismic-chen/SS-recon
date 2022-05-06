% May 2, 2022, Yunfeng Chen, Global Seismology Group, Zhejiang University
% Basic workflow for SS precursor processing and imaging
clear; close all; clc;
addpath ./utils//MatSAC/
addpath ./utils/
addpath ./ss/
javaaddpath ./utils/FMI/lib/FMI.jar
addpath ./utils/FMI/matTaup
addpath ./utils/m_map/
addpath ./utils/open-source/
%% Read in data and calculate bounce point (midpoint)
datadir = './data/*';
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
% save 'ss.mat' 'ss';
%% load the data if necessary
% load ss.mat;
figure;
plot([ss.bplon],[ss.bplat],'.')
k=2;
figure;
plot(ss(k).t,ss(k).d/max(ss(k).d)); hold on;
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
%% apply cross-correlation to all traces
nt=length(t);
for k=1:length(ss)
    ss(k).d=ss(k).d(1:nt);
end
din = [ss.d]; % flatten the tensor
N=5; % number of iteration for cross-correlation measurments
t=0:nt-1;
times = repmat(t(:),1,size(din,2));
t0=ones(1,size(din,2))*900; % SS arrival
xwin=[-100 100]; % cross-correlation window
maxlag=50; % maximum time lag
is_plot=0; % flag controls plotting
dout = ss_align_v2(din,times,N,t0,xwin,maxlag,is_plot);
for k=1:length(ss)
    ss(k).d=dout(:,k);
end
%% binning
% dx=5; dy=5; dh=2;
dx=2.5; dy=2.5; dh=2;
xmin=110; ymin=20; hmin=100;
xmax=160; ymax=60; hmax=170;
lonlim=[xmin xmax];
latlim=[ymin ymax];
% define grid center
x = xmin+dx/2:dx:xmax;
y = ymin+dy/2:dy:ymax;
h = hmin+dh/2:dh:hmax;
t = ss(1).t;
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
    d_filt=d_filt/max(abs(d_filt));
    d1(:,i,j,k)=d1(:,i,j,k)+d_filt(1:nt);
end
% nomalization
for i=1:nx
    for j=1:ny
        for k=1:nh
            if fold_map(i,j,k)>0
               d1(:,i,j,k)=d1(:,i,j,k)/fold_map(i,j,k); 
            end
        end
    end
end
% plot fold map
fold_map_xy=sum(fold_map,3);
% figure;
% set(gcf,'Position',[100 100 1000 800],'color','w')
% imagesc(x,y,fold_map_xy'); hold on;
% cm=colormap('gray');
% colormap(flipud(cm));
% colorbar;
% axis tight;
% plot([ss.bplon],[ss.bplat],'b.'); 
% xlabel('Longitude (deg)');
% ylabel('Latitude (deg)');
% colorbar;
% set(gca,'fontsize',14)
% axis equal;
d4d=d1;
save 'd4d_ss.mat' d4d h t x y fold_map
figure
imagesc(squeeze(d4d(:,10,10,:)))
%% plot fold map
figure;
set(gcf,'Position',[100 100 1600 800],'color','w')
subplot(121)
m_proj('lambert','long', lonlim, 'lat', latlim); hold on;
set(gcf,'color','w')
[X,Y]=meshgrid(x,y);
hs=m_scatter([ss.bplon],[ss.bplat],5,'b','filled');
alpha(hs,0.5);
m_gshhs('i','line','color','k','linewidth',1)
m_gshhs('lb2','line','color','k')
m_grid('linewidth',2,'tickdir','out',...
    'xaxisloc','bottom','yaxisloc','left','fontsize',24);
text(-0.12,0.98,'(a)','Units','normalized','FontSize',32)

subplot(122)
m_proj('lambert','long', lonlim, 'lat', latlim); hold on;
[X,Y]=meshgrid(x,y);
hh=m_pcolor(X,Y,fold_map_xy');
set(hh,'edgecolor','none')
cm=colormap('gray');
colormap(flipud(cm));
caxis([0 300])
m_gshhs('i','line','color','k','linewidth',1)
m_gshhs('lb2','line','color','k')
m_grid('linewidth',2,'tickdir','out',...
    'xaxisloc','bottom','yaxisloc','left','fontsize',24);
hh=colorbar('h');
set(hh,'fontsize',24);
set(hh,'Position',[0.6 0.1256 0.3 0.0250])
xlabel(hh,'Count');
text(-0.12,0.98,'(b)','Units','normalized','FontSize',32)
% title('Fold map','fontsize',24)
%% stack the CMP bin 
for n=1:nh
    times=taupTime('ak135',10,'SS,S^410S,S^660S','deg',h(n));
    indices = find(strcmp({times.phaseName},'S^660S'));
    t660(n)=times(indices(1)).time;
    indices = find(strcmp({times.phaseName},'S^410S'));
    t410(n)=times(indices(1)).time;
    indices = find(strcmp({times.phaseName},'SS'));
    tss(n)=times(indices(1)).time;
end

d2d = squeeze(mean(mean(d1,3),2)); % simple averaging
W = any(d1);    % obtain the non-zero trace
w = squeeze(sum(sum(W,3),2));  % calcualte the weight
d2d_w = squeeze(sum(sum(d1,3),2))*diag(1./w); % weighted averaging

% find the time of SS phase and set it to 0 time
[~,index] = max(sum(d2d,2));
tshift = t(index);
t=t-tshift;
ntraces = squeeze(sum(sum(fold_map,2),1));
%% plot results from two averaging approaches
figure;
subplot(511)
bar(h,ntraces)
subplot(5,1,2:5)
set(gcf,'Position',[0 0 1000 1000],'Color','w')
wigb(d2d,10,h,t)
plot(h,t660-tss,'--r')
plot(h,t410-tss,'--r')
axis xy
ylim([-500 100])
ylabel('Time (s)')
xlabel('Distance (deg)')
set(gca,'fontsize',14)

figure;
subplot(511)
bar(h,ntraces)
subplot(5,1,2:5)
set(gcf,'Position',[0 0 1000 1000],'Color','w')
wigb(d2d_w,10,h,t)
plot(h,t660-tss,'--r')
plot(h,t410-tss,'--r')
axis xy
ylim([-500 100])
ylabel('Time (s)')
xlabel('Distance (deg)')
set(gca,'fontsize',14)
%% conduct NMO correction
% calcualte travel time table
dist = 95:5:170;
depth = 0:5:1000; 
[tt, f]=ss_tt_table(dist,depth);
d1_nmo=zeros(size(d1));
h0=135;
is_plot=false;
for i=1:nx
    for j=1:ny
        for k=1:nh
            din = d1(:,i,j,k);
            if any(d)
                [dout,t410_ref,t660_ref] = ss_nmo_v2(din,t,h(k),h0,f,is_plot);
                d1_nmo(:,i,j,k)=dout;
            end
        end
    end
end
%% perform stacking for each CMP gather
d3d_nmo = zeros(nt,nx,ny);
d3d = zeros(nt,nx,ny);
for i=1:nx
    for j=1:ny
        % move-out corrected cmp
        d_cmp = squeeze(d1_nmo(:,i,j,:));
        nstack = sum(any(d_cmp));
        if nstack>0
            d_stack = sum(d_cmp,2)/nstack;
            d_stack = d_stack/rms(d_stack);
            d3d_nmo(:,i,j)=d_stack;
        end
        % non-move-out corrected cmp
        d_cmp = squeeze(d1(:,i,j,:));
        nstack = sum(any(d_cmp));
        if nstack>0
            d_stack = sum(d_cmp,2)/nstack;
            d_stack = d_stack/rms(d_stack);
            d3d(:,i,j)=d_stack;
        end
    end
end
%% plot all stacked CMP gather
figure;
set(gcf,'Position',[100 100 1600 500],'color','w')
% wigb(reshape(d3d_nmo,nt,nx*ny),20,1:nx*ny,t)
imagesc(1:nx*ny,t,reshape(d3d_nmo,nt,nx*ny)); hold on;
colormap(seismic(3));
caxis([-0.3 0.3])
plot(1:nx*ny,t410_ref*ones(1,nx*ny),'--r','linewidth',2)
plot(1:nx*ny,t660_ref*ones(1,nx*ny),'--r','linewidth',2)
axis xy
ylim([-300 100])
ylabel('Time (s)')
xlabel('CMP#')
title('NMO corrected post-stacked data')
set(gca,'fontsize',14)
figure;
set(gcf,'Position',[100 100 1600 500],'color','w')
% wigb(reshape(d3d,nt,nx*ny),20,1:nx*ny,t)
imagesc(1:nx*ny,t,reshape(d3d,nt,nx*ny)); hold on;
colormap(seismic(3));
caxis([-0.3 0.3])
plot(1:nx*ny,t410_ref*ones(1,nx*ny),'--r','linewidth',2)
plot(1:nx*ny,t660_ref*ones(1,nx*ny),'--r','linewidth',2)
axis xy
ylim([-300 100])
ylabel('Time (s)')
xlabel('CMP#')
title('Raw post-stacked data')
set(gca,'fontsize',14)
%% SSA
save 'd3d_ss.mat' d3d h t x y 
% addpath ~/MATLAB/yangkang/open-source/
keep = t>=-300 & t<=-100;
din = d3d_nmo(keep,:,:);
mask = ones(size(din));
[nt,nx,ny] = size(din);
for ix=1:nx
    for iy=1:ny
        a = sum(din(:,ix,iy));
        if a == 0
            mask(:,ix,iy) = 0;
        else
            din(:,ix,iy)=din(:,ix,iy)/rms(din(:,ix,iy));
        end
    end
end
dt=1;
flow = 1/75.;
fhigh= 1/15.;
% flow = 0.001;
% fhigh = 0.45;
N=5; Niter=20; mode=1;verb=1;
coef=(Niter-(1:Niter))/(Niter-1); %linearly decreasing
dout=fxymssa_recon(din,mask,flow,fhigh,dt,N,Niter,eps,verb,mode,coef);
% K=5;
% dout=fxydmssa_recon(din,mask,flow,fhigh,dt,N,K,Niter,eps,verb,mode,coef);
disp('SSA Denoising and Interpolation Complete!');
%% plot denosing result
figure;
set(gcf,'Position',[100 100 1600 500],'color','w')
% wigb(reshape(din,nt,nx*ny),1,1:nx*ny,t(keep))
imagesc(1:nx*ny,t(keep),reshape(din,nt,nx*ny)); hold on;
colormap(seismic(3));
cmax=2*rms(din(:));
caxis([-cmax cmax])
plot(1:nx*ny,t410_ref*ones(1,nx*ny),'--r','linewidth',2)
plot(1:nx*ny,t660_ref*ones(1,nx*ny),'--r','linewidth',2)
axis xy
ylim([-300 -100])
ylabel('Time (s)')
xlabel('CMP#')
title('Raw post-stacked data')
set(gca,'fontsize',14)

figure;
set(gcf,'Position',[100 100 1600 500],'color','w')
% wigb(reshape(dout,nt,nx*ny),1,1:nx*ny,t(keep))
imagesc(1:nx*ny,t(keep),reshape(dout,nt,nx*ny)); hold on;
colormap(seismic(3));
caxis([-cmax cmax])
plot(1:nx*ny,t410_ref*ones(1,nx*ny),'--r','linewidth',2)
plot(1:nx*ny,t660_ref*ones(1,nx*ny),'--r','linewidth',2)
axis xy
ylim([-300 -100])
ylabel('Time (s)')
xlabel('CMP#')
title('Denoised post-stacked data')
set(gca,'fontsize',14)
%% time to depth conversion
z0=0:1:1000; 
h0=135;
x0=h0*ones(size(z0));
t0=f(x0,z0);
ti=t(keep);
nz=length(z0);
d3d_depth_ssa=zeros(nz,nx,ny);
d3d_depth=zeros(nz,nx,ny);
for i = 1:nx
    for j = 1:ny
        d=din(:,i,j);
        dtmp=interp1(ti,d,t0,'linear',0);
        d3d_depth(:,i,j)=dtmp;
        % denoising result
        d=dout(:,i,j);
        dtmp=interp1(ti,d,t0,'linear',0);
        d3d_depth_ssa(:,i,j)=dtmp;
    end
end
%% find the maximum amplitude in the given depth intervals
depth0=z0;
drange1=[390,430];
drange2=[640,680];

amp410=zeros(nx,ny);
amp660=zeros(nx,ny);
d410=zeros(nx,ny);
d660=zeros(nx,ny);

amp410_ssa=zeros(nx,ny);
amp660_ssa=zeros(nx,ny);
d410_ssa=zeros(nx,ny);
d660_ssa=zeros(nx,ny);
for i=1:nx
    for j=1:ny
        d0=squeeze(d3d_depth(:,i,j));
        keep1=find(depth0>=drange1(1) & depth0<=drange1(2));
        keep2=find(depth0>=drange2(1) & depth0<=drange2(2));
        % find the maximum amplitude
        [amp410(i,j),i1]=max(d0(keep1));
        [amp660(i,j),i2]=max(d0(keep2));
        ind410=keep1(i1);
        ind660=keep2(i2);
        d410(i,j)=depth0(ind410);
        d660(i,j)=depth0(ind660);
        
        d0=squeeze(d3d_depth_ssa(:,i,j));
        keep1=find(depth0>=drange1(1) & depth0<=drange1(2));
        keep2=find(depth0>=drange2(1) & depth0<=drange2(2));
        % find the maximum amplitude
        [amp410_ssa(i,j),i1]=max(d0(keep1));
        [amp660_ssa(i,j),i2]=max(d0(keep2));
        ind410=keep1(i1);
        ind660=keep2(i2);
        d410_ssa(i,j)=depth0(ind410);
        d660_ssa(i,j)=depth0(ind660);
    end
end
d410=d410';
d660=d660';
d410_ssa=d410_ssa';
d660_ssa=d660_ssa';
%% plot 410
lonlim=[xmin xmax];
latlim=[ymin ymax];

[Xgrid,Ygrid]=meshgrid(x,y);
% smooth the results
ngrid=2;
K = (1/ngrid^2)*ones(ngrid,ngrid);
d410_smooth = conv2(d410,K,'same');
d410_ssa_smooth = conv2(d410_ssa,K,'same');

Vgrid=d410_smooth;
Vgrid_ssa=d410_ssa_smooth;
F410=scatteredInterpolant(Xgrid(:),Ygrid(:),Vgrid(:),'natural','none');
F410_ssa=scatteredInterpolant(Xgrid(:),Ygrid(:),Vgrid_ssa(:),'natural','none');

xi=lonlim(1):0.2:lonlim(2);
yi=latlim(1):0.2:latlim(2);
[XI,YI]=meshgrid(xi,yi);
d410_smooth=F410(XI,YI);
d410_ssa_smooth=F410_ssa(XI,YI);

figure;
set(gcf,'Position',[100 100 1600 800],'color','w')
subplot(121)
m_proj('lambert','long', lonlim, 'lat', latlim); hold on;
set(gcf,'color','w')
h=m_pcolor(XI,YI,d410_smooth);
set(h,'edgecolor','none')
cm=colormap(jet);
colormap(flipud(cm));
caxis([390 440])

m_gshhs('i','line','color','k','linewidth',1)
m_gshhs('lb2','line','color','k')
m_grid('linewidth',2,'tickdir','out',...
    'xaxisloc','bottom','yaxisloc','left','fontsize',24);
hh=colorbar('h');
set(hh,'fontsize',24);
set(hh,'Position',[0.3065 0.1256-0.045 0.4220 0.0250])
xlabel(hh,'Depth (km)');
title('MTZ 410 depth ?raw)','fontsize',24)
text(-0.1,0.98,'(a)','Units','normalized','FontSize',32)

subplot(122)
m_proj('lambert','long', lonlim, 'lat', latlim); hold on;
set(gcf,'color','w')
h=m_pcolor(XI,YI,d410_ssa_smooth);
set(h,'edgecolor','none')
cm=colormap(jet);
colormap(flipud(cm));
caxis([390 440])

m_gshhs('i','line','color','k','linewidth',1)
m_gshhs('lb2','line','color','k')
m_grid('linewidth',2,'tickdir','out',...
    'xaxisloc','bottom','yaxisloc','left','fontsize',24);
hh=colorbar('h');
set(hh,'fontsize',24);
set(hh,'Position',[0.3065 0.1256-0.045 0.4220 0.0250])
xlabel(hh,'Depth (km)');
title('MTZ 410 depth (denoised)','fontsize',24)
text(-0.1,0.98,'(b)','Units','normalized','FontSize',32)
%% plot 660
K = (1/ngrid^2)*ones(ngrid,ngrid);
d660_smooth = conv2(d660,K,'same');
d660_ssa_smooth = conv2(d660_ssa,K,'same');
Vgrid=d660_smooth;
Vgrid_ssa=d660_ssa_smooth;
F660=scatteredInterpolant(Xgrid(:),Ygrid(:),Vgrid(:),'natural','none');
F660_ssa=scatteredInterpolant(Xgrid(:),Ygrid(:),Vgrid_ssa(:),'natural','none');

xi=lonlim(1):0.2:lonlim(2);
yi=latlim(1):0.2:latlim(2);
[XI,YI]=meshgrid(xi,yi);
d660_smooth=F660(XI,YI);
d660_ssa_smooth=F660_ssa(XI,YI);

figure;
set(gcf,'Position',[100 100 1600 800],'Color','w')
subplot(121)
m_proj('lambert','long', lonlim, 'lat', latlim); hold on;
set(gcf,'color','w')
h=m_pcolor(XI,YI,d660_smooth);
set(h,'edgecolor','none')
cm=colormap(jet);
colormap(flipud(cm));
caxis([630 690])

m_gshhs('i','line','color','k','linewidth',1)
m_gshhs('lb2','line','color','k')
m_grid('linewidth',2,'tickdir','out',...
    'xaxisloc','bottom','yaxisloc','left','fontsize',24);
hh=colorbar('h');
set(hh,'fontsize',24);
set(hh,'Position',[0.3065 0.1256-0.045 0.4220 0.0250])
xlabel(hh,'Depth (km)');
title('MTZ 660 depth (raw)','fontsize',24)
text(-0.1,0.98,'(a)','Units','normalized','FontSize',32)

subplot(122)
m_proj('lambert','long', lonlim, 'lat', latlim); hold on;
set(gcf,'color','w')
h=m_pcolor(XI,YI,d660_ssa_smooth);
set(h,'edgecolor','none')
cm=colormap(jet);
colormap(flipud(cm));
caxis([630 690])

m_gshhs('i','line','color','k','linewidth',1)
m_gshhs('lb2','line','color','k')
m_grid('linewidth',2,'tickdir','out',...
    'xaxisloc','bottom','yaxisloc','left','fontsize',24);
hh=colorbar('h');
set(hh,'fontsize',24);
set(hh,'Position',[0.3065 0.1256-0.045 0.4220 0.0250])
xlabel(hh,'Depth (km)');
title('MTZ 660 depth (denoised)','fontsize',24)
text(-0.1,0.98,'(b)','Units','normalized','FontSize',32)
%% plot thickness
thi=d660_smooth-d410_smooth;
thi_ssa=d660_ssa_smooth-d410_ssa_smooth;
figure;
set(gcf,'Position',[100 100 1600 800],'Color','w')
subplot(121)
m_proj('lambert','long', lonlim, 'lat', latlim); hold on;
set(gcf,'color','w')
h=m_pcolor(XI,YI,thi);
set(h,'edgecolor','none')
cm=colormap(jet);
colormap(flipud(cm));
caxis([220 270])

m_gshhs('i','line','color','k','linewidth',1)
m_gshhs('lb2','line','color','k')
m_grid('linewidth',2,'tickdir','out',...
    'xaxisloc','bottom','yaxisloc','left','fontsize',24);
hh=colorbar('h');
set(hh,'fontsize',24);
set(hh,'Position',[0.3065 0.1256-0.045 0.4220 0.0250])
xlabel(hh,'Thickness (km)');
title('MTZ thickness (raw)','fontsize',24)
text(-0.1,0.98,'(a)','Units','normalized','FontSize',32)

subplot(122)
m_proj('lambert','long', lonlim, 'lat', latlim); hold on;
set(gcf,'color','w')
h=m_pcolor(XI,YI,thi_ssa);
set(h,'edgecolor','none')
cm=colormap(jet);
colormap(flipud(cm));
caxis([220 270])

m_gshhs('i','line','color','k','linewidth',1)
m_gshhs('lb2','line','color','k')
m_grid('linewidth',2,'tickdir','out',...
    'xaxisloc','bottom','yaxisloc','left','fontsize',24);
hh=colorbar('h');
set(hh,'fontsize',24);
set(hh,'Position',[0.3065 0.1256-0.045 0.4220 0.0250])
xlabel(hh,'Thickness (km)');
title('MTZ thickness (denoised)','fontsize',24)
text(-0.1,0.98,'(b)','Units','normalized','FontSize',32)