% Demonstration script for 
% simultaneous seismic data denoising and reconstruction via
% damped multichannel singular spectrum analysis
%
%  Copyright (C) 2015 The University of Texas at Austin
%  Copyright (C) 2015 Yangkang Chen
%
%  This program is free software: you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published
%  by the Free Software Foundation, either version 3 of the License, or
%  any later version.
%
%  This program is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details: http://www.gnu.org/licenses/
%
%  Reference:   
%  [1] Huang, W., R. Wang, Y. Chen, H. Li, and S. Gan, 2016, Damped multichannel singular spectrum analysis for 3D random noise attenuation, Geophysics, 81, V261-V270.
%  [2] Chen, Y., D. Zhang, Z. Jin, X. Chen, S. Zu, W. Huang, and S. Gan, 2016, Simultaneous denoising and reconstruction of 5D seismic data via damped rank-reduction method, Geophysical Journal International, 206, 1695-1717.
%  [3] Chen, Y., W. Huang, D. Zhang, W. Chen, 2016, An open-source matlab code package for improved rank-reduction 3D seismic data denoising and reconstruction, Computers & Geosciences, 95, 59-66.
%				
%
clc;clear;close all;

%% generate 3D synthetic data
a1=zeros(300,20);
[n,m]=size(a1);
a3=a1;
a4=a1;

k=0;
a=0.1;
b=1;
for t=-0.055:0.002:0.055
    k=k+1;
    b1(k)=(1-2*(pi*30*t).^2).*exp(-(pi*30*t).^2);
    b2(k)=(1-2*(pi*40*t).^2).*exp(-(pi*40*t).^2);
    b3(k)=(1-2*(pi*40*t).^2).*exp(-(pi*40*t).^2);
    b4(k)=(1-2*(pi*30*t).^2).*exp(-(pi*30*t).^2);
end
for i=1:m
  t1(i)=round(140);
  t3(i)=round(-6*i+180);
  t4(i)=round(6*i+10);
  a1(t1(i):t1(i)+k-1,i)=b1; 
  a3(t3(i):t3(i)+k-1,i)=b1; 
  a4(t4(i):t4(i)+k-1,i)=b1;
end

temp=a1(1:300,:)+a3(1:300,:)+a4(1:300,:);
for j=1:20
    a4=zeros(300,20);
    for i=1:m
  t4(i)=round(6*i+10+3*j); 
  a4(t4(i):t4(i)+k-1,i)=b1;
  
  t1(i)=round(140-2*j);
  a1(t1(i):t1(i)+k-1,i)=b1;
    end
    shot(:,:,j)=a1(1:300,:)+a3(1:300,:)+a4(1:300,:);
end
plane3d=shot;
d=plane3d/max(max(max(plane3d)));

%% without noise
dn=d;

%% decimate
[nt,nx,ny]=size(d);
ratio=0.5;
mask=genmask(reshape(d,nt,nx*ny),ratio,'c',201415);
mask=reshape(mask,nt,nx,ny);
d0=dn.*mask;

%% reconstruct (without denoising)
flow=0;fhigh=125;dt=0.004;N=3;Niter=10;mode=0;verb=1;
d1=fxymssa_recon(d0,mask,flow,fhigh,dt,N,Niter,eps,verb,mode);

% 2D quick comparison (clean,noisy,observed,reconstructed using MSSA)
figure;imagesc([d(:,:,9),d0(:,:,9),d1(:,:,9)]);caxis([-0.5,0.5]);

%% simultaneous denoising and reconstruction
% adding noise
randn('state',201314);
var=0.2;
dn=d+var*randn(size(d));
d0=dn.*mask;

% using MSSA
flow=0;fhigh=250;dt=0.002;N=3;Niter=10;mode=1;verb=1;
a=(Niter-(1:Niter))/(Niter-1); %linearly decreasing
d1=fxymssa_recon(d0,mask,flow,fhigh,dt,N,Niter,eps,verb,mode,a);

% 2D quick comparison
figure;imagesc([d(:,:,9),d0(:,:,9),d1(:,:,9)]);caxis([-0.5,0.5]);

% using DMSSA
flow=0;fhigh=250;dt=0.002;N=3;Niter=10;mode=1;verb=1;K=2;
a=(Niter-(1:Niter))/(Niter-1); %linearly decreasing
d2=fxydmssa_recon(d0,mask,flow,fhigh,dt,N,K,Niter,eps,verb,mode,a);

% 2D quick comparison
figure;imagesc([d(:,:,9),d0(:,:,9),d2(:,:,9)]);caxis([-0.5,0.5]);

%% calculate Signal-to-noise Ratio (SNR)
snrcyk(d,d0,2) %observed data
snrcyk(d,d1,2) %MSSA method
snrcyk(d,d2,2) %DMSSA method

%SNR results (might be slightly different for different PC platform)
%d0: -5.9853
%d1: 1.4503
%d2: 4.1965

%% save to binary files for generating final figures in Madagascar
 fid1=fopen('syn_clean.bin','w'); 		%Figure 1a
 fid2=fopen('syn_noisy.bin','w');		%Figure 1b
 fid3=fopen('syn_obs.bin','w');         %Figure 1c
 fid4=fopen('syn_mssa.bin','w');		%Figure 1d
 fid5=fopen('syn_dmssa.bin','w');		%Figure 1e
 
 fwrite(fid1,reshape(d,300,20*20),'float');
 fwrite(fid2,reshape(dn,300,20*20),'float');
 fwrite(fid3,reshape(d0,300,20*20),'float');
 fwrite(fid4,reshape(d1,300,20*20),'float');
 fwrite(fid5,reshape(d2,300,20*20),'float');