% process stacks of images for fit to polynomial (at least 5 dz's)
% calls fitonepixelALL (least-squares fit to polynomial of I(z) for one
% pixel), TIE.
% Laura Waller, Jan 2010, MIT
tic
ordnum=1;        %pick orders to make for higher order TIE
load testInten.mat;
alpha=1;
%solver settings
zpaditer=2048;      %zero padding for FFT or iteration # for MGD
solvertype='FFT';   %'FFT' for fft and 'MGD' for multigrid
s=1;                %increase this to downsample for speed in testing
zstart=1;
zend=size(Ividmeas,3);
zstep=1;
zarray=z(zstart:zstep:zend)';

%_________Higher order TIE - fit to polynomials___________________
dIdz=zeros(size(I,1),size(I,2),20);
RMSE=zeros(size(I,1),size(I,2),20);
Ampl=zeros(size(I,1),size(I,2),20);
pxx=1;
for pxx=pxx:size(I,1)
    for pxy=1:size(I,2)
        [dI,A,RM]=fitpixel1to20(squeeze(I(pxx,pxy,:)),zarray);
        pxx
        dIdz(pxx,pxy,:)=dI;
        Ampl(pxx,pxy,:)=A;
        RMSE(pxx,pxy,:)=RM;
    end
    clc;toc
    pxx
    %save parwaydone.mat
end
%compute phase for all orders
phi=zeros(size(dIdz));
del2psi=(-2*pi/lambda)*dIdz(:,:,ordnum);
del2psi=del2psi-mean(mean(del2psi));
%solve poisson
phi(:,:,ordnum)=-poissonFFT(del2psi,zpaditer)*ps^2*alpha;
phiHOTIE=-phi(:,:,ordnum)-mean2(-phi(:,:,ordnum));
ampHOTIE=Ampl(:,:,ordnum);

%compute errors
errorHOTIE=mean2(abs((phiHOTIE-phi0).^2))

%__________________DISPLAY-PHASE________________
imagesc(-phiHOTIE);axis image;axis off;colorbar 
title(sprintf('%dth order',ordnum));
colormap gray