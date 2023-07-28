%% --Define axes and return propagating orders for field and intensity

%Parameters
% NA=1.4;
% N=568;
% ps = 4.4;lambda = 193;


k0 = NA*2*pi/lambda;
center_pixel = N/2+1; %Only for on-axis
[Xtemp,Ytemp] = meshgrid((1:N)-center_pixel,(1:N)-center_pixel);
kx_mat = Xtemp/N/ps*2*pi;kx_vec = kx_mat(1,:);
ky_mat = Ytemp/N/ps*2*pi; ky_vec = ky_mat(:,1);
prop_orders = double((k0^2*ones(N)-kx_mat.^2-ky_mat.^2)>=0);
kz_mat = sqrt((k0^2*ones(N)-kx_mat.^2-ky_mat.^2).*prop_orders);
xaxis = (1:N)*ps/1000; %in micron

% figure(100);imagesc(kx_vec,ky_vec,kz_mat);title('kz mat');xlabel('spatial frequency (radian,spans 2*pi/ps)')


%Find propagating orders for field
[row,col] = find(kz_mat); %For non zero kz_mat values
rowrange = min(row):max(row);
colrange = min(col):max(col);

%Find propagating orders for intensity (twice the size of the pupil at the
%same sampling
rowrI = min(center_pixel + 2*(rowrange - center_pixel)):max(center_pixel + 2*(rowrange - center_pixel));
colrI = min(center_pixel + 2*(colrange - center_pixel)):max(center_pixel + 2*(colrange - center_pixel));

% rng = 40;
% rowrI = center_pixel-rng:center_pixel+rng;
% colrI = center_pixel-rng:center_pixel+rng;

%Find angle vector within the range of propagating orders
ky_angle = atan(ky_mat(rowrange,center_pixel)./kz_mat(rowrange,center_pixel))*180/pi;
kx_angle = atan(kx_mat(center_pixel,colrange)./kz_mat(center_pixel,colrange))*180/pi;

%Find tan of angle vector within the range of intensity orders
% ky_k = ky_mat(rowrI,center_pixel)./k0;
% kx_k = kx_mat(center_pixel,colrI)./k0;