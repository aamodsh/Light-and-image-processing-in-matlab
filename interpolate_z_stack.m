function Iout = interpolate_z_stack(Iin,interpfactor)
%function Iout = interpolate_z_stack(Iin,interpfactor)
ny = size(Iin,1);
nx = size(Iin,2);
nz = size(Iin,3);

%% Using 3D interpolation
X = 1:nx;
Y= (1:ny)';
Z = 1:nz;

Xq = 1:nx;
Yq = (1:ny)';
Zq = linspace(1,nz,nz*interpfactor);
Iout = interp3(X,Y,Z,Iin,Xq,Yq,Zq,'spline');


%% Using interpmat pointwise

% npts0 = size(Iin,3);
% polydeg = npts0;
% zin= 1:npts0;
% npts = interpfactor*npts0;
% Iout = zeros(size(Iin,1), size(Iin,2),npts);

% for i = 1:(nx)
%     for j = 1:ny    
%         Iout(i,j,:) = interp_mat(squeeze(Iin(i,j,:)),zin,npts,polydeg);  
%         i %print i
%         
%         j %print j
%     end
% end



end