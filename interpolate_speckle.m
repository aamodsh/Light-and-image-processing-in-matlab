%% read image stack
%creates the 3d stack

%%-read strong speckle---%
% load('matlab_workspace_diffuser.mat', 'mat_interp_fine_zflip');
for i = 0:10
    mat_w(:,:,i+1) = double(imread(['/Users/aamod/Dropbox/Speckle_data_11_20/im' num2str(i) 'in.bmp']));
end

mat = register_im_stack(mat_w);

%%-read weak speckle ---%
for j = 1:9
    im(:,:,j) = imread(['/Users/aamod/Dropbox/Cusgadmin Documents symlink/[Research]/Waller/2017 speckle/SHARP data for Stuart/SHARP through-focus/IMO113715-140110-0015-000' num2str(j) '.png']);
end

im(:,:,10) = imread(['/Users/aamod/Dropbox/Cusgadmin Documents symlink/[Research]/Waller/2017 speckle/SHARP data for Stuart/SHARP through-focus/IMO113715-140110-0015-0010.png']);
mat = register_im_stack(im);

%% Interpolate in z and plot contours
mat_cont = mat_w;
% mat_cont = im;
%--z interp factor
factor = 7;
%--for weak
% yrange = 650:950;
% xrange = yrange+100;
%--for strong
yrange = 350:650;
xrange = yrange+100;

mat_interp = interpolate_z_stack(double(mat_cont(yrange,xrange,:)),factor);
zmax = size(mat_interp,3);
imagesc(mat_interp(:,:,1));


%Flip and tile in the z coordinate
mat_interp_flip = cat(3,mat_interp(:,:,2:(zmax-1)),flip(mat_interp,3));
mat_interp_flip = real(sqrt(mat_interp_flip));
ps = 5;

while 1
    for j = 1:size(mat_interp_flip,3)
        %--contour through-focus
        %         contour(sqrt(mat_interp_flip(:,:,j)),20, 'Linewidth',1.5);
        %--fft through-focus
        imagesc(rmvDC((abs(fftshift(fft2(sqrt(mat_interp_flip(:,:,j))))))));
        %--mesh plot
        %         mesh(ps*(1:size(mat_interp_flip,2)),ps*(1:size(mat_interp_flip,1)),(mat_interp_flip(:,:,j)<10));
        %         axis([0 ps*size(mat_interp_flip,1) 0 ps*size(mat_interp,2) 0 4000]);
        %           view([-40 70]); %[v,a] = view (default [-40 75]);
        %          caxis([0 3000]);
        
        
        title(['speckle image# ' num2str(j)]);
        xlabel('um');
        ylabel('um');
        zoom(2);
        waitforbuttonpress;
        %         caxis([5 6])
        %         f(:,:,j) = getframe;
        %         pause(0.01);
    end
end

%important to consider the field at the boundaries, derivative
%singularities as the attractors, going to zero or to infinity

%% save movie from frames
saveMovie(f,'sharp_speckle_ptc_weak.gif',10,'k');

%% plot through 1D
xpt = 45;
ypt = 58;

% h1 = figure(1);
% contour(mat_interp(:,:,40));
% makedatatip(ypt,xpt);

figure(3);
plot(squeeze(sqrt(mat_interp(ypt,xpt,:))));

soundsc(squeeze(sqrt(mat_interp(ypt,xpt,:))),80);
%% Interpolate through x-y-z
yrange = 200:600;
xrange = 200:600;

% yrange = 1:size(mat,1);
% xrange = yrange;
% % %z interpolation
factor = 2;
mat_interp = (interpolate_z_stack(double(mat(yrange,xrange,:)),factor));
[size(mat)]
%xy interpolations
factor = 3;
sizex = size(mat_interp,1);sizey = size(mat_interp,2);
[x,y] = meshgrid(linspace(1,sizex,sizex),linspace(1,sizey,sizey));
[xq,yq] = meshgrid(linspace(1,sizex,factor*sizex),linspace(1,sizey,factor*sizey));

clear mat_interp_fine im ;
mat_interp_fine = zeros(size(mat_interp,1)*factor,size(mat_interp,2)*factor,size(mat_interp,3));
for j = 1:size(mat_interp,3)
    mat_interp_fine(:,:,j) = interp2(y',x',mat_interp(:,:,j),yq',xq');
end
size(mat_interp_fine)


%% Plot mesh plots through z
%z flip and tile for continuity
zmax = size(mat_interp_fine,3);
% mat_interp_fine_zflip = zeros(size(mat_interp_fine,1),size(mat_interp_fine,2),2*size(interp_mat,3));
mat_interp_fine_zflip = cat(3,mat_interp_fine(:,:,2:(zmax-1)),flip(mat_interp_fine,3));
scanstep = 1400;
% figure(4);
xstart = 330;
xrange = xstart:(xstart+21);
yrange = xstart:(xstart+21);
ps = 5;
while 1
    for j = 1:size(mat_interp_fine_zflip,3)
        
        % Contour plot
%         contour(mat_interp_fine_zflip(yrange,xrange,j),40);%zoom(5);
        %       --MESH PLOT
        mesh(ps*(1:size(xrange(:))),ps*(1:size(yrange(:))),real(sqrt((mat_interp_fine_zflip(yrange,xrange,j)))));
        mesh(real(sqrt((mat_interp_fine_zflip(yrange,xrange,j)))));
        caxis([2 3]);
%       axis([0 ps*size(xrange,2) 0 ps*size(yrange,2) 0 60]);
        view([40 120]);
        title(['speckle image# ' num2str(j)]);
        xlabel('um');
        ylabel('um');
        %       axis off;
        %       axis image;
        f(:,:,j) = getframe;
        pause(0.01);
        %         waitforbuttonpress
    end
end
%%
saveMovie(f,'ptc_single_speckle.gif',15,'k');

%% plot  zero contours
mat_interp_flip = real(sqrt(mat_interp_fine(:,:,:)));
while 1
    for j = 1:size(mat_interp_flip,3)
        contour(squeeze(mat_interp_flip(:,:,j)),'Linewidth',2);
        %         caxis([0 3000]);
        title(['speckle image# ' num2str(j)]);
        xlabel('um');
        ylabel('um');
        %         zoom(2);
        %         waitforbuttonpress;
        % f(:,:,i) = getframe;
        pause(0.1);
    end
end

%% x-z through y
mat_interp_flip = ((mat_interp_fine(:,:,:)));
while 1
    for j = 1:size(mat_interp_flip,3)
        contour(squeeze(mat_interp_flip(:,j,:)), 'Linewidth',2);
        %         caxis([0 3000]);
        title(['speckle image# ' num2str(j)]);
        xlabel('um');
        ylabel('um');
        zoom(1);
        waitforbuttonpress;
        % f(:,:,i) = getframe;
        %pause(0.1);
    end
end

%%
% contour(squeeze(mat_interp_flip(:,160,:)),30, 'Linewidth',2);
imagesc(squeeze(mat_interp_flip(:,:,:)));
% zoom(-11);

















