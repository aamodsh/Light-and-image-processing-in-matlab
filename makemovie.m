clear im_stack;
%Aberrations
% im_stack = Hmat;
% % im_stack = squeeze(F2Ioutmat(:,:,1,:));
% zstack = -kx_shift_vec/k/NA;

%Ptch dataset
    for i = 1:120;
    im_stack(:,:,i) = log(abs(diffractogram(img_tmp{i}(1:300,1:300))));
    end
    zstack = 1:120;%-kx_shift_vec/k/NA;

%%
fig = figure;
% set(fig, 'Position', [0 0 100 100])

for i = 1:(size(im_stack,3)-3)
figure(fig);
% subplot(2,1,1);
imagesc(abs(im_stack(100:900,100:900,i)).^2-abs(im_stack(1,1,i)).^2);
% title([ num2str(zstack(i)) ' NA']);
% caxis([0 1.2]);
axis off;
set(gcf,'color','w');

% colormap gray;
% subplot(2,1,2)
% plot(abs(im_stack(251,:,i)),'Linewidth',2);
% colorbar;
% caxis([-1.2*pi 1.2*pi]);

m(:,:,i) = getframe(gcf);

waitforbuttonpress;
end

saveMovie(m,'ptc2.gif',10,'k');