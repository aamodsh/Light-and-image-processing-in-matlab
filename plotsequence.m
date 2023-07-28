function plotsequence(Imstack,pauselen,coloraxis)
%function plotsequence(Imstack,pauselen)
%Plots image sequence at gaps of pauselen

% figure;
for i = 1:size(Imstack,3)
    imagesc(Imstack(:,:,i));
    axis image;
    colorbar;
    title(num2str(i));
    if nargin>2
        caxis(coloraxis);
    end
    pause(pauselen);
end