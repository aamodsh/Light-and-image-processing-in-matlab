function mat_interp_fine = interpolate_x_y(mat_interp,factor)
sizex = size(mat_interp,1);sizey = size(mat_interp,2);
[x,y] = meshgrid(linspace(1,sizex,sizex),linspace(1,sizey,sizey));
[xq,yq] = meshgrid(linspace(1,sizex,factor*sizex),linspace(1,sizey,factor*sizey));

clear mat_interp_fine im ;
mat_interp_fine = zeros(size(mat_interp,1)*factor,size(mat_interp,2)*factor,size(mat_interp,3));
for j = 1:size(mat_interp,3)
    mat_interp_fine(:,:,j) = interp2(y',x',mat_interp(:,:,j),yq',xq');
end

end