function [outmat, best_focus]  = interp_mat(inmat,z_in,npts,polydeg)
%function [outmat, best_focus]  = interp_mat(inmat,z_in,npts,polydeg)
%Interpolate inmat rowwise to npts with polynomial of degree polydeg

zinterp = linspace(min(z_in),max(z_in),npts);

for  i =  1:size(inmat,2)
    %--Generate interpolated matrix
    cd_temp = inmat(:,i);
    outmat(i,:) = polyval(polyfit(z_in,cd_temp',polydeg),zinterp);
    
    %--Find peaks of curves
    [a,b] = max(outmat(i,:));
%   [a,b] = min(y(i,:)); %could be min
    best_focus(i) = zinterp(b);

end
end