function dPhi_dz =  inverse_tpe(Efield,ps,lam, z)
%First term
[gx, gy] = gradient(angle(Efield));
gradPhi_sq = gx.^2 + gy.^2;

%Second term
Eamp = abs(Efield);
[gx, gy] = gradient(abs(Efield));
[g2x,dummy] = gradient(gx);
[dummy,g2y ] = gradient(gy);
grad2_Eamp = g2x+g2y;

dPhi_dz = grad2_Eamp./Eamp.^2 + gradPhi_sq/2;

end
