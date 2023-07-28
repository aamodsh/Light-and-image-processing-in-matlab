function lhs_estimate  = inverse_tie(I0,phase_in,ps)
% function lhs_estimate  = inverse_tie(I0,phase_in,ps)
%Recover dI/dz estimate (div(IgradPhi))
[gradPhix,gradPhiy] = gradient(phase_in);
I_gradPhix = (I0).*gradPhix/ps;
I_gradPhiy = (I0).*gradPhiy/ps;
%Take x gradient of IgradPhi_X and y gradient of the I_gradPhiy
[grad2x,dummy] = gradient(I_gradPhix);[dummy,grad2y] = gradient(I_gradPhiy);
lhs_estimate = (grad2x+grad2y)/ps;