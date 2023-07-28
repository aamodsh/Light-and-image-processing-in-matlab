%% Function to calculate Fourier transform
%function P = diffractogram(im,yrange,xrange)

function P = diffractogram(im)

    % Take FT of image 
    % Remove DC (equivalent to enforcing mean == 0)
    [ft, ~] = rmvDC(fft2c(im));
    % Compute 2D periodogram
      P = (ft);
%     P = Pfull(yrange,xrange);
    
end
    
    