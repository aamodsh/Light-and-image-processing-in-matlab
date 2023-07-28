function outmat = fft_convolve(E1,E2)
fftmacros;
Ef = F2(E1).*((F2(E2)));
outmat = iF2(Ef);
end