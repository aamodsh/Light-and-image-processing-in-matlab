function outmat = fft_correlate(E1,E2)
fftmacros;
Ef = F2(E1).*(conj(F2(E2)));
outmat = iF2(Ef);
end