function outmat = dcmask(N)
%N is even and 2D
outmat = ones(N,N);
outmat(N/2+1,N/2+1) = 0;
end