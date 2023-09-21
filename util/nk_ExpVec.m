function y = nk_ExpVec(a,b,n)
gf = 10^(1/(n - 1));
nx=1:1:n;
y = rescale(1*gf.^(nx-1),a,b);
