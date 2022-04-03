function Sj_masked = nk_SparseMaskObserv(Sj, p1)

m = size(Sj,2);
vec = zeros(1,m);
k = round(m*p1);
idx = randperm(m);
vec(idx(1:k)) = rand(k,1);
Sj_masked = Sj.*vec;