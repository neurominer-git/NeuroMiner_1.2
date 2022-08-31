function [dY, mY, sY] = discretize(Y, astart, astep, aend, mY, sY)
% =========================================================================
% FORMAT [dY, IN] = nk_PerfDiscretizeObj(Y, IN)
% =========================================================================
%
% Discretizes Y columnwise to mean +/- alpha*std, where alpha is a
% vector of discrete values. By default alpha is set to 0:0.5:4.
% _________________________________
% (c) Nikolaos Koutsouleris 01/2022
    
if ~exist("mY","var")
    mY = mean(Y); sY = std(Y);
end

[m,n] = size(Y);
dY=zeros(m,n);

if nargin < 2,
    astart = 0; astep  = 0.5; aend = 4;
end

alphas = astart:astep:aend;

tmY = repmat(mY,m,1); tsY = repmat(sY,m,1);

for j=1:numel(alphas)-1
    dY( Y > ( tmY + alphas(j) .* tsY) ) = alphas(j+1);
    dY( Y < ( tmY - alphas(j) .* tsY) ) = -alphas(j+1);
end

dY( Y > tmY + alphas(end) .* tsY ) = aend;
dY( Y < tmY - alphas(end) .* tsY ) = -(aend);
