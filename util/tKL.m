function delta = tKL(x,y)
%% total Kullback Leibler divergence (tKL) delta_f(x,y)
% f is the convex and differentiable generator function
% f(x) is the negative entropy;
% Input: 
%     x: x in n-simplex;
%     y: y in n-simplex;
% Output: 
%     delta: the tKL between x and y.
% Signature 
%   Author: Meizhu Liu
%   E-Mail: mliu@cise.ufl.edu
%   Date  : December 28 2010
% References:
%   Baba C. Vemuri, Meizhu Liu, Shun-Ichi Amari and Frank Nielsen, 
%   Total Bregman Divergence and its Applications to DTI Analysis, 
%   IEEE Transactions on Medical Imaging (TMI'10), 2010.
% 
%   Meizhu Liu, Baba C. Vemuri, Shun-Ichi Amari and Frank Nielsen, 
%   Total Bregman Divergence and its Applications to Shape Retrieval, 
%   IEEE Conference on Computer Vision and Pattern Recognition (CVPR’10),
%   2010.
n = size(x,1);
n2 = size(y,1);
if n2~=n
    fprintf('The dimension of x and y should be the same');
    exit;
end
idx = find(y==0);
z = y;
z(idx) = eps;
numerator = 0;
denominator = 0;
for i = 1:n
    numerator = numerator + x(i)*log(x(i)/z(i));
    denominator = denominator + y(i)*(1+log(z(i)))^2;
end
delta = numerator/sqrt(1+denominator);