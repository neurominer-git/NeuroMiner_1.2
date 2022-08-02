% SB1_KERNELFUNCTION    Compute kernel functions for the RVM model
%
%       K = SB1_KERNELFUNCTION(X1,X2,KERNEL,LENGTH)
%
% OUTPUT ARGUMENTS:
%
%       K       N1 x N2 kernel matrix.
% 
% INPUT ARGUMENTS:
%
%       X1      N1 x d data matrix
%       X2      N2 x d data matrix
%       KERNEL  Kernel type: currently one of
%               'rbf'         Gaussian
%               'laplace'       Laplacian
%               'poly'         Polynomial of order 'N'
%               'hpoly'        Homogeneous Polynomial of order 'N'
%               'spline'        Linear spline [Vapnik et al]
%               'cauchy'        Cauchy (heavy tailed) in distance
%               'cubic'         Cube of distance
%               'r'             Distance
%               'tps'           'Thin-plate' spline
%               'bubble'        Neighbourhood indicator
%       LENGTH  Input length scale
% 
%
% Copyright 2009 :: Michael E. Tipping
%
% This file is part of the SPARSEBAYES baseline implementation (V1.10)
%
% Contact the author: m a i l [at] m i k e t i p p i n g . c o m
%
function K = SB1_KernelFunction(X1,X2,kernel_,lengthScale)

[N1 ]		= size(X1);
[N2 ]		= size(X2);

if isempty(lengthScale), lengthScale=1;end

switch lower(kernel_)
 
 case {'lin','linear',' -t 0', 'lin_kernel'}
    
     K = X2*X1';
     
 case 'expo'
     
     eta	= 1/lengthScale^2;
     K = exp(-eta*sqrt(distSqrd(X1,X2)));
     
 case {'rbf', 'gauss', 'gaussian'}
    
     eta	= 1/lengthScale^2;
     K	= exp(-eta*distSqrd(X1,X2));
  
 case 'tps'
     
     eta	= 1/lengthScale^2;
     r2	= eta*distSqrd(X1,X2);
     K	= 0.5 * r2.*log(r2+(r2==0));
  
 case 'cauchy'
    
     eta	= 1/lengthScale^2;
     r2	= eta*distSqrd(X1,X2);
     K	= 1./(1+r2);
  
 case 'cubic'
    
     eta	= 1/lengthScale^2;
     r2	= eta*distSqrd(X1,X2);
     K	= r2.*sqrt(r2);
  
 case 'r'
    
     eta	= 1/lengthScale^2;
     K	= sqrt(eta)*sqrt(distSqrd(X1,X2));
  
 case 'bubble'
    
     eta	= 1/lengthScale^2;
     K	= eta*distSqrd(X1,X2);
     K	= K<1;
  
 case 'laplace'

     eta	= 1/lengthScale^2;
     K	= exp(-sqrt(eta*distSqrd(X1,X2)));
  
 case {'poly','polynomial','polyn'}
     
     eta	= 1/lengthScale(1)^2;
     K	= (eta*X1*X2' + lengthScale(2)).^ lengthScale(3);

 case {'hpoly', 'hpolyn'}
     
     eta	= 1/lengthScale(1)^2;
     K	= (eta*X1*X2').^ lengthScale(2);

 case 'spline'

    K	= 1;
    X1	= X1/lengthScale;
    X2	= X2/lengthScale;
    for i=1:d
        XX		= X1(:,i)*X2(:,i)';
        Xx1		= X1(:,i)*ones(1,N2);
        Xx2		= ones(N1,1)*X2(:,i)';
        minXX	= min(Xx1,Xx2);
        K	= K .* (1 + XX + XX.*minXX-(Xx1+Xx2)/2.*(minXX.^2) + (minXX.^3)/3);
    end
  
 otherwise
    error('Unrecognised kernel function type: %s', kernel_)
end 

%%
%% Support function: squared distance
%%
function D2 = distSqrd(X,Y)
nx	= size(X,1);
ny	= size(Y,1);
D2	= sum(X.^2,2)*ones(1,ny) + ones(nx,1)*sum(Y.^2,2)' - 2*(X*Y');
