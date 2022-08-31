% authors: Ioannis Psorakis psorakis@gmail.com, Theodoros Damoulas theo@dcs.gla.ac.uk
%    Copyright NCR, 2009 
%
%    Use of this code is subject to the Terms and Conditions of the Research License Agreement,
%    agreed to when this code was downloaded, and a copy of which is available at 
%    http://www.dcs.gla.ac.uk/inference/pMKL/Terms_and_Conditions.html.

function K = build_standarized_kernel(X,Y,type,Param)

N1 = size(X,1);
N2 = size(Y,1);

K_XY = SB1_KernelFunction(X, Y, type, Param);
K_XX = SB1_KernelFunction(X, X, type, Param);
K_YY = SB1_KernelFunction(Y, Y, type, Param);

% old kernel functions 
% K_XY = kernel_func_RVM(X,Y,type,param);
% K_XX = kernel_func_RVM(X,X,type,param);
% K_YY = kernel_func_RVM(Y,Y,type,param);

KXs = repmat(diag(K_XX),1,N2);
KYs = repmat(diag(K_YY)',N1,1);
KXYs = sqrt(KXs .* KYs);

K = K_XY ./ KXYs;