% authors: Ioannis Psorakis psorakis@gmail.com, Theodoros Damoulas theo@dcs.gla.ac.uk
%    Copyright NCR, 2009 
%
%    Use of this code is subject to the Terms and Conditions of the Research License Agreement,
%    agreed to when this code was downloaded, and a copy of which is available at 
%    http://www.dcs.gla.ac.uk/inference/pMKL/Terms_and_Conditions.html.

function K = kernel_func_RVM(X,Y,type,param)
%disp('---Kernel function---')

if strcmp(type,'linear') || strcmp(type,'lin') || strcmp(type,' -t 0')
   
    K=X*Y'; 
    
elseif strcmp(type,'gaussian') || strcmp(type,'rbf') || strcmp(type,'gauss')
    N = size(X,1);
    M = size(Y,1);
    D = size(X,2);
    
    if ~exist('param','var')
        param = 1/D;
    end
    
    if size(param)==ones(1,2)
        param = repmat(param,1,D);
    end
    
    T=diag(-param);% The diagonal matrix DxD with D parameters
    K=exp(repmat(diag(X*T*X'),1,M)+repmat(diag(Y*T*Y')',N,1)-2*X*T*Y');
  
elseif strcmp(type,'polynomial') || strcmp(type,'polyN') || strcmp(type,'poly')
  
    if ~exist('param','var')
        param = 2;
    end
    
    prod = X*Y';
    K = (prod + ones(size(prod))).^param;
end
