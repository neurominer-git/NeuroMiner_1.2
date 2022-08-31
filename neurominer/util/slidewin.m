function [ tr, x ] = slidewin(X, E, winsize)

nX = numel(X);
if winsize > nX, error('Window size must be smaller than X'); end
[~,sI] = sort(X, 'ascend');
X = X(sI);
E = E(sI);

for i = 1 : nX - winsize
   
    tr(i) = sum(E(i:winsize+i)==1)/winsize;
    x(i) = mean(X(i:winsize+i));
end