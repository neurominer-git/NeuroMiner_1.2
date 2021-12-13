function [ res ] = ECE( labels, predictions )
% Expected calibration error
    % Check 
    labels(labels==2)=0;
    labels(labels==-1)=0;
    ordered = sortrows([predictions,labels]);
    N = size(ordered,1);
    rest = mod(N,10);
    S = 0;
    B = min(N,10);
    for i = 1 : B
        if (i <= rest)
            group = ordered((i-1) * ceil(N / 10) + 1 : i * ceil(N / 10),:);
        else
            group = ordered(rest + (i-1) * fix(N / 10) + 1 : rest + i * fix(N / 10),:);
        end
        n = size(group,1);
        observed = mean(group(:,2));
        expected = mean(group(:,1));
        S(i,1) = abs(expected-observed);
        W(i,1) = n/N;
    end
    res = S'*W;

end
