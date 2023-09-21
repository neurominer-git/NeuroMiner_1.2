function [C, cAUC] = nk_CreatPredFeat(AUC_upper, AUC_lower, L, verbose)

cAUC = 0;
c1_x0_rng = linspace(0, .25, 5); c1n = numel(c1_x0_rng);
c2_x0_rng = linspace(-0.25, 0, 5); c2n = numel(c2_x0_rng);
w_st = 0.5; w_end = 1.5;
c_w_rng = linspace(w_st, w_end, 5); cwn = numel(c_w_rng);
c1_caseN = sum(L==1);
c2_caseN = sum(L==-1);
max_iter = 1000; cnt=1;

while cAUC > AUC_upper || cAUC < AUC_lower
    
    c1_x0_I = randi(c1n);
    c2_x0_I = randi(c2n);
    c1_w_I = randi(cwn);
    c2_w_I = randi(cwn);
    
    c1 = GaussDist(c1_x0_rng(c1_x0_I), c_w_rng(c1_w_I), c1_caseN);
    c2 = GaussDist(c2_x0_rng(c2_x0_I), c_w_rng(c2_w_I), c2_caseN);
    C = [c1; c2];
    
    cAUC = AUC(L,C); 
    cnt=cnt+1;
    if cnt == max_iter
        c1_x0_rng = c1_x0_rng + 0.05;
        c2_x0_rng = c2_x0_rng - 0.05;
        w_end = w_end - 0.1; 
        if w_end == w_st, w_st = w_st-0.05; end
        c_w_rng = linspace(w_st, w_end, 10);
        cnt=1;
    end
end

if verbose
    figure;histogram(c1,10,'Normalization','probability'); hold on; histogram(c2,10,'Normalization','probability')
end