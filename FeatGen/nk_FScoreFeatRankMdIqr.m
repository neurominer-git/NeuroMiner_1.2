function D = nk_FScoreFeatRankMdIqr(Y, L, N)

% Y     :       Data
% L     :       Target Labelss

ix = unique(L);

if numel(ix) == 2

    indP = L==ix(1); indM = L==ix(2);
    YP = Y(indP,:); YM = Y(indM,:);
    
    % Mean of Positive / Negative Labels
    mP = nm_nanmedian(YP); mM = nm_nanmedian(YM);
    
    % Standard Deviations of Positive / Negative Labels
    sP = iqr(YP); sM = iqr(YM);

    if exist('N','var') && ~isempty(N)
        
        ND = zeros(size(Y,2),size(N,2));
        % Compute nuisance-insensitive Fscore(s)
        for i = 1:size(N,2)
            nx = unique(N(:,i));
            indNP = N(:,i)==nx(1); indNM = N(:,i)==nx(2);
            mNP = nm_nanmedian(Y(indNP,:)); mNM = nm_nanmedian(Y(indNM,:));
            sNP = iqr(Y(indNP,:)); sNM = iqr(Y(indNM,:));
            tND = (mNP-mNM).^2 ./ (sNP + sNM) ; 
            ND(:,i) = tND';   
        end
    end

    % Compute Fscore and optionally weight down nuisance voxels
    D = (mP-mM).^2 ./ (sP + sM); D = D'; 
    
    if exist('ND','var')
        D = D ./ nm_nansum(ND,2);
    end
else
    error('F-Score feature ranking works only for binary problems! Change pre-processing group mode to binary.')
end