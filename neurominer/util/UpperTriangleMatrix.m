function [Mxx, T, ax] = UpperTriangleMatrix(M, T, L, CL, FL)

if ~exist('T','var') || isempty(T), [~,T] = histcounts(M, 10) ;end
nT = numel(T); 

% Prepare P value matrix
%pval_mat(pval_mat==0)=realmin; pval_mat = -log10(pval_mat);
Mx = M;
Mx ( M < T(1) ) = NaN;
for j = 1:nT-1, Mx ( M>= T(j) & M<T(j+1) ) = j; end

Mxx = nT + 2.*eye(size(Mx,1)); 
Mxx( Mxx == nT )=NaN;
Ixx = itril(size(Mx,1),-1);
Jxx = itriu(size(Mx,1),1);
Mxx(Jxx) = Mx(Jxx);  
Mxx(Ixx) = nT+1;
Mxx(Mxx > nT+2) = nT-1;

% Display matrix
if nargin>2
    figure; ax=axes;
    imagesc(ax, Mxx, 'AlphaData', ~isnan(Mxx));
    ax.Position = [0.2 0.2 0.7 0.7];
    cmap1 = jet(nT-1);
    cmap2 = flipud(gray(3));
    cmap = [cmap1;cmap2];
    colormap(ax,cmap);
    ax.CLim =[ 1 nT+3 ];
    
    if exist('L','var') && ~isempty(L), 
        ax.XTickLabel = L; 
        ax.XTickLabelRotation=45; 
        ax.YTickLabel = L; 
    end
    CLabels = cellstr(num2str(T','%1.2f')); %CLabels{end} = sprintf('>%1.2f',T(end-1));
    cbar1 = colorbar(ax,'Ticks', 1:nT, 'TickLabels', CLabels); cbar1.Limits=[1 nT];
    if exist('CL','var') && ~isempty(CL), 
        cbar1.Label.String=CL; 
        cbar1.Label.FontWeight='bold';cbar1.Label.FontSize = 12;
    end
    
    if exist('FL','var') && ~isempty(FL), 
        ax.Title.String=FL; 
        ax.Title.FontWeight='bold';ax.Title.FontSize = 12;
    end
else
    ax=[];
end