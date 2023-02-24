function CM = cv_graph_constructionKLS(A, method, brainmask, atlas)
% KL divergence method following Kong et al. 2014
%
% INPUT: nii images (already smoothed, it's another preprocessing step of
% course which will always be done first) 
%
% Steps: 
%   1. nk_WriteVol --> save image to disk (.nii); compute network; delete
%       this temp.nii; loop through images 
%   within Python script
%       1. extract voxel values of each region
%       2. compute probability density functions of the gray matter volumes 
%           for each region of the provided atlas 
%       3. compute symmetric KL divergence between each pair of probability
%           density functions
%       4. bring data in the right format (flat connectivity matrix, upper
%           triangle); dimensions: n people x (n/2-n) edges; either safe
%           this vector to csv-file or add to existing one 
% 
% OUTPUT: matrix format, dimensions: n people x (n/2-n) edges 
% 
% --> adapted Python script from my thesis
% potential issue: I also used R scripts for the computation ...

% 1. extract voxel values per region 
S.Vm                         = spm_vol(brainmask);
[S.dims, S.indvol, ~, S.vox] = nk_ReadMaskIndVol(S.Vm, []);

image_for_size = char(brainmask);

if ~exist('atlas','var') || isempty(atlas)
    maskVec = nk_Vol2Vec(brainmask);
else 
    ROI_matrix = resize_img_useTemp_imcalc(char(atlas),image_for_size); %from JuSpace Toolbox
    maskVec = reshape(ROI_matrix,size(ROI_matrix,1)*size(ROI_matrix,2)*size(ROI_matrix,3),1)';
end 

maskVec = int64(maskVec(maskVec > 0));
n_rois = length(unique(maskVec));
rois = unique(maskVec);

voxelValues_struct = struct;
pdf_struct = struct;
sp = 512; % can be adjusted! perhaps in config
for i = 1:size(A,1)
    pp_values = A(i,:)
    ppname = sprintf('pp_%d',i);
    for j = 1:n_rois
        roi = rois(j); 
        voxel_values = pp_values(maskVec == roi);
        roiname = sprintf('roi_%d',j);
        voxelValues_struct.(roiname) = voxel_values
        [bw, density, grid, cdf_values] = kde(voxel_values, n = sp);
        pdf_struct.(ppname).(roiname) = cdf_values;
    end
end


n_edges = (n_rois*n_rois-n_rois)/2; % upper triangle of the connectivity matrix (without diagonal)
kls_df = zeros(size(A,1),n_edges);

for i = 1:size(A,1)
    con_mat = diag(ones(n_rois,1)); % values on diagonal are one (since the pdfs are equal)
    ppname = sprintf('pp_%d',i);
    pp_rois = pdf_struct.(ppname);
    for j = 1:n_rois
        for k = j+1:n_rois % upper triangle of square matrix excluding diagonal = all values where column index > row index
            roiname1 = sprintf('roi_%d',j);
            roiname2 = sprintf('roi_%d',k);
            pdf_j = pp_rois.(roiname1);
            pdf_k = pp_rois.(roiname2);
            kls = symmetric_kldiv(pdf_j, pdf_k);
            con_mat(j,k) = kls;
            con_mat(k,j) = kls; % because it's symmetric (actually arbitrary since lower half of matrix will be removed
        end
    end
    % remove lower triangle + diagonal and flatten matrix 
    % upper_tri = triu(con_mat,1); % without diagonal
    con_mat_t = con_mat.';
    rem_idx = (1:size(con_mat_t,1)).' > (1:size(con_mat_t,2));
    edges_vec = con_mat_t(rem_idx);
    kls_df(i,:) = edges_vec;
end
CM = kls_df;
end

function sym_kldiv = symmetric_kldiv(p,q)
    % computes a symmetric version of the KL divergence (which is not
    % symmetric) 
% 
%     for i = 1:length(p)
%         sym_kldiv = sym_kldiv + (p(i) * log2(p(i)/q(i)) + q(i) * log2(q(i)/p(i))); 
%     end

    P = reshape(p,[numel(p),1]);
    Q = reshape(q,[numel(p),1]); % p und q have to have the same n of bins (length)
    sym_kldiv = sum( P .* log2(P./Q) + Q .* log2(Q./P),'omitnan');
end