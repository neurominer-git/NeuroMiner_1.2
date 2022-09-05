function [vROI, tROI, yROI, yT] = ...
    nk_ROImapper(Pimg, Paddimg, ROIimg, ROIlist, MNIlist_index, MNIlist_labels, options)
% =========================================================================
% function [vROI, tROI, yROI, yT] = ...
%                      nk_aal_volume_roi(Pimg, Paddimg, ...
%                                        ROIimg, ROIlist, ...
%                                        MNIlist_index, MNIlist_labels, ...
%                                        options)
% =========================================================================
% The tool maps statistical image to neuroanatomical atlases and optionally 
% compares them to univariate measures with the chosen parcellation scheme.
%
% Inputs:
% Pimg :            paths to statistical images <char array>
% Paddimg :         (optional) paths to additional statistical 
%                   images to be co-analyzed with Pimg <char array>
% ROIimg :          path to parcellation image <char>
% ROIlist :         path to csv or excel file with parcellation info, which
%                   should contain the the id column indexing parcellations 
%                   in the parcellation image <char>
% MNIlist_index :   The column number in ROIlist which contains the IDs
%                   <double>
% MNIlist_labels :  The column number with the anatomical labels belonging 
%                   to the IDs <double>
% options :         a structure containing the following runtime options
%                   for the script
%   vxcheck :           check voxel dimensions across statistical images
%   thresh :            threshold applied to statistical images
%                       (scalar or 2-element vector)
%   typthresh :         thresholding operation (e.g. '>' or >=')
%   minext :            minimum percentage of voxels fulfilling minext 
%                       cutoff in ROI
%   saving :            save parcellation results to CSV/excel file 
%   name :              if saving=1, name of the CSV/excel file
%   ROIsel :            Pre-select specific ROIs from parcellation image as
%                       identified by their id
%   pattern_names :     Names of staistical images to be used in output
%   data :              optional, struct array of additional subject-level 
%                       data to be analyzed with univariate stats. Can be 
%                       either one matrix to which the different
%                       statistical maps are projected, or one matrix per
%                       statistical image file. This enables comparisons
%                       between univariate and multivariate statistical
%                       metrics which could be important to identify
%                       strengths of either approach over the other.
%       data.Y          The matrix with the subject-level data
%       data.label      Group index (at the moment only classification
%                       framework supported)
%       data.normalize  boolean flag to z-normalize entire matrix 
%       data.reflabel   if normalize = 1, index vector to cases used to
%                       normalize the entire matrix before computing
%                       univariate stats
%   compute_voxelstats : if data is provided, decided whether to analyze
%                       it.
%
% Outputs:
% vROI :            the parcellation struct array
% tROI :            the parecellation results table
% yROI :            the univariate analysis struct array
% yT :              the univariate results table
% *************************************************************************
% (c) Nikolaos Koutsouleris, 08/2022

if ~exist("ROIimg","var") || isempty(ROIimg)
    ROIimg = spm_select(1,'image','Select labeled atlas image');
end
if ~exist("ROIlist","var") || isempty(ROIlist)
    ROIlist = spm_select(1,'csv','Select parcellation label list');
end
if ~exist("MNIlist_index","var") || isempty(MNIlist_index)
    MNIlist_index = spm_input('Define the name of the table column containing the ROI indices', 0,'s');
end
if ~exist("MNIlist_labels","var") || isempty(MNIlist_labels)
    MNIlist_labels = spm_input('Define the name of the table column containing the ROI labels', 0,'s');
end

V_MNI = spm_vol(ROIimg);
V_dim = V_MNI.dim(3);

T = readtable(ROIlist, 'Delimiter', ';');
fprintf('\nPrinting first table row in %s', ROIlist); T(1,:)
ROIindex = T.(MNIlist_index);
ROInom_L = T.(MNIlist_labels);
for n=1:numel(ROIindex)
    ROI(n).ID = ROIindex(n); ROI(n).Nom_L = ROInom_L{n};
end

nROI = size(ROI,2);
if ~exist("Pimg","var") || isempty(Pimg)
    Pimg = spm_select(Inf,'image','Select cluster mask image(s)'); 
end
nP = size(Pimg,1);
if ~exist("Paddimg","var") || isempty(Paddimg)
    %addfl = spm_input('Do you want to read in an additional images for ROI analysis',0,'yes|no',[1 0]);
    addfl = 0;
    for i=1:nP
        if addfl
            Paddimg{i} = spm_select(Inf,'image',['Select additional images for cluster mask ' num2str(i)]);
        else
            Paddimg{i} = [];
        end
    end
end

patternnames = cellstr([repmat('P', nP, 1) num2str((1:nP)')]);
ROIsel = true(1,nROI);
compute_voxelstats = false;

if ~exist("options","var") || isempty(options)
    fprintf('\nDefine options for ROI parecellation analysis');
    thresh = spm_input('Cluster threshold (val or -/+)',0,'e',0);
    vxcheck = spm_input('Check image equality across signatures ?',0,'yes|no',[1,0],2);
    if ~isempty(thresh) && numel(thresh) == 2
        typthresh = spm_input('Type of threshold',0,'m','>|>=|<|<=', 1:4,1);
    else
        typthresh = spm_input('Type of threshold',0,'m','>|>=|<|<=|==', 1:5 ,1);
    end

    minext = spm_input('Minimum percentage of ROI',0,'e',5);
    saving = spm_input('Save results to disk?',0,'yes|no',[1 0],1);
    
    if saving
        px = fileparts(deblank(Pimg(i,:)));
        name = fullfile(px, 'ROI-Parcellation-Analysis.csv');
    end
    
else
    if isfield(options,'pattern_names') && numel(options.pattern_names) == nP
        patternnames = options.pattern_names;
    end
    vxcheck = options.vxcheck;
    thresh = options.thresh;
    typthresh = str2num(options.typthresh);
    minext = options.minext;
    saving = options.saving;
    name = options.name;
    if isfield(options,'ROIsel') && numel(options.ROIsel) == nROI
        ROIsel = options.ROIsel;
    end
    if isfield(options,'data') && isfield(options,'compute_voxelstats')  && options.compute_voxelstats
        compute_voxelstats = options.compute_voxelstats;
    end
end
 
for h=1:nP
    
    V{h} = spm_vol(deblank(Pimg(h,:)));
    Vx{h} = spm_vol(Paddimg{h});
    if ~isempty(Paddimg{h})
        nx(h) = size(Paddimg{h},1);
    else
        nx(h) = 0;
    end

    if nx(h)>0 && vxcheck == 1
        if any(any(diff(cat(1,V.dim,cat(1,Vx.dim)),1,1),1)&[1,1,1])
            error('images don''t all have same dimensions'), end
        if any(any(any(diff(cat(3,V.mat,cat(3,Vx.mat)),1,3),3)))
            error('images don''t all have same orientation & voxel size'), end
    end

    if nx(h) > 0
        for i=1:nx(h)
            [~,namx,ex] = fileparts(deblank(Paddimg{h}(i,:)));
            imgname{h}{i} = spm_input(['Image descriptor #' num2str(i) ' (' namx ex ')'],0,'s'); 
        end
        slice_valxcat{h} = cell(nx(h) , nROI);
    end

    M{h}   = V{h}.mat;

end

slice_valcat = cell(nROI, nP);
slice_ROIcat = zeros(nROI, 1);

for i = 1:nROI % Loop through all ROIs

    fprintf('\nROI %g/%g: %s',i,nROI, ROI(i).Nom_L );
    %if ~ROIsel(i), fprintf(' ...skip.'); end

    for j = 1:V_dim % Loop through all dimensions
        
        % Read ROI slice
        Mi = spm_matrix([0 0 j 0 0 0 1 1 1]);
        ROIimg = spm_slice_vol(V_MNI, Mi, V_MNI.dim(1:2),1);
        
        % compute ROI mask
        msk = round(ROIimg) == ROI(i).ID;

        % sum-up number of voxel in ROI across V_dim
        slice_ROIcat(i) = slice_ROIcat(i) + sum(msk(:));

        for h = 1:nP % Loop through all target images  

            % Read target image slice
            M1  = V_MNI.mat \ M{h} \ Mi;
            img = spm_slice_vol( V{h}, M1, V_MNI.dim(1:2), 1 );
            
            % Compute index of voxels to be analyzed as defined in thresh
            imgind = return_imgind(typthresh, thresh, img, msk);
            
            % Now limit the voxel space to the ROI and the voxels returned
            % by thresholding operation
            slice_val = img(imgind).*msk(imgind);
            
            % Concatenate voxels into vector
            slice_valcat{i, h} = [slice_valcat{i, h}; slice_val];

            if ~isempty(slice_val)
                for l=1:nx(h) % Loop through additional images
                    try
                        M1  = Vx{l}.mat\M{h}\Mi;
                        imgx = spm_slice_vol(Vx{l},M1,V_MNI.dim(1:2),[1 0]);
                        slice_valx = imgx(imgind).*msk(imgind);
                        slice_valxcat{h, l, i} = [slice_valxcat{h, l,i};slice_valx];
                    catch
                        fprintf('-');
                    end
                end
            end
        end
    end
    
    fprintf('\t%g => ',slice_ROIcat(i));

    % Compute stats
    for h = 1:nP
        if h>1, fprintf(', '); end
        fprintf('%g',sum(slice_valcat{i, h}(:)~=0))
        vROI(i).ref.id = ROI(i).ID;
        vROI(i).ref.name = ROI(i).Nom_L;
        vROI(i).ref.nvoxK(h) = length(slice_valcat{i, h});
        vROI(i).ref.volK(h) = length(slice_valcat{i, h}) * abs(det(V_MNI.mat(1:3,1:3)));
        vROI(i).ref.nvoxROI = slice_ROIcat(i);
        vROI(i).ref.volROI = slice_ROIcat(i) * abs(det(V_MNI.mat(1:3,1:3)));
        vROI(i).ref.percROI(h) = vROI(i).ref.nvoxK(h)*(100/slice_ROIcat(i));
        if ~isempty(slice_valcat{i,h})
            vl = slice_valcat{i,h};
            for l=1:nx(h), vlx(l) = slice_valxcat{h,l,i}; end
        else
            vl = nan; vlx = nan(1,nx(h));
        end
        vROI(i).ref.mean(h) = mean(vl);
        vROI(i).ref.std(h) = std(vl);
        vROI(i).ref.min(h) = min(vl);
        vROI(i).ref.max(h) = max(vl);
        for l=1:nx(h)
            vROI(i).img(h,l).mean = mean(vlx(l));
            vROI(i).img(h,l).std = std(vlx(l));
            vROI(i).img(h,l).min = min(vlx(l));
            vROI(i).img(h,l).max = max(vlx(l));
        end
    end
end

if saving && nargout > 1
    %fid = fopen(name,'w');
    fprintf(1,'\nWriting ROI data to file %s:',name)
    % Creating header line
    H = {'AnatomicalRegion', 'ROIvox', 'ROImm3'};
  
    for h=1:nP
        H = [H   sprintf('K_%s[vox]', patternnames{h}) ...
                 sprintf('K_%s[mm3]', patternnames{h}) ...
                 sprintf('K_ROI_%s[%%]', patternnames{h}) ...
                 sprintf('Mean_%s', patternnames{h}) ...
                 sprintf('SD_%s', patternnames{h}) ...
                 sprintf('Min_%s', patternnames{h}) ...
                 sprintf('Max_%s', patternnames{h}) ];
        for l=1:nx(h) 
            H = [H sprintf('Mean_%s_%g', patternnames{h},l) ...
                   sprintf('SD_%s_%g', patternnames{h},l) ...
                   sprintf('Min_%_s%g', patternnames{h},l) ...
                   sprintf('Max_%_s%g', patternnames{h},l) ];
        end
    end
    M = cell(nROI, numel(H)); delvec = [];
    for i=1:nROI      
        if (sum(vROI(i).ref.nvoxK==0)==nP || sum(vROI(i).ref.percROI<minext)==nP) || ~ROIsel(i)
            delvec = [delvec i];
            vROI(i).ref.sel = false(1,nP);
            continue; 
        end     
        M(i,1:3) = {ROI(i).Nom_L, vROI(i).ref.nvoxROI, vROI(i).ref.volROI};
        cnt=4;
        for h = 1:nP
           if vROI(i).ref.percROI(h)>=minext     
               M(i, cnt:cnt+6) = { vROI(i).ref.nvoxK(h), ...
                                vROI(i).ref.volK(h), ...
                                vROI(i).ref.percROI(h), ...
                                vROI(i).ref.mean(h), ...
                                vROI(i).ref.std(h), ...
                                vROI(i).ref.min(h), ...
                                vROI(i).ref.max(h) };
               cnt=cnt+7;
               for l=1:nx(h)
                    M(i, cnt : cnt+3) = { vROI(i).img(h,l).mean, ...
                                      vROI(i).img(h,l).std, ...
                                      vROI(i).img(h,l).min, ...
                                      vROI(i).img(h,l).max };
                    cnt=cnt+4;
               end
               fl=true;
           else
               M(i, cnt:cnt+6) = {nan}; cnt=cnt+7; 
               for l=1:nx(h), M(i, cnt : cnt+3) = {nan}; cnt=cnt+4; end
               fl=false;
           end
           vROI(i).ref.sel(h) = fl;
        end
    end
    M(delvec,:) = [];
    tROI = cell2table(M, 'VariableNames', H);
    writetable(tROI, fullfile(options.name, options.pattern_name), 'FileType', 'spreadsheet', 'Sheet', 'Vol2ROIparc-Main');
end

if nargout > 2 && isfield(options,'data') && ( numel(options.data) == nP || numel(options.data) == 1 )
    
    yROI=struct('data',[],'mean',[],'std',[],'voxROI',[], ...
        "F", [], "T", [], "Pval", [], "Pval_fdr", [], ...
        "minF", [], "maxF", [], "minT", [], "maxT", [], ...
        "minPval", [], "maxPval", [], "minPval_fdr", [], "maxPval_fdr", []);

    nD = numel(options.data); 
    yT = cell(1,nP);
    
    for i=1:nP

        if nD>1, idx= i; else, idx=1; end
        voxROI=[];
        
        % Read brainmask image from file
        if ischar(options.data(idx).brainmask)
            V_brainmask = spm_vol(options.data(idx).brainmask); 
        else
            V_brainmask = options.data(idx).brainmask; 
        end
    
        yROIvec = [];
        
        for j = 1:V_brainmask.dim(3) % Loop through all dimensions of the brainmask image
            
            % Read ROI slice and reslice it into brainmask image
            Mi = spm_matrix([0 0 j 0 0 0 1 1 1]);
            M1 = V_brainmask.mat\ V_MNI.mat \ Mi;

            % Read ROI slice
            ROIimg = spm_slice_vol(V_MNI, M1, V_brainmask.dim(1:2),1);
            ROIimg = round(ROIimg);
            
            % Read brainmask slice
            Mskimg = spm_slice_vol(V_brainmask, Mi, V_brainmask.dim(1:2),1);
            ind0 = Mskimg>0; icnt=0;

            % Read target image slice
            M1  = V_brainmask.mat \ V{idx}.mat \ Mi;
            Timg = spm_slice_vol( V{idx}, M1, V_brainmask.dim(1:2),1 );
            Timgind = return_imgind(typthresh, thresh, Timg);

            % Set un-selected ROIs in slice to 0
            for n = 1:nROI
                if ~any(vROI(n).ref.sel), 
                    ROIimg(ROIimg == vROI(n).ref.id) = 0; 
                else
                    icnt=icnt+1;
                end
            end
            
            % Set voxels in ROIs which are 0 in Timg to 0
            ROIimg(~Timgind) = 0;

            % Add pruned slice to vector (ind0 = msk > 0)
            yROIvec = [ yROIvec; ROIimg(ind0(:)) ];
        end

        % Now compute univariate stats (classification setting only at the moment)
        uL = unique(options.data(idx).label); nL = numel(uL);
        IN.X = double(nk_MakeDummyVariables(options.data(idx).label));
        Y=[]; 
        
        if compute_voxelstats
            voxMaxF = []; voxMaxT = []; voxMinF = []; voxMinT = []; voxMaxPval = [];
        end

        % Loop through ROIs and compute mean across voxels in ROI
        for n = 1:nROI
            if any(vROI(n).ref.sel) && ROIsel(n)
                if vROI(n).ref.sel(i)
                    % Find columns in matrix that belong to current ROI
                    Ix = yROIvec == vROI(n).ref.id;
                    % Compute mean across columns belonging to current ROI
                    nY = options.data(idx).Y(:, Ix);
                    Y = [ Y mean(nY,2) ];
                    if compute_voxelstats 
                        if sum(Ix)
                             nvoxROI = compute_ANOVA([],nY, IN);
                             voxROI{n} = nvoxROI;
                             voxMaxF = [ voxMaxF nvoxROI.maxF ];
                             voxMaxT = [ voxMaxT nvoxROI.maxT ];
                             voxMinF = [ voxMinF nvoxROI.minF ];
                             voxMinT = [ voxMinT nvoxROI.minT ];
                             voxMaxPval = [ voxMaxPval nvoxROI.maxPval ];
                        else
                            voxMaxF = [ voxMaxF nan ];
                            voxMaxT = [ voxMaxT nan ];
                            voxMinF = [ voxMinF nan ];
                            voxMinT = [ voxMinT nan ];
                            voxMaxPval = [ voxMaxPval nan ];
                        end
                    end
                else
                    % if ROI not selected add nan column to Y                    
                    Y = [ Y nan(size(options.data(idx).Y,1),1) ];
                    if compute_voxelstats 
                        voxMaxF = [ voxMaxF nan ];
                        voxMaxT = [ voxMaxT nan ];
                        voxMinF = [ voxMinF nan ];
                        voxMinT = [ voxMinT nan ];
                        voxMaxPval = [ voxMaxPval nan ];
                    end
                end
            end
        end

        if isfield(options.data(idx),'normalize') && options.data(idx).normalize
            if isfield(options.data(idx),'')
                Y = Y./nm_nanmean(Y(options.data(idx).reflabel==1,:));
            else
                Y = Y./nm_nanmean(Y);
            end
        end
        
        yROI(i).data = Y;
        yROI(i).mean = zeros(nL, size(Y,2));
        yROI(i).std = zeros(nL, size(Y,2));
        for q = 1 : nL
            yROI(i).mean(q,:) = mean(Y(options.data(idx).label== uL(q),:));
            yROI(i).std(q,:) = std(Y(options.data(idx).label== uL(q),:));
        end
        % Compute ANOVA over mean ROI values
        yROI(i) = compute_ANOVA(yROI(i), Y, IN);
        yy = [yROI(i).mean' yROI(i).std' yROI(i).F yROI(i).Pval yROI(i).Pval_fdr ];
        if compute_voxelstats
            yROI(i).voxMaxF = voxMaxF;
            yROI(i).voxMaxT = voxMaxT;
            yROI(i).voxMinF = voxMinF;
            yROI(i).voxMinT = voxMinT;
            yROI(i).voxMaxPval = voxMaxPval;
            yy = [yy yROI(i).voxMaxF' yROI(i).voxMinF' yROI(i).voxMaxT' yROI(i).voxMinT' yROI(i).voxMaxPval' ];
        end
        MeanH = cell(1,nL); StdH = cell(1,nL);
        for q = 1:nL
            MeanH{q} = sprintf('Mean_%s_Label%g ', patternnames{i}, uL(q)); 
            StdH{q} =  sprintf('Std_%s_Label%g ', patternnames{i}, uL(q)); 
        end
        H = ['Anatomical Region' MeanH StdH sprintf('F_%s', patternnames{i}) sprintf('P_%s', patternnames{i}) sprintf('Pfdr_%s', patternnames{i})];
        if compute_voxelstats
            H = [H  sprintf('voxMax-F_%s', patternnames{i}) ...
                    sprintf('voxMin-F_%s', patternnames{i}) ...
                    sprintf('voxMax-T_%s', patternnames{i}) ...
                    sprintf('voxMin-T_%s', patternnames{i}) ...
                         sprintf('voxMax-P_%s', patternnames{i})];
        end
        yT{i} = cell2table([tROI{:,1} num2cell(yy)],'VariableNames', H);
        writetable(yT{i}, options.name, 'FileType','spreadsheet','Sheet',sprintf('Vol2ROIparc_%s',patternnames{i}));
    end
else
    yROI=[];
    yT=[];
end

function imgind = return_imgind(typthresh, thresh, img, msk)

if ~exist('msk','var') || isempty(msk)
    msk = true(size(img));
end

if length(thresh) > 1
    switch typthresh
        case 1
            imgind = (img < thresh(1) | img > thresh(2)) & msk == 1;
        case 2
            imgind = (img <= thresh(1) | img >= thresh(2)) & msk == 1;
        case 3
            imgind = (img > thresh(1) | img < thresh(2)) & msk == 1;
        case 4
            imgind = (img >= thresh(1) | img <= thresh(2)) & msk == 1;
    end
else
    switch typthresh
        case 1
            imgind = img < thresh & msk == 1 ;
        case 2
            imgind = img <= thresh & msk == 1 ;
        case 3
            imgind = img > thresh & msk == 1 ;
        case 4
            imgind = img >= thresh & msk == 1 ;
        case 5
            imgind = img == thresh & msk == 1 ;
    end
end

function yROI = compute_ANOVA(yROI, Y, IN)

IN = nk_PerfANOVAObj(Y, IN);
yROI.F = IN.F; 
yROI.T = sqrt(IN.F); 
yROI.Pval = IN.p; 
[~,~,~,yROI.Pval_fdr] = fdr_bh(yROI.Pval);
yROI.Pval ( yROI.Pval==0 ) = realmin;
yROI.Pval_fdr(yROI.Pval_fdr==0 ) = realmin;
yROI.Pval = -log10(yROI.Pval); 
yROI.Pval_fdr = -log10(yROI.Pval_fdr);
yROI.minF = min(yROI.F);
yROI.maxF = max(yROI.F);
yROI.minT = sqrt(yROI.minF);
yROI.maxT = sqrt(yROI.maxF);
yROI.minPval = min(yROI.Pval);
yROI.maxPval = max(yROI.Pval);
yROI.minPval_fdr = min(yROI.Pval_fdr);
yROI.maxPval_fdr = max(yROI.Pval_fdr);
