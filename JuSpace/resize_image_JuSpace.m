function image_resized = resize_image_JuSpace(image_to_resize,image_for_size)

Vi = {image_for_size;image_to_resize};
Vo = image_for_size;

f = 'i2';

flags.dmtx = 0;
flags.mask = 0;
flags.interp = 0;
flags.dtype = 4;

if ~isstruct(Vi), Vi = spm_vol(char(Vi)); end

%% part spm_imcalc
if isfield(flags,'dmtx'),   dmtx   = flags.dmtx;   else dmtx   = []; end
if isfield(flags,'mask'),   mask   = flags.mask;   else mask   = []; end
if isfield(flags,'interp'), interp = flags.interp; else interp = []; end
if isfield(flags,'dtype'),  dtype  = flags.dtype;  else dtype  = []; end

% if ischar(Vo)
%     [p, n, e] = spm_fileparts(Vo);
%     Vo = struct('fname',   fullfile(p, [n, e]),...
%                 'dim',     Vi(1).dim(1:3),...
%                 'dt',      [dtype spm_platform('bigend')],...
%                 'pinfo',   [Inf Inf Inf]',...
%                 'mat',     Vi(1).mat,...
%                 'n',       1,...
%                 'descrip', 'spm - algebra');
% end

Vo = spm_vol(Vo);

n = numel(Vi);
Y = zeros(Vo.dim(1:3));

for p = 1:Vo.dim(3)
    B = spm_matrix([0 0 -p 0 0 0 1 1 1]);

    if dmtx, X = zeros(n,prod(Vo.dim(1:2))); end
    for i = 1:n
        M = inv(B * inv(Vo.mat) * Vi(i).mat);
        d = spm_slice_vol(Vi(i), M, Vo.dim(1:2), [interp,NaN]);
        if dmtx, X(i,:) = d(:)'; else eval(['i',num2str(i),'=d;']); end
    end
    
    eval(['Yp = ' f ';']);


    Y(:,:,p) = reshape(Yp,Vo.dim(1:2));
end

% Y(isnan(Y)) = 0;

image_resized = Y;

% V = Vo;
% 
% %% part spm_write_vol
% 
% use_offset = false;
% 
% dim = [size(Y) 1 1 1];
% 
% if ~isfield(V,'pinfo')
%     V.pinfo = [1;0;0];
%     rescal  = true;
% elseif ~all(isfinite(V.pinfo(1:2))) || V.pinfo(1) == 0
%     V.pinfo(1:2) = [1;0];
%     rescal  = true;
% else
%     rescal  = false;
% end
% 
% clearvars -except Y V use_offset rescal dim
% 
% if rescal
%     dt           = V.dt(1);
%     s            = find(dt == [2 4 8 256 512 768]);
%     if isempty(s)
%         V.pinfo(1:2) = [1;0];
%     else
%         dmnmx        = [0 -2^15 -2^31 -2^7 0 0 ; 2^8-1 2^15-1 2^31-1 2^7-1 2^16-1 2^32-1];
%         dmnmx        = dmnmx(:,s);
%         mxs          = zeros(dim(3),1)+NaN;
%         mns          = zeros(dim(3),1)+NaN;
% 
%         for p=1:dim(3)
%             tmp      = double(Y(:,:,p));
%             tmp      = tmp(isfinite(tmp));
%             if ~isempty(tmp)
%                 mxs(p) = max(tmp);
%                 mns(p) = min(tmp);
%             end
%         end
% 
%         mx = max(mxs(isfinite(mxs)));
%         mn = min(mns(isfinite(mns)));
%         if isempty(mx), mx = 0; end
%         if isempty(mn), mn = 0; end
%         if mx ~= mn
%             if use_offset
%                 V.pinfo(1,1) = (mx-mn)/(dmnmx(2)-dmnmx(1));
%                 V.pinfo(2,1) = (dmnmx(2)*mn-dmnmx(1)*mx)/(dmnmx(2)-dmnmx(1));
%             else
%                 if dmnmx(1) < 0
%                     V.pinfo(1) = max(mx/dmnmx(2),mn/dmnmx(1));
%                 else
%                     V.pinfo(1) = mx/dmnmx(2);
%                 end
%                 V.pinfo(2) = 0;
%             end
%         else
%             V.pinfo(1,1) = mx/dmnmx(2);
%             V.pinfo(2,1) = 0;
%         end
%     end
% end
% 
% 
% %% part spm_create_vol
% 
% clearvars -except V Y
% 
% for i=1:numel(V)
%     v = create_vol(V(i));
%     
%     f = fieldnames(v);
%     for j=1:numel(f)
%         V(i).(f{j}) = v.(f{j});
%     end
% end
% 
% dat = Y;
% 
% %% part spm_write_plane
% 
% clearvars -except V dat
% 
% n = ':';
% 
% n1 = num2cell(V.n);
% n  = {n n1{:}};
% 
% S      = struct('type','()','subs',{{':',':',n{:}}});
% V.private.dat = subsasgn(V.private.dat,S,dat);
% 
% %% part spm_read_vols
% 
% image_resized = spm_read_vols(V);
% 
% %% part spm_vol
% 
% n = [];
% 
% N = V.private;
% 
% % V = struct('fname',   {},...
% %                'dim',     {},...
% %                'dt',      {},...
% %                'pinfo',   {},...
% %                'mat',     {},...
% %                'n',       {},...
% %                'descrip', {},...
% %                'private', {});
% 
% n  = [n 1 1];
% n  = n(1:2);
% dm = [N.dat.dim 1 1 1 1];
% 
% dt = struct(N.dat);
% dt = double([dt.dtype dt.be]);
% 
% off = (n(1)-1+dm(4)*(n(2)-1))*ceil(spm_type(dt(1),'bits')*dm(1)*dm(2)/8)*dm(3) + ...
%     N.dat.offset;
% 
% V.pinfo = [];
% 
% V.pinfo = [N.dat.scl_slope N.dat.scl_inter off]';
% 
% 
% %% part spm_read_vols
% 
% clearvars -except V
% 
% mask = 0;
% 
% n = numel(V);
% Y = zeros([V(1).dim(1:3),n]);
% 
% for i=1:n, for p=1:V(1).dim(3)
%     Y(:,:,p,i) = spm_slice_vol(V(i),spm_matrix([0 0 p]),V(i).dim(1:2),0);
% end, end
% 
% if n==1, Y=Y(:,:,:,1); end
% 
% 
% 
% 
% %% part spm_create_vol
% 
% 
% %==========================================================================
% %-function V = create_vol(V)
% %==========================================================================
% 
% function V = create_vol(V)
%     
% if ~isfield(V,'n')
%     V.n = [1 1];
% else
%     V.n = [V.n(:)' 1 1];
%     V.n =  V.n(1:2);
% end
% 
% if ~isfield(V,'dt')
%     V.dt = [spm_type('float64') spm_platform('bigend')];
% end
% dt{1} = spm_type(V.dt(1));
% 
% if V.dt(2), dt{2} = 'BE'; else dt{2} = 'LE'; end
% 
% if ~isfield(V,'pinfo'), V.pinfo      = [Inf Inf 0]'; end
% if size(V.pinfo,1)==2,  V.pinfo(3,:) = 0;            end
% 
% V.fname  = deblank(V.fname);
% ext      = spm_file(V.fname,'ext');
% 
% switch ext
% case {'img'}
%     minoff = 0;
% case {'nii'}
%     minoff = 352; % or 544 for NIfTI-2
% end
% bits   = spm_type(V.dt(1),'bits');
% minoff = minoff + ceil(prod(V.dim(1:2))*bits/8)*V.dim(3)*(V.n(1)-1+V.n(2)-1);
% V.pinfo(3,1) = max(V.pinfo(3,:),minoff);
% 
% if ~isfield(V,'descrip'), V.descrip = '';     end
% if ~isfield(V,'private'), V.private = struct; end
% 
% dim    = [V.dim(1:3) V.n];
% dat    = file_array(V.fname,dim,[dt{1} '-' dt{2}],0,V.pinfo(1),V.pinfo(2));
% N      = nifti;
% N.dat  = dat;
% N.mat  = V.mat;
% N.mat0 = V.mat;
% N.mat_intent  = 'Aligned';
% N.mat0_intent = 'Aligned';
% N.descrip = V.descrip;
% %try, N.timing = V.private.timing; end
% 
% % try
% %     N0  = nifti(V.fname);
% % 
% %     % Just overwrite if both are single volume files.
% %     tmp = [N0.dat.dim ones(1,5)];
% %     if prod(tmp(4:end))==1 && prod(dim(4:end))==1
% %         N0 = [];
% %     end
% % catch
%     N0  = [];
% % end
% 
% if ~isempty(N0)
% 
%     % If the dimensions differ, then there is the potential for things to go badly wrong.
%     tmp = [N0.dat.dim ones(1,5)];
% 
%     N.dat.dim = [dim(1:3) max(dim(4:5),tmp(4:5))];
% 
%     if V.n(1)==1
% 
%         % Ensure volumes 2..N have the original matrix
%         nt = size(N.dat,4);
%         if nt>1 && sum(sum((N0.mat-V.mat).^2))>1e-8
%             M0 = N0.mat;
%             if ~isfield(N0.extras,'mat')
%                 N0.extras.mat = zeros([4 4 nt]);
%             else
%                 if size(N0.extras.mat,4)<nt
%                     N0.extras.mat(:,:,nt) = zeros(4);
%                 end
%             end
%             for i=2:nt
%                 if sum(sum(N0.extras.mat(:,:,i).^2))==0
%                     N0.extras.mat(:,:,i) = M0;
%                 end
%             end
%             N.extras.mat = N0.extras.mat;
%         end
% 
%         N0.mat = V.mat;
%         if strcmp(N0.mat0_intent,'Aligned'), N.mat0 = V.mat; end
%         if ~isempty(N.extras) && isstruct(N.extras) && isfield(N.extras,'mat') &&...
%             size(N.extras.mat,3)>=1
%             N.extras.mat(:,:,V.n(1)) = V.mat;
%         end
%     else
%         if sum(sum((N0.mat-V.mat).^2))>1e-8
%             N.extras.mat(:,:,V.n(1)) = V.mat;
%         end
%     end
% 
%     if ~isempty(N0.extras) && isstruct(N0.extras) && isfield(N0.extras,'mat')
%         N0.extras.mat(:,:,V.n(1)) = N.mat;
%         N.extras                  = N0.extras;
%     end
%     if sum((V.mat(:)-N0.mat(:)).^2) > 1e-4
%         N.extras.mat(:,:,V.n(1)) = V.mat;
%     end
% end
% 
% if isfield(N.extras,'mat')
%     M0 = N.mat;
%     for i=1:size(N.extras.mat,3)
%         if sum((M0-N.extras.mat(:,:,i)).^2) < 1e-8
%             N.extras.mat(:,:,i) = 0;
%         end
%     end
%     if sum(N.extras.mat(:).^2) < 1e-8*size(N.extras.mat,3)
%         N.extras = [];
%         if ~isempty(N0) && ~isempty(N0.extras) && isstruct(N0.extras) && isfield(N0.extras,'mat')
%             spm_unlink(spm_file(N.dat.fname,'ext','mat'));
%         end
%     end
% end
% 
% V.private = N;

