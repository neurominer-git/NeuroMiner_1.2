function res = JuSpace_noGUI_2D(list1,atlas,options,PET_list,image_for_size,dir_save)


addpath(genpath('/opt/matlab/JuSpace_v1.3/'));

list2 = [];

if ~ischar(image_for_size)
   image_for_size = char(image_for_size); 
end

% list1: cell with filenames
% list2 (options

% JuSpace directory
dir_tool = fileparts(which('JuSpace'));

% atlas file
% atlas = fullfile(dir_tool,'atlas',atlas);

% PET directory and selected files
dir_PET = fullfile(dir_tool,'PETatlas/');
aa = dir(dir_PET);
names_PET = {aa(3:end).name}';
files_PET_all = strcat(repmat(dir_PET, size(aa(3:end),1),1), names_PET);
filesPET = files_PET_all(PET_list,1);


for i = 1:size(names_PET,1)
   tt = regexp(names_PET{i},'_','Split');
   Rec_list{i}=tt{1};
end

PET_names = Rec_list(PET_list);

opt_comp = options(1);
opt_ana = options(2);
opt_perm = options(3);
opt_auto = options(4);
opt_perm_spatial = options(5);

file_part = '';
switch options(1)
    case 1
        file_part = [file_part '_ESb'];
    case 2
        file_part = [file_part '_ESw'];
    case 3
        file_part = [file_part '_mList1'];
    case 4
        file_part = [file_part '_List1Each'];
    case 5
        file_part = [file_part '_indZ'];
    case 6
        file_part = [file_part '_pwDiff'];
    case 7
        file_part = [file_part '_looList1'];
    case 8
        file_part = [file_part '_List1EachGroupTest'];
end

switch options(2)
    case 1
        file_part = [file_part '_Spearman'];
    case 2
        file_part = [file_part '_Pearson'];
    case 3
        file_part = [file_part '_multReg'];
end

if options(3) == 1 
    file_part = [file_part '_withExactP'];
end

if options(5) == 1
    file_part = [file_part '_withExactSpatialP'];
end

% global files_PET
% aa = get(handles.PETlist,'Value');
% str_PET = get(handles.PETlist,'String');
% str_PET = str_PET(aa);
% for i = 1:length(str_PET);
%    tt = regexp(str_PET{i},'_','Split');
%    Rec_list{i}=tt{1};
% end
% filesPET = files_PET(aa);
% atlas = char(get(handles.atlaslist,'String'));

time_now = datestr(datetime('now'),'ddmmmyyyy_HHMMSS');

% file_save = fullfile(dir_save,['Results_',name_save file_part '_' time_now '.mat']);

if options(1)<4
%     image_save = fullfile(dir_save,[name_save file_part '_' time_now '.nii']);
    [res,p_all,stats,data,D1,D2, data_PET,Resh,T1] = compute_DomainGauges_2D(list1,list2,filesPET,atlas,options,image_for_size,dir_save,image_save);
else
    [res,p_all,stats,data,D1,D2, data_PET,Resh,T1] = compute_DomainGauges_2D(list1,list2,filesPET,atlas,options,image_for_size,dir_save);
end

% opt_for_perm = [1,2,5,6];
% opt_for_spat_perm = [3,4,8];
% 
% if options(3)==1 && ismember(options(1),opt_for_perm)% && options(2)~=3
%     disp('Computing exact p-value');
%    [p_exact,dist_r] = compute_exact_pvalue(D1,D2,data_PET,res,Nperm,options,T1);
%    Resh(:,end+1) = [{'p_exact'}; num2cell_my(p_exact')];
% end
% 
% if options(5)==1 && ismember(options(1),opt_for_spat_perm)
%     disp('Computing exact spatial p-value')
%     [p_exact,dist_r] = compute_exact_spatial_pvalue(D1,data_PET,atlas,res,Nperm,options,filesPET, T1);
%     Resh(:,end+1) = [{'p_exact'}; num2cell_my(p_exact')];
% end

% 
% switch options(1)
%     case 1
%         ff = 'Effect Size';
%     case 2
%         ff = 'Effect Size';
%     case 3
%         ff = 'Mean list';
%     case 4
%         ff = 'List 1 image';
%     case 5
%         ff = 'Z-score file';
%     case 6
%         ff = 'Delta';
%     case 7
%         ff = 'Z-score loo';
%     case 8
%         ff = 'List 1 all against null';
% end
% 
% 
% if exist('dist_r','var')
%     save(file_save,'res','p_all','stats','data','Resh','D1','D2','data_PET','list1','list2','filesPET','atlas','options','dist_r','p_exact','Nperm');
% else
%     save(file_save,'res','p_all','stats','data','Resh','D1','D2','data_PET','list1','list2','filesPET','atlas','options');
% end
% 
% res = stats.res_ind;
% ind_plot =[];
% for i = 1:size(res,1)
%     for j = 1:size(res,2)
%        ind_plot(end+1,1) = j;
%     end
% end
% 
% vals_plot= res; %cell2num_my(Resh(2:end,3));
% for i = 1:length(filesPET)
%    all_i = vals_plot(:,i);%vals_plot(ind_plot==i);
%    
%    size_y_xi = sum(ind_plot==i);
%    y_dist = abs(vals_plot(ind_plot==i) - mean(vals_plot(ind_plot==i)));
%    y_dist_inv = abs(y_dist -max(y_dist))./max(abs(y_dist -max(y_dist)));
%    if isnan(y_dist_inv)
%        y_dist_inv = rand(size(y_dist_inv));
%    end
%    
%    x_n_n(ind_plot==i) = ind_plot(ind_plot==i) + 1.*(rand(size_y_xi,1)-0.5).*y_dist_inv;
%    all_ii = all_i(~isinf(all_i));
%    m_x(i) = mean(all_ii);
%    std_x(i,1) = 1.95996.*std(all_ii)./sqrt(length(all_ii));
% 
% end
% 
% % plot specs
% marker_color = [0 0.7 1];
% bar_color = generate_colors_nice_my(length(filesPET));
% dot_edge = [0 0.4 1];
% error_color = [0.2 0.3 0.4];
% error_thickness = 3;
% marker_size = 20;
% %_______
% 
% %-----------------------
% % generate bar plots
% %-----------------------
% % hold('off');
% if length(filesPET)<3
%     h1 = figure('Position',[400 400 500 500],'visible','off');
% else
%     h1 = figure('Position',[400 400 length(filesPET).*120 500],'visible','off');
% end
% 
% h2 = bar(1:length(filesPET),diag(m_x),0.9,'stacked');
% %colormap(bar_color)
% % bar(handles.plot_bars, 1:length(filesPET),m_x,0.9,'EdgeColor',bar_edge);
% hold('on');
% 
% ax = gca;
% set(ax, 'XTick', 1:size(PET_names,2));
% set(ax,'XTickLabel',PET_names,'FontSize',16);
% vv = vals_plot';
% vv2 = vv(:);
% 
% for i = 1:length(filesPET)
%     set(h2(i),'facecolor',bar_color(i,:));
%     if opt_comp>=4
%         plot(x_n_n(ind_plot==i),vv2(ind_plot==i),'o','MarkerSize',8,'MarkerFaceColor',bar_color(i,:),'MarkerEdgeColor',error_color,'LineWidth',1);
%     end
% end
% if opt_comp>=4
%         h2 = errorbar(1:length(filesPET),m_x,std_x,'.');
%         title('Error bars represent 95% confidence interval','FontSize',12);
%         set(h2,'LineWidth',error_thickness,'Color',error_color,'MarkerSize',marker_size,'MarkerEdgeColor',error_color);
% end
% 
% 
% 
% % if opt_ana<3
% %    ylim([-1 1]);
% % end
% 
% 
% % t = findall(gcf,'Type','Axes');
% % pp = FDR(p_all,0.05);
% % for i = 1:length(p_all)
% %     if p_all(i)<=pp
% %         if opt_comp<4
% %             text(t(end),i,max(m_x)+0.2,'*','FontSize',20);
% %     end
% % end
% 
% switch opt_ana
%     case 1
%         ylabel(ax,'Fisher''s z (Spearman rho)','fontweight','bold');
%     case 2
%         ylabel(ax,'Fisher''s z (Pearson r)','fontweight','bold');
%     case 3
%         ylabel(ax,'Beta coeffiecient','fontweight','bold');
% end
% savefig(h1,fullfile(dir_save,['Bar_' name_save file_part '_' time_now '.fig']));
% print(h1,fullfile(dir_save,['Bar_' name_save file_part '_' time_now '.png']),'-dpng','-r300');
% %-----------------------
% % end generate bar plots
% %-----------------------
% 
% %-----------------------
% % generate scatter plots (only for less than 10 lines)
% %-----------------------
% % if [size(data,1).*size(data_PET,1)]<10
% 
% if opt_ana<3
% %     n_color = 1./(n_draw+1);
%     n_color = 1./(size(data,1)+1);
%     n_rows = ceil(size(data_PET,1)./3);
%    
%   
%         if size(data_PET,1)==1
%             ind_col = 1;
%         elseif size(data_PET,1)==2
%             ind_col = 2;
%         else
%             ind_col = 3;
%         end
%         
%     
%     if ind_col==1
%         hh = figure('Position',[400 0 n_rows.*500 n_rows.*500],'visible','off');
%         subplot(n_rows,size(data_PET,1),1);
%     else
% 
%         hh = figure('Position',[400 0 900 n_rows.*300],'visible','off');
%         subplot(n_rows,ind_col,1);
%     end
%     
%     n_draw = 31;
%     for j = 1:size(data_PET,1)
%         nn = 1;
%         colors_all = [];
%         color_i = [1 0.5 0];
% %         legend_all = {};
% 
%          subplot(n_rows,ind_col,j);
% 
% %         if opt_ana == 1
% %             xlabel(['Relative rank ' Rec_list{j}],'FontSize',12);
% %         else
%             xlabel(['Data ' PET_names{j}],'FontSize',16);
% %         end
%         hold('on');
%             for i = 1:size(data,1)
%                     nn = nn + 1;
%                     if nn < n_draw+1
%                         x = data_PET(j,:);
%                         y = data(i,:)';
% %                         if opt_ana == 1
% %                             [temp,x]  = ismember(x,unique(x));
% %                             [temp,y]  = ismember(y,unique(y));
% %                         end
%                         plot(x,y,'o','MarkerSize',6,'MarkerFaceColor',color_i,'MarkerEdgeColor','black');
%                         colors_all = [colors_all; color_i];
%                         color_i(1) = color_i(1)-n_color;
%                         color_i(3) = color_i(3)+n_color;
% %                         legend_all{end+1} = [ff num2str(i) ' with ' Rec_list{j}];
%                     end
% 
%                 end
% 
% %         if size(data,1)==1
% %             legend(legend_all,'AutoUpdate','off','Location','southoutside');
% %         else
% %         end
%         h1 = lsline;
%         set(h1,'LineWidth',3);
%         for i = 1:length(h1)
%            h1(i).Color = colors_all(i,:); 
%         end
% 
%     end
%     
%     if size(data_PET,1)==1
%         ylabel('Data Modality','FontSize',16);
% %         xlabel('PET','FontSize',20);
%     else
% 
%         h1=subplot(n_rows,ind_col,1);
%         h2=subplot(n_rows,ind_col,ind_col);
%         h3=subplot(n_rows,ind_col,n_rows.*ind_col - (ind_col-1));
% 
% %         h4=subplot(n_rows,3,n_rows.*3);
%         p1=get(h1,'position');
%         p2=get(h2,'position');
%         p3=get(h3,'position');
% %         p4=get(h4,'position');
%         height=p3(2) + (p1(2)-p3(2)) + p1(4)./2;
%         width=p3(1)+ (p2(1)-p3(1)+ p3(3)./2);
%         h5=axes('position',[p3(1) p3(2).*(1 + 1./n_rows) width height],'visible','off'); 
%         h5.XLabel.Visible='on';
%         h5.YLabel.Visible='on';
%         set(gcf, 'CurrentAxes', h5)
% %         axes(h5);
% %         if opt_ana == 1
% %             ylabel('Relative rank modality','FontSize',20);
% %         else
%         ylabel('Data Modality','FontSize',20);
% %         end
% %         xlabel('PET','FontSize',20);
%     end
%     
%     
% %     suplabel('Data PET','x','FontSize',20);
% %     suplabel('Data Modality','y','FontSize',20);
%     savefig(hh,fullfile(dir_save,['Line_' name_save file_part '_' time_now '.fig']));
%     print(hh,fullfile(dir_save,['Line_' name_save file_part '_' time_now '.png']),'-dpng','-r300');
% 
% %     close all
    
end
