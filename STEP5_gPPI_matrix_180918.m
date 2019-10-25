clear all; close all;

% 2018-Sep-18 Yun-An Huang
% using 18ROIS from Laura's study with all subjects, FWE 0-05, Voxel 20, brain mask, 12mm.

% 2018-Jun-23 Yun-An Huang
% using 19ROIS from Laura's study with all subjects, FWE 0-05, Voxel 20.

% 2018-May-28 Yun-An Huang
% using ROIs from Laura's  study

% 2017-Nov-03 Yun-An Huang
% use the 51 ROIs.

% 2017-Sep-04 Yun-An Huang
% use the 62 symmetric ROIs.

% 2017-Jul-31 Yun-An Huang
% correct the contrast name

% 27-Jun-2017 Yun-An Huang
% applied the gPPI analysis to the healthy group of temporal lobectomy
% study.
% analyze 61 ROIs only. 4 overlapped ROI were removed.


% 24-Mar-2017 Yun-An Huang
% applied the permutation test on the common area of general emotion
% pathway.


% 12-Jan-2017 Yun-An Huang
% calculating the common area of general emotion pahway.


% 14-Nov-2016 Yun-An Huang 
% print out the connections 

% 12-Oct-2016 Yun-An Hunag
% 1. use the weighted matrix instead of the binary matrix
% 2. remove the overlap ROI, ROI5: Cerebellum-Crus1-L-1, ROI17: vmPFC-R-1,
% ROI34: Ant-Insula-R,  ROI63: pdACC-R 
% 3. during the dijkstra track the sink could be hub or receiver

% 8-June-2016 Yun-An Huang
% draw the gPPI matrix with percentage scale

% 7-June-2016 Yun-An Huang 
% generating the connection view of significant nodes in BrainNetViewer.

% 23-Mar-2016 Yun-An Huang
% generating the node and edge used for BrainNetViewer.
% 

% 18-Feb-2016 Yun-An Huang
% the Density_threshold_string{density_temp} means the end density
% threshold. but not the fix matrix threshold ex. the range is 0.02-0.08,
% and the fix matrix threshold is 0.05.
% d range: 0.02----0.05------0.08
% d matrix at 0.05

% 15-Feb-2016 Yun-An Huang
% corrected the permutation test result. 
% the significant level increase from > to >=.

% 12-Feb-2016 Yun-An Huang
% corrected the figure of "significant nodes across contrast 2-5 from
% 2adjacency"

% 11-Feb-2016 Yun-An Huang
% modified the significant node type figure

% 21-Aug-2015 Yun-An Huang
% 8 cluster by anatomy, manually.

% 20-Aug-2015 Yun-An Huang
% 1. compute the conjunction of 2 adjacency matrix
% 2. add the reduced gPPI matrix with connected node to significant node type 
% 3. show the causal index with varied density only if the nodes are significant
% 4. show the significant matrix with ROIs are significant.
% 5. separate the contrast 1 and 2-5 in the plot of density and threshold.

% 31-Jul-2015 Yun-An Huang
% 1. drawing the significant matrix.


% 30-Jul-2015 Yun-An Huang
% 1. printing the significant node.

% 28-Jul-2015 Yun-An Huang
% 1. output the image as high resolution png or eps file

% 16-Jul-2015 Yun-An Huang
% 1. build the whole mask of ROIs

% 23-Jun-2015 Yun-An Huang
% 1. computing the 3D brain surface.
% 2. plot the path map on the 2D contour.

% 22-Jun-2015 Yun-An Huang
% computing the conjunction of primary contrast and the difference contrast (the contrast from anova) 
%

% 18-Jun-2015 Yun-An Huang
% 1. draw gPPI matrix with density thresholded
% 2. draw causal array with density thresholded

% 15-Jun-2015 Yun-An Haung
% correct the significant node type with the threshold at density level.

% 12-Jun-2015 Yun-An Huang
% 1. calculating the track from source to sink of gPPI_matrix_sig_node by dijkstra algorithm (breadth first search)
% 2. reduce the gPPI matrix based on tracks results
% 3. drawing the track graph.

% 10-Jun-2015 Yun-An Huang
% 1. calculate the gPPI matrix with significant node type.

% 03-Jun-2015 Yun-An Huang
% 1. calculate the causal index array with varied density threshold
% 2. test by permutation test
% 3. collect source, intermediate, and sink.

% 01-Jun-2015 Yun-An Huang
% 1. compute the density function.
% 2. computing causal index array
% 3. testing by permutation.

% 27-May-2015 Yun-An Huang
% add the anova contrast

% 15-May-2015 Yun-An Huang
% 1. correct the contrast (con) of emotion vs. neutral. since the contrast
% of individual level is set "1 1 1 1 -4", if would be four times large
% than the contrast "1/4 1/4 1/4 1/4 -1".

% 12-May-2015 Yun-An Huang
% 1. substract all emotion from each emotion.

% 01-Apr-2015 Yun-An Huang
% separate the analysis of positive and negative gPPI

% 27-Mar-2015 Yun-An
% the DD, DO, DI, DT is replace by causal index.
% causal index included outdegree, indegree, totaldegree, causalflow.
%

% 17-Feb-2015 Yun-An Huang
% the positive and negative modulation is considered.

% 15-Jan-2015 Yun-An Huang
% the script is used to calculate the contrast of gPPI matrix.
% emotion_vs_neutral
% angry_vs_neutral
% happy_vs_neutral
% fear_vs_neutral
% sad_vs_neutral
%
% moreover, the similarity matrix was calculated.

% 04-Dec-2014 Yun-An Huang
% this is the script by sorting the ROIs according to their pathlength

% 25-Nov-2014 Yun-An Huang
% this script is used to draw the result of group tested gPPI.

% define parameter

ROI_num = 18;
ROI_num_fine = 18; % the exactly analyze ROI 
sbj_num= 12;  % the number of the subjects.


radius_ROI = 6; %6mm

Index_num = 4; % 1=outdegree, 2=indegree, 3= totaldegree, 4=flow.

permute_test_num = 2000; 

Output_contrast_name{1,1} = 'emotion_vs_neutral';
Output_contrast_name{2,1} = 'angry_vs_neutral';
Output_contrast_name{3,1} = 'fear_vs_neutral';
Output_contrast_name{4,1} = 'sad_vs_neutral';
Output_contrast_name{5,1} = 'happy_vs_neutral';

Output_contrast_name{6,1} = 'one_way_anova_f';
Output_contrast_name{7,1} = 'one_way_anova_Angry_t';
Output_contrast_name{8,1} = 'one_way_anova_Fear_t';
Output_contrast_name{9,1} = 'one_way_anova_Sad_t';
Output_contrast_name{10,1} = 'one_way_anova_Happy_t';

Primary_contrast_num = 5; % the contrast 1-5


Sign_num = 2; % Sing=1 for positive PPI, Sign = 2 for negative PPI.

Sign_name{1} = 'positive';
Sign_name{2} = 'negative';


outpath = 'output_data6_gPPI_group_Laura18ROI/'; % "output_data2_gPPI_group_space_65ROI/"1
datafolder = 'output_data2_gPPI_Laura18ROI';


% %%%%%%%%%%%%%%% the output directory %%%%%%
if exist(outpath)
   
    
else
    mkdir(outpath);

end



% %%%%%%%%%%%%%%%%%%%%%%%%% Setting the log file %%%%%%%%%%%%%%%%%%%%%%%%%
diary ([outpath,'/log_draw ',datestr(now,'yyyy-mm-dd HH-MM'),'.txt']);
diary on
warning off 
display('Version: gPPI draw sort 2015.Apr.01');
display('programer: Yun-An Huang');
display(' ');
display(['start at   ',datestr(clock)]);
display(' ');
display(['the output directory: ', outpath]);
display(' ');

% load ROI list
ROIlist_filename = fullfile('ROI_Laura_wholeSB_FWE_05_v20_12mm_brainmask.txt');
ROIname_filename= fullfile('ROIlistname_Laura18ROI_20180918.txt'); % modified 2017-Nov-10


% %%%%% load ROI information %%%%

fid=fopen(ROIlist_filename,'r','n','UTF-8');
ROI_str=fread(fid,inf,'*char');

peak_str_str = strfind(ROI_str','=(');
peak_str_end = strfind(ROI_str',');');

peak_region = zeros(3,ROI_num);
for itemp = 1:ROI_num % retrieve the ROI only from emotional correlation area.
    peak_str =  ROI_str((peak_str_str(itemp)+2):(peak_str_end(itemp)-1));
    peak_temp = strsplit((peak_str'),',');
    peak_region(1,itemp) = str2num(peak_temp{1});
    peak_region(2,itemp) = str2num(peak_temp{2});
    peak_region(3,itemp) = str2num(peak_temp{3});
    
end

% %%%%%%%%%%%%%%%%%% load ROI name

fid=fopen(ROIname_filename,'r','n','UTF-8');
readstr_temp = fread(fid,inf,'*char');

ROI_name=strsplit(readstr_temp','\n');

% %%%%%%%%%%%%% build mask for all ROI %%%%%%%%%%%%%%%%%%

mask_path = fullfile(outpath,'mask'); 
if exist(mask_path)
else
    mkdir(mask_path);
end

VI=spm_vol(fullfile(datafolder,'ROI_01','results','spmT_emotion_vs_neutral.nii'));
img_temp = spm_read_vols(VI);

for rtemp = 1:ROI_num
    mask_new=zeros(VI.dim);
    
    if exist(fullfile(mask_path,['ROI_mask_',num2str(rtemp,'%02d'),'.nii']))
    else
        for xtemp = 1:VI.dim(1)
            for ytemp = 1:VI.dim(2)
                for ztemp = 1:VI.dim(3)
                    peak_temp=VI.mat*[xtemp; ytemp; ztemp; 1];
                    
                    if ((peak_region(:,rtemp)-peak_temp(1:3))'*(peak_region(:,rtemp)-peak_temp(1:3)) <= radius_ROI^2 ) && ~isnan(img_temp(xtemp,ytemp,ztemp))
                    
                        mask_new(xtemp,ytemp,ztemp)=1;
                    end
                    
                end
            end
        end
        
        VO=VI;
        VO.fname = fullfile(mask_path,['ROI_mask_',num2str(rtemp,'%02d'),'.nii']);
        spm_write_vol(VO,mask_new);
                
    end
end

% check the overlap of ROIs

overlap_matrix = zeros(ROI_num,ROI_num);
distance_matrix = zeros(ROI_num,ROI_num);
for r1_temp  = 1:ROI_num
    for r2_temp=r1_temp:ROI_num

        if r1_temp~=r2_temp
        
            dis_temp=(peak_region(:,r1_temp)-peak_region(:,r2_temp))'*(peak_region(:,r1_temp)-peak_region(:,r2_temp))
            if dis_temp <= (2*radius_ROI)^2
                overlap_matrix(r1_temp,r2_temp)=1;
                overlap_matrix(r2_temp,r1_temp)=1;
                distance_matrix(r1_temp,r2_temp)=sqrt((peak_region(:,r1_temp)-peak_region(:,r2_temp))'*(peak_region(:,r1_temp)-peak_region(:,r2_temp)));
                distance_matrix(r2_temp,r1_temp)=sqrt((peak_region(:,r1_temp)-peak_region(:,r2_temp))'*(peak_region(:,r1_temp)-peak_region(:,r2_temp)));
            else
                overlap_matrix(r1_temp,r2_temp)=0;
                overlap_matrix(r2_temp,r1_temp)=0;
            end
        
        end
    
    end
end

idx=find(overlap_matrix);
[idx_raw, idx_column] = ind2sub(size(overlap_matrix),idx)
distance_matrix(idx)

% the new 8 cluster by anatomy manually.

outperm = [1:ROI_num]; % the ROIs is sort in order
cluster_label = [ones(1,ROI_num) ]; % 2017-Sep-04 % 6 cluster only
cluster_label_num = 1;
cluster_color = hsv(cluster_label_num);
% cluster_color(2:4,:)=[1.0000 0.6000 0; 0.9  0.9 0;  0 0.9000 0];
edge_color = cluster_color;
edge_color(find(edge_color>1))=1;

% %%%%%%%% build the whole mask of ROIs %%%%%%
% 16-Jul-2015
%

if exist(fullfile(mask_path,['ROI_mask_all.nii']))

    display(['ROI_maks_all.nii is existed.']);
else
    
    % load parameter from a file    
    VI=spm_vol(fullfile(mask_path,['ROI_mask_',num2str(1,'%02d'),'.nii']));
    img_temp = spm_read_vols(VI);
    
    % create empty space
    mask_all = zeros(VI.dim);
    mask_all_ROI = zeros(VI.dim);
    
    for rtemp = 1:ROI_num
        
        VI_ROI = spm_vol(fullfile(mask_path,['ROI_mask_',num2str(outperm(rtemp),'%02d'),'.nii']));
        VI_ROI_img = spm_read_vols(VI_ROI);
        
        ind_temp = find(VI_ROI_img>0);
        
        mask_all(ind_temp) = cluster_label(rtemp);
        mask_all_ROI(ind_temp) = rtemp;
    end

    VO=VI;
    VO.fname = fullfile(mask_path,['ROI_mask_all.nii']);
    spm_write_vol(VO,mask_all);
    
    VO=VI;
    VO.fname = fullfile(mask_path,['ROI_mask_all_roi.nii']);
    spm_write_vol(VO,mask_all_ROI);
    
    
end


if exist(fullfile(mask_path,['ROI_mask_all_axis_icbm152.nii']))

    display(['ROI_mask_all_axis_icbm152.nii is existed.']);
else
   

    %load the template
    VI_T1=spm_vol(fullfile('mni_icbm152_t1_tal_nlin_asym_09a.nii'));
    img_T1 = spm_read_vols(VI_T1);

    
    % interpolate the img_ROI to img_T1

    [T1_X,T1_Y,T1_Z]= meshgrid(1:VI_T1.dim(1),1:VI_T1.dim(2),1:VI_T1.dim(3)); % the mesh dimension of X and Y is change, if dim(1)=m,dim(2)=n then size(T1_X) = [n,m]
    T1_X_re=reshape(T1_X, 1, length(T1_X(:)));
    T1_Y_re=reshape(T1_Y, 1, length(T1_Y(:)));
    T1_Z_re=reshape(T1_Z, 1, length(T1_Z(:)));
    VI_T1_coordinate_mm = VI_T1.mat*[T1_X_re; T1_Y_re; T1_Z_re; ones(1,length(T1_X_re(:)))];
    T1_X_mm = reshape(VI_T1_coordinate_mm(1,:),VI_T1.dim(2), VI_T1.dim(1),VI_T1.dim(3));
    T1_Y_mm = reshape(VI_T1_coordinate_mm(2,:),VI_T1.dim(2), VI_T1.dim(1),VI_T1.dim(3));
    T1_Z_mm = reshape(VI_T1_coordinate_mm(3,:),VI_T1.dim(2), VI_T1.dim(1),VI_T1.dim(3));

    img_ROI_new_axis = zeros(VI_T1.dim);
    img_ROI_new_axis_ROI = zeros(VI_T1.dim);
    
    for rtemp = 1:ROI_num
       
        % load ROIs
        VI_ROI = spm_vol(fullfile(mask_path,['ROI_mask_',num2str(outperm(rtemp),'%02d'),'.nii']));
        img_ROI = spm_read_vols(VI_ROI);
        
        [ROI_X,ROI_Y,ROI_Z]= meshgrid(1:VI_ROI.dim(1),1:VI_ROI.dim(2),1:VI_ROI.dim(3));
        ROI_X_re=reshape(ROI_X, 1, length(ROI_X(:)));
        ROI_Y_re=reshape(ROI_Y, 1, length(ROI_Y(:)));
        ROI_Z_re=reshape(ROI_Z, 1, length(ROI_Z(:)));
        VI_ROI_coordinate_mm = VI_ROI.mat*[ROI_X_re; ROI_Y_re; ROI_Z_re; ones(1,length(ROI_X_re(:)))];
        ROI_X_mm = reshape(VI_ROI_coordinate_mm(1,:),VI_ROI.dim(2), VI_ROI.dim(1),VI_ROI.dim(3)); 
        ROI_Y_mm = reshape(VI_ROI_coordinate_mm(2,:),VI_ROI.dim(2), VI_ROI.dim(1),VI_ROI.dim(3));
        ROI_Z_mm = reshape(VI_ROI_coordinate_mm(3,:),VI_ROI.dim(2), VI_ROI.dim(1),VI_ROI.dim(3));


        img_ROI_interp = interp3(ROI_X_mm,ROI_Y_mm,ROI_Z_mm,permute(img_ROI,[2 1 3]),T1_X_mm,T1_Y_mm,T1_Z_mm,'linear',0);
        img_ROI_interp_permute = permute(img_ROI_interp,[2 1 3]);
        ind_temp = find(img_ROI_interp_permute>0);
        img_ROI_new_axis(ind_temp) = cluster_label(rtemp);
        img_ROI_new_axis_ROI(ind_temp) = rtemp;
    end


    VO=VI_T1;
    VO.fname = fullfile(mask_path,'ROI_mask_all_axis_icbm152.nii');
    spm_write_vol(VO,img_ROI_new_axis);
    
    
    VO=VI_T1;
    VO.fname = fullfile(mask_path,'ROI_mask_all_axis_roi_icbm152.nii');
    spm_write_vol(VO,img_ROI_new_axis_ROI);
end


%% Computing the gPPI matrix

display(['Computing the gPPI matrix: started at ', datestr(clock)]);


for ctemp = 1:size(Output_contrast_name,1)
        
    gPPI_matrix = zeros(ROI_num_fine,ROI_num_fine,Sign_num); % the 3rd dimension is sign.
    gPPI_matrix_con = zeros(ROI_num_fine,ROI_num_fine,Sign_num); % the 3rd dimension is sign.
    survive_matrix = zeros(ROI_num_fine,ROI_num_fine,Sign_num); % the 3rd dimension is sign.
    display('drawing contrast ', Output_contrast_name{ctemp,1});

    % construct the gPPI matrix.
    if exist(fullfile(outpath,['gPPI_matrix_',Output_contrast_name{ctemp,1},'.mat']))
        load(fullfile(outpath,['gPPI_matrix_',Output_contrast_name{ctemp,1},'.mat']));
    else
        for r1temp = 1:ROI_num_fine

            % the load file was according to the new order from the clustering result
            if ctemp ==6 % the f test of anova
                
                VI=spm_vol(fullfile(datafolder,['ROI_',num2str(outperm(r1temp),'%02d')],'results',['spmF_',Output_contrast_name{ctemp,1},'.nii']));
                c_img = spm_read_vols(VI);

                VI_con=spm_vol(fullfile(datafolder,['ROI_',num2str(outperm(r1temp),'%02d')],'results',['ess_',Output_contrast_name{ctemp,1},'.nii']));
                con_img = spm_read_vols(VI_con);
                
            else % for the t-test
                VI=spm_vol(fullfile(datafolder,['ROI_',num2str(outperm(r1temp),'%02d')],'results',['spmT_',Output_contrast_name{ctemp,1},'.nii']));
                c_img = spm_read_vols(VI);

                VI_con=spm_vol(fullfile(datafolder,['ROI_',num2str(outperm(r1temp),'%02d')],'results',['con_',Output_contrast_name{ctemp,1},'.nii']));
                con_img = spm_read_vols(VI_con);
            end
                
            for r2temp = 1:ROI_num_fine
                % the load file was according to the new order from clustering result
                Vmask = spm_vol(fullfile(mask_path,['ROI_mask_',num2str(outperm(r2temp),'%02d'),'.nii']));
                mask_img = spm_read_vols(Vmask);



                ROI_T = c_img.*mask_img;
                ROI_con = con_img.*mask_img;

                ROI_con(isnan(ROI_T))=[];
                ROI_T(isnan(ROI_T))=[];



                % 2015-Mar-30
                % separate the positive and negative connection

                % positive
                if isempty(find(ROI_T>0))

                    gPPI_matrix(r1temp,r2temp,1) = 0;
                    gPPI_matrix_con(r1temp,r2temp,1) = 0;

                else


                    gPPI_matrix(r1temp,r2temp,1) = max(ROI_T(:));

                    ind_temp = find(ROI_T==max(ROI_T(:)));

                    if ctemp==1

                        gPPI_matrix_con(r1temp,r2temp,1) = max(ROI_con(ind_temp))/4 ; % for the correction of contrast. 2015-May-15
                    else
                        gPPI_matrix_con(r1temp,r2temp,1) = max(ROI_con(ind_temp)) ;
                    end
                end

                % negative
                if isempty(find(ROI_T<0))
                    gPPI_matrix(r1temp,r2temp,2) = 0;
                    gPPI_matrix_con(r1temp,r2temp,2) = 0;
                else
                    gPPI_matrix(r1temp,r2temp,2) = min(ROI_T(:));

                    ind_temp = find(ROI_T==min(ROI_T(:)));

                    if ctemp==1

                        gPPI_matrix_con(r1temp,r2temp,2) = min(ROI_con(ind_temp))/4; % for the correction of contrast. 2015-May-15
                    else
                        gPPI_matrix_con(r1temp,r2temp,2) = min(ROI_con(ind_temp));
                    end
                end

            end

        end

        save(fullfile(outpath,['gPPI_matrix_',Output_contrast_name{ctemp,1},'.mat']),'gPPI_matrix','gPPI_matrix_con');

    end
end

clear gPPI_matrix gPPI_matrix_con

display(['Computing the gPPI matrix: completed at ', datestr(clock)]);



%% Computing the weighted gPPI matrix 
% 2016-OCT-17 Yun-An Huang

display(['Computing the weighted gPPI matrix : started at ', datestr(clock)]);


for ctemp = 1:5 %the one sample t-test
% for ctemp = 1




   load(fullfile(outpath,['gPPI_matrix_',Output_contrast_name{ctemp,1},'.mat'])); % 'gPPI_matrix','gPPI_matrix_con' is loaded
   gPPI_matrix_temp=size(gPPI_matrix(:,:,1));
   DOF = sbj_num-1;% degree of freedom 
   for itemp = 1:size(gPPI_matrix(:,:,1),1)
        for jtemp = 1:size(gPPI_matrix(:,:,1),2)
   
            gPPI_matrix_temp(itemp,jtemp)= spm_t2z(gPPI_matrix(itemp,jtemp,1),DOF);
            
        end
   end

   gPPI_matrix_weighted = zeros(size(gPPI_matrix(:,:,1)));
   gPPI_matrix_weighted = (2*normcdf(gPPI_matrix_temp,0,1)-1).^4;
        
    
    save(fullfile(outpath,['gPPI_matrix_',Output_contrast_name{ctemp,1},'_weighted.mat']),'gPPI_matrix_weighted');
end

        


% clear gPPI_matrix gPPI_matrix_con gPPI_matrix_weighted

display(['Computing the weighted gPPI matrix: completed at ', datestr(clock)]);


% %% Computing the conjunction of gPPI matrix of contrast 7-10
% % 2015-Aug-06 Yun-An Huang
% 
display(['Computing the conjunction of gPPI matrix of contrast 7-10: started at ', datestr(clock)]);

gPPI_matrix_conjunction_7_10 = zeros(ROI_num_fine,ROI_num_fine,Sign_num);
gPPI_matrix_conjunction_7_10_temp = zeros(ROI_num_fine,ROI_num_fine,Sign_num);


for ctemp = 7:10
    
    load(fullfile(outpath,['gPPI_matrix_',Output_contrast_name{ctemp,1},'.mat']));

        
    ind = find(abs(gPPI_matrix)>=abs(gPPI_matrix_conjunction_7_10_temp));
       
    gPPI_matrix_conjunction_7_10(ind)=ctemp-6;
    gPPI_matrix_conjunction_7_10_temp(ind) = gPPI_matrix(ind); 
        
        
end

gPPI_array_conjunction_7_10 = zeros(ROI_num_fine,2,Sign_num); % the second column is output or input

for r_temp = 1:ROI_num_fine

    for sign_temp = 1:Sign_num
        
        max_contrast_num = [0 0];
    
        for ctemp = 1:4

            % for output
            ind = find(gPPI_matrix_conjunction_7_10(r_temp,:,sign_temp)==ctemp);
            
            if length(ind)>=max_contrast_num(1) 
            
                max_contrast_num(1)=length(ind);
                gPPI_array_conjunction_7_10(r_temp,1,sign_temp)=ctemp;
            end
            
            % for input
            ind = find(gPPI_matrix_conjunction_7_10(:,r_temp,sign_temp)==ctemp);
            
            if length(ind)>=max_contrast_num(2) 
            
                max_contrast_num(2)=length(ind);
                gPPI_array_conjunction_7_10(r_temp,2,sign_temp)=ctemp;
            end

        end
    
    end
end

save(fullfile(outpath,['gPPI_matrix_conjunction_7_10.mat']),'gPPI_matrix_conjunction_7_10','gPPI_array_conjunction_7_10');

display(['Computing the conjunction of gPPI matrix of contrast 7-10: completed at ', datestr(clock)]);



%% computing the causal index array of by weighted gPPI matrix function.
% and test causal index array by permutation test
% after the permutation test, collect the source, intermediate and sink.
% source: significant high in outdegree and significant directionality in positive
% intermediate: significant high in outdegree and significant high in indegree
% sink: significant high in indegree and significant directionality in negative
% 
% 2015-Jun-03 Yun-An Huang
%

display(['Computing the causal index array of weighted gPPI matrix: started at ', datestr(clock)]);

% define


% permute_test_num = 20; 

causal_array_all_contrast_weighted_arr_permute = zeros(ROI_num_fine, Index_num, size(Output_contrast_name,1),Sign_num,permute_test_num); % store the distribution of causal array of permutation test
causal_array_all_contrast_weighted_arr_sig = zeros(ROI_num_fine, Index_num, size(Output_contrast_name,1),Sign_num,2); % store the significance of causal array of permutation test
causal_array_all_contrast_weighted_arr = zeros(ROI_num_fine, Index_num, size(Output_contrast_name,1),Sign_num); % store the causal array of permutation test
causal_array_all_contrast_weighted_arr_node = zeros(ROI_num_fine, size(Output_contrast_name,1),Sign_num,length(Alpha_arr)); % store the nodes type, the last dimension 1, 0.05 ; 2, 0.25; 3. 0.05/65
% causal_array_all_contrast_weighted_arr_sig_node=zeros(ROI_num_fine,size(Output_contrast_name,1),Sign_num,length(Alpha_arr),length(Density_threshold_index)); % store the significant node type
% causal_array_all_contrast_weighted_arr_threshold=zeros(size(Output_contrast_name,1),Sign_num,length(Density_thre_arr_fine)); % store the threshold of the density


for ctemp = 1:5 % for one sample t-test only
    display(Output_contrast_name{ctemp})
    for sign_temp = 1:Sign_num
        % load gPPI matrix
        load(fullfile(outpath,['gPPI_matrix_',Output_contrast_name{ctemp,1},'_weighted.mat'])); % gPPI_matrix_weighted is loaded.

       % calculate causal index

            causal_array_all_contrast_weighted_arr(:,1,ctemp,sign_temp) = sum(gPPI_matrix_weighted,2); % outdegree
            causal_array_all_contrast_weighted_arr(:,2,ctemp,sign_temp) = sum(gPPI_matrix_weighted,1); % indegree    
            causal_array_all_contrast_weighted_arr(:,3,ctemp,sign_temp) = causal_array_all_contrast_weighted_arr(:,1,ctemp,sign_temp)+causal_array_all_contrast_weighted_arr(:,2,ctemp,sign_temp); % totaldegree    
            causal_array_all_contrast_weighted_arr(:,4,ctemp,sign_temp) = causal_array_all_contrast_weighted_arr(:,1,ctemp,sign_temp)-causal_array_all_contrast_weighted_arr(:,2,ctemp,sign_temp); % flow    

            % permutation test

            for p_temp = 1:permute_test_num
                original_array = reshape(gPPI_matrix_weighted,ROI_num_fine*ROI_num_fine,1);
                random_array = rand(ROI_num_fine*ROI_num_fine,1);
                [sort_temp IX]=sort(random_array);

                gPPI_matrix_weighted_permute = reshape(original_array(IX),ROI_num_fine,ROI_num_fine);

                causal_array_all_contrast_weighted_arr_permute(:,1,ctemp,sign_temp,p_temp) = sum(gPPI_matrix_weighted_permute,2); % outdegree
                causal_array_all_contrast_weighted_arr_permute(:,2,ctemp,sign_temp,p_temp) = sum(gPPI_matrix_weighted_permute,1); % indegree    
                causal_array_all_contrast_weighted_arr_permute(:,3,ctemp,sign_temp,p_temp) = causal_array_all_contrast_weighted_arr_permute(:,1,ctemp,sign_temp,p_temp)+causal_array_all_contrast_weighted_arr_permute(:,2,ctemp,sign_temp,p_temp); % totaldegree    
                causal_array_all_contrast_weighted_arr_permute(:,4,ctemp,sign_temp,p_temp) = causal_array_all_contrast_weighted_arr_permute(:,1,ctemp,sign_temp,p_temp)-causal_array_all_contrast_weighted_arr_permute(:,2,ctemp,sign_temp,p_temp); % flow    

            end

             % calculate the distribution of permutation.
             
                for rtemp = 1:ROI_num_fine
                    for index_temp = 1:Index_num

                        % for causal degree
                        sig_temp = causal_array_all_contrast_weighted_arr_permute(rtemp,index_temp,ctemp,sign_temp,:) >= causal_array_all_contrast_weighted_arr(rtemp,index_temp,ctemp,sign_temp);
                        sig_level = sum(sig_temp(:))/length(sig_temp(:));

                        sig_temp_l = causal_array_all_contrast_weighted_arr_permute(rtemp,index_temp,ctemp,sign_temp,:) <= causal_array_all_contrast_weighted_arr(rtemp,index_temp,ctemp,sign_temp);
                        sig_level_l = sum(sig_temp_l(:))/length(sig_temp_l(:));



                        if sig_level < sig_level_l
                            causal_array_all_contrast_weighted_arr_sig(rtemp,index_temp,ctemp,sign_temp,1) = sig_level;
                            causal_array_all_contrast_weighted_arr_sig(rtemp,index_temp,ctemp,sign_temp,2) = 1;

                        else
                            causal_array_all_contrast_weighted_arr_sig(rtemp,index_temp,ctemp,sign_temp,1) = sig_level_l;
                            causal_array_all_contrast_weighted_arr_sig(rtemp,index_temp,ctemp,sign_temp,2) = -1;

                        end

                    end


                end

            
            % seperate the nodes type
            for rtemp = 1:ROI_num_fine

                for alpha_temp = 1:length(Alpha_arr)
                    
                   alpha_t=Alpha_arr(alpha_temp);               


%                     % 2016-Feb-15 new
%                     
                    if causal_array_all_contrast_weighted_arr_sig(rtemp,1,ctemp,sign_temp,1)<alpha_t && causal_array_all_contrast_weighted_arr_sig(rtemp,1,ctemp,sign_temp,2)==1 && causal_array_all_contrast_weighted_arr_sig(rtemp,4,ctemp,sign_temp,1)<alpha_t/2 && causal_array_all_contrast_weighted_arr_sig(rtemp,4,ctemp,sign_temp,2)==1
                    % the outdegree significant higher and directionality significant positive.

                        causal_array_all_contrast_weighted_arr_node(rtemp,ctemp,sign_temp,alpha_temp) = 1;

                    elseif causal_array_all_contrast_weighted_arr_sig(rtemp,2,ctemp,sign_temp,1)<alpha_t && causal_array_all_contrast_weighted_arr_sig(rtemp,2,ctemp,sign_temp,2)==1 && causal_array_all_contrast_weighted_arr_sig(rtemp,4,ctemp,sign_temp,1)<alpha_t/2 && causal_array_all_contrast_weighted_arr_sig(rtemp,4,ctemp,sign_temp,2)==-1
                     % the indegree significant higher and directionality significant negative.   

                        causal_array_all_contrast_weighted_arr_node(rtemp,ctemp,sign_temp,alpha_temp) = 3;
                        
                    elseif causal_array_all_contrast_weighted_arr_sig(rtemp,3,ctemp,sign_temp,1)< alpha_t && causal_array_all_contrast_weighted_arr_sig(rtemp,3,ctemp,sign_temp,2)==1 
                     % the outdegree, indegree and totaldegree significant higher and not
                     % significant in directionality.
                        causal_array_all_contrast_weighted_arr_node(rtemp,ctemp,sign_temp,alpha_temp) = 2;
                       
                    end
                end

            end
            
  
            
        
        
        
    end
end

save(fullfile(outpath,'causal_array_all_contrast_weighted_arr.mat'),'causal_array_all_contrast_weighted_arr','causal_array_all_contrast_weighted_arr_sig','causal_array_all_contrast_weighted_arr_node');

display(['Computing the causal index array of weighted gPPI matrix: completed at ', datestr(clock)]);



%% drawing weighted gPPI matrix 

% 2016-OCT-19 Yun-An Huang
%

display(['Drawing the weighted gPPI matrix: started at ', datestr(clock)]);

load(fullfile(outpath,'causal_array_all_contrast_weighted_arr.mat')); % 'causal_array_all_contrast_weighted_arr','causal_array_all_contrast_weighted_arr_sig','causal_array_all_contrast_weighted_arr_node') is loaded

significant_temp=zeros(ROI_num_fine,Index_num, size(Output_contrast_name,1),Sign_num,2); % the last dimension 1 for p<.05 , 2 for bafferoni correction.

% outdegree, indegree, totaldegree
significant_temp(:,1:3,:,:,1) = causal_array_all_contrast_weighted_arr_sig(:,1:3,:,:,1) < 0.05 & causal_array_all_contrast_weighted_arr_sig(:,1:3,:,:,2) == 1 ; % the significant level of 2 sample t-test with permutation test, one tail <0.05.
significant_temp(:,1:3,:,:,2) = causal_array_all_contrast_weighted_arr_sig(:,1:3,:,:,1) < 0.05/ROI_num_fine & causal_array_all_contrast_weighted_arr_sig(:,1:3,:,:,2) == 1 ; % one tail <0.05

% flow 
significant_temp(:,4,:,:,1) = causal_array_all_contrast_weighted_arr_sig(:,4,:,:,1) < 0.025 ; % the significant level of 2 sample t-test with permutation test, two tail <0.025.
significant_temp(:,4,:,:,2) = causal_array_all_contrast_weighted_arr_sig(:,4,:,:,1) < 0.025/ROI_num_fine ; % two tail <0.025

for ctemp = 1 % for the general emotion only
    
    load(fullfile(outpath,['gPPI_matrix_',Output_contrast_name{ctemp,1},'_weighted.mat'])); % 'gPPI_matrix_weighted' is loaded;
    
    
    
    for sign_temp = 1 % positive only
        
        T_range = [-1 1];
        C_coefficients = 100;
        cmap = [jet((T_range(2)-T_range(1))*C_coefficients+1); cluster_color;1 1 1;]; % merge two colormap
%         cmap_cir = [jet(T_range(2)-T_range(1)+1); cluster_color;1 1 1;]; % merge two colormap
        % the circle color map
        reso_num = 200;
        colormap_matrix = zeros(reso_num,reso_num);
        color_radius1 = reso_num/2*0.9;
        color_radius2 = reso_num/2*0.4;

        centor_cor= [reso_num/2 reso_num/2];
        for itemp = 1:reso_num
            for jtemp = 1:reso_num

                radius_temp=sqrt((itemp-centor_cor(1))^2+(jtemp-centor_cor(2))^2);
                angle_temp = angle((itemp-centor_cor(1))+i*(jtemp-centor_cor(2)));
                if angle_temp <0
                    angle_temp = angle_temp+2*pi;
                end
                if radius_temp<color_radius1 && radius_temp>color_radius2 && angle_temp >=0 && angle_temp<=(3/2*pi) 

                    if sign_temp
                        colormap_matrix(itemp,jtemp)=(1-angle_temp/(3/2*pi))*T_range(2)*C_coefficients;
                    else
                        colormap_matrix(itemp,jtemp)=(1-angle_temp/(3/2*pi))*T_range(1)*C_coefficients;
                    end
                else
                    
                    
                    colormap_matrix(itemp,jtemp)=T_range(2)*C_coefficients+cluster_label_num+1;
                end
            end
        end
        
            % combine causal degree index together.
                     
            
            pic = figure(3);
            clf
            set(gcf,'InvertHardCopy','off');
            set(gcf,'Color',[1,1,1]);
            set(gcf,'Position',[200 200 500 500]);
            
            yRange=[0 ROI_num_fine];
            subplot('position',[0.15 0.86 0.65 0.14]) % indegree bar plot
            for btemp = 1:cluster_label_num
                index_temp=find(cluster_label==btemp);
                h=bar(index_temp,causal_array_all_contrast_weighted_arr(index_temp,2,ctemp,sign_temp),'facecolor',cluster_color(btemp,:),'edgecolor',edge_color(btemp,:));ylim(yRange);xlim([0.5 ROI_num_fine+0.5]);axis off; 
                hold on
            end
            
            x=1:ROI_num_fine;sig_temp = significant_temp(:,2,ctemp,sign_temp,1); x(sig_temp==0)=[];cluster_temp=cluster_label;cluster_temp(sig_temp==0)=[];sig_temp(sig_temp==0)=[];
            for ptemp = 1:length(sig_temp)
                plot(x(ptemp),sig_temp(ptemp)*yRange(2)*0.9,'O','markerfacecolor',cluster_color(cluster_temp(ptemp),:),'color',[0 0 0],'markersize',2.5);            
            end
            hold off

            subplot('position',[0 0.2 0.14 0.65]) % outdegree bar plot
            for btemp = 1:cluster_label_num 
                index_temp=find(cluster_label==btemp);
                h=barh(index_temp,causal_array_all_contrast_weighted_arr(index_temp,1,ctemp,sign_temp),'facecolor',cluster_color(btemp,:),'edgecolor',edge_color(btemp,:));xlim(yRange);ylim([0.5 ROI_num_fine+0.5]);axis off; 
                hold on
            end
            
            x=1:ROI_num_fine;sig_temp = significant_temp(:,1,ctemp,sign_temp,1); x(sig_temp==0)=[];cluster_temp=cluster_label;cluster_temp(sig_temp==0)=[];sig_temp(sig_temp==0)=[];
            for ptemp = 1:length(sig_temp)
                plot(sig_temp(ptemp)*yRange(2)*0.9,x(ptemp),'O','markerfacecolor',cluster_color(cluster_temp(ptemp),:),'color',[0 0 0],'markersize',2.5);     
            end
            
            set(get(h,'parent'),'XDir','reverse')
            set(get(h,'parent'),'YDir','reverse')
            hold off

            subplot('position',[0.15 0.2 0.65 0.65]) % gppi matrix
            imagesc(gPPI_matrix_weighted(:,:)*C_coefficients,[T_range(1)*C_coefficients T_range(2)*C_coefficients+cluster_label_num+1]); colormap(cmap);
            set(gca,'xtick',linspace(0.5,ROI_num_fine+0.5, ROI_num_fine+1),'ytick',linspace(0.5,ROI_num_fine+0.5, ROI_num_fine+1));
            set(gca,'xticklabel',[],'yticklabel',[]);
            set(gca,'xgrid','on','ygrid','on','gridlinestyle','-','xcolor',[0.4 0.8 0.4],'ycolor',[0.4 0.8 0.4]);
            saveas(pic, fullfile(outpath,['gPPI_matrixnbar_',Output_contrast_name{ctemp,1},'_',Sign_name{sign_temp},'_weighted.jpg']));

            subplot('position',[0.81 0.04 0.14 0.14]) % the half circle color bar.
            imagesc(colormap_matrix,[T_range(1)*C_coefficients T_range(2)*C_coefficients+cluster_label_num+1]); 
            if ctemp ==6
                text(reso_num*0,reso_num*0.85,{'Weighted'},'fontsize',6,'FontWeight','bold');
            else
                text(reso_num*0,reso_num*0.85,{'Weighted'},'fontsize',6,'FontWeight','bold');
            end
            if sign_temp
                text(reso_num*0.1,reso_num*0.4,{'w=0'},'fontsize',6,'color','k');
                text(reso_num*0.5,reso_num*0.85,{num2str(T_range(2))},'fontsize',6,'color','w');
            else
                text(reso_num*0.1,reso_num*0.4,{'w=0'},'fontsize',6,'color','k');
                text(reso_num*0.5,reso_num*0.85,{num2str(T_range(1))},'fontsize',6,'color','w');
            end
            axis off;

%             subplot('position',[0.805 0.2 0.02 0.65]) % square color bar 
%             imagesc((T_range(2):-1:T_range(1))',[T_range(1) T_range(2)+cluster_label_num]); axis off; colormap(cmap);
% %             text([0.3 0.5],[-20 1020],{num2str(T_range(2)),num2str(T_range(1))});
%             if ctemp == 6
%                 text(1,length(T_range(1):T_range(2))*0.5,'F value','HorizontalAlignment','center','VerticalAlignment','middle','Rotation',90,'color',[1 0.2 1],'fontsize',8);
%                 text(1,length(T_range(1):T_range(2))*0.1,['F = ',num2str(T_range(2))],'HorizontalAlignment','center','VerticalAlignment','middle','Rotation',90,'color',[1 1 1],'fontsize',8);
%                 text(1,length(T_range(1):T_range(2))*0.9,['F = ',num2str(T_range(1))],'HorizontalAlignment','center','VerticalAlignment','middle','Rotation',90,'color',[1 1 1],'fontsize',8);
%             else
%                 text(1,length(T_range(1):T_range(2))*0.5,'T value','HorizontalAlignment','center','VerticalAlignment','middle','Rotation',90,'color',[1 0.2 1],'fontsize',8);
%                 text(1,length(T_range(1):T_range(2))*0.1,['T = ',num2str(T_range(2))],'HorizontalAlignment','center','VerticalAlignment','middle','Rotation',90,'color',[1 1 1],'fontsize',8);
%                 text(1,length(T_range(1):T_range(2))*0.9,['T = ',num2str(T_range(1))],'HorizontalAlignment','center','VerticalAlignment','middle','Rotation',90,'color',[1 1 1],'fontsize',8);
%             end

            % plot cluster
            ax3=subplot('Position',[0.805,0.2 ,0.02, 0.65]);
            imagesc(cluster_label(1:ROI_num_fine)'+T_range(2)*C_coefficients,[T_range(1)*C_coefficients T_range(2)*C_coefficients+cluster_label_num+1]); colormap(cmap); % for two colormap, cluster_label should above T_range
            text(1.7*ones(1,ROI_num_fine),1:ROI_num_fine,{ROI_name{outperm(1:ROI_num_fine)}},'fontsize',4);
            set(ax3,'ytick',linspace(0.5,ROI_num_fine+0.5, ROI_num_fine+1),'xtick',[]);
            set(ax3,'xticklabel',[],'yticklabel',[]);
            set(ax3,'ygrid','on','gridlinestyle','-','ycolor',[0.8 0.8 0.8],'xcolor',[0.8 0.8 0.8]);
            
            % plot cluster
            ax3=subplot('Position',[0.15,0.170 ,0.65, 0.02]);
            imagesc(cluster_label(1:ROI_num_fine)+T_range(2)*C_coefficients,[T_range(1)*C_coefficients T_range(2)*C_coefficients+cluster_label_num+1]); colormap(cmap);
            text(1:ROI_num_fine,1.7*ones(1,ROI_num_fine),{ROI_name{outperm(1:ROI_num_fine)}},'HorizontalAlignment','Right','Rotation',90,'fontsize',5);
            set(ax3,'xtick',linspace(0.5,ROI_num_fine+0.5, ROI_num_fine+1),'ytick',[]);
            set(ax3,'xticklabel',[],'yticklabel',[]);
            set(ax3,'xgrid','on','gridlinestyle','-','xcolor',[0.8 0.8 0.8],'ycolor',[0.8 0.8 0.8]');
            
            
%             subplot('position',[0.856 0.2 0.144 0.65]) % source label
%             text(3*ones(1,ROI_num),ROI_num-0.4:-1:0,{ROI_name{outperm}});
%             axis off;ylim([0 65]);xlim([0 150])
% 
%             subplot('position',[0.15 0 0.65 0.14]) % sink label
%             axis off;xlim([0 65]);ylim([0 150]);
%             text(0.5:ROI_num, 150*ones(1,ROI_num),{ROI_name{outperm}},'HorizontalAlignment','Right','Rotation',90);

            saveas(pic, fullfile(outpath,['gPPI_matrixnbar_',Output_contrast_name{ctemp,1},'_',Sign_name{sign_temp},'_weighted_label.jpg']));
            

            saveas(pic, fullfile(outpath,['gPPI_matrixnbar_',Output_contrast_name{ctemp,1},'_',Sign_name{sign_temp},'_weighted_label.eps']),'epsc');
            export_fig(pic, '-dpng','-r400', fullfile(outpath,['gPPI_matrixnbar_',Output_contrast_name{ctemp,1},'_',Sign_name{sign_temp},'_weighted_label.png']));


        
    end
    
end

 


display(['Drawing the weighted gPPI matrix: completed at ', datestr(clock)]);
