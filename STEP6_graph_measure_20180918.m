clear all; close all;

% 2018-Sep-18 Yun-An Huang
% using 18ROIS from Laura's study with all subjects, FWE 0-05, Voxel 20, brain mask, 12mm.


% 2018-Jun-23 Yun-An Huang
% using 19ROIS from Laura's study with all subjects, FWE 0-05, Voxel 20.


% 2018-May-28 Yun-An Huang
% using ROIs from Laura's study

% Yun-An Huang 2018-March-01
% Calculating the graphical measures for the group analysis.


data_folder=fullfile('output_data6_gPPI_group_Laura18ROI');

contrast_num = 5;
Contrast_name{1,1} = 'emotion_vs_neutral';
Contrast_name{2,1} = 'angry_vs_neutral';
Contrast_name{3,1} = 'fear_vs_neutral';
Contrast_name{4,1} = 'sad_vs_neutral';
Contrast_name{5,1} = 'happy_vs_neutral';

random_num = 200;
graphical_measures_Healthy = {};
graphical_measures_Healthy_normalized = {};
for ctemp = 1:contrast_num;

    
    filename = ['gPPI_matrix_',Contrast_name{ctemp},'_weighted.mat']; % load the gPPI_matrix_weighted
    load(fullfile(data_folder,filename));
    
  
    gm = gPPI_weighted_graph_measures(gPPI_matrix_weighted);
    graphical_measures_Healthy.contrast{ctemp}.measures = gm;

         
    gm_nor = gPPI_weighted_graph_measures_normalize(gPPI_matrix_weighted,random_num);
    graphical_measures_Healthy_normalized.contrast{ctemp}.measures = gm_nor;
  
      
    
end

save(fullfile(data_folder,'graphical_measures_Healthy_new.mat'),'graphical_measures_Healthy','graphical_measures_Healthy_normalized');