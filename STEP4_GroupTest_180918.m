clear all;close all;

% 2018-Sep-18 Yun-An Huang
% using 18ROIS from Laura's study with all subjects, FWE 0-05, Voxel 20, brain mask, 12mm.


% 2018-Jun-23 Yun-An Huang
% using 19ROIS from Laura's study with all subjects, FWE 0-05, Voxel 20.



% 2018-May-28 Yun-An Huang
% for 23 ROIs from laura's study

% 2017-Aug-25 Yun-An Huang
% for sym 62 ROIs

% 2017-Jul-31 Yun-An Huang
% correct the contrast name

% 2017-Jun-26 Yun-An Huang
% applied the gPPI analysis to the healthy group of temporal lobectomy
% study.


% 2015-May-27 Yun-An Huang
% add the one way anova test of the contrast 2-5


% 2014/11/24 Yun-an Huang
% this script is used to examine the PPI results on the second level.
% including the one sample t-test:
% 1. emo vs neutral
% 2. angry vs neutral
% 3. happy vs neutral
% 4. fear vs neutral
% 5. sad vs neutral
% 
% also, perform the paired t test:
%
% 1. angry vs happy 
% 2. angry vs fear 
% 3. angry vs sad 
% 4. happy vs fear
% 5. happy vs sad
% 6. fear  vs sad
% 
% also, perform the ANOVA to examine:
% if there are significance difference between: angry, happy, fear and sad.
%

 
%% define parameter

ROI_num=18;
group_test_num = 12; % 2015-May-27
path_top = pwd;
data_folder = 'data_Laura18ROI';

subject_num = 12;
Snames{1} = 'EriClae_09022015';
Snames{2} = 'FraBar_10122015';
Snames{3} = 'GeeHen_08052015';
Snames{4} = 'GodVan_04052015';
Snames{5} = 'JohCha';
Snames{6} = 'MarDer_13052015';
Snames{7} = 'MarHil_13072015';
Snames{8} = 'MyrWel';
Snames{9} = 'PatVan_09072015';
Snames{10} = 'TinVan_12022015';
Snames{11} = 'VANMAR_27042015';
Snames{12} = 'VicDeb_16072015';

Contrast_name{1,1} = 'emotion_vs_neutral';
Contrast_name{2,1} = 'angry_vs_neutral';
Contrast_name{3,1} = 'fear_vs_neutral';
Contrast_name{4,1} = 'sad_vs_neutral';
Contrast_name{5,1} = 'happy_vs_neutral';
Contrast_name{6,1} = 'angry_vs_fear';
Contrast_name{6,2} = 'fear_vs_angry';
Contrast_name{7,1} = 'angry_vs_sad';
Contrast_name{7,2} = 'sad_vs_angry';
Contrast_name{8,1} = 'angry_vs_happy';
Contrast_name{8,2} = 'happy_vs_angry';
Contrast_name{9,1} = 'fear_vs_sad';
Contrast_name{9,2} = 'sad_vs_fear';
Contrast_name{10,1} = 'fear_vs_happy';
Contrast_name{10,2} = 'happy_vs_fear';
Contrast_name{11,1} = 'sad_vs_happy';
Contrast_name{11,2} = 'happy_vs_sad';
Contrast_name{12,1} = 'one_way_anova_f';
Contrast_name{12,2} = 'one_way_anova_Angry_t';
Contrast_name{12,3} = 'one_way_anova_Fear_t';
Contrast_name{12,4} = 'one_way_anova_Sad_t';
Contrast_name{12,5} = 'one_way_anova_Happy_t';


%% create output directory
outfolder = 'output_data2_gPPI_Laura18ROI';
if exist(outfolder,'dir')
    for jtemp = 1:ROI_num
        dirtemp =strcat(outfolder,'/ROI_',num2str(jtemp,'%02d'));    
        if exist(dirtemp) 

        else
            mkdir(dirtemp);
        end
    end
else
    mkdir(outfolder);
    for jtemp = 1:ROI_num
        dirtemp =strcat(outfolder,'/ROI_',num2str(jtemp,'%02d'));    
        mkdir(dirtemp);
    end

end



%%%%%%%%%%%%%%%%%%%%%%%%%% Setting the log file %%%%%%%%%%%%%%%%%%%%%%%%%
diary ([outfolder,'/log ',datestr(now,'yyyy-mm-dd HH-MM'),'.txt']);
diary on
warning off 
display('Version: gPPI group test 2014.11.24');
display('programer: Yun-An Huang');
display(' ');
display(['start at   ',datestr(clock)]);
display(' ');
display(['the output directory: ', outfolder]);
display(' ');


%% computing second level analysis
% initialize spm
spm('Defaults','fMRI');
spm_jobman('initcfg');


for rtemp = 1:ROI_num
    
    cd(fullfile(outfolder,['ROI_',num2str(rtemp,'%02d')]));

    display(['Process ROI',num2str(rtemp)]);
    
    for gtemp = 1: group_test_num
        % create directory of each contrast     
        if exist('results')
        else
            mkdir('results')
        end
    
        group_test_dir = fullfile(path_top,outfolder,['ROI_',num2str(rtemp,'%02d')],['group_test_',num2str(gtemp,'%02d')]);
        if exist(group_test_dir)
        else
            mkdir(group_test_dir);
        end
    
        if exist(fullfile(group_test_dir,'SPM.mat'))
            delete(fullfile(group_test_dir,'SPM.mat'));
        end
        clear matlabbatch

        matlabbatch{1}.spm.stats.factorial_design.dir = cellstr(group_test_dir);
        matlabbatch{1}.spm.stats.factorial_design.masking.im = 0;
        matlabbatch{1}.spm.stats.factorial_design.masking.im = 0;
        matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};


        %% build up 2nd level model
        switch gtemp
            case 1 % one sample t: emotional - neutral
                f = spm_select('ExtFPList', fullfile(path_top,data_folder,['ROI_',num2str(rtemp,'%02d')]), strcat('con_00',num2str(gtemp)),1);
                matlabbatch{1}.spm.stats.factorial_design.des.t1.scans =  cellstr(f);

            case 2 % one sample t: angry - neutral
                f = spm_select('ExtFPList', fullfile(path_top,data_folder,['ROI_',num2str(rtemp,'%02d')]), strcat('con_00',num2str(gtemp)),1);
                matlabbatch{1}.spm.stats.factorial_design.des.t1.scans =  cellstr(f);

            case 3 % one sample t: fear - neutral
                f = spm_select('ExtFPList', fullfile(path_top,data_folder,['ROI_',num2str(rtemp,'%02d')]), strcat('con_00',num2str(gtemp)),1);
                matlabbatch{1}.spm.stats.factorial_design.des.t1.scans =  cellstr(f);

            case 4 % one sample t: sad - neutral
                f = spm_select('ExtFPList', fullfile(path_top,data_folder,['ROI_',num2str(rtemp,'%02d')]), strcat('con_00',num2str(gtemp)),1);
                matlabbatch{1}.spm.stats.factorial_design.des.t1.scans =  cellstr(f);

            case 5 % one sample t: happy - neutral
                f = spm_select('ExtFPList', fullfile(path_top,data_folder,['ROI_',num2str(rtemp,'%02d')]), strcat('con_00',num2str(gtemp)),1);
                matlabbatch{1}.spm.stats.factorial_design.des.t1.scans =  cellstr(f);

            case 6 % two sample t: angry-fear
                for stemp = 1:subject_num
                    f1 = spm_select('ExtFPList', fullfile(path_top,data_folder,['ROI_',num2str(rtemp,'%02d')]), strcat('con_00',num2str(2),'_',Snames(stemp)),1);
                    f2 = spm_select('ExtFPList', fullfile(path_top,data_folder,['ROI_',num2str(rtemp,'%02d')]), strcat('con_00',num2str(3),'_',Snames(stemp)),1);
                    matlabbatch{1}.spm.stats.factorial_design.des.pt.pair(stemp).scans =  cellstr(strvcat(f1,f2));
                end
                
            case 7 % two sample t: angry-sad
                for stemp = 1:subject_num
                    f1 = spm_select('ExtFPList', fullfile(path_top,data_folder,['ROI_',num2str(rtemp,'%02d')]), strcat('con_00',num2str(2),'_',Snames(stemp)),1);
                    f2 = spm_select('ExtFPList', fullfile(path_top,data_folder,['ROI_',num2str(rtemp,'%02d')]), strcat('con_00',num2str(4),'_',Snames(stemp)),1);
                    matlabbatch{1}.spm.stats.factorial_design.des.pt.pair(stemp).scans =  cellstr(strvcat(f1,f2));
                end
                
            case 8 % two sample t: angry-happy
                for stemp = 1:subject_num
                    f1 = spm_select('ExtFPList', fullfile(path_top,data_folder,['ROI_',num2str(rtemp,'%02d')]), strcat('con_00',num2str(2),'_',Snames(stemp)),1);
                    f2 = spm_select('ExtFPList', fullfile(path_top,data_folder,['ROI_',num2str(rtemp,'%02d')]), strcat('con_00',num2str(5),'_',Snames(stemp)),1);
                    matlabbatch{1}.spm.stats.factorial_design.des.pt.pair(stemp).scans =  cellstr(strvcat(f1,f2));
                end
                
            case 9 % two sample t: fear-sad
                for stemp = 1:subject_num
                    f1 = spm_select('ExtFPList', fullfile(path_top,data_folder,['ROI_',num2str(rtemp,'%02d')]), strcat('con_00',num2str(3),'_',Snames(stemp)),1);
                    f2 = spm_select('ExtFPList', fullfile(path_top,data_folder,['ROI_',num2str(rtemp,'%02d')]), strcat('con_00',num2str(4),'_',Snames(stemp)),1);
                    matlabbatch{1}.spm.stats.factorial_design.des.pt.pair(stemp).scans =  cellstr(strvcat(f1,f2));
                end
                
            case 10 % two sample t: fear-happy
                for stemp = 1:subject_num
                    f1 = spm_select('ExtFPList', fullfile(path_top,data_folder,['ROI_',num2str(rtemp,'%02d')]), strcat('con_00',num2str(3),'_',Snames(stemp)),1);
                    f2 = spm_select('ExtFPList', fullfile(path_top,data_folder,['ROI_',num2str(rtemp,'%02d')]), strcat('con_00',num2str(5),'_',Snames(stemp)),1);
                    matlabbatch{1}.spm.stats.factorial_design.des.pt.pair(stemp).scans =  cellstr(strvcat(f1,f2));
                end
            case 11 % two sample t: sad-happy
                for stemp = 1:subject_num
                    f1 = spm_select('ExtFPList', fullfile(path_top,data_folder,['ROI_',num2str(rtemp,'%02d')]), strcat('con_00',num2str(4),'_',Snames(stemp)),1);
                    f2 = spm_select('ExtFPList', fullfile(path_top,data_folder,['ROI_',num2str(rtemp,'%02d')]), strcat('con_00',num2str(5),'_',Snames(stemp)),1);
                    matlabbatch{1}.spm.stats.factorial_design.des.pt.pair(stemp).scans =  cellstr(strvcat(f1,f2));
                end
    
            case 12 % one way anova test
                for stemp = 1:subject_num
                
                    f = spm_select('ExtFPList', fullfile(path_top,data_folder,['ROI_',num2str(rtemp,'%02d')]), strcat('con_00[',num2str([2:5]),']_',Snames(stemp)),1); % select contrast 2:5
                    matlabbatch{1}.spm.stats.factorial_design.des.anovaw.fsubject(stemp).scans = cellstr(f);
                    matlabbatch{1}.spm.stats.factorial_design.des.anovaw.fsubject(stemp).conds = [2:5];
                end                
                
        end
        
        %% model estimate
        matlabbatch{2}.spm.stats.fmri_est.spmmat = cellstr(fullfile(group_test_dir,'SPM.mat'));
    
        %% inference and contrast define
        matlabbatch{3}.spm.stats.con.spmmat = cellstr(fullfile(group_test_dir,'SPM.mat'));
        if gtemp<=5 % one sample t
            matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = Contrast_name{gtemp,1};
            matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = [1];
        elseif gtemp <=11 % pair t
            matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = Contrast_name{gtemp,1};
            matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = [1 -1];
            matlabbatch{3}.spm.stats.con.consess{2}.tcon.name = Contrast_name{gtemp,2};
            matlabbatch{3}.spm.stats.con.consess{2}.tcon.weights = [-1 1];
        else % one way anova test
            matlabbatch{3}.spm.stats.con.consess{1}.fcon.name = Contrast_name{gtemp,1};
            matlabbatch{3}.spm.stats.con.consess{1}.fcon.weights = [ones(4,4)*(-1/3)+eye(4,4)*(1+1/3) zeros(4,subject_num)];
            matlabbatch{3}.spm.stats.con.consess{2}.tcon.name = Contrast_name{gtemp,2};
            matlabbatch{3}.spm.stats.con.consess{2}.tcon.weights = [1 -1/3 -1/3 -1/3 zeros(1,subject_num)];
            matlabbatch{3}.spm.stats.con.consess{3}.tcon.name = Contrast_name{gtemp,3};
            matlabbatch{3}.spm.stats.con.consess{3}.tcon.weights = [-1/3 1 -1/3 -1/3 zeros(1,subject_num)];
            matlabbatch{3}.spm.stats.con.consess{4}.tcon.name = Contrast_name{gtemp,4};
            matlabbatch{3}.spm.stats.con.consess{4}.tcon.weights = [-1/3 -1/3 1 -1/3 zeros(1,subject_num)];
            matlabbatch{3}.spm.stats.con.consess{5}.tcon.name = Contrast_name{gtemp,5};
            matlabbatch{3}.spm.stats.con.consess{5}.tcon.weights = [-1/3 -1/3 -1/3 1 zeros(1,subject_num)];
            
            
        end
        
%          %% results
%         matlabbatch{4}.spm.stats.results.spmmat = cellstr(fullfile(group_test_dir,'SPM.mat'));
%         matlabbatch{4}.spm.stats.results.conspec.contrasts = 1;
%         matlabbatch{4}.spm.stats.results.conspec.threshdesc = 'none';
%         matlabbatch{4}.spm.stats.results.conspec.thresh = 0.01;
%         matlabbatch{4}.spm.stats.results.conspec.extent = 3;
%         matlabbatch{4}.spm.stats.results.print = false;
        
        spm_jobman('run',matlabbatch);
        
        if gtemp <=5
            copyfile(fullfile(group_test_dir,'spmT_0001.nii'),fullfile(path_top,outfolder,['ROI_',num2str(rtemp,'%02d')],'results',['spmT_',Contrast_name{gtemp,1},'.nii']));
            copyfile(fullfile(group_test_dir,'con_0001.nii'),fullfile(path_top,outfolder,['ROI_',num2str(rtemp,'%02d')],'results',['con_',Contrast_name{gtemp,1},'.nii']));
        elseif gtemp <=11
            copyfile(fullfile(group_test_dir,'spmT_0001.nii'),fullfile(path_top,outfolder,['ROI_',num2str(rtemp,'%02d')],'results',['spmT_',Contrast_name{gtemp,1},'.nii']));
            copyfile(fullfile(group_test_dir,'spmT_0002.nii'),fullfile(path_top,outfolder,['ROI_',num2str(rtemp,'%02d')],'results',['spmT_',Contrast_name{gtemp,2},'.nii']));
            
            copyfile(fullfile(group_test_dir,'con_0001.nii'),fullfile(path_top,outfolder,['ROI_',num2str(rtemp,'%02d')],'results',['con_',Contrast_name{gtemp,1},'.nii']));
            copyfile(fullfile(group_test_dir,'con_0002.nii'),fullfile(path_top,outfolder,['ROI_',num2str(rtemp,'%02d')],'results',['con_',Contrast_name{gtemp,2},'.nii']));
        else % one way anova test
            % 2015-May-27
            copyfile(fullfile(group_test_dir,'spmF_0001.nii'),fullfile(path_top,outfolder,['ROI_',num2str(rtemp,'%02d')],'results',['spmF_',Contrast_name{gtemp,1},'.nii']));
            copyfile(fullfile(group_test_dir,'ess_0001.nii'),fullfile(path_top,outfolder,['ROI_',num2str(rtemp,'%02d')],'results',['ess_',Contrast_name{gtemp,1},'.nii']));
            
            copyfile(fullfile(group_test_dir,'spmT_0002.nii'),fullfile(path_top,outfolder,['ROI_',num2str(rtemp,'%02d')],'results',['spmT_',Contrast_name{gtemp,2},'.nii']));
            copyfile(fullfile(group_test_dir,'con_0002.nii'),fullfile(path_top,outfolder,['ROI_',num2str(rtemp,'%02d')],'results',['con_',Contrast_name{gtemp,2},'.nii']));
            
            copyfile(fullfile(group_test_dir,'spmT_0003.nii'),fullfile(path_top,outfolder,['ROI_',num2str(rtemp,'%02d')],'results',['spmT_',Contrast_name{gtemp,3},'.nii']));
            copyfile(fullfile(group_test_dir,'con_0003.nii'),fullfile(path_top,outfolder,['ROI_',num2str(rtemp,'%02d')],'results',['con_',Contrast_name{gtemp,3},'.nii']));
            
            copyfile(fullfile(group_test_dir,'spmT_0004.nii'),fullfile(path_top,outfolder,['ROI_',num2str(rtemp,'%02d')],'results',['spmT_',Contrast_name{gtemp,4},'.nii']));
            copyfile(fullfile(group_test_dir,'con_0004.nii'),fullfile(path_top,outfolder,['ROI_',num2str(rtemp,'%02d')],'results',['con_',Contrast_name{gtemp,4},'.nii']));
            
            copyfile(fullfile(group_test_dir,'spmT_0005.nii'),fullfile(path_top,outfolder,['ROI_',num2str(rtemp,'%02d')],'results',['spmT_',Contrast_name{gtemp,5},'.nii']));
            copyfile(fullfile(group_test_dir,'con_0005.nii'),fullfile(path_top,outfolder,['ROI_',num2str(rtemp,'%02d')],'results',['con_',Contrast_name{gtemp,5},'.nii']));
            
            
            
        end
    end
    cd '/'
    cd(path_top);
end

display(['end at   ',datestr(clock)]);
diary off;