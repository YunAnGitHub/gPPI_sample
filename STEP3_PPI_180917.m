
clear all; close all;

% 2018-Sep-18 Yun-An Huang
% using 18ROIS from Laura's study with all subjects, FWE 0-05, Voxel 20, brain mask, 12mm.


% 2018-Jun-21 Yun-An Huang
% using 19ROIS from Laura's study with all subjects, FWE 0-05, Voxel 20.

% 2017-08-18 Yun-An Huang
% change the ROI.

% 2017-06-23 Yun-An Huang
% to analyze the heathly control data
% analyze 61 ROIs only. 4 overlapped ROI were removed.

% 2017-May-31 Yun-An Huang
% to run the gPPI analysis of ATL patients
% please run the preporcessing and GLM part in advanced.


% 2017-Apr-12 Yun-An Huang
% to analyze the data of ATL patients by using GLM.
% the data should be preprocessed by the other procedure.
% after testing by GLM, we can analyze the gPPI for each subject.


% 2016-Dec-21 Yun-An Huang
% to analyze the data of ATL patients with gPPI matrix
% 


% 2016-Dec-20 Yun-An Huang
% use the preprocess method of SPM12
% slice timing, realign, coregister, segment, normalize, smooth



% 2015-Jan-08 Yun-An Huang
% remove the mask in the first level GLM analysis, the mask were conduct
% the missing result in the second level analysis in cerebellum.

% 19/11/14 Yun-An Huang
% this script is used to analysis the general PPI results of all subjects data
% which is finish the first level analysis.
% input: 
%
% output:
%
% the script is modified the PPI example provided by Guillaume Flandin & 
% Darren Gitelmanfrom Wellcome Trust Centre for neuroimaging. the original
% message is listed below:

%/
% This batch script analyses the Attention to Visual Motion fMRI dataset
% available from the SPM website using PPI:
%   http://www.fil.ion.ucl.ac.uk/spm/data/attention/
% as described in the SPM manual:
%   http://www.fil.ion.ucl.ac.uk/spm/doc/manual.pdf#Chap:data:ppi
%__________________________________________________________________________
% Copyright (C) 2014 Wellcome Trust Centre for Neuroimaging
% Guillaume Flandin & Darren Gitelman
% $Id: ppi_spm12_batch.m 11 2014-09-29 18:54:09Z guillaume $
%/



% define parameter:
% subject name
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
sbj_arr = 1:4;

% parameters
ROI_num=18;
Sess_num=4;
Contrast_num=5; % there are totally 5 contrast, including emotion-neutral, angry-neutral, happy-neutral, fear-neutral and sad-neutral
Rep_num=217; % the repetition
TR = 2;
Condition_num=8; % there are totally 8 conditions, including angry, happy, fear, sad, control, neutral, baseline, press   
nslice = 36;
smooth_kernel = [6 6 6];
%outpath
outpath = 'z_group_PPI/data_Laura18ROI';
folder_str = 'gPPI_Laura18ROI';
if exist(outpath)
    for jtemp = 1:ROI_num
        dirtemp =strcat(outpath,'/ROI_',num2str(jtemp,'%02d'));    
        if exist(dirtemp) 
            
        else
            mkdir(dirtemp);
        end
    end
else
    mkdir(outpath);
    for jtemp = 1:ROI_num
        dirtemp =strcat(outpath,'/ROI_',num2str(jtemp,'%02d'));    
        mkdir(dirtemp); 
        
    end
end

GLM_path = ['GLM_Laura18ROI'];

%%%%%%%%%%%%%%%%%%%%%%%%%% Setting the log file %%%%%%%%%%%%%%%%%%%%%%%%%
diary ([outpath,'/log ',datestr(now,'yyyy-mm-dd HH-MM'),'.txt']);
diary on
warning off 
display('Version: gPPI 2016.12.09');
display('programer: Yun-An Huang');
display(' ');
display(['start at   ',datestr(clock)]);
display(' ');
display(['the output directory: ', outpath]);
display(' ');

% load ROI list
ROIlist_filename = fullfile('ROI_Laura_wholeSB_FWE_05_v20_12mm_brainmask.txt');


% load ROI information

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

% recursive processing each sujects data.
path_top=pwd;

%initialize spm
spm('Defaults','fMRI');
spm_jobman('initcfg');

% %
% 
%% running gPPI analysis
%initialize spm
spm('Defaults','fMRI');
spm_jobman('initcfg');
for stemp = sbj_arr
% for stemp = 1:subject_num
   
    display(['gPPI: processing ',Snames{stemp},' data files....']);
   %% the preparation, including create directory and copy file 
    % create analysis folder
    ana_dir = strcat(path_top,'/',Snames{stemp},'_analysis_gPPI');
    if exist(ana_dir)
        dirtemp= strcat(ana_dir,'/functional');
        if exist(dirtemp)
            if exist(strcat(dirtemp,'/image'))
                
            else
                mkdir(strcat(dirtemp,'/image'));
                mkdir(strcat(ana_dir,'/functional/image/sess1'));
                mkdir(strcat(ana_dir,'/functional/image/sess2'));
                mkdir(strcat(ana_dir,'/functional/image/sess3'));
                mkdir(strcat(ana_dir,'/functional/image/sess4'));
                mkdir(strcat(ana_dir,'/functional/image/T1'));
            end
        else
            mkdir(dirtemp);
            mkdir(strcat(dirtemp,'/image'));
            mkdir(strcat(ana_dir,'/functional/image/sess1'));
            mkdir(strcat(ana_dir,'/functional/image/sess2'));
            mkdir(strcat(ana_dir,'/functional/image/sess3'));
            mkdir(strcat(ana_dir,'/functional/image/sess4'));
            mkdir(strcat(ana_dir,'/functional/image/T1'));
        end
        
        if exist(fullfile(ana_dir,folder_str))
        else
            mkdir(fullfile(ana_dir,folder_str));
        end
    else
        mkdir(ana_dir);
        mkdir(strcat(ana_dir,'/functional'));
        mkdir(strcat(ana_dir,'/functional/image'));
        mkdir(strcat(ana_dir,'/functional/image/sess1'));
        mkdir(strcat(ana_dir,'/functional/image/sess2'));
        mkdir(strcat(ana_dir,'/functional/image/sess3'));
        mkdir(strcat(ana_dir,'/functional/image/sess4'));
        mkdir(strcat(ana_dir,'/functional/image/T1'));
        mkdir(fullfile(ana_dir,folder_str));
        
    end
    
    cd(ana_dir);
    
    display('extracting ROI....');
    %% extract VOI   
    % VOI: EXTRACTING TIME SERIES: 
    batch_index=1;
    clear matlabbatch

    for rtemp = 1:ROI_num
        for itemp = 1:Sess_num % 4 sessions
%             if exist(fullfile(GLM_path,['VOI_ROI',num2str(rtemp,'%02d'),'_',folder_str,'_',num2str(itemp),'.mat'])); 
                % if VOI file existed then skipp ththe process, this will save a lot of time if some VOIs has beeb extracted. otherwise it will take a long time to process all the data.
                % however, it will be a problem, if you change the extract procedure below without deleting the existing file.
%                 display(strcat('VOI file: ',fullfile(GLM_path,['VOI_ROI',num2str(rtemp,'%02d'),'_',folder_str,'_',num2str(itemp),'.mat']), ' is existed. the VOI file has not been re-calculated. if you would like to process the VOI file agian, please deleted the existing file'));
%             else
                
                    %==========================================================================
                    
                    matlabbatch{batch_index}.spm.util.voi.spmmat = cellstr(fullfile(GLM_path,'SPM.mat'));
                    matlabbatch{batch_index}.spm.util.voi.adjust = 1;
                    matlabbatch{batch_index}.spm.util.voi.session = itemp;
                    matlabbatch{batch_index}.spm.util.voi.name = strcat('ROI',num2str(rtemp,'%02d'),'_',folder_str);
                    matlabbatch{batch_index}.spm.util.voi.roi{1}.spm.spmmat = {''};
                    matlabbatch{batch_index}.spm.util.voi.roi{1}.spm.contrast = 6;
                    matlabbatch{batch_index}.spm.util.voi.roi{1}.spm.threshdesc = 'none';
                    matlabbatch{batch_index}.spm.util.voi.roi{1}.spm.thresh = 0.05;
                    matlabbatch{batch_index}.spm.util.voi.roi{1}.spm.extent = 0;
                    matlabbatch{batch_index}.spm.util.voi.roi{2}.sphere.centre = peak_region(:,rtemp)';
                    matlabbatch{batch_index}.spm.util.voi.roi{2}.sphere.radius = 6;
                    matlabbatch{batch_index}.spm.util.voi.roi{2}.sphere.move.local.spm = 1;
                    matlabbatch{batch_index}.spm.util.voi.expression = 'i1&i2';
                    batch_index = batch_index+1;
%             end
        end
    end
    if batch_index<2
    else
        spm_jobman('run',matlabbatch);
    end
    cd /;
    cd(ana_dir);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % PSYCHO-PHYSIOLOGIC INTERACTION
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    
    %% create PPI variable
    display('creating PPI structure.....');
    clear matlabbatch
    batch_index=1;
    for rtemp = 1:ROI_num
        for itemp = 1:Sess_num
            for ctemp = 1:Condition_num
                
%                 if exist(fullfile(GLM_path,['PPI_ROI',num2str(rtemp,'%02d'),'_',folder_str,'_Sess',num2str(itemp),'_x_Condition',num2str(ctemp),'.mat'])); 
                % if PPI file existed then skipp ththe process, this will save a lot of time if some VOIs has beeb extracted. otherwise it will take a long time to process all the data.
                % however, it will be a problem, if you change the extract procedure below without deleting the existing file.
%                     display(strcat('PPI file: ',fullfile(GLM_path,['PPI_ROI',num2str(rtemp,'%02d'),'_',folder_str,'_Sess',num2str(itemp),'_x_Condition',num2str(ctemp),'.mat']), ' is existed. the PPI file has not been re-calculated. if you would like to process the PPI file agian, please deleted the existing file'));
%                     
%                 else
                    
                    % GENERATE PPI STRUCTURE
                    %==========================================================================
                    matlabbatch{batch_index}.spm.stats.ppi.spmmat = cellstr(fullfile(GLM_path,'SPM.mat'));
                    matlabbatch{batch_index}.spm.stats.ppi.type.ppi.voi = cellstr(fullfile(GLM_path,strcat('VOI_ROI',num2str(rtemp,'%02d'),'_',folder_str,'_',num2str(itemp),'.mat')));

                    matlabbatch{batch_index}.spm.stats.ppi.type.ppi.u = [ctemp 1 1;];
                    
                    matlabbatch{batch_index}.spm.stats.ppi.name = strcat('ROI',num2str(rtemp,'%02d'),'_',folder_str,'_Sess',num2str(itemp),'_x_Condition',num2str(ctemp));
                    matlabbatch{batch_index}.spm.stats.ppi.disp = 0;

                    batch_index=batch_index+1;
%                 end
                
            end
        end
    end
    if batch_index<2 % to check whether this process have been done or not.
    else
        spm_jobman('run',matlabbatch);
    end
    cd /;
    cd(ana_dir);

    %% PPI GLM analysis
    display('PPI GLM analysis.....');
  
    for rtemp = 1:ROI_num
%         for ctemp = 1:Contrast_num
            clear matlabbatch
            % MODEL SPECIFICATION
            %==========================================================================

            % Directory
            %--------------------------------------------------------------------------
            PPI_dir = fullfile(ana_dir,folder_str,strcat('ROI',num2str(rtemp,'%02d')));
            if exist(fullfile(PPI_dir))
            else
                mkdir(PPI_dir);
            end
            
            if exist(fullfile(PPI_dir,'SPM.mat')) % run PPI again.
                delete(fullfile(PPI_dir,'SPM.mat'));
            end
            
            matlabbatch{1}.spm.stats.fmri_spec.dir = cellstr(PPI_dir);

            % Timing
            %--------------------------------------------------------------------------
            matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'secs';
            matlabbatch{1}.spm.stats.fmri_spec.timing.RT = TR;

            % mask 
            matlabbatch{1}.spm.stats.fmri_spec.mthresh=0;
            
            % Session
            %--------------------------------------------------------------------------
            
            for itemp = 1:Sess_num
                
                % load img file
                f = spm_select('ExtFPList', fullfile(ana_dir,'functional','image',['sess',num2str(itemp)]), '^swra',1:1000);
                matlabbatch{1}.spm.stats.fmri_spec.sess(itemp).scans = cellstr(f);
                
                R=zeros(size(f,1),Condition_num*2+1+6);% Condition_num*PPI, PPI.Y, Condition_num*PPI.P, 6 regressor
                for ctemp=1:Condition_num
                    % Regressors
                    % build regressor ( interaction, physiology,psychology,
                    % R1-R6), totally Condition_num*2+1+6 regressors in gPPI.

                    clear PPI; % to remove the previous PPI data
                    load(fullfile(ana_dir,GLM_path,['PPI_ROI',num2str(rtemp,'%02d'),'_',folder_str,'_Sess',num2str(itemp),'_x_Condition',num2str(ctemp),'.mat']));
                    clear SPM; % to remove the previous SPM data
                    load(fullfile(ana_dir,GLM_path,'SPM.mat'));
                    
                    R(:,ctemp)=PPI.ppi;
                    R(:,Condition_num+1)=PPI.Y;
                    R(:,Condition_num+1+ctemp)=PPI.P;
                end     
                R(:,(Condition_num*2+1+1):(Condition_num*2+1+6)) = [SPM.Sess(itemp).C.C ]; %PPI.ppi: interaction, PPI.Y: physiology, PPI.P: psychology.
                save(fullfile(PPI_dir,['regressor_ROI',num2str(rtemp,'%02d'),'_',folder_str,'_Sess',num2str(itemp),'.mat']),'R');


                matlabbatch{1}.spm.stats.fmri_spec.sess(itemp).multi_reg = {...
                fullfile(PPI_dir,['regressor_ROI',num2str(rtemp,'%02d'),'_',folder_str,'_Sess',num2str(itemp),'.mat'])};
            
                
            end
            
            % MODEL ESTIMATION
            %==========================================================================
            matlabbatch{2}.spm.stats.fmri_est.spmmat = cellstr(fullfile(PPI_dir,'SPM.mat'));

            % INFERENCE
            %==========================================================================
            matlabbatch{3}.spm.stats.con.spmmat = cellstr(fullfile(PPI_dir,'SPM.mat'));
            
            for ctemp=1:Contrast_num
                matlabbatch{3}.spm.stats.con.consess{ctemp}.tcon.name = ['PPI-Interaction-Contrast',num2str(ctemp)]; 
                % contrast 1: emotion-neutral, contrast 2: angry-neutral,
                % contrast 3: happy-neutral, contrast 4: fear-neutral.
                % contrast 5: sad-neutral
                
                con_array = zeros(1,(Condition_num*2+1+6)*Sess_num+Sess_num);
                switch ctemp
                
                    case 1 % emotional-neutral
                
                    
                        con_array(1:(Condition_num*2+1+6):((Condition_num*2+1+6)*Sess_num))=1;
                        con_array(2:(Condition_num*2+1+6):((Condition_num*2+1+6)*Sess_num))=1;
                        con_array(3:(Condition_num*2+1+6):((Condition_num*2+1+6)*Sess_num))=1;
                        con_array(4:(Condition_num*2+1+6):((Condition_num*2+1+6)*Sess_num))=1;
                        con_array(6:(Condition_num*2+1+6):((Condition_num*2+1+6)*Sess_num))=-4;
                        
                    case 2 % angry-neutral
                        con_array(1:(Condition_num*2+1+6):((Condition_num*2+1+6)*Sess_num))=1;
                        con_array(6:(Condition_num*2+1+6):((Condition_num*2+1+6)*Sess_num))=-1;
                        
                
                    case 3 % fear-neutral
                        con_array(2:(Condition_num*2+1+6):((Condition_num*2+1+6)*Sess_num))=1;
                        con_array(6:(Condition_num*2+1+6):((Condition_num*2+1+6)*Sess_num))=-1;
                    
                    case 4 % sad-neutral
                        con_array(3:(Condition_num*2+1+6):((Condition_num*2+1+6)*Sess_num))=1;
                        con_array(6:(Condition_num*2+1+6):((Condition_num*2+1+6)*Sess_num))=-1;
                        
                    case 5 % happy-neutral
                        con_array(4:(Condition_num*2+1+6):((Condition_num*2+1+6)*Sess_num))=1;
                        con_array(6:(Condition_num*2+1+6):((Condition_num*2+1+6)*Sess_num))=-1;
                        
                       

                end
                matlabbatch{3}.spm.stats.con.consess{ctemp}.tcon.weights = con_array;
    
            end
            
            % RESULTS
            %==========================================================================
%             matlabbatch{4}.spm.stats.results.spmmat = cellstr(fullfile(PPI_dir,'SPM.mat'));
%             matlabbatch{4}.spm.stats.results.conspec.contrasts = 1;
%             matlabbatch{4}.spm.stats.results.conspec.threshdesc = 'none';
%             matlabbatch{4}.spm.stats.results.conspec.thresh = 0.01;
%             matlabbatch{4}.spm.stats.results.conspec.extent = 3;
%             matlabbatch{4}.spm.stats.results.print = false;

            spm_jobman('run',matlabbatch);

            cd /;
            cd(ana_dir); 
            
            copyfile(fullfile(PPI_dir,'spmT_0001.nii'),fullfile(path_top,outpath,['ROI_',num2str(rtemp,'%02d')],['spmT_001_',Snames{stemp},'.nii']));
            copyfile(fullfile(PPI_dir,'con_0001.nii'),fullfile(path_top,outpath,['ROI_',num2str(rtemp,'%02d')],['con_001_',Snames{stemp},'.nii']));
            copyfile(fullfile(PPI_dir,'spmT_0002.nii'),fullfile(path_top,outpath,['ROI_',num2str(rtemp,'%02d')],['spmT_002_',Snames{stemp},'.nii']));
            copyfile(fullfile(PPI_dir,'con_0002.nii'),fullfile(path_top,outpath,['ROI_',num2str(rtemp,'%02d')],['con_002_',Snames{stemp},'.nii']));
            copyfile(fullfile(PPI_dir,'spmT_0003.nii'),fullfile(path_top,outpath,['ROI_',num2str(rtemp,'%02d')],['spmT_003_',Snames{stemp},'.nii']));
            copyfile(fullfile(PPI_dir,'con_0003.nii'),fullfile(path_top,outpath,['ROI_',num2str(rtemp,'%02d')],['con_003_',Snames{stemp},'.nii']));
            copyfile(fullfile(PPI_dir,'spmT_0004.nii'),fullfile(path_top,outpath,['ROI_',num2str(rtemp,'%02d')],['spmT_004_',Snames{stemp},'.nii']));
            copyfile(fullfile(PPI_dir,'con_0004.nii'),fullfile(path_top,outpath,['ROI_',num2str(rtemp,'%02d')],['con_004_',Snames{stemp},'.nii']));
            copyfile(fullfile(PPI_dir,'spmT_0005.nii'),fullfile(path_top,outpath,['ROI_',num2str(rtemp,'%02d')],['spmT_005_',Snames{stemp},'.nii']));
            copyfile(fullfile(PPI_dir,'con_0005.nii'),fullfile(path_top,outpath,['ROI_',num2str(rtemp,'%02d')],['con_005_',Snames{stemp},'.nii']));
        
%         end
    end
    cd /;
    cd(path_top);
end
cd /;
cd(path_top);

display(['end at   ',datestr(clock)]);
diary off;

