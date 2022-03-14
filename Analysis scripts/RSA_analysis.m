% Build RSA matrices

clear all
sublist = [1:11 13 15:18 20 22:31 34:36 38];

ROIs = ['rapc';'lapc';'rppc';'lppc';'rofc';'lofc';'ramg';'lamg';'rhip'; 'lhip'];

sessions = {'PrecEPI1'; 'PostEPI3'; 'P2EPI4'};

homed = pwd;

for roiLoop = 1:5
    
    curROI = ROIs(2*roiLoop-1,:);
    curROI2 = ROIs(2*roiLoop,:);
    
    eval(['RSA_', curROI, '_AllAveCorrMatrix = [];'])
    eval(['RSA_', curROI, '_AllAveCorrMatrixZ = [];'])
    
    
    for sub = sublist
            
        curSub = num2str(sub);      
        
        % loop through all three sessions
        for sessionLoop = 1:3
            
            curSession = sessions{sessionLoop,:};
                        
            % loop through each condition
            for myCond = [1 2 3 4 5]
                
                curCond = num2str(myCond);
                
                betaM = [];
                eval(['load D:\lilab\OCF_Analyses\SVM_Processing_no_cap_new32\sub',curSub,'\beta_matrices_LSA\sub',curSub,'_',curSession,'_',curROI,'_cond',curCond,'.mat betaMatrix']);
                betaM = [betaM betaMatrix];
                
                eval(['load D:\lilab\OCF_Analyses\SVM_Processing_no_cap_new32\sub',curSub,'\beta_matrices_LSA\sub',curSub,'_',curSession,'_',curROI2,'_cond',curCond,'.mat betaMatrix']);
                betaM = [betaM betaMatrix];
                
                switch myCond
                    case 1
                        
                        eval(['PV1_' num2str(sessionLoop) ' = nanmean(betaM,1);']);
                    case 2
                        
                        eval(['PV2_' num2str(sessionLoop) ' = nanmean(betaM,1);']);
                    case 3
                        
                        eval(['PV3_' num2str(sessionLoop) ' = nanmean(betaM,1);']);
                    case 4
                        
                        eval(['PV4_' num2str(sessionLoop) ' = nanmean(betaM,1);']);
                    case 5
                        
                        eval(['PV5_' num2str(sessionLoop) ' = nanmean(betaM,1);']);
                end                
                
            end
                        
        end
        
        
        %RSA correlation matrix arranged as Pre CSsafe -> CSthreat, Post
        %CSsafe -> CSthreat, Post2 CSsafe -> CSthreat
        switch rem(sub,2)
            case(0)
                IVs = [PV5_1' PV4_1' PV3_1' PV2_1' PV1_1' PV5_2' PV4_2'  PV3_2' PV2_2' PV1_2' PV5_3' PV4_3' PV3_3' PV2_3' PV1_3'];
            case(1)
                IVs = [PV1_1' PV2_1' PV3_1' PV4_1' PV5_1' PV1_2' PV2_2' PV3_2' PV4_2' PV5_2' PV1_3' PV2_3' PV3_3' PV4_3' PV5_3'];
                
        end
        
        corr_matrix = corrcoef(IVs);
        corr_matrix_z = .5.*log((1+corr_matrix)./(1-corr_matrix));
        dissimilarity = 1 - corr_matrix;
        dissimilarity_z = -1.*corr_matrix_z;
        
    eval(['RSA_', curROI, '_AllAveCorrMatrix = cat(3,RSA_', curROI, '_AllAveCorrMatrix, corr_matrix);'])
    eval(['RSA_', curROI, '_AllAveCorrMatrixZ = cat(3,RSA_', curROI, '_AllAveCorrMatrixZ, corr_matrix_z);'])
        
        
    end
    
    eval(['RSA_', curROI, '_AllAveRDM = 1 - mean(RSA_', curROI,'_AllAveCorrMatrix,3);']);
    
    eval(['RSA_', curROI, '_PostPre_RDM = RSA_', curROI, '_AllAveRDM(6:10,6:10) - RSA_', curROI, '_AllAveRDM(1:5,1:5) ']);
    
    eval(['RSA_', curROI, '_Post2Pre_RDM = RSA_', curROI, '_AllAveRDM(11:15,11:15) - RSA_', curROI, '_AllAveRDM(1:5,1:5) ']);
    
    eval(['save D:\lilab\OCF_Analyses\RSA_l', curROI, '_31SubAveCorrMatrix_CSm2CSp.mat RSA_', curROI, '_AllAveCorrMatrix RSA_', curROI, '_AllAveCorrMatrixZ RSA_', curROI, '_AllAveRDM RSA_', curROI, '_PostPre_RDM RSA_', curROI, '_Post2Pre_RDM']);
    
end
%% sort individual distance data

ROIs = ['apc';'ppc';'ofc';'amg';'hip'];
PDIIndex = zeros(31,15);

for roiLoop =  1:size(ROIs,1)
    curROI = ROIs(roiLoop,:);
    eval(['load D:\lilab\OCF_Analyses\RSA_lr', curROI, '_31SubAveCorrMatrix_CSm2CSp.mat']);
    
    eval(['Pre_', curROI,'_All = zeros(31,10);']);
    eval(['Post_', curROI, '_All = zeros(31,10);']);
    eval(['Post2_', curROI, '_All = zeros(31,10);']);

    for s = 1:31
        eval([' dissimilarity_z = -1*RSA_r' ,curROI, '_AllAveCorrMatrixZ(:,:,s);']);
        eval([' corr_matrix_z = RSA_r' ,curROI, '_AllAveCorrMatrixZ(:,:,s);']);
        
        %%Calculate pattern differentation index (PDI)
        PDIPre = dissimilarity_z(1,2)-dissimilarity_z(2,3)-dissimilarity_z(3,4)+dissimilarity_z(4,5);
        PDIPost = dissimilarity_z(6,7)-dissimilarity_z(7,8)-dissimilarity_z(8,9)+dissimilarity_z(9,10);
        PDIPost2 = dissimilarity_z(11,12)-dissimilarity_z(12,13)-dissimilarity_z(13,14)+dissimilarity_z(14,15);
       
        PDIIndex(s,1+(roiLoop-1)*3) = PDIPre;
        PDIIndex(s,2+(roiLoop-1)*3) = PDIPost;
        PDIIndex(s,3+(roiLoop-1)*3) = PDIPost2;
   
    end
end

%% Plot RDM

load \\psy-wlhost1\lilab\OCF_Analyses\MVPA_Pipeline_scripts\APC_rdm_c1_map.mat
ROIs = ['apc';'ppc';'ofc';'amg';'hip'];

for roiLoop = 1:size(ROIs,1)
    curROI = ROIs(roiLoop,:);
    eval(['load \\psy-wlhost1\lilab\OCF_Analyses\RSA_lr', curROI, '_31SubAveCorrMatrix_CSm2CSp.mat']);

    figure;
    eval(['imagesc(RSA_r' curROI '_AllAveRDM, [0.2 0.5])']);

    colormap(c1)
    colorbar
    eval(['title(''RSA ', curROI, ' RDM PrePostPost2 CSmGS3GS2GS1CSp'')']);
    eval(['saveas(gcf,''\\psy-wlhost1\lilab\OCF_Analyses\MVPA_Pipeline_scripts\RSA_RDM_CSm2CSp_', curROI, '_31subs.eps'');']);

end
