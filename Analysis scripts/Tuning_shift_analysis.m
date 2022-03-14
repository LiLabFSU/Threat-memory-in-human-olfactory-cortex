%% Calculate Mutual Information for each voxel based on
%% ALL trials

clear all
sublist = [1:11 13 15:18 20 22:31 34:36 38];%31subs

ROIs = ['rapc';'lapc';'rppc';'lppc';'rofc';'lofc';'ramg';'lamg';'rhip'; 'lhip'];
ROIlist = ['apc';'ppc';'ofc';'amg';'hip'];
sessions = {'PrecEPI1'; 'PostEPI3'; 'P2EPI4'};
conds = [1 2 3 4 5];
k = 5;

for r = 1:5
    region = ROIlist(r,:);
    eval(['MI_',region,' = cell(31,3);']); %mutual information
end

%bin sizes pre-calculated for each sub/condition
allsubbins = [9	9	9	9	9
    9	9	9	9	9
    8	8	8	8	8
    8	8	8	8	8
    9	8	8	8	8
    10	11	11	10	10
    9	9	8	9	8
    9	9	9	9	9
    8	8	9	8	8
    8	8	9	8	8
    8	9	10	9	10
    8	8	8	8	8
    9	9	8	8	8
    8	8	9	8	8
    9	8	9	8	8
    8	9	8	8	8
    11	11	11	12	11
    9	8	9	8	8
    8	8	8	8	8
    8	8	8	9	9
    8	9	8	8	8
    8	8	9	8	8
    9	8	8	8	8
    9	9	9	9	9
    9	9	10	9	10
    8	9	8	8	8
    15	15	15	14	14
    8	8	8	8	8
    9	8	8	8	8
    9	9	10	10	10
    9	9	9	9	9
    
    ];

subcount = 0;
for sub = sublist
    subcount = subcount+1;
    
    for roi =  1:2:9 %collapsing left and right
        roiname = ROIs(roi,2:end);
        
        for session = 1:size(sessions,1)
            
            % load in each matrix for each condition and assign name
            
            for cond = [1 2 3 4 5]
                
                betaM = [];
                
                eval(['load  \\psy-wlhost1\lilab\OCF_Analyses\SVM_Processing_no_cap_new32\sub',num2str(sub),'\beta_matrices_LSA\sub',num2str(sub),'_',sessions{session,:},'_',ROIs(roi,:),'_cond',num2str(cond),'.mat;']);
                betaM = [betaM betaMatrix];
                
                eval(['load  \\psy-wlhost1\lilab\OCF_Analyses\SVM_Processing_no_cap_new32\sub',num2str(sub),'\beta_matrices_LSA\sub',num2str(sub),'_',sessions{session,:},'_',ROIs(roi+1,:),'_cond',num2str(cond),'.mat;']);
                betaM = [betaM betaMatrix];
                
                % remove the mean beta value from each trial
                eval(['cond',num2str(cond),' = betaM - repmat(mean(betaM'')'',1,size(betaM,2));']);
            end
            
            if sub==28 && session==2 %fill in the missing trials with zeros (since it will be averaged anyway)
                %  k = 14;
                cond1 = [cond1; nan(1,size(cond1,2))];
                cond2 = [cond2; nan(1,size(cond2,2))];
                cond3 = [cond3; nan(1,size(cond3,2))];
            end
            
            %             %for each voxel, normalize across all trials (i.e. across
            %             %conditions)
            allcondBeta = [cond1; cond2; cond3; cond4; cond5];
            allcondnorm = zeros(size(allcondBeta,1), size(allcondBeta,2));
            
            %          normalize each voxel (z-score)
            for voxel = 1:size(cond1,2)
                normVox = (allcondBeta(:,voxel) - repmat(nanmean(allcondBeta(:,voxel)),size(allcondBeta,1),1))/nanstd(allcondBeta(:,voxel));
                allcondnorm(:,voxel) = normVox;
            end
            
            %vectorize allcondnorm (75 trials X N voxels)%Serences et al 2009: the
            %bins were defined based on the range of responses across all
            %voxels in a given ROI in a scanning session
            allcondv = reshape(allcondnorm, 1, size(allcondnorm,1)*size(allcondnorm,2));
            
            %find the range of values (max, min)
            maxR = ceil(max(allcondv));
            minR = floor(min(allcondv));
            
            binsize = allsubbins(subcount,session); %e.g. ceil(1+log2(75)) = 8;
            inc = (maxR - minR)/binsize;
            
            allvoxMI = zeros(size(allcondnorm,2),1); %store MI calculated for each voxel
            sigvoxMI = zeros(size(allcondnorm,2),1);
            
            %for each voxel and each discrete response window (1sd), find
            %the P(R) and P(R|S)
            for v = 1:size(allcondnorm,2)
                
                voxR = allcondnorm(:,v); %trial-wise repsonse for each single voxel
                %                 binsize = ceil((max(voxR)-min(voxR))/(2*iqr(voxR)*75^(-1/3))); %Freedman-Diaconis' rule for bin width
                %                 inc = (max(voxR)-min(voxR))/binsize;
                %                allbinsizes = [allbinsizes; binsize];
                
                PR = zeros(binsize,1);%P(R): holds probability for each discrete response window
                cm = 0;
                for m = minR:inc:(maxR-inc)
                    cm = cm+1;
                    PR(cm) = sum(voxR >=m & voxR < m+1)/size(voxR,1); %cm: number of bins
                end
                
                PSR = zeros(binsize,5); %P(R|S): conditioningal probability for each response window given each odor cond
                PS = zeros(5,1)+0.2;
                cm = 0;
                for m = minR:inc:(maxR-inc)
                    cm = cm+1;
                    for sm = 1:5
                        voxRS = voxR((sm-1)*15+1:sm*15);
                        PSR(cm,sm) = sum(voxRS >=m & voxRS < m+1)/size(voxRS,1);%cm: number of bins; sm: number of conditions
                    end
                end
                
                Htotal = -nansum(log2(PR).*(PR)); %skipping the nans (and log2(0) problem) is equivalent to treating them as 0
                Hnoise = 0;
                for sm = 1:5
                    psrs = PSR(:,sm);
                    Hnoise = Hnoise - PS(sm)*nansum(log2(psrs).*(psrs));
                end
                
                voxMI = Htotal - Hnoise;
                allvoxMI(v) = voxMI;
                
            end
            
            
            %             %saving results
            eval(['MI_', roiname,'{subcount,session} = allvoxMI;']);

        end
    end
end

cd('\\psy-wlhost1\lilab\OCF_Analyses');
save OCF_MI_newbinsize_allnorm_31subs.mat MI_apc MI_ppc MI_ofc MI_amg MI_hip

%% Identify preferred odor (out of 5 conditions) for each voxel

sublist = [1:11 13 15:18 20 22:31 34:36 38];%31subs

ROIs = ['rapc';'lapc';'rppc';'lppc';'rofc';'lofc';'ramg';'lamg';'rhip'; 'lhip'];
ROIlist = ['apc';'ppc';'ofc';'amg';'hip'];
sessions = {'PrecEPI1'; 'PostEPI3'; 'P2EPI4'};
conds = [1 2 3 4 5];
k = 5;

%Ogroup_roiname: subNo X 3 sessions cell, each cell containing Ogroup
%assignment
%OvoxC_roiname: subNo X 3sessions (X5 conditions)
for r = 1:5
    region = ROIlist(r,:);
    eval(['Ogroup_',region,' = cell(31,3);']); %preferred odor based on highest freq out of 15 training
    eval(['OBeta_',region,' = cell(31,3);']); %mean beta for each odor
    eval(['OvoxC_',region,' = zeros(31,15);']); %percentage of voxels counted for each preferred odor at pre post post2
end


subcount = 0;
for sub = sublist
    subcount = subcount+1;
    
    for roi = 1:2:9 %collapsing left and right
        roiname = ROIs(roi,2:end);
        
        for session = 1:size(sessions,1)
            
            for cond = [1 2 3 4 5]
                
                betaM = [];
                
                eval(['load  \\psy-wlhost1\lilab\OCF_Analyses\SVM_Processing_no_cap_new32\sub',num2str(sub),'\beta_matrices_LSA\sub',num2str(sub),'_',sessions{session,:},'_',ROIs(roi,:),'_cond',num2str(cond),'.mat;']);
                betaM = [betaM betaMatrix];
                
                eval(['load  \\psy-wlhost1\lilab\OCF_Analyses\SVM_Processing_no_cap_new32\sub',num2str(sub),'\beta_matrices_LSA\sub',num2str(sub),'_',sessions{session,:},'_',ROIs(roi+1,:),'_cond',num2str(cond),'.mat;']);
                betaM = [betaM betaMatrix];
                
                % remove the mean beta value from each trial
                eval(['cond',num2str(cond),' = betaM - repmat(mean(betaM'')'',1,size(betaM,2));']);
            end
            
            if sub==28 && session==2 %fill in the missing trials with zeros (since it will be averaged anyway)
                %  k = 14;
                cond1 = [cond1; nan(1,size(cond1,2))];
                cond2 = [cond2; nan(1,size(cond2,2))];
                cond3 = [cond3; nan(1,size(cond3,2))];
            end
            
            %             %for each voxel, normalize across all trials (i.e. across
            %             %conditions)
            allcondBeta = [cond1; cond2; cond3; cond4; cond5];
            allcondnorm = zeros(size(allcondBeta,1), size(allcondBeta,2));
            
            %          normalize each voxel (z-score)
            for voxel = 1:size(cond1,2)
                normVox = (allcondBeta(:,voxel) - repmat(nanmean(allcondBeta(:,voxel)),size(allcondBeta,1),1))/nanstd(allcondBeta(:,voxel));
                allcondnorm(:,voxel) = normVox;
            end
            
            cond1beta = nanmean(allcondnorm(1:15,:));
            cond2beta = nanmean(allcondnorm(16:30,:));
            cond3beta = nanmean(allcondnorm(31:45,:));
            cond4beta = nanmean(allcondnorm(46:60,:));
            cond5beta = nanmean(allcondnorm(61:75,:));
            
            OResponse = [cond1beta; cond2beta; cond3beta; cond4beta; cond5beta];%holding odor histogram for each voxel; each row for each odor group
            %categorize voxels according to its group assignment
            
            [val, Ogroup] = max(OResponse); %Ogroup indicates which group each voxel belongs to - this is the hard grouping that remained the same through Post/Post2
            %Ogroup = allids_hg;
            
            for o = 1:5
                eval(['o' num2str(o) 'voxC = sum(Ogroup ==o)/size(Ogroup,2);']);
            end
            
            eval(['Ogroup_', roiname,'{subcount,session} = Ogroup;']);
            eval(['OBeta_', roiname,'{subcount,session} = OResponse;']);
            
            if rem(sub,2)==0
                eval(['OvoxC_', roiname,'(subcount,1+(session-1)*5:5+(session-1)*5) = [o1voxC o2voxC o3voxC o4voxC o5voxC];']);
            else
                eval(['OvoxC_', roiname,'(subcount,1+(session-1)*5:5+(session-1)*5) = [o5voxC o4voxC o3voxC o2voxC o1voxC];']);
            end
            
            clear OResponse
            
        end
    end
end

save \\psy-wlhost1\lilab\OCF_Analyses\OCF_allnorm_CondMean_Ogroup_31subs.mat Ogroup_apc Ogroup_ppc Ogroup_ofc Ogroup_amg Ogroup_hip OvoxC_apc OvoxC_ppc OvoxC_amg OvoxC_ofc OvoxC_hip OBeta_apc OBeta_ppc OBeta_ofc OBeta_amg OBeta_hip


%% Find percentages of voxels that shifted tuning after conditioning
%clear all

load \\psy-wlhost1\lilab\OCF_Analyses\Tuningfunctions\OCF_allnorm_CondMean_Ogroup_31subs.mat
load \\psy-wlhost1\lilab\OCF_Analyses\Tuningfunctions\OCF_MI_newbinsize_allnorm_31subs.mat

ROIlist = ['apc';'ppc';'ofc';'amg';'hip'];
sublist = [1:11 13 15:18 20 22:31 34:36 38];%31subs

for r = 1:5
    roiname = ROIlist(r,:);
    eval(['VoxRt_',roiname,' = cell(31,2);']);
    eval(['VoxRt_totalper_',roiname,' = cell(31,2);']);
    eval(['VoxRtBeta_',roiname,' = cell(31,3);']);
    
    subcount = 0;
    
    for s = 1:11%sublist
        
        subcount = subcount+1;
        
        voxloc_post = zeros(5,5); %row: original (pre) voxel odor preference; column: percentage retained for post/post2 O1-O5
        voxloc_post2 = zeros(5,5);
        
        voxloc_post_pertotal = zeros(5,5); %pertotal: percentage of total voxels
        voxloc_post2_pertotal = zeros(5,5);
        
        beta_post = zeros(5,5);
        beta_post2 = zeros(5,5);
        
        eval(['ROItotalvx = size(Ogroup_',roiname,'{subcount,1},2);']);
        
        MIcutoff = .10; %threshold below which to be discarded
        
        %load pre MI and Ogroup
        eval(['MIpre = MI_',roiname,'{subcount,1};']); %grab MI info
        eval(['Ogrouppre = Ogroup_',roiname,'{subcount,1}'';']); %grab Group info
        %discard low MI voxels
        MIOGpre = [(1:ROItotalvx)' Ogrouppre MIpre];
        Ogrouppre_s = MIOGpre(MIpre > quantile(MIpre, MIcutoff),:);
        
        %load post MI and Ogroup-double trimming
        eval(['MIpost = MI_',roiname,'{subcount,2};']); %grab MI info
        MIpost = MIpost(Ogrouppre_s(:,1));
        eval(['Ogrouppost = Ogroup_',roiname,'{subcount,2}'';']); %grab Group info
        Ogrouppost = Ogrouppost(Ogrouppre_s(:,1));
        %discard low MI voxels
        MIOGpost = [Ogrouppre_s(:,1) Ogrouppost MIpost];
        Ogrouppost_s = MIOGpost(MIpost > quantile(MIpost, MIcutoff),:);
        
        %load post2 MI and Ogroup
        eval(['MIpost2 = MI_',roiname,'{subcount,3};']); %grab MI info
        MIpost2 = MIpost2(Ogrouppre_s(:,1));
        eval(['Ogrouppost2 = Ogroup_',roiname,'{subcount,3}'';']); %grab Group info
        Ogrouppost2 = Ogrouppost2(Ogrouppre_s(:,1));
        %discard low MI voxels
        MIOGpost2 = [Ogrouppre_s(:,1) Ogrouppost2 MIpost2];
        Ogrouppost2_s = MIOGpost2(MIpost2 > quantile(MIpost2, MIcutoff),:);
        
        for ogrp = 1:5
            
            %look at preferred voxels at pre
            eval(['lo' num2str(ogrp) 'pre = Ogrouppre_s(Ogrouppre_s(:,2) == ogrp, 1);']);
            eval(['lo' num2str(ogrp) 'pre_perct = size(lo' num2str(ogrp) 'pre,1)/ROItotalvx;']);
            
            for cond = 1:5
                eval(['lo' num2str(cond) 'post = Ogrouppost_s(Ogrouppost_s(:,2) == cond, 1);']);
                eval(['lo' num2str(cond) 'post2 = Ogrouppost2_s(Ogrouppost2_s(:,2) == cond, 1);']);
                
                eval(['voxloc_post(ogrp,cond) = size(intersect(lo' num2str(ogrp) 'pre, lo' num2str(cond) 'post),1)/size(lo' num2str(ogrp) 'pre,1);']);
                eval(['voxloc_post2(ogrp,cond) = size(intersect(lo' num2str(ogrp) 'pre, lo' num2str(cond) 'post2),1)/size(lo' num2str(ogrp) 'pre,1);']);
                
                eval(['voxloc_post_pertotal(ogrp,cond) = size(intersect(lo' num2str(ogrp) 'pre, lo' num2str(cond) 'post),1)/ROItotalvx;']); %abs percentage of total voxels of the ROI
                eval(['voxloc_post2_pertotal(ogrp,cond) = size(intersect(lo' num2str(ogrp) 'pre, lo' num2str(cond) 'post2),1)/ROItotalvx;']);
                
                eval(['beta_pre(ogrp,cond) = mean(OBeta_',roiname,'{subcount,1}(cond,lo' num2str(ogrp) 'pre));']);
                eval(['beta_post(ogrp,cond) = mean(OBeta_',roiname,'{subcount,2}(cond,lo' num2str(ogrp) 'pre));']);
                eval(['beta_post2(ogrp,cond) = mean(OBeta_',roiname,'{subcount,3}(cond,lo' num2str(ogrp) 'pre));']);
                
            end
            
        end
        
        voxloc_pre_pertotal = [lo1pre_perct lo2pre_perct lo3pre_perct lo4pre_perct lo5pre_perct];
        
        if rem(s,2)==1 %rotate matrix to CS+-CS- gradient for odd sub
            voxloc_post = flipud(fliplr(voxloc_post));
            voxloc_post2 = flipud(fliplr(voxloc_post2));
            
            voxloc_pre_pertotal = fliplr(voxloc_pre_pertotal);
            voxloc_post_pertotal = flipud(fliplr(voxloc_post_pertotal));
            voxloc_post2_pertotal = flipud(fliplr(voxloc_post2_pertotal));
            
            beta_pre = flipud(fliplr(beta_post));
            beta_post = flipud(fliplr(beta_post));
            beta_post2 = flipud(fliplr(beta_post2));
        end
        
        eval(['VoxRt_',roiname,'{subcount,1} = voxloc_post;']); %post vox assignment change (retuning)
        eval(['VoxRt_',roiname,'{subcount,2} = voxloc_post2;']); %post2 vox assignment change (retuning)
        
        eval(['VoxRt_totalper_',roiname,'{subcount,1} = voxloc_pre_pertotal;']); 
        eval(['VoxRt_totalper_',roiname,'{subcount,2} = voxloc_post_pertotal;']); 
        eval(['VoxRt_totalper_',roiname,'{subcount,3} = voxloc_post2_pertotal;']); 
        
        eval(['VoxRtBeta_',roiname,'{subcount,1} = beta_pre;']); 
        eval(['VoxRtBeta_',roiname,'{subcount,2} = beta_post;']); 
        eval(['VoxRtBeta_',roiname,'{subcount,3} = beta_post2;']); 
        
    end
    
end

save \\psy-wlhost1\lilab\OCF_Analyses\Tuningfunctions\OCF_allnorm_Condmean_Ogrp_retuneVox_PrePostMI10_31subs.mat

%% Retuning - get the % of Vox remain retuned to its
%% pre-conditioning preference
load \\psy-wlhost1\lilab\OCF_Analyses\Tuningfunctions\OCF_allnorm_Condmean_Ogrp_retuneVox_PrePostMI10_31subs.mat

ROIlist = ['apc';'ppc';'ofc';'amg';'hip'];
GS1PostPost2 = cell(5,1); %i.e. nCSt retuning at post/post2
GS3PostPost2 = cell(5,1); %i.e. nCSs retuning at post/post2

for r = 1:5
    roiname = ROIlist(r,:);
    GS1retune = []; 
    GS3retune = [];
    
    for s = 1:31
        eval(['GS1retune = [GS1retune; VoxRt_',roiname,'{s,1}(2,:) VoxRt_',roiname,'{s,2}(2,:)];']);
        eval(['GS3retune = [GS3retune; VoxRt_',roiname,'{s,1}(4,:) VoxRt_',roiname,'{s,2}(4,:)];']);
        
    end
    GS1PostPost2{r,1} = GS1retune;
    GS3PostPost2{r,1} = GS3retune;
    
end

