%Behavior Analzyer
clear all

SubInit = ['BM';'KM';'NG';'MW';'KP';'JM';'BC';'MK';'DD';'KS';'EM';'AK';'BC';'AC';'EE';'MD';'KP';'HO';'CH';'KJ';'XZ';'AJ';'MS';'MS';'RJ';'KC';'MD';'QH';'LX';'JC';'AA';'TP';'AE';'MS';'AO';'AL';'ZZ';'KL'];

%% Check no responses - 1st session all blocks; check odor recognition accuracy during
%% conditioning

NorespTally = [];
CondAcc = [];
SUDS = [];

allsubs = 1:38;

for s = allsubs
    initials = SubInit(s,:);
    eval(['load OCF_precond_sub' num2str(s) '_' initials '.mat rtypes distressed_preCond1;']);
    
    response = rtypes(:,3);
    Nonan1 = sum(isnan(response)); %sum the No. of no response
    
    eval(['load OCF_cond_sub' num2str(s) '_' initials '.mat rtypes distressed_Cond;']);
    
    response = rtypes(:,3);
    Nonan2 = sum(isnan(response));
    rtypes = rtypes(~isnan(response),:); %exclude no response
    
    odor1 = rtypes(rtypes(:,2)==1,3); %find responses regarding OdorA
    odor5 = rtypes(rtypes(:,2)==5,3); %find responses regarding OdorB
    accOdorA = sum(odor1==7)/length(odor1)*100;
    accOdorB = sum(odor5==9)/length(odor5)*100;
    
    eval(['load OCF_postcond_sub' num2str(s) '_' initials '.mat rtypes distressed_postCond1;']);
    
    response = rtypes(:,3);
    Nonan3 = sum(isnan(response));
    
    NorespTally = [NorespTally; s Nonan1 Nonan2 Nonan3];
    CondAcc = [CondAcc; s accOdorA accOdorB];
    SUDS = [SUDS; s distressed_preCond1 distressed_Cond distressed_postCond1];
    
end

%% Check no responses - 2st session all blocks;

NorespTally = zeros(38,3)+999;
SUDS = zeros(38,3)+999;

%allsubs = [1:13 15:18 20 22:31 33:36 38]; %for postcond2
allsubs = [1:7 9:13 15:18 20 22:31 33:36 38]; %for olocalizer, excluding sub8

for s = allsubs
    initials = SubInit(s,:);
    eval(['load OCF_postcond2_sub' num2str(s) '_' initials '.mat rtypes distressed_postCond1;']);
    
    response = rtypes(:,3);
    Nonan4 = sum(isnan(response));
    
    eval(['load OCF_localizer_sub' num2str(s) '_' initials '.mat resp_all distressed_preLocalizer;']);
    
    response = resp_all(:,2);
    Nonan5 = sum(isnan(response));
    
    
    NorespTally(s,:) = [s Nonan4 Nonan5];
    SUDS(s,:) = [s distressed_postCond1 distressed_preLocalizer];
end

%% Check StimRs
allsubs = [1:7 9:13 15:18 20 22:31 33:36 38]; %for olocalizer, excluding sub8
tally = [];

for s = allsubs
    initials = SubInit(s,:);
    eval(['load OCF_postcond_sub' num2str(s) '_' initials '.mat StimR;']);
    
    S1 = StimR;
    
    eval(['load OCF_postcond2_sub' num2str(s) '_' initials '.mat StimR;']);
    
    S2 = StimR;
    Sdiff = S1-S2;
    
    tally = [tally; s Sdiff'];
end


%% Precond - %CS+ response; CS+ GS1 GS2 GS3 CS-
%Odd No. Sub: B as CS+; Even No. Sub: A as CS+
AllpreCS_acc = [];
AllpreCS_rt = [];%average RTs for all/correct odorA response

AllpreodorA = [];
AllodorA_rt = [];
AllpreodorB = [];
AllodorB_rt = [];

allsubsE = 2:2:38;
for s = allsubsE
    initials = SubInit(s,:);
    eval(['load OCF_precond_sub' num2str(s) '_' initials '.mat rtypes;']);
    
    response = rtypes(:,3);
    rtypes = rtypes(~isnan(response),:); %exclude no response
    
    odor1 = rtypes(rtypes(:,2)==1, 3);
    odor2 = rtypes(rtypes(:,2)==2, 3);
    odor3 = rtypes(rtypes(:,2)==3, 3);
    odor4 = rtypes(rtypes(:,2)==4, 3);
    odor5 = rtypes(rtypes(:,2)==5, 3);
    
    odor1A = length(find(odor1 == 7)) /length(odor1)*100; % proportion of odor A judgment
    odor2A = length(find(odor2 == 7)) /length(odor2)*100; % proportion of odor A judgment
    odor3A = length(find(odor3 == 7)) /length(odor3)*100; % proportion of odor A judgment
    odor4A = length(find(odor4 == 7)) /length(odor4)*100; % proportion of odor A judgment
    odor5A = length(find(odor5 == 7)) /length(odor5)*100; % proportion of odor A judgment
    
    allrt = rtypes(:,4);
    rtindex = (allrt > mean(allrt) - 3*std(allrt)) & (allrt < mean(allrt)+3*std(allrt));
    rtypes_rt = rtypes(rtindex,:);
    
    odor1rt = [rtypes_rt(rtypes_rt(:,2)==1, 4)];
    odor2rt = [rtypes_rt(rtypes_rt(:,2)==2, 4)];
    odor3rt = [rtypes_rt(rtypes_rt(:,2)==3, 4)];
    odor4rt = [rtypes_rt(rtypes_rt(:,2)==4, 4)];
    odor5rt = [rtypes_rt(rtypes_rt(:,2)==5, 4)];
    
    
    AllpreodorA = [AllpreodorA; s odor1A odor2A odor3A odor4A odor5A];
    AllodorA_rt = [AllodorA_rt; s mean(odor1rt) mean(odor2rt) mean(odor3rt) mean(odor4rt) mean(odor5rt)];
    
end

allsubsO = 1:2:37;
for s = allsubsO
    initials = SubInit(s,:);
    eval(['load OCF_precond_sub' num2str(s) '_' initials '.mat rtypes;']);
    
    response = rtypes(:,3);
    rtypes = rtypes(~isnan(response),:); %exclude no response
    
    odor1 = rtypes(rtypes(:,2)==1, 3);
    odor2 = rtypes(rtypes(:,2)==2, 3);
    odor3 = rtypes(rtypes(:,2)==3, 3);
    odor4 = rtypes(rtypes(:,2)==4, 3);
    odor5 = rtypes(rtypes(:,2)==5, 3);
    
    odor1B = length(find(odor1 == 9)) /length(odor1)*100; % proportion of odor A judgment
    odor2B = length(find(odor2 == 9)) /length(odor2)*100; % proportion of odor A judgment
    odor3B = length(find(odor3 == 9)) /length(odor3)*100; % proportion of odor A judgment
    odor4B = length(find(odor4 == 9)) /length(odor4)*100; % proportion of odor A judgment
    odor5B = length(find(odor5 == 9)) /length(odor5)*100; % proportion of odor A judgment
    
    allrt = rtypes(:,4);
    rtindex = (allrt > mean(allrt) - 3*std(allrt)) & (allrt < mean(allrt)+3*std(allrt));
    rtypes_rt = rtypes(rtindex,:);
    
    odor1rt = [rtypes_rt(rtypes_rt(:,2)==1, 4)];
    odor2rt = [rtypes_rt(rtypes_rt(:,2)==2, 4)];
    odor3rt = [rtypes_rt(rtypes_rt(:,2)==3, 4)];
    odor4rt = [rtypes_rt(rtypes_rt(:,2)==4, 4)];
    odor5rt = [rtypes_rt(rtypes_rt(:,2)==5, 4)];
    
    
    AllpreodorB = [AllpreodorB; s odor5B odor4B odor3B odor2B odor1B];
    AllodorB_rt = [AllodorB_rt; s mean(odor5rt) mean(odor4rt) mean(odor3rt) mean(odor2rt) mean(odor1rt)];
    
end

AllpreCS_acc = [AllpreodorA; AllpreodorB];
AllpreCS_rt = [AllodorA_rt; AllodorB_rt];



%% Postcond - %CS+ response; CS+ GS1 GS2 GS3 CS-
AllpostCS_acc = [];
AllpostCS_rt = [];

AllpostodorA = [];
AllodorA_rt = [];%average RTs for all/correct odorA response

allsubs = 2:2:38; %OdorA CS+
%allsubs = 1:2:37; %OdorB CS+

for s = allsubs
    
    initials = SubInit(s,:);
    eval(['load OCF_postcond_sub' num2str(s) '_' initials '.mat rtypes;']);
    
    if s==1 || s == 2 || s ==3 || s==4 || s==5
        %remove CS+ reinforcement trials
        rtypes(rtypes(:,1)==3|rtypes(:,1)==28|rtypes(:,1)==38|rtypes(:,1)==59|rtypes(:,1)==76,:) = [];
    else
        %remove CS+ reinforcement trials
        rtypes(rtypes(:,1)==3|rtypes(:,1)==14|rtypes(:,1)==38|rtypes(:,1)==59|rtypes(:,1)==76,:) = [];
    end
    
    response = rtypes(:,3);
    rtypes = rtypes(~isnan(response),:); %exclude no response
    
    odor1 = rtypes(rtypes(:,2)==1, 3);
    odor2 = rtypes(rtypes(:,2)==2, 3);
    odor3 = rtypes(rtypes(:,2)==3, 3);
    odor4 = rtypes(rtypes(:,2)==4, 3);
    odor5 = rtypes(rtypes(:,2)==5, 3);
    
    odor1A = length(find(odor1 == 7)) /length(odor1)*100; % proportion of odor A judgment
    odor2A = length(find(odor2 == 7)) /length(odor2)*100; % proportion of odor A judgment
    odor3A = length(find(odor3 == 7)) /length(odor3)*100; % proportion of odor A judgment
    odor4A = length(find(odor4 == 7)) /length(odor4)*100; % proportion of odor A judgment
    odor5A = length(find(odor5 == 7)) /length(odor5)*100; % proportion of odor A judgment
    
    allrt = rtypes(:,4);
    rtindex = (allrt > mean(allrt) - 3*std(allrt)) & (allrt < mean(allrt)+3*std(allrt));
    rtypes_rt = rtypes(rtindex,:);
    
    odor1rt = [rtypes_rt(rtypes_rt(:,2)==1, 4)];
    odor2rt = [rtypes_rt(rtypes_rt(:,2)==2, 4)];
    odor3rt = [rtypes_rt(rtypes_rt(:,2)==3, 4)];
    odor4rt = [rtypes_rt(rtypes_rt(:,2)==4, 4)];
    odor5rt = [rtypes_rt(rtypes_rt(:,2)==5, 4)];
    
    AllpostodorA = [AllpostodorA; s odor1A odor2A odor3A odor4A odor5A];
    AllodorA_rt = [AllodorA_rt; s mean(odor1rt) mean(odor2rt) mean(odor3rt) mean(odor4rt) mean(odor5rt)];
    
end


AllpostodorB = [];
AllodorB_rt = [];%average RTs for all/correct odorA response

%allsubs = 2:2:38; %OdorA CS+
allsubs = 1:2:37; %OdorB CS+

for s = allsubs
    
    initials = SubInit(s,:);
    eval(['load OCF_postcond_sub' num2str(s) '_' initials '.mat rtypes;']);
    
    if s==1 || s == 2 || s ==3 || s==4 || s==5
        %remove CS+ reinforcement trials
        rtypes(rtypes(:,1)==3|rtypes(:,1)==28|rtypes(:,1)==38|rtypes(:,1)==59|rtypes(:,1)==76,:) = [];
    else
        %remove CS+ reinforcement trials
        rtypes(rtypes(:,1)==3|rtypes(:,1)==14|rtypes(:,1)==38|rtypes(:,1)==59|rtypes(:,1)==76,:) = [];
    end
    
    response = rtypes(:,3);
    rtypes = rtypes(~isnan(response),:); %exclude no response
    
    odor1 = rtypes(rtypes(:,2)==1, 3);
    odor2 = rtypes(rtypes(:,2)==2, 3);
    odor3 = rtypes(rtypes(:,2)==3, 3);
    odor4 = rtypes(rtypes(:,2)==4, 3);
    odor5 = rtypes(rtypes(:,2)==5, 3);
    
    odor1B = length(find(odor1 == 9)) /length(odor1)*100; % proportion of odor A judgment
    odor2B = length(find(odor2 == 9)) /length(odor2)*100; % proportion of odor A judgment
    odor3B = length(find(odor3 == 9)) /length(odor3)*100; % proportion of odor A judgment
    odor4B = length(find(odor4 == 9)) /length(odor4)*100; % proportion of odor A judgment
    odor5B = length(find(odor5 == 9)) /length(odor5)*100; % proportion of odor A judgment
    
    allrt = rtypes(:,4);
    rtindex = (allrt > mean(allrt) - 3*std(allrt)) & (allrt < mean(allrt)+3*std(allrt));
    rtypes_rt = rtypes(rtindex,:);
    
    odor1rt = [rtypes_rt(rtypes_rt(:,2)==1, 4)];
    odor2rt = [rtypes_rt(rtypes_rt(:,2)==2, 4)];
    odor3rt = [rtypes_rt(rtypes_rt(:,2)==3, 4)];
    odor4rt = [rtypes_rt(rtypes_rt(:,2)==4, 4)];
    odor5rt = [rtypes_rt(rtypes_rt(:,2)==5, 4)];
    
    AllpostodorB = [AllpostodorB; s odor5B odor4B odor3B odor2B odor1B];
    AllodorB_rt = [AllodorB_rt; s mean(odor5rt) mean(odor4rt) mean(odor3rt) mean(odor2rt) mean(odor1rt)];
end

AllpostCS_acc = [AllpostodorA; AllpostodorB];
AllpostCS_rt = [AllodorA_rt; AllodorB_rt];


%% Postcond2 - Sub2

csp_trials = [3,14,38,59,76];

AllpostodorA = zeros(19,6)+999;
AllodorA_rt = zeros(19,6)+999;%average RTs for all/correct odorA response

%allsubs = [4 6 8 10 12 16 18 20 22 24 26 28 30 34 36 38]; %OdorA CS+; Sub2 needs to be treated separately for it has 2 parts
s = 2; %OdorB CS+

initials = SubInit(s,:);
eval(['load OCF_postcond2_sub' num2str(s) '_' initials '_combined.mat rtypes1 rtypes2;']);

rtypes = [rtypes1(1:22,:); rtypes2];

%remove CS+ reinforcement trials
rtypes(rtypes(:,1)==3|rtypes(:,1)==14|rtypes(:,1)==38|rtypes(:,1)==59|rtypes(:,1)==76,:) = [];

response = rtypes(:,3);
rtypes = rtypes(~isnan(response),:); %exclude no response

odor1 = rtypes(rtypes(:,2)==1, 3);
odor2 = rtypes(rtypes(:,2)==2, 3);
odor3 = rtypes(rtypes(:,2)==3, 3);
odor4 = rtypes(rtypes(:,2)==4, 3);
odor5 = rtypes(rtypes(:,2)==5, 3);

odor1A = length(find(odor1 == 7)) /length(odor1)*100; % proportion of odor A judgment
odor2A = length(find(odor2 == 7)) /length(odor2)*100; % proportion of odor A judgment
odor3A = length(find(odor3 == 7)) /length(odor3)*100; % proportion of odor A judgment
odor4A = length(find(odor4 == 7)) /length(odor4)*100; % proportion of odor A judgment
odor5A = length(find(odor5 == 7)) /length(odor5)*100; % proportion of odor A judgment

allrt = rtypes(:,4);
rtindex = (allrt > mean(allrt) - 3*std(allrt)) & (allrt < mean(allrt)+3*std(allrt));
rtypes_rt = rtypes(rtindex,:);

odor1rt = [rtypes_rt(rtypes_rt(:,2)==1, 4)];
odor2rt = [rtypes_rt(rtypes_rt(:,2)==2, 4)];
odor3rt = [rtypes_rt(rtypes_rt(:,2)==3, 4)];
odor4rt = [rtypes_rt(rtypes_rt(:,2)==4, 4)];
odor5rt = [rtypes_rt(rtypes_rt(:,2)==5, 4)];


AllpostodorA(s/2,:) = [s odor1A odor2A odor3A odor4A odor5A];
AllodorA_rt(s/2,:) = [s mean(odor1rt) mean(odor2rt) mean(odor3rt) mean(odor4rt) mean(odor5rt)];

%% Postcond2 - %CS+ response; CS+ GS1 GS2 GS3 CS-
Allpost2CS_acc = [];
Allpost2CS_rt = [];

Allpost2odorA = [];
AllodorA_rt = [];%average RTs for all/correct odorA response

Allpost2odorB = [];
AllodorB_rt = [];%average RTs for all/correct odorA response

allsubs = [2 4 6 8 10 12 16 18 20 22 24 26 28 30 34 36 38];
%allsubs = 1:2:37; %OdorB CS+

for s = allsubs
    
    initials = SubInit(s,:);
    if s == 2
        eval(['load OCF_postcond2_sub' num2str(s) '_' initials '_combined.mat rtypes1 rtypes2;']);
        rtypes = [rtypes1(1:22,:); rtypes2];
    else
        eval(['load OCF_postcond2_sub' num2str(s) '_' initials '.mat rtypes;']);
    end
    
    
    %remove CS+ reinforcement trials
    rtypes(rtypes(:,1)==3|rtypes(:,1)==14|rtypes(:,1)==38|rtypes(:,1)==59|rtypes(:,1)==76,:) = [];
    
    
    response = rtypes(:,3);
    rtypes = rtypes(~isnan(response),:); %exclude no response
    
    odor1 = rtypes(rtypes(:,2)==1, 3);
    odor2 = rtypes(rtypes(:,2)==2, 3);
    odor3 = rtypes(rtypes(:,2)==3, 3);
    odor4 = rtypes(rtypes(:,2)==4, 3);
    odor5 = rtypes(rtypes(:,2)==5, 3);
    
    odor1A = length(find(odor1 == 7)) /length(odor1)*100; % proportion of odor A judgment
    odor2A = length(find(odor2 == 7)) /length(odor2)*100; % proportion of odor A judgment
    odor3A = length(find(odor3 == 7)) /length(odor3)*100; % proportion of odor A judgment
    odor4A = length(find(odor4 == 7)) /length(odor4)*100; % proportion of odor A judgment
    odor5A = length(find(odor5 == 7)) /length(odor5)*100; % proportion of odor A judgment
    
    allrt = rtypes(:,4);
    rtindex = (allrt > mean(allrt) - 3*std(allrt)) & (allrt < mean(allrt)+3*std(allrt));
    rtypes_rt = rtypes(rtindex,:);
    
    odor1rt = [rtypes_rt(rtypes_rt(:,2)==1, 4)];
    odor2rt = [rtypes_rt(rtypes_rt(:,2)==2, 4)];
    odor3rt = [rtypes_rt(rtypes_rt(:,2)==3, 4)];
    odor4rt = [rtypes_rt(rtypes_rt(:,2)==4, 4)];
    odor5rt = [rtypes_rt(rtypes_rt(:,2)==5, 4)];
    
    Allpost2odorA = [Allpost2odorA; s odor1A odor2A odor3A odor4A odor5A];
    AllodorA_rt = [AllodorA_rt; s mean(odor1rt) mean(odor2rt) mean(odor3rt) mean(odor4rt) mean(odor5rt)];
    
end


%allsubs = 2:2:38; %OdorA CS+
allsubs = [1 3 5 7 9 11 13 15 17 23 25 27 29 31 33 35]; %OdorB CS+

for s = allsubs
    
    initials = SubInit(s,:);
    eval(['load OCF_postcond2_sub' num2str(s) '_' initials '.mat rtypes;']);
    
    
    %remove CS+ reinforcement trials
    rtypes(rtypes(:,1)==3|rtypes(:,1)==14|rtypes(:,1)==38|rtypes(:,1)==59|rtypes(:,1)==76,:) = [];
    
    
    response = rtypes(:,3);
    rtypes = rtypes(~isnan(response),:); %exclude no response
    
    odor1 = rtypes(rtypes(:,2)==1, 3);
    odor2 = rtypes(rtypes(:,2)==2, 3);
    odor3 = rtypes(rtypes(:,2)==3, 3);
    odor4 = rtypes(rtypes(:,2)==4, 3);
    odor5 = rtypes(rtypes(:,2)==5, 3);
    
    odor1B = length(find(odor1 == 9)) /length(odor1)*100; % proportion of odor A judgment
    odor2B = length(find(odor2 == 9)) /length(odor2)*100; % proportion of odor A judgment
    odor3B = length(find(odor3 == 9)) /length(odor3)*100; % proportion of odor A judgment
    odor4B = length(find(odor4 == 9)) /length(odor4)*100; % proportion of odor A judgment
    odor5B = length(find(odor5 == 9)) /length(odor5)*100; % proportion of odor A judgment
    
    allrt = rtypes(:,4);
    rtindex = (allrt > mean(allrt) - 3*std(allrt)) & (allrt < mean(allrt)+3*std(allrt));
    rtypes_rt = rtypes(rtindex,:);
    
    odor1rt = [rtypes_rt(rtypes_rt(:,2)==1, 4)];
    odor2rt = [rtypes_rt(rtypes_rt(:,2)==2, 4)];
    odor3rt = [rtypes_rt(rtypes_rt(:,2)==3, 4)];
    odor4rt = [rtypes_rt(rtypes_rt(:,2)==4, 4)];
    odor5rt = [rtypes_rt(rtypes_rt(:,2)==5, 4)];
    
    Allpost2odorB = [Allpost2odorB; s odor5B odor4B odor3B odor2B odor1B];
    AllodorB_rt = [AllodorB_rt; s mean(odor5rt) mean(odor4rt) mean(odor3rt) mean(odor2rt) mean(odor1rt)];
end

Allpost2CS_acc = [Allpost2odorA; Allpost2odorB];
Allpost2CS_rt = [AllodorA_rt; AllodorB_rt];


%% Risk Rating Post

%allsubs = [2:9 11 12 15 17 18 20 22 23 26:30 33 35 37 38];
allsubs = [2:18 20:36 38];
AllRisk = zeros(38,6)+999;

for s = allsubs
    initials = SubInit(s,:);
    eval(['load OCF_risk_sub' num2str(s) '_' initials '.mat risk_sorted;']);
    
    O1 = mean(risk_sorted(1:3,2));
    O2 = mean(risk_sorted(4:6,2));
    O3 = mean(risk_sorted(7:9,2));
    O4 = mean(risk_sorted(10:12,2));
    O5 = mean(risk_sorted(13:15,2));
    
    if rem(s,2)==0
        AllRisk(s,:) = [s O1 O2 O3 O4 O5];
    else
        AllRisk(s,:) = [s O5 O4 O3 O2 O1];
    end
    
end

%% Risk Rating Post2

%allsubs = [2:9 11 12 15 17 18 20 22 23 26:30 33 35 37 38];
allsubs = [1:13 15:18 20 22:31 33:36 38];
AllRisk2 = zeros(38,6)+999;

for s = allsubs
    initials = SubInit(s,:);
    eval(['load OCF_risk2_sub' num2str(s) '_' initials '.mat risk_sorted;']);
    
    O1 = mean(risk_sorted(1:3,2));
    O2 = mean(risk_sorted(4:6,2));
    O3 = mean(risk_sorted(7:9,2));
    O4 = mean(risk_sorted(10:12,2));
    O5 = mean(risk_sorted(13:15,2));
    
    if rem(s,2)==0
        AllRisk2(s,:) = [s O1 O2 O3 O4 O5];
    else
        AllRisk2(s,:) = [s O5 O4 O3 O2 O1];
    end
    
end