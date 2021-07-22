%%Quantatively determine good events for statistical study 2019/09/18
tic
close all
clear
warning('off','MATLAB:polyfit:RepeatedPointsOrRescale')
cd '/Users/andrewvu/Library/Mobile Documents/com~apple~CloudDocs/Research/Analysis/Criterion'
%% Set up
SSthreshold = 1;
durationForAverage = 6;
edgeNRatio = 0; %1.25;
%Structure For Storing Values
criterionStructure = struct(...
    'eventNumber', [],...
    'coreDuration',[],...
    'coreSpeed',[],...
    'coreSpeedx',[],...
    'upSpeed',[],...
    'downSpeed', [], ...
    'vRatio', [], ...
    'vxRatio', [], ...
    'vxMinRatio', [], ...
    'vDifference',[],...
    'DeflectionAngle',[],...
    'coreDensity', [], ...
    'coreMinDensity',[],...
    'downDensity', [], ...
    'upDensity', [],...
    'nRatio', [],...
    'nminRatio',[],...
    'eTpara', [],...
    'eTperp', [],...
    'criteriaPassed',[],...
    'vCorr', [], ...
    'vxCorr', [], ...
    'tempparaCorr', [], ...
    'tempperpCorr', [], ...
    'tempCorr', [],...
    'BmagCorr',[],...
    'l_nBCorr',[],...
    'l_nvCorr',[],...
    'l_vBCorr',[],...
    'l_nTCorr',[],...
    'l_meanB',[],...
    'l_maxB',[],...
    'l_deltaB',[],...
    'l_maxBoverAmbientB',[],...
    'l_meanN',[],...
    'l_maxN',[],...
    'l_deltaN',[],...
    'l_maxNoverAmbientN',[],...
    'l_meanV',[],...
    'l_maxV',[],...
    'l_deltaV',[],...
    'l_maxVoverAmbientV',[],...
    't_nBCorr',[],...
    't_nvCorr',[],...
    't_vBCorr',[],...
    't_nTCorr',[],...
    't_meanB',[],...
    't_maxB',[],...
    't_deltaB',[],...
    't_maxBoverAmbientB',[],...
    't_meanN',[],...
    't_maxN',[],...
    't_deltaN',[],...
    't_maxNoverAmbientN',[],...
    't_meanV',[],...
    't_maxV',[],...
    't_deltaV',[],...
    't_maxVoverAmbientV',[]);
%% Chosens
AtypicalEvents = [4,10,13,23,26,27,40,44,45,46,49,59,64,69,72,75,88,89,91,96,102,107,108,109,110,112,113,114,119,120,126,128,129,131,137,140,141,144,145,146,147,149,150,151,163,173,177,180];
TypicalEvents = [2,3,5,6,8,11,12,14,16,17,18,20,21,24,29,31,32,35,37,38,41,42,43,48,50,52,54,55,58,60,65,67,68,76,79,81,82,83,84,85,93,94,98,100,105,111,116,117,121,122,123,124,127,132,133,135,136,138,142,153,154,157,158,159,160,161,162,164,165,166,167,168,169,170];
TypicalEvents = [2;3;5;6;7;8;9;11;12;14;15;16;17;18;20;21;22;24;25;28;29;30;31;32;34;35;36;37;38;39;41;42;43;47;48;50;51;52;54;55;56;57;58;60;61;65;66;67;68;76;77;78;79;80;81;82;83;84;85;87;93;94;97;98;99;100;104;105;106;111;116;117;121;122;123;124;127;130;132;133;135;136;138;142;152;153;154;157;158;159;160;161;162;164;165;166;167;168;169;170];
Event_number = [2:174,179:207];
%TypicalEvents7[2;3;5;6;7;8;9;11;12;14;15;16;17;18;20;21;22;24;25;28;29;30;31;32;34;35;36;37;38;39;41;42;43;47;48;50;51;52;54;55;56;57;58;60;61;65;66;67;68;76;77;78;79;80;81;82;83;84;85;87;93;94;97;98;99;100;105;106;111;116;117;121;122;123;124;127;130;132;133;135;136;138;142;143;152;153;154;157;158;159;160;161;162;164;165;167;168;169;170;179;183;184;185;187;188;189;191;199;200;201;202;203;204;205]
%% Loop
for i=104%Event_number
        if sum(i==AtypicalEvents) == 1
            continue
        end
%     if sum(i==TypicalEvents) ~= 1
%         continue
%     end
    i
    %Get Event Specifications
    [Event_Type, ~,~,~,...
        event_start,event_end,...
        left_OuterEdge, left_InnerEdge,...
        right_InnerEdge, right_OuterEdge,...
        ~, ~,...
        ~, ~,...
        ~, ~,...
        ~, ~] = get_eventTimes(i);
    
    eventDataFileName = strcat('MMS1_Data_EventNumber_',num2str(i),'.mat');
    eventDataDirectory =  '/Users/andrewvu/data/Event Data/';
    eventDataDirectoryFileName = strcat(eventDataDirectory,eventDataFileName);

    
    if exist(eventDataDirectoryFileName,'file') ~= 2
        
        %Load MEC Data
        [mms1_mec_timedata_raw, mms1_mec_rdata_raw] = load_mec(event_start,event_end,1,'brst');
        [mms2_mec_timedata_raw, mms2_mec_rdata_raw] = load_mec(event_start,event_end,2,'brst');
        [mms3_mec_timedata_raw, mms3_mec_rdata_raw] = load_mec(event_start,event_end,3,'brst');
        [mms4_mec_timedata_raw, mms4_mec_rdata_raw] = load_mec(event_start,event_end,4,'brst');
        
        %load FGM data
        [mms1_fgm_timedata_raw, mms1_fgm_bdata_raw, ~, ~] = load_fgm(event_start,event_end,1,'brst');
        [mms2_fgm_timedata_raw, mms2_fgm_bdata_raw, ~, ~] = load_fgm(event_start,event_end,2,'brst');
        [mms3_fgm_timedata_raw, mms3_fgm_bdata_raw, ~, ~] = load_fgm(event_start,event_end,3,'brst');
        [mms4_fgm_timedata_raw, mms4_fgm_bdata_raw, ~, ~] = load_fgm(event_start,event_end,4,'brst');
        
        [mms1_fgm_timedata_srvy, mms1_fgm_bdata_srvy] = load_fgm(event_start,event_end,1,'srvy'); %For Sliding Window
        
        %Load FPI_e
        [fpi_e_timedata,fpi_e_ndata,fpi_e_vdata,fpi_e_tparadata,fpi_e_tperpdata,...
            fpi_e_edata,fpi_e_espectdata,fpi_e_pressdata] = load_fpi(event_start,event_end,1,'brst','e');
        %Load FPI_i
        [fpi_i_timedata,fpi_i_ndata,fpi_i_vdata,fpi_i_tparadata,fpi_i_tperpdata,...
            fpi_i_edata,fpi_i_espectdata,fpi_i_pressdata] = load_fpi(event_start,event_end,1,'brst','i');
        
        %Load FPI Dist
        [time,phi_vector,theta_vector,energy_vector,dist] = load_dist(event_start,event_end,1,'brst','i');

        
        
        save(eventDataDirectoryFileName,...
            'mms1_mec_timedata_raw',...
            'mms1_mec_rdata_raw',...
            'mms2_mec_timedata_raw',...
            'mms2_mec_rdata_raw',...
            'mms3_mec_timedata_raw',...
            'mms3_mec_rdata_raw',...
            'mms4_mec_timedata_raw',...
            'mms4_mec_rdata_raw',...
            'mms1_fgm_timedata_raw',...
            'mms1_fgm_bdata_raw',...
            'mms2_fgm_timedata_raw',...
            'mms2_fgm_bdata_raw',...
            'mms3_fgm_timedata_raw',...
            'mms3_fgm_bdata_raw',...
            'mms4_fgm_timedata_raw',...
            'mms4_fgm_bdata_raw',...
            'mms1_fgm_timedata_srvy',...
            'mms1_fgm_bdata_srvy',...
            'fpi_e_timedata',...
            'fpi_e_ndata',...
            'fpi_e_vdata',...
            'fpi_e_tparadata',...
            'fpi_e_tperpdata',...
            'fpi_e_edata',...
            'fpi_e_espectdata',...
            'fpi_e_pressdata',...
            'fpi_i_timedata',...
            'fpi_i_ndata',...
            'fpi_i_vdata',...
            'fpi_i_tparadata',...
            'fpi_i_tperpdata',...
            'fpi_i_edata',...
            'fpi_i_espectdata',...
            'fpi_i_pressdata',...
            'time',...
            'phi_vector',...
            'theta_vector',...
            'energy_vector',...
            'dist')
    else
        load(eventDataDirectoryFileName) %#ok<LOAD>
    end
    % % % %
% % % %     %Load MEC Data
% % % %     % [mms1_mec_timedata_raw, mms1_mec_rdata_raw] = load_mec(event_start,1,data_type);
% % % %     % [mms2_mec_timedata_raw, mms2_mec_rdata_raw] = load_mec(event_start,2,data_type);
% % % %     % [mms3_mec_timedata_raw, mms3_mec_rdata_raw] = load_mec(event_start,3,data_type);
% % % %     % [mms4_mec_timedata_raw, mms4_mec_rdata_raw] = load_mec(event_start,4,data_type);
% % % %     
% % % %     %load FGM data
% % % %     [mms1_fgm_timedata_raw, mms1_fgm_bdata_raw, ~, ~] = load_fgm(event_start,event_end,1,'brst');
% % % %     % [mms2_fgm_timedata_raw, mms2_fgm_bdata_raw, ~, ~] = load_fgm(event_start,event_end,2,'brst');
% % % %     % [mms3_fgm_timedata_raw, mms3_fgm_bdata_raw, ~, ~] = load_fgm(event_start,event_end,3,'brst');
% % % %     % [mms4_fgm_timedata_raw, mms4_fgm_bdata_raw, ~, ~] = load_fgm(event_start,event_end,4,'brst');
% % % %     
% % % %     % [mms1_fgm_timedata_srvy, mms1_fgm_bdata_srvy] = load_fgm(event_start,event_end,1,'srvy'); %For Sliding Window
% % % %     
% % % %     %Load FPI_e
% % % %     [fpi_e_timedata,fpi_e_ndata,~,fpi_e_tparadata,fpi_e_tperpdata,...
% % % %         ~,~,~] = load_fpi(event_start,event_end,1,'brst','e');
% % % %     
% % % %     %Load FPI_i
% % % %     [fpi_i_timedata,fpi_i_ndata,fpi_i_vdata,fpi_i_tparadata,fpi_i_tperpdata,...
% % % %         fpi_i_edata,fpi_i_espectdata] = load_fpi(event_start,event_end,1,'brst','i');
    %% Duration Criteria
    Duration = findDuration(left_InnerEdge,right_InnerEdge);
    %% Speed Criteria
    %Calculate mean Core Speed Deflection
    fpi_i_vdata(:,4) = vecnorm(fpi_i_vdata,2,2);
    [meanCoreSpeed] = calculate_meanDepletion(fpi_i_timedata,fpi_i_vdata(:,4),left_InnerEdge,right_InnerEdge,SSthreshold);
    [meanCoreSpeedx] = calculate_meanDepletion(fpi_i_timedata,fpi_i_vdata(:,1),left_InnerEdge,right_InnerEdge,SSthreshold);
    [~,minCoreSpeedx]=calculate_minCore(fpi_i_timedata,abs(fpi_i_vdata(:,1)),left_InnerEdge,right_InnerEdge);
    %Calculate Ambient Solar Wind Speed
    [downstreamSpeed,upstreamSpeed,~,~] = calculate_updownSW(event_start,event_end,i);
    close all
    %Calculate Ratio of of much speed has been deflected
    deflectionRatio = min(meanCoreSpeed./[downstreamSpeed upstreamSpeed],[],'all');
    deflectionxRatio = min(abs(meanCoreSpeedx)./[downstreamSpeed upstreamSpeed],[],'all');
    
    %Calculate magnitude difference
    speedDiff = max([downstreamSpeed upstreamSpeed]-abs(meanCoreSpeedx));
    
    %Calculate angle of deflection
    [dangle,dratio] = calculate_vDeflection(fpi_i_timedata,fpi_i_vdata,event_start,left_OuterEdge,left_InnerEdge,right_InnerEdge,right_OuterEdge,event_end);
    deflectionxminRatio = min(abs(minCoreSpeedx)./[downstreamSpeed upstreamSpeed],[],'all'); 
    %% Density Criteria
    %Calculate mean core Density without substructures (SSThreshold)
    [meanCoreDensity] = calculate_meanDepletion(fpi_i_timedata,fpi_i_ndata,left_InnerEdge,right_InnerEdge,SSthreshold);
    [~,minCoreDensity]=calculate_minCore(fpi_i_timedata,fpi_i_ndata,left_InnerEdge,right_InnerEdge);

    %Calculate ambient densities
    [downstreamDensity,upstreamDensity] = calculate_prepostAverages(fpi_i_timedata,fpi_i_ndata,event_start,event_end,'i',durationForAverage);
    
    %Calculate Ratio of how much density is depleted within Core
    depletionRatio = min(meanCoreDensity./[downstreamDensity upstreamDensity]);
    depletionminRatio = min(minCoreDensity./[downstreamDensity upstreamDensity]);
    %% Electron Temperature Heating
    [meanCoreETperp] = calculate_meanDepletion(fpi_e_timedata,fpi_e_tperpdata,left_InnerEdge,right_InnerEdge,0);
    [meanCoreETpara] = calculate_meanDepletion(fpi_e_timedata,fpi_e_tparadata,left_InnerEdge,right_InnerEdge,0);
    %% Correlations
    %Interpolate FGM to FPI
    [~,fgm_bdata] = interpxyz(mms1_fgm_timedata_raw,mms1_fgm_bdata_raw(:,1:3),fpi_i_timedata);
    [~,fpi_e_ndata] = interpxyz(fpi_e_timedata,fpi_e_ndata,fpi_i_timedata);
    
    [criterionStructure.vCorr(i),criterionStructure.vxCorr(i),...
        criterionStructure.tempparaCorr(i),criterionStructure.tempperpCorr(i),...
        criterionStructure.tempCorr(i),...
        criterionStructure.BmagCorr(i)] = calculate_correlation(left_InnerEdge,right_InnerEdge,...
        fpi_i_timedata,fpi_i_ndata,fpi_i_vdata(:,1:3),fpi_i_tparadata,fpi_i_tperpdata,fgm_bdata);
    
    %Event Edges Analysis
    [criterionStructure.l_nBCorr(i),criterionStructure.l_nvCorr(i),criterionStructure.l_vBCorr(i),criterionStructure.l_nTCorr(i),...
        criterionStructure.l_meanB(i),criterionStructure.l_maxB(i),criterionStructure.l_deltaB(i),criterionStructure.l_maxBoverAmbientB(i),...
        criterionStructure.l_meanN(i),criterionStructure.l_maxN(i),criterionStructure.l_deltaN(i),l_maxNoverAmbientN_i,...
        criterionStructure.l_meanV(i),criterionStructure.l_maxV(i),criterionStructure.l_deltaV(i),criterionStructure.l_maxVoverAmbientV(i)]...
        = calculate_edgesAnalysis(event_start,event_end,left_OuterEdge,left_InnerEdge,fpi_i_timedata,fpi_i_ndata,fpi_i_vdata(:,1:3),fpi_i_tparadata,fpi_i_tperpdata,fgm_bdata,durationForAverage,'leading');
    
    [~,~,~,~,...
        ~,~,~,~,...
        ~,~,~,l_maxNoverAmbientN_e,...
        ~,~,~,~]...
        = calculate_edgesAnalysis(event_start,event_end,left_OuterEdge,left_InnerEdge,fpi_i_timedata,fpi_e_ndata,fpi_i_vdata(:,1:3),fpi_i_tparadata,fpi_i_tperpdata,fgm_bdata,durationForAverage,'leading');
    
    criterionStructure.l_maxNoverAmbientN(i)= max([l_maxNoverAmbientN_e l_maxNoverAmbientN_i]);
    
    [criterionStructure.t_nBCorr(i),criterionStructure.t_nvCorr(i),criterionStructure.t_vBCorr(i),criterionStructure.t_nTCorr(i),...
        criterionStructure.t_meanB(i),criterionStructure.t_maxB(i),criterionStructure.t_deltaB(i),criterionStructure.t_maxBoverAmbientB(i),...
        criterionStructure.t_meanN(i),criterionStructure.t_maxN(i),criterionStructure.t_deltaN(i),t_maxNoverAmbientN_i,...
        criterionStructure.t_meanV(i),criterionStructure.t_maxV(i),criterionStructure.t_deltaV(i),criterionStructure.t_maxVoverAmbientV(i)]...
        = calculate_edgesAnalysis(event_start,event_end,right_InnerEdge,right_OuterEdge,fpi_i_timedata,fpi_i_ndata,fpi_i_vdata(:,1:3),fpi_i_tparadata,fpi_i_tperpdata,fgm_bdata,durationForAverage,'trailing');
    
    [~,~,~,~,...
        ~,~,~,~,...
        ~,~,~,t_maxNoverAmbientN_e,...
        ~,~,~,~]...
        = calculate_edgesAnalysis(event_start,event_end,right_InnerEdge,right_OuterEdge,fpi_i_timedata,fpi_e_ndata,fpi_i_vdata(:,1:3),fpi_i_tparadata,fpi_i_tperpdata,fgm_bdata,durationForAverage,'trailing');
    
    criterionStructure.t_maxNoverAmbientN(i)= max([t_maxNoverAmbientN_e t_maxNoverAmbientN_i]);
    
    
    %% Save Values to Structure
    
    criterionStructure.eventNumber(i) =  i;
    criterionStructure.coreDuration(i) = Duration;
    
    criterionStructure.coreSpeed(i) =  meanCoreSpeed;
    criterionStructure.coreSpeedx(i) =  meanCoreSpeedx;
    criterionStructure.downSpeed(i) =  downstreamSpeed;
    criterionStructure.upSpeed(i) =  upstreamSpeed;
    criterionStructure.vRatio(i) =  deflectionRatio;
    criterionStructure.vxRatio(i) =  deflectionxRatio;
    criterionStructure.vxMinRatio(i)= deflectionxminRatio;
    criterionStructure.DeflectionAngle(i)= dangle;

    criterionStructure.vDifference(i) =  speedDiff;
    
    criterionStructure.coreDensity(i) =  meanCoreDensity;
    criterionStructure.coreMinDensity(i) =  minCoreDensity;
    criterionStructure.downDensity(i) =  downstreamDensity;
    criterionStructure.upDensity(i) =  upstreamDensity;
    criterionStructure.nRatio(i) =  depletionRatio;
    criterionStructure.nminRatio(i) =  depletionminRatio;
    
    criterionStructure.eTperp(i) =  meanCoreETperp;
    criterionStructure.eTpara(i) =  meanCoreETpara;
    %38 events for 2.0 maxNoverAmbientN
    %67 events for 1.5, 70
    %Criterion Application
    
%     if ((meanCoreETperp >= 15 || meanCoreETpara >= 15) &&...
%             speedDiff >= 135 && depletionRatio <= 0.75 &&...
%             criterionStructure.l_maxNoverAmbientN(i) >= edgeNRatio &&...
%             criterionStructure.t_maxNoverAmbientN(i) >= edgeNRatio)
%         criterionStructure.criteriaPassed(i) = i;
%     end
    
    if ((meanCoreETperp >= 15 || meanCoreETpara >= 15) &&...
            deflectionxminRatio <= 0.5 && depletionRatio <= 0.75)
        criterionStructure.criteriaPassed(i) = i;
    end
    
%         criterionStructure.criteriaPassed= [];
%         for k=1:length(criterionStructure.eventNumber)
%             if (criterionStructure.eTperp(k) >= 15 || criterionStructure.eTpara(k) >= 15) &&  criterionStructure.vDifference(k) >= 140 && criterionStructure.nRatio(k) <= 0.75 && criterionStructure.l_maxNoverAmbientN(k) >= 1.25 && criterionStructure.t_maxNoverAmbientN(k) >= edgeNRatio
%                 criterionStructure.criteriaPassed(k) = criterionStructure.eventNumber(k);
%             end
%         end
        
        %61 25 to 30s, 1.25
%     TypicalEvents(~ismember(TypicalEvents,nonzeros(criterionStructure.criteriaPassed)))
end
size(nonzeros(criterionStructure.criteriaPassed))

%Get rid of nonzeros
S = fieldnames(criterionStructure);
for j=1:length(S)
    criterionStructure.(string(S(j))) = nonzeros(criterionStructure.(string(S(j))));
end

% Plotting
plot_histograms({'Core Duration'},0:5:150,criterionStructure.coreDuration,'All');
plot_histograms({'Core Density Ratio'},0:0.1:2.5,criterionStructure.nRatio,'All');
plot_histograms({'Core Min Density Ratio'},0:0.1:2.5,criterionStructure.nminRatio,'All');
plot_histograms({'Core Speed Ratio'},0:0.05:1.1,criterionStructure.vRatio,'All');
plot_histograms({'Core Vx Ratio'},0:0.05:1.1,criterionStructure.vxRatio, 'All');
plot_histograms({'Core Vx Min Ratio'},0:0.05:1.1,criterionStructure.vxMinRatio, 'All');
plot_histograms({'Velocity Deflection Angle'},0:10:100,criterionStructure.DeflectionAngle, 'All');
plot_histograms({'Speed Difference'},0:50:800,criterionStructure.vDifference,'All');
plot_histograms({'Electron Parallel Temperature'},0:2.5:50,criterionStructure.eTpara,'All');
plot_histograms({'Electron Perpendicular Temperature'},0:2.5:50,criterionStructure.eTperp,'All');


plot_histograms({'V Correlation Coefficient'},-1:0.2:1,criterionStructure.vCorr,'All');
plot_histograms({'Vx Correlation Coefficient'},-1:0.2:1,criterionStructure.vxCorr, 'All');
plot_histograms({'Tpara Correlation Coefficient'},-1:0.2:1,criterionStructure.tempparaCorr, 'All');
plot_histograms({'Tperp Correlation Coefficient'},-1:0.2:1,criterionStructure.tempperpCorr, 'All');
plot_histograms({'T Correlation Coefficient'},-1:0.2:1,criterionStructure.tempCorr, 'All');
plot_histograms({'B Correlation Coefficient'},-1:0.2:1,criterionStructure.BmagCorr, 'All');

plot_histograms({'Event Edges n B Correlation'},-1:0.2:1,criterionStructure.t_nBCorr, 'Trailing Edge',...
    criterionStructure.l_nBCorr, 'Leading Edge'    );
plot_histograms({'Event Edges n V Correlation'},-1:0.2:1,criterionStructure.t_nvCorr, 'Trailing Edge',...
    criterionStructure.l_nvCorr, 'Leading Edge'    );
plot_histograms({'Event Edges V B Correlation'},-1:0.2:1,criterionStructure.t_vBCorr, 'Trailing Edge',...
    criterionStructure.l_vBCorr, 'Leading Edge'    );
plot_histograms({'Event Edges n T Correlation'},-1:0.2:1,criterionStructure.t_nTCorr, 'Trailing Edge',...
    criterionStructure.l_nTCorr, 'Leading Edge'    );
plot_histograms({'Event Edges max B over Ambient B'},0:1:10,criterionStructure.t_maxBoverAmbientB, 'Trailing Edge',...
    criterionStructure.l_maxBoverAmbientB, 'Leading Edge'    );
plot_histograms({'Event Edges max N over Ambient N'},0:1.5:10,criterionStructure.t_maxNoverAmbientN, 'Trailing Edge',...
    criterionStructure.l_maxNoverAmbientN, 'Leading Edge'    );
plot_histograms({'Event Edges max V over Ambient V'},0.8:0.1:2,criterionStructure.t_maxVoverAmbientV, 'Trailing Edge',...
    criterionStructure.l_maxVoverAmbientV, 'Leading Edge'    );

toc
%% Plot the histogram, up to 3 data points
function [] = plot_histograms(parameter,bin_size,data1,label1,data2,label2,data3)
    screenSize = get(0,'ScreenSize');
    figure('Position',[0 screenSize(4) 350 350])
    %     co = [0 0 1;
    %         0 1 0;
    %         1 0 0;
    %         0 0 0];
    %     set(gcf,'defaultAxesColorOrder',co)
    set(gcf,'color','w');
    
    if nargin == 4
        data2=[];
        label2=[];
    end
    
    if size(bin_size) == 1
        
        %     data1=rmoutliers(data1,'ThresholdFactor',25);
        %     data2=rmoutliers(data2,'ThresholdFactor',25);
        if min([data1;data2])  >= 0 && max([data1;data2]) <= 1
            binRange = floor(min([data1;data2])):bin_size:ceil(max([data1;data2]));
        elseif min([data1;data2])  >= -1 && max([data1;data2]) <= 1
            binRange = floor(min([data1;data2])):bin_size:ceil(max([data1;data2]));
        elseif min([data1;data2])  >= 0 && max([data1;data2]) <= 75
            binRange = floor(min([data1;data2])):bin_size:ceil(max([data1;data2]));
        elseif min([data1;data2])  >= 0 && max([data1;data2]) <= 180
            binRange = 0:bin_size:180;
        elseif  min([data1;data2]) >= -90 && max([data1;data2]) <= 90
            binRange = -90:bin_size:90;
        else
            binRange = floor(min([data1;data2])):bin_size:ceil(max([data1;data2]));
        end
        custom_bins=0;
    else
        binRange=bin_size;
        custom_bins = 1;
        bin_size = binRange(3)-binRange(2);
    end
    
    %Normalization
    if nargin == 7
        C = histcounts(data3,[-Inf binRange Inf]);
        if sum(contains(parameter,'ClosestDistance')) > 0
            %2 minutes per data point, per hour is 30 times.
            C = C./30;
        end
        %         C = C./sum(C);
        Cerror = C.^(-1/2);
        Cerror(Cerror == inf) = NaN;
        Cerror(Cerror == 0) = NaN;
        if sum(contains(parameter,'MMS')) > 0
            yAxisLabel = '# of Events normalized to MMS observations';
        else
            yAxisLabel = '# of Events normalized to WIND observations';
        end
    else
        C=1;
        Cerror = 1;
        yAxisLabel = '# of Events / Total # Events';
    end
    
    %weight = C./sum(C);
    if nargin ~= 4 %Wewighting from Spacecraft Data Over All Time Range
        A = histcounts(data1,[-Inf binRange Inf]);
        B = histcounts(data2,[-Inf binRange Inf]);
        
        A_normalized = A./C;
        B_normalized =B./C;
        
        if length(C) == 1
            A_normalized = A./sum(A);
            B_normalized =B./sum(B);
            Aerror = (A).^(-1/2)./sum(A);
            Berror = (B).^(-1/2)./sum(B);
        else
            Aerror = A_normalized.*(A.^(-1) + C.^(-1)).^(1/2);
            Berror = B_normalized.*(B.^(-1) + C.^(-1)).^(1/2);
            
        end
        
        if A(1) == 0 && B(1) == 0
            A = A(2:end);
            A_normalized = A_normalized(2:end);
            
            B = B(2:end);
            B_normalized = B_normalized(2:end);
            
            if length(C) ~= 1
                C = C(2:end);
                Cerror = Cerror(2:end);
                Aerror = Aerror(2:end);
                Berror = Berror(2:end);
            else
                Aerror = Aerror(2:end);
                Berror = Berror(2:end);
            end
        else
            binRange = [binRange(1) - bin_size, binRange];
        end
        xlabels = [binRange(2)-(binRange(3)-binRange(2)),binRange(2:end-1),binRange(end-1)+(binRange(3)-binRange(2))];
        
        bar(xlabels+bin_size/2,A_normalized',1.0,'edgecolor','none'); hold on
        [XB,YB] = stairs(xlabels,B_normalized');
        XB = [XB;XB(end)+bin_size];
        YB = [YB;YB(end)];
        plot(XB,YB,'linewidth',3.5);
        
        %         stairs(xlabels,B_normalized','linewidth',3.5);
        
        if A(end) ~= 0 || B(end) ~= 0
            xticks([binRange binRange(end)+(binRange(3)-binRange(2))])
            xlim([binRange(1) binRange(end)+(binRange(3)-binRange(2))])
        else
            xticks(binRange)
            xlim([binRange(1) binRange(end)])
        end
        
        %Labels
        if custom_bins == 1 && xlabels(1) ~= 0 && (max(data1) >= binRange(end) || max(data2) >= binRange(end) )
            xlabelCell = xticklabels;
            xlabelCell(1) = strcat('\leq',xlabelCell(1));
            xlabelCell(end) = strcat('\geq',xlabelCell(end));
            xticklabels(xlabelCell);
        elseif custom_bins == 1 && (max(data1) >= binRange(end) || max(data2) >= binRange(end) )
            xlabelCell = xticklabels;
            xlabelCell(end) = strcat('\geq',xlabelCell(end));
            xticklabels(xlabelCell);
        end
        
        [maxA,Ia] = max(A_normalized(A_normalized~=inf));
        [maxB,Ib] = max(B_normalized(B_normalized~=inf));
        if maxA > maxB
            I=Ia;
        else
            I = Ib;
        end
        
        if I > length(binRange)/2
            legend(label1,label2,'FontSize',9,'Location','Northwest')
        else
            legend(label1,label2,'FontSize',9,'Location','Northeast')
        end
        
        hold on
        
        
        errorbar(xlabels+8*bin_size/16,A_normalized,Aerror,'LineStyle','none','Color','black','LineWidth',1.5,'HandleVisibility','off')
        
        errorbar(xlabels+8*bin_size/16,B_normalized,Berror,'LineStyle','none','Color','black','LineWidth',1.5,'HandleVisibility','off')
        
        
        hold off
    else
        A = histcounts(data1,[binRange Inf]);
        A_normalized = A./sum(A);
        Aerror = A.^(-1/2)./sum(A);
        xlabels = [binRange(2)-(binRange(3)-binRange(2)),binRange(2:end-1),binRange(end-1)+(binRange(3)-binRange(2))];
        bar(xlabels+bin_size/2,A_normalized',1,'edgecolor','none'); hold on
        errorbar(xlabels+8*bin_size/16,A_normalized,Aerror,'LineStyle','none','Color','black','LineWidth',1.5,'HandleVisibility','off')
        xticks(binRange)
        
        %Labels
        if custom_bins == 1 && xlabels(1) ~= 0
            xlabelCell = xticklabels;
            xlabelCell(1) = strcat('\leq',xlabelCell(1));
            xlabelCell(end) = strcat('\geq',xlabelCell(end));
            xticklabels(xlabelCell);
        elseif custom_bins == 1
            xlabelCell = xticklabels;
            xlabelCell(end) = strcat('\geq',xlabelCell(end));
            xticklabels(xlabelCell);
        end
        
        %Legend
        %Location of Max
        [~,I] = max(A(A~=inf));
        if I > length(binRange)/2
            legend(label1,'FontSize',9,'Location','Northwest')
        else
            legend(label1,'FontSize',9,'Location','Northeast')
        end
        hold on
        hold off
        
    end
    ylimits = ylim;
    if ylimits(1) < 0
        ylim([0 ylimits(2)])
        %         xlim([binRange(1) binRange(end)])
    end
    
    fileName = strcat(cell2mat(parameter(1)),'_',label1);
    colormap('winter');
    xlabel(parameter,'FontSize',11)
    ylabel(yAxisLabel,'FontSize',11)
    title(strcat('f vs.', {' '}, parameter(1)),'FontSize',12,'FontWeight', 'normal')
    set(gca,'XMinorTick','off','TickDir','out','YMinorTick','on','linewidth',2)
    
    print(gcf,'-dpng','-r300', '-loose', strcat(fileName));
    savefig(strcat(fileName));
    
end
function [maximum,minimum] = calculate_minCore(fpi_timedata,fpi_ndata,left_InnerEdge,right_InnerEdge)
    %Calculate max or min.
    formatIn='yyyy-mm-dd HH:MM:SS.FFF';
    
    %Core Time limits
    coreStart = datenum(left_InnerEdge,formatIn);
    coreEnd = datenum(right_InnerEdge,formatIn);
    
    coreStart_index = find(fpi_timedata > coreStart, 1);
    coreEnd_index = find(fpi_timedata > coreEnd, 1);
    
    fpi_coreNdata = fpi_ndata(coreStart_index:coreEnd_index,1);
    
    minimum = min(fpi_coreNdata);
    maximum = max(fpi_coreNdata);
end

function [dangle,dratio] = calculate_vDeflection(fpi_timedata,fpi_vdata,event_start,left_OuterEdge,left_InnerEdge,right_InnerEdge,right_OuterEdge,event_end)
    
    %fpi_vdata(:,4) = vecnorm(fpi_vdata(:,1:3),2,2);
    
        %Calculate max or min.
    formatIn='yyyy-mm-dd HH:MM:SS.FFF';
    
    % Time limits
    coreStart = datenum(left_InnerEdge,formatIn);
    coreStart_index = find(fpi_timedata > coreStart, 1);
    coreEnd = datenum(right_InnerEdge,formatIn);
    coreEnd_index = find(fpi_timedata > coreEnd, 1);
    
    leading_start = datenum(left_OuterEdge,formatIn);
    leading_start_index = find(fpi_timedata > leading_start, 1);
    trailing_end = datenum(right_OuterEdge,formatIn);
    trailing_end_index = find(fpi_timedata > trailing_end, 1);
    
    eventStart = datenum(event_start,formatIn);
    eventStart_index = find(fpi_timedata > eventStart, 1);
    eventEnd = datenum(event_end,formatIn);
    eventEnd_index = find(fpi_timedata > eventEnd, 1);
    
    coreV = fpi_vdata(coreStart_index:coreEnd_index,:);
    [~,minCoreIndex] = min(abs(coreV(:,1)));
    mincoreV = coreV(minCoreIndex,1:3);
    
    preVmean = mean(fpi_vdata(eventStart_index:leading_start_index,1:3));
    postVmean = mean(fpi_vdata(trailing_end_index:eventEnd_index,1:3));
    
    preVminVangle = angle(mincoreV,preVmean);
    postVminVangle = angle(mincoreV,postVmean);
    
    dangle = max([preVminVangle,postVminVangle]);
    
    dratio = min(mincoreV(:,1)./([preVmean(:,1),postVmean(:,1)]),[],'all');
    
end