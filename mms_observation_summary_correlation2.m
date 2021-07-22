%Matlab CDF Plotter for MMS
%Andrew Vu 9/19/18
%data=spdfcdfread(filename);
%datainfo=spdfcdfinfo(filename);
% clear; close all
warning('off','MATLAB:polyfit:RepeatedPointsOrRescale')
cd '~/Library/Mobile Documents/com~apple~CloudDocs/Research/Analysis/Correlation2/'
mms_directory = '/Users/andrewvu/data/mms/';
num_plots = 9;
data_type = 'brst';
plot_gap=1.25;
probe_num = '1';

ssSTDs = 2.00;
minSSDataPoints = 5; %Changed up from 4, 9172018
durationForAverage = 4;
Events = 2:174;
Events = [2;3;5;6;7;8;9;11;12;14;15;16;17;18;20;21;22;24;25;28;29;30;31;32;34;35;36;37;38;39;41;42;43;47;48;50;51;52;54;55;56;57;58;60;61;65;66;67;68;76;77;78;79;80;81;82;83;84;85;87;93;94;97;98;99;100;104;105;106;111;116;117;121;122;123;124;127;130;132;133;135;136;138;142;152;153;154;157;158;159;160;161;162;164;165;166;167;168;169;170];

% Without Data
%% Create Structures
eventStructure = struct('eventNumber', zeros(size(Events)),...
    'densityCV',zeros(size(Events)),...
    'meanDensity',zeros(size(Events)),...
    'relativeMaxDensity',zeros(size(Events)),...
    'Shear',zeros(size(Events)),...
    'vCorr', zeros(size(Events)), ...
    'vxCorr', zeros(size(Events)), ...
    'tempparaCorr', zeros(size(Events)), ...
    'tempperpCorr', zeros(size(Events)), ...
    'tempCorr', zeros(size(Events)),...
    'BmagCorr',zeros(size(Events)),...
    'massFluxCV',zeros(size(Events)),...
    'l_nBCorr',zeros(size(Events)),...
    'l_nvCorr',zeros(size(Events)),...
    'l_vBCorr',zeros(size(Events)),...
    'l_nTCorr',zeros(size(Events)),...
    'l_meanB',zeros(size(Events)),...
    'l_maxB',zeros(size(Events)),...
    'l_deltaB',zeros(size(Events)),...
    'l_maxBoverAmbientB',zeros(size(Events)),...
    'l_meanN',zeros(size(Events)),...
    'l_maxN',zeros(size(Events)),...
    'l_deltaN',zeros(size(Events)),...
    'l_maxNoverAmbientN',zeros(size(Events)),...
    'l_meanV',zeros(size(Events)),...
    'l_maxV',zeros(size(Events)),...
    'l_deltaV',zeros(size(Events)),...
    'l_maxVoverAmbientV',zeros(size(Events)),...
    't_nBCorr',zeros(size(Events)),...
    't_nvCorr',zeros(size(Events)),...
    't_vBCorr',zeros(size(Events)),...
    't_nTCorr',zeros(size(Events)),...
    't_meanB',zeros(size(Events)),...
    't_maxB',zeros(size(Events)),...
    't_deltaB',zeros(size(Events)),...
    't_maxBoverAmbientB',zeros(size(Events)),...
    't_meanN',zeros(size(Events)),...
    't_maxN',zeros(size(Events)),...
    't_deltaN',zeros(size(Events)),...
    't_maxNoverAmbientN',zeros(size(Events)),...
    't_meanV',zeros(size(Events)),...
    't_maxV',zeros(size(Events)),...
    't_deltaV',zeros(size(Events)),...
    't_maxVoverAmbientV',zeros(size(Events)));

eventSub = struct('eventNumber', zeros(size(Events)),...
    'densityCV',zeros(size(Events)),...
    'meanDensity',zeros(size(Events)),...
    'relativeMaxDensity',zeros(size(Events)),...
    'Shear',zeros(size(Events)),...
    'vCorr', zeros(size(Events)), ...
    'vxCorr', zeros(size(Events)), ...
    'tempparaCorr', zeros(size(Events)), ...
    'tempperpCorr', zeros(size(Events)), ...
    'tempCorr', zeros(size(Events)),...
    'BmagCorr',zeros(size(Events)),...
    'massFluxCV',zeros(size(Events)),...
    'l_nBCorr',zeros(size(Events)),...
    'l_nvCorr',zeros(size(Events)),...
    'l_vBCorr',zeros(size(Events)),...
    'l_nTCorr',zeros(size(Events)),...
    'l_meanB',zeros(size(Events)),...
    'l_maxB',zeros(size(Events)),...
    'l_deltaB',zeros(size(Events)),...
    'l_maxBoverAmbientB',zeros(size(Events)),...
    'l_meanN',zeros(size(Events)),...
    'l_maxN',zeros(size(Events)),...
    'l_deltaN',zeros(size(Events)),...
    'l_maxNoverAmbientN',zeros(size(Events)),...
    'l_meanV',zeros(size(Events)),...
    'l_maxV',zeros(size(Events)),...
    'l_deltaV',zeros(size(Events)),...
    'l_maxVoverAmbientV',zeros(size(Events)),...
    't_nBCorr',zeros(size(Events)),...
    't_nvCorr',zeros(size(Events)),...
    't_vBCorr',zeros(size(Events)),...
    't_nTCorr',zeros(size(Events)),...
    't_meanB',zeros(size(Events)),...
    't_maxB',zeros(size(Events)),...
    't_deltaB',zeros(size(Events)),...
    't_maxBoverAmbientB',zeros(size(Events)),...
    't_meanN',zeros(size(Events)),...
    't_maxN',zeros(size(Events)),...
    't_deltaN',zeros(size(Events)),...
    't_maxNoverAmbientN',zeros(size(Events)),...
    't_meanV',zeros(size(Events)),...
    't_maxV',zeros(size(Events)),...
    't_deltaV',zeros(size(Events)),...
    't_maxVoverAmbientV',zeros(size(Events)));

eventNoSub = struct('eventNumber', zeros(size(Events)),...
    'densityCV',zeros(size(Events)),...
    'meanDensity',zeros(size(Events)),...
    'relativeMaxDensity',zeros(size(Events)),...
    'Shear',zeros(size(Events)),...
    'vCorr', zeros(size(Events)), ...
    'vxCorr', zeros(size(Events)), ...
    'tempparaCorr', zeros(size(Events)), ...
    'tempperpCorr', zeros(size(Events)), ...
    'tempCorr', zeros(size(Events)),...
    'BmagCorr',zeros(size(Events)),...
    'massFluxCV',zeros(size(Events)),...
    'l_nBCorr',zeros(size(Events)),...
    'l_nvCorr',zeros(size(Events)),...
    'l_vBCorr',zeros(size(Events)),...
    'l_nTCorr',zeros(size(Events)),...
    'l_meanB',zeros(size(Events)),...
    'l_maxB',zeros(size(Events)),...
    'l_deltaB',zeros(size(Events)),...
    'l_maxBoverAmbientB',zeros(size(Events)),...
    'l_meanN',zeros(size(Events)),...
    'l_maxN',zeros(size(Events)),...
    'l_deltaN',zeros(size(Events)),...
    'l_maxNoverAmbientN',zeros(size(Events)),...
    'l_meanV',zeros(size(Events)),...
    'l_maxV',zeros(size(Events)),...
    'l_deltaV',zeros(size(Events)),...
    'l_maxVoverAmbientV',zeros(size(Events)),...
    't_nBCorr',zeros(size(Events)),...
    't_nvCorr',zeros(size(Events)),...
    't_vBCorr',zeros(size(Events)),...
    't_nTCorr',zeros(size(Events)),...
    't_meanB',zeros(size(Events)),...
    't_maxB',zeros(size(Events)),...
    't_deltaB',zeros(size(Events)),...
    't_maxBoverAmbientB',zeros(size(Events)),...
    't_meanN',zeros(size(Events)),...
    't_maxN',zeros(size(Events)),...
    't_deltaN',zeros(size(Events)),...
    't_maxNoverAmbientN',zeros(size(Events)),...
    't_meanV',zeros(size(Events)),...
    't_maxV',zeros(size(Events)),...
    't_deltaV',zeros(size(Events)),...
    't_maxVoverAmbientV',zeros(size(Events)));

%% Loop for each Event
for i=1:length(Events)
    
    close all
    figure('Position',[1 1 650 850])
    co = [0 0 1;
        0 1 0;
        1 0 0;
        0 0 0];
    set(gcf,'defaultAxesColorOrder',co)
    set(gcf,'color','w');
    
    Event_number = Events(i)
    
    [Event_Type, ~,threshold_std,~,...
        event_start,event_end,...
        left_OuterEdge, left_InnerEdge,...
        right_InnerEdge, right_OuterEdge,...
        leading_leftmost_date, leading_rightmost_date,...
        trailing_leftmost_date, trailing_rightmost_date,...
        ~, ~,...
        ~, ~] = get_eventTimes(Event_number);
    
    
    formatIn='yyyy-mm-dd HH:MM:SS.FFF';
    tstart = datenum(event_start,formatIn);
    tend = datenum(event_end,formatIn);
    
    %% Load Data
    %%%%%%%%%%%%%%%%%%%%%%%%%%%Position Data%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %load mec data
    [mec_timedata,mec_r_gsedata] = load_mec(event_start,probe_num,'srvy');
    [mec_timedata,mec_r_gsedata,~,~] = crop(mec_timedata,mec_r_gsedata,event_start,event_end);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%Load Data%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %load fgm data
    [fgm_timedata, fgm_bdata,~,~]= load_fgm(event_start,event_end,probe_num,data_type);
    [fgm_timedata_srvy, fgm_bdata_srvy,~,~]= load_fgm(event_start,event_end,probe_num,'srvy');
    
    %Load FPI_e
    [fpi_e_timedata,fpi_e_ndata,fpi_e_vdata,fpi_e_tparadata,fpi_e_tperpdata,...
        fpi_e_edata,fpi_e_espectdata] = load_fpi(event_start,event_end,probe_num,data_type,'e');
    %Load FPI_i
    [fpi_i_timedata,fpi_i_ndata,fpi_i_vdata,fpi_i_tparadata,fpi_i_tperpdata,...
        fpi_i_edata,fpi_i_espectdata] = load_fpi(event_start,event_end,probe_num,data_type,'i');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %Interpolate FGM to FPI
    [~,fgm_bdata] = interpxyz(fgm_timedata,fgm_bdata(:,1:3),fpi_i_timedata);
    
    
    
    %% binary for sigmaSigma threshold
    [substructurePresent,~,fpi_coreNdata,~,~,~,~] = calculate_densityStats(fpi_i_timedata,fpi_i_ndata,left_InnerEdge,right_InnerEdge,ssSTDs,minSSDataPoints);
    [~,coreDensity_CV] = calculate_coreDensitySTD(left_InnerEdge,right_InnerEdge,fpi_i_timedata,fpi_i_ndata);
    [~,coreMassFlux_CV] = calculate_coreDensitySTD(left_InnerEdge,right_InnerEdge,fpi_i_timedata,fpi_i_ndata.*vecnorm(fpi_i_vdata,2,2));
    
    %% Storing Values into Structures
    
    if threshold_std == 0
        %Current Sheet Normal Calculation Plot, manually.
        [~,B_pre,B_post] = manualCurrentSheet(event_start,event_end,leading_leftmost_date,leading_rightmost_date,...
            trailing_leftmost_date,trailing_rightmost_date,fgm_timedata_srvy,fgm_bdata_srvy);
    else
        %Calculate Current Sheet
        [~,B_pre,B_post] = calculateCurrentSheet(event_start,event_end,fgm_timedata_srvy,fgm_bdata_srvy,threshold_std);
    end
    
    substructurePresent
    %     if contains(Event_Type,'HFA')
    eventStructure.Shear(i) = angle(B_pre,B_post);
    if substructurePresent == 0
        eventStructure.eventNumber(i) = Event_number;
        eventNoSub.eventNumber(i) = Event_number;
        
        eventStructure.densityCV(i) = coreDensity_CV;
        eventNoSub.densityCV(i) = coreDensity_CV;
        
        eventStructure.massFluxCV(i) = coreMassFlux_CV;
        eventNoSub.massFluxCV(i) = coreMassFlux_CV;
        
        eventStructure.meanDensity(i) = mean(fpi_coreNdata);
        eventStructure.relativeMaxDensity(i) = abs(mean(fpi_coreNdata)-max(fpi_coreNdata))/mean(fpi_coreNdata);
        eventNoSub.meanDensity(i) = mean(fpi_coreNdata);
        eventNoSub.relativeMaxDensity(i) = abs(mean(fpi_coreNdata)-max(fpi_coreNdata))/mean(fpi_coreNdata);
        eventNoSub.Shear(i) = angle(B_pre,B_post);
        
        Substructure = 0;
        
    elseif substructurePresent == 1
        eventStructure.eventNumber(i) = Event_number;
        eventSub.eventNumber(i) = Event_number;
        
        eventStructure.densityCV(i) = coreDensity_CV;
        eventSub.densityCV(i) = coreDensity_CV;
        
        eventStructure.massFluxCV(i) = coreMassFlux_CV;
        eventSub.massFluxCV(i) = coreMassFlux_CV;
        
        eventStructure.meanDensity(i) = mean(fpi_coreNdata);
        eventStructure.relativeMaxDensity(i) = abs(mean(fpi_coreNdata)-max(fpi_coreNdata))/mean(fpi_coreNdata);
        eventSub.meanDensity(i) = mean(fpi_coreNdata);
        eventSub.relativeMaxDensity(i) = abs(mean(fpi_coreNdata)-max(fpi_coreNdata))/mean(fpi_coreNdata);
        eventSub.Shear(i) = angle(B_pre,B_post);
        Substructure = 1;
        
    end
    
    %     else
    %         %If it is not an HFA, mve on.
    %         continue
    %     end
    
    
    %Correlation Coefficients
    [eventStructure.vCorr(i),eventStructure.vxCorr(i),...
        eventStructure.tempparaCorr(i),eventStructure.tempperpCorr(i),...
        eventStructure.tempCorr(i),...
        eventStructure.BmagCorr(i)] = calculate_correlation(left_InnerEdge,right_InnerEdge,...
        fpi_i_timedata,fpi_i_ndata,fpi_i_vdata,fpi_i_tparadata,fpi_i_tperpdata,fgm_bdata);
    
    %Event Edges Analysis
    [eventStructure.l_nBCorr(i),eventStructure.l_nvCorr(i),eventStructure.l_vBCorr(i),eventStructure.l_nTCorr(i),...
        eventStructure.l_meanB(i),eventStructure.l_maxB(i),eventStructure.l_deltaB(i),eventStructure.l_maxBoverAmbientB(i),...
        eventStructure.l_meanN(i),eventStructure.l_maxN(i),eventStructure.l_deltaN(i),eventStructure.l_maxNoverAmbientN(i),...
        eventStructure.l_meanV(i),eventStructure.l_maxV(i),eventStructure.l_deltaV(i),eventStructure.l_maxVoverAmbientV(i)]...
        = calculate_edgesAnalysis(event_start,event_end,left_OuterEdge,left_InnerEdge,fpi_i_timedata,fpi_i_ndata,fpi_i_vdata,fpi_i_tparadata,fpi_i_tperpdata,fgm_bdata,durationForAverage,'leading');
    
    [eventStructure.t_nBCorr(i),eventStructure.t_nvCorr(i),eventStructure.t_vBCorr(i),eventStructure.t_nTCorr(i),...
        eventStructure.t_meanB(i),eventStructure.t_maxB(i),eventStructure.t_deltaB(i),eventStructure.t_maxBoverAmbientB(i),...
        eventStructure.t_meanN(i),eventStructure.t_maxN(i),eventStructure.t_deltaN(i),eventStructure.t_maxNoverAmbientN(i),...
        eventStructure.t_meanV(i),eventStructure.t_maxV(i),eventStructure.t_deltaV(i),eventStructure.t_maxVoverAmbientV(i)]...
        = calculate_edgesAnalysis(event_start,event_end,right_InnerEdge,right_OuterEdge,fpi_i_timedata,fpi_i_ndata,fpi_i_vdata,fpi_i_tparadata,fpi_i_tperpdata,fgm_bdata,durationForAverage,'trailing');
    
    
    if substructurePresent == 0
        S = fieldnames(eventStructure);
        for j=1:length(S)
            eventNoSub.(string(S(j)))(i) = eventStructure.(string(S(j)))(i);
        end
        
        
    elseif substructurePresent == 1
        S = fieldnames(eventStructure);
        for j=1:length(S)
            eventSub.(string(S(j)))(i) = eventStructure.(string(S(j)))(i);
        end
    end
    
    
end

S = fieldnames(eventStructure);
for j=1:length(S)
    A = eventStructure.(string(S(j)));
    A(A==0) = [];
    eventStructure.(string(S(j)))= A;
    
    B = eventNoSub.(string(S(j)));
    B(B==0) = [];
    eventNoSub.(string(S(j))) = B;
    
    C  = eventSub.(string(S(j)));
    C(C==0) = [];
    eventSub.(string(S(j))) = C;
    
end
%% Plotting
% %Get rid of nonzeros because it isn't an HFA
% eventStructure.densitySTD(eventStructure.densitySTD==0) = [];
% eventStructure.meanDensity(eventStructure.meanDensity==0) = [];
% eventStructure.relativeMaxDensity(eventStructure.relativeMaxDensity==0) = [];
% eventStructure.Shear(eventStructure.Shear==0) = [];
% eventStructure.vCorr(eventStructure.vCorr==0) = [];
% eventStructure.vxCorr(eventStructure.vxCorr==0) = [];
% eventStructure.tempparaCorr(eventStructure.tempparaCorr==0) = [];
% eventStructure.tempperpCorr(eventStructure.tempperpCorr==0) = [];
% eventStructure.tempCorr(eventStructure.tempCorr==0) = [];
% eventStructure.BmagCorr(eventStructure.BmagCorr==0) = [];
%
% eventNoSub.densitySTD(eventNoSub.densitySTD==0) = [];
% eventNoSub.meanDensity(eventNoSub.meanDensity==0) = [];
% eventNoSub.relativeMaxDensity(eventNoSub.relativeMaxDensity==0) = [];
% eventNoSub.Shear(eventNoSub.Shear==0) = [];
% eventNoSub.vCorr(eventNoSub.vCorr==0) = [];
% eventNoSub.vxCorr(eventNoSub.vxCorr==0) = [];
% eventNoSub.tempparaCorr(eventNoSub.tempparaCorr==0) = [];
% eventNoSub.tempperpCorr(eventNoSub.tempperpCorr==0) = [];
% eventNoSub.tempCorr(eventNoSub.tempCorr==0) = [];
% eventNoSub.BmagCorr(eventNoSub.BmagCorr==0) = [];
%
% eventSub.densitySTD(eventSub.densitySTD==0) = [];
% eventSub.meanDensity(eventSub.meanDensity==0) = [];
% eventSub.relativeMaxDensity(eventSub.relativeMaxDensity==0) = [];
% eventSub.Shear(eventSub.Shear==0) = [];
% eventSub.vCorr(eventSub.vCorr==0) = [];
% eventSub.vxCorr(eventSub.vxCorr==0) = [];
% eventSub.tempparaCorr(eventSub.tempparaCorr==0) = [];
% eventSub.tempperpCorr(eventSub.tempperpCorr==0) = [];
% eventSub.tempCorr(eventSub.tempCorr==0) = [];
% eventSub.BmagCorr(eventSub.BmagCorr==0) = [];
save('correlationStructures',...
    'eventStructure','eventSub','eventNoSub');
%% Core Plots based on CV
%All Events
plot_histograms({'Core Density CV'},0.2,eventStructure.densityCV,'All');
plot_histograms({'Core Mass Flux CV'},0:0.1:1.2,eventStructure.massFluxCV,'All');
plot_histograms({'V Correlation Coefficient'},-1:0.2:1,eventStructure.vCorr,'All');
plot_histograms({'Vx Correlation Coefficient'},-1:0.2:1,eventStructure.vxCorr, 'All');
plot_histograms({'Tpara Correlation Coefficient'},-1:0.2:1,eventStructure.tempparaCorr, 'All');
plot_histograms({'Tperp Correlation Coefficient'},-1:0.2:1,eventStructure.tempperpCorr, 'All');
plot_histograms({'T Correlation Coefficient'},-1:0.2:1,eventStructure.tempCorr, 'All');
plot_histograms({'B Correlation Coefficient'},-1:0.2:1,eventStructure.BmagCorr, 'All');

%Classifying by CV
indicesForVaried = eventStructure.densityCV > 0.35;

plot_histograms({'V Correlation Coefficient'},-1:0.25:1,eventStructure.vCorr(~indicesForVaried),'Low CV',...
    eventStructure.vCorr(indicesForVaried),'High CV' );
plot_histograms({'Vx Correlation Coefficient'},-1:0.25:1,eventStructure.vxCorr(~indicesForVaried), 'Low CV',...
    eventStructure.vxCorr(indicesForVaried),'High CV' );
plot_histograms({'Tpara Correlation Coefficient'},-1:0.25:1,eventStructure.tempparaCorr(~indicesForVaried), 'Low CV',...
    eventStructure.tempparaCorr(indicesForVaried),'High CV' );
plot_histograms({'Tperp Correlation Coefficient'},-1:0.25:1,eventStructure.tempperpCorr(~indicesForVaried), 'Low CV',...
    eventStructure.tempperpCorr(indicesForVaried),'High CV' );
plot_histograms({'T Correlation Coefficient'},-1:0.25:1,eventStructure.tempCorr(~indicesForVaried), 'Low CV',...
    eventStructure.tempCorr(indicesForVaried),'High CV' );
plot_histograms({'B Correlation Coefficient'},-1:0.25:1,eventStructure.BmagCorr(~indicesForVaried), 'Low CV',...
    eventStructure.BmagCorr(indicesForVaried),'High CV'    );


%Classifying by CV
indicesForVaried = eventStructure.massFluxCV > 0.3;

plot_histograms({'Mass Flux CV'},0:0.2:2,eventStructure.massFluxCV(~indicesForVaried),'Low CV',...
    eventStructure.massFluxCV(indicesForVaried),'High CV' );
plot_histograms({'Mass Flux and V Correlation Coefficient'},-1:0.2:1,eventStructure.vCorr(~indicesForVaried),'Low CV',...
    eventStructure.vCorr(indicesForVaried),'High CV' );
plot_histograms({'Mass Flux and Vx Correlation Coefficient'},-1:0.2:1,eventStructure.vxCorr(~indicesForVaried), 'Low CV',...
    eventStructure.vxCorr(indicesForVaried),'High CV' );
plot_histograms({'Mass Flux and Tpara Correlation Coefficient'},-1:0.2:1,eventStructure.tempparaCorr(~indicesForVaried), 'Low CV',...
    eventStructure.tempparaCorr(indicesForVaried),'High CV' );
plot_histograms({'Mass Flux and Tperp Correlation Coefficient'},-1:0.2:1,eventStructure.tempperpCorr(~indicesForVaried), 'Low CV',...
    eventStructure.tempperpCorr(indicesForVaried),'High CV' );
plot_histograms({'Mass Flux and T Correlation Coefficient'},-1:0.2:1,eventStructure.tempCorr(~indicesForVaried), 'Low CV',...
    eventStructure.tempCorr(indicesForVaried),'High CV' );
plot_histograms({'Mass Flux and B Correlation Coefficient'},-1:0.2:1,eventStructure.BmagCorr(~indicesForVaried), 'Low CV',...
    eventStructure.BmagCorr(indicesForVaried),'High CV'    );

%Classify by SS or NS
plot_histograms({'V Correlation Coefficient'},-1:0.4:1,eventSub.vCorr,'Substructure',...
    eventNoSub.vCorr,'No Substructure');
plot_histograms({'Vx Correlation Coefficient'},-1:0.4:1,eventSub.vxCorr, 'Substructure',...
    eventNoSub.vxCorr,'No Substructure');
plot_histograms({'Tpara Correlation Coefficient'},-1:0.25:1,eventSub.tempparaCorr, 'Substructure',...
    eventNoSub.tempparaCorr,'No Substructure' );
plot_histograms({'Tperp Correlation Coefficient'},-1:0.25:1,eventSub.tempperpCorr, 'Substructure',...
    eventNoSub.tempperpCorr,'No Substructure' );
plot_histograms({'T Correlation Coefficient'},-1:0.25:1,eventSub.tempCorr, 'Substructure',...
    eventNoSub.tempCorr,'No Substructure');
plot_histograms({'B Correlation Coefficient'},-1:0.4:1,eventSub.BmagCorr, 'Substructure',...
    eventNoSub.BmagCorr,'No Substructure');

%% Edges Plots -
plot_histograms({'Leading Edge n B Correlation'},-1:0.4:1,...
    eventSub.l_nBCorr, 'Structure',...
    eventNoSub.l_nBCorr, 'No Structure');
plot_histograms({'Leading Edge n V Correlation'},-1:0.4:1,...
    eventSub.l_nvCorr, 'Structure',...
    eventNoSub.l_nBCorr, 'No Structure');
plot_histograms({'Leading Edge V B Correlation'},-1:0.4:1,...
    eventSub.l_vBCorr, 'Structure',...
    eventNoSub.l_vBCorr, 'No Structure');
plot_histograms({'Leading Edge n T Correlation'},-1:0.4:1,...
    eventSub.l_nTCorr, 'Structure',...
    eventNoSub.l_nTCorr, 'No Structure');
plot_histograms({'Leading Edge max B over Ambient B'},0:0.5:5,...
    eventSub.l_maxBoverAmbientB, 'Structure',...
    eventNoSub.l_maxBoverAmbientB, 'No Structure');
plot_histograms({'Leading Edge max N over Ambient N'},0:0.5:5.5,...
    eventSub.l_maxNoverAmbientN, 'Structure',...
    eventNoSub.l_maxNoverAmbientN, 'No Structure');
plot_histograms({'Leading Edge max V over Ambient V'},0.5:0.2:1.5,...
    eventSub.l_maxVoverAmbientV, 'Structure',...
    eventNoSub.l_maxVoverAmbientV, 'No Structure');
%
plot_histograms({'Trailing Edge n B Correlation'},-1:0.4:1,...
    eventSub.t_nBCorr, 'Structure',...
    eventNoSub.t_nBCorr, 'No Structure');
plot_histograms({'Trailing Edge n V Correlation'},-1:0.4:1,...
    eventSub.t_nvCorr, 'Structure',...
    eventNoSub.t_nvCorr, 'No Structure');
plot_histograms({'Trailing Edge V B Correlation'},-1:0.4:1,...
    eventSub.t_vBCorr, 'Structure',...
    eventNoSub.t_vBCorr, 'No Structure');
plot_histograms({'Trailing Edge n T Correlation'},-1:0.4:1,...
    eventSub.t_nTCorr, 'Structure',...
    eventNoSub.t_nTCorr, 'No Structure');
plot_histograms({'Trailing Edge max B over Ambient B'},0:1:7,...
    eventSub.t_maxBoverAmbientB, 'Structure',...
    eventNoSub.t_maxBoverAmbientB, 'No Structure');
plot_histograms({'Trailing Edge max N over Ambient N'},0:1:7,...
    eventSub.t_maxNoverAmbientN, 'Structure',...
    eventNoSub.t_maxNoverAmbientN, 'No Structure');
plot_histograms({'Trailing Edge max V over Ambient V'},0.5:0.2:1.5,...
    eventSub.t_maxVoverAmbientV, 'Structure',...
    eventNoSub.t_maxVoverAmbientV, 'No Structure');




plot_histograms({'Event Edges n B Correlation'},-1:0.2:1,eventStructure.t_nBCorr, 'Trailing Edge',...
    eventStructure.l_nBCorr, 'Leading Edge'    );
plot_histograms({'Event Edges n V Correlation'},-1:0.2:1,eventStructure.t_nvCorr, 'Trailing Edge',...
    eventStructure.l_nvCorr, 'Leading Edge'    );
plot_histograms({'Event Edges V B Correlation'},-1:0.2:1,eventStructure.t_vBCorr, 'Trailing Edge',...
    eventStructure.l_vBCorr, 'Leading Edge'    );
plot_histograms({'Event Edges n T Correlation'},-1:0.2:1,eventStructure.t_nTCorr, 'Trailing Edge',...
    eventStructure.l_nTCorr, 'Leading Edge'    );
plot_histograms({'Event Edges max B over Ambient B'},0:1:10,eventStructure.t_maxBoverAmbientB, 'Trailing Edge',...
    eventStructure.l_maxBoverAmbientB, 'Leading Edge'    );
plot_histograms({'Event Edges max N over Ambient N'},0:1:10,eventStructure.t_maxNoverAmbientN, 'Trailing Edge',...
    eventStructure.l_maxNoverAmbientN, 'Leading Edge'    );
plot_histograms({'Event Edges max V over Ambient V'},0.8:0.1:2,eventStructure.t_maxVoverAmbientV, 'Trailing Edge',...
    eventStructure.l_maxVoverAmbientV, 'Leading Edge'    );


% %plot Histograms for all events
% plot_scatter({'Shear Angle';'V Correlation Coefficient'},eventStructure.Shear,eventStructure.vCorr, 'All');
% plot_scatter({'Shear Angle';'V_x Correlation Coefficient'},eventStructure.Shear,eventStructure.vxCorr, 'All');
% plot_scatter({'Shear Angle';'T_{para} Correlation Coefficient'},eventStructure.Shear,eventStructure.tempparaCorr, 'All');
% plot_scatter({'Shear Angle';'T_{perp} Correlation Coefficient'},eventStructure.Shear,eventStructure.tempperpCorr, 'All');
% plot_scatter({'Shear Angle';'T Correlation Coefficient'},eventStructure.Shear,eventStructure.tempCorr, 'All');
% plot_scatter({'Shear Angle';'B Correlation Coefficient'},eventStructure.Shear,eventStructure.BmagCorr, 'All');
%
% %plot Histograms for SS and noSS
% plot_histograms({'V Correlation Coefficient'},-1:0.2:1,eventSub.vCorr,'Sub',eventNoSub.vCorr,'NoSub');
% plot_histograms({'Vx Correlation Coefficient'},-1:0.2:1,eventSub.vxCorr, 'Sub',eventNoSub.vxCorr,'NoSub');
% plot_histograms({'Tpara Correlation Coefficient'},-1:0.2:1,eventSub.tempparaCorr, 'Sub',eventNoSub.tempparaCorr,'NoSub');
% plot_histograms({'Tperp Correlation Coefficient'},-1:0.2:1,eventSub.tempperpCorr, 'Sub',eventNoSub.tempperpCorr,'NoSub');
% plot_histograms({'T Correlation Coefficient'},-1:0.2:1,eventSub.tempCorr, 'Sub',eventNoSub.tempCorr,'NoSub');
% plot_histograms({'B Correlation Coefficient'},-1:0.2:1,eventSub.BmagCorr, 'Sub',eventNoSub.BmagCorr,'NoSub');
%
% %plot Histograms for SS and noSS
% plot_scatter({'Shear Angle';'V Correlation Coefficient'},eventSub.Shear,eventSub.vCorr, 'Sub',eventNoSub.Shear,eventNoSub.vCorr, 'NoSub');
% plot_scatter({'Shear Angle';'V_x Correlation Coefficient'},eventSub.Shear,eventSub.vxCorr,'Sub',eventNoSub.Shear,eventNoSub.vxCorr, 'NoSub');
% plot_scatter({'Shear Angle';'T_{para} Correlation Coefficient'},eventSub.Shear,eventSub.tempparaCorr,'Sub',eventNoSub.Shear,eventNoSub.tempparaCorr, 'NoSub');
% plot_scatter({'Shear Angle';'T_{perp} Correlation Coefficient'},eventSub.Shear,eventSub.tempperpCorr, 'Sub',eventNoSub.Shear,eventNoSub.tempperpCorr, 'NoSub');
% plot_scatter({'Shear Angle';'T Correlation Coefficient'},eventSub.Shear,eventSub.tempCorr, 'Sub',eventNoSub.Shear,eventNoSub.tempCorr, 'NoSub');
% plot_scatter({'Shear Angle';'B Correlation Coefficient'},eventSub.Shear,eventSub.BmagCorr, 'Sub',eventNoSub.Shear,eventNoSub.BmagCorr, 'NoSub');
%
% %plot Histograms for SS and noSS
% plot_scatter({'densitySTD';'V Correlation Coefficient'},eventSub.densityCV,eventSub.vCorr, 'Sub',eventNoSub.densityCV,eventNoSub.vCorr, 'NoSub');
% plot_scatter({'densitySTD';'V_x Correlation Coefficient'},eventSub.densityCV,eventSub.vxCorr,'Sub',eventNoSub.densityCV,eventNoSub.vxCorr, 'NoSub');
% plot_scatter({'densitySTD';'T_{para} Correlation Coefficient'},eventSub.densityCV,eventSub.tempparaCorr,'Sub',eventNoSub.densityCV,eventNoSub.tempparaCorr, 'NoSub');
% plot_scatter({'densitySTD';'T_{perp} Correlation Coefficient'},eventSub.densityCV,eventSub.tempperpCorr, 'Sub',eventNoSub.densityCV,eventNoSub.tempperpCorr, 'NoSub');
% plot_scatter({'densitySTD';'T Correlation Coefficient'},eventSub.densityCV,eventSub.tempCorr, 'Sub',eventNoSub.densityCV,eventNoSub.tempCorr, 'NoSub');
% plot_scatter({'densitySTD';'B Correlation Coefficient'},eventSub.densityCV,eventSub.BmagCorr, 'Sub',eventNoSub.densityCV,eventNoSub.BmagCorr, 'NoSub');

% %plot Histograms for SS and noSS
% plot_histograms({'Mean Density'},0:1:30,eventStructure.meanDensity,'All');
% plot_histograms({'Mean Density'},0:4:40,eventSub.meanDensity,'Sub',eventNoSub.meanDensity,'NoSub');
% plot_histograms({'Relative Peak Density'},0:0.5:8,eventSub.relativeMaxDensity,'Sub',eventNoSub.relativeMaxDensity,'NoSub');
% plot_histogram({'Relative Peak Density'},20,eventSub.relativeMaxDensity,'Sub',eventNoSub.relativeMaxDensity,'NoSub');
% plot_scatter({'densitySTD';'Mean Density'},eventSub.densitySTD,eventSub.meanDensity, 'Sub',eventNoSub.densitySTD,eventNoSub.meanDensity, 'NoSub');
%% functions
function plot_correlationHistograms(corrCoeffs,parameter,bins)
    screenSize = get(0,'ScreenSize');
    figure('Position',[0 screenSize(4) 600 400])
    co = [0 0 1;
        0 1 0;
        1 0 0;
        0 0 0];
    set(gcf,'defaultAxesColorOrder',co)
    set(gcf,'color','w');
    
    histogram(corrCoeffs,bins)
    colormap('winter');
    ylabel('f','FontSize',14)
    xlabel({'Correlation Coefficient'},'FontSize',14)
    title(strcat('f vs.', {' '}, parameter),'FontSize',16,'FontWeight', 'normal')
    set(gca,'XMinorTick','off','TickDir','out','YMinorTick','on','linewidth',2)
    
    print(gcf,'-dpng','-r300',strcat('CorrcoeffHistograms_',parameter));
end

function [] = plot_correlationScatters(shear,corrCoeffs,parameter)
    screenSize = get(0,'ScreenSize');
    figure('Position',[0 screenSize(4) 600 400])
    co = [0 0 1;
        0 1 0;
        1 0 0;
        0 0 0];
    set(gcf,'defaultAxesColorOrder',co)
    set(gcf,'color','w');
    
    scatter(shear,corrCoeffs,'filled')
    
    
    xlabel('Magnetic Shear','FontSize',14)
    ylabel('Correlation Coefficient','FontSize',14)
    ylim([-1.0 1.0])
    title(strcat('Magnetic Shear vs.', { ' ' }, parameter),'FontSize',16,'FontWeight', 'normal')
    set(gca,'XMinorTick','on','TickDir','out','YMinorTick','on','linewidth',2)
    box on
    
    
    
    print(gcf,'-dpng','-r300',strcat('CorrcoeffScatters_',parameter));
    
end


function [] = plot_histogram(parameter,bin_size,data1,label1,data2,label2,data3,label3)
    screenSize = get(0,'ScreenSize');
    figure('Position',[0 screenSize(4) 600 400])
    co = [0 0 1;
        0 1 0;
        1 0 0;
        0 0 0];
    set(gcf,'defaultAxesColorOrder',co)
    set(gcf,'color','w');
    
    
    
    if nargin == 4
        binEdges = linspace(floor(min(data1)),ceil(max(data1)),bin_size);
        histogram(data1,binEdges,'FaceColor','r','FaceAlpha',0.25,'LineWidth',2,'Displayname',label1)
        fileName = strcat(cell2mat(parameter(1)),'_',label1);
    elseif nargin == 6
        binEdges = linspace(floor(min([data1,data2])),ceil(max([data1,data2])),bin_size);
        histogram(data1,binEdges,'FaceColor','r','FaceAlpha',0.25,'LineWidth',2,'Displayname',label1)
        hold on
        histogram(data2,binEdges,'FaceColor','y','FaceAlpha',0.25,'LineWidth',2,'Displayname',label2)
        fileName = strcat(cell2mat(parameter(1)),'_',label1,'_',label2);
    elseif nargin == 8
        binEdges = linspace(floor(min([data1,data2,data3])),ceil(max([data1,data2,data3])),bin_size);
        histogram(data1,binEdges,'FaceColor','r','FaceAlpha',0.25,'LineWidth',2,'Displayname',label1)
        hold on
        histogram(data2,binEdges,'FaceColor','y','FaceAlpha',0.25,'LineWidth',2,'Displayname',label2)
        histogram(data3,binEdges,'FaceColor','c','FaceAlpha',0.25,'LineWidth',2,'Displayname',label3)
        fileName = strcat(cell2mat(parameter(1)),'_',label1,'_',label2,'_',label3);
    end
    
    colormap('winter');
    xlabel(parameter,'FontSize',14)
    if max(abs(binEdges)) > 1.0
        xlim([0 max(binEdges)])
    else
        xlim([-1 1])
    end
    ylabel({'# of Events'},'FontSize',14)
    title(strcat('f vs.', {' '}, parameter(1)),'FontSize',16,'FontWeight', 'normal')
    set(gca,'XMinorTick','off','TickDir','out','YMinorTick','on','linewidth',2)
    
    legend
    
    print(gcf,'-depsc2', '-loose', strcat(fileName,'.eps'));
    
    
end


%Plot the scatter, up to 3 data points
function [] = plot_scatter(parameter,data1a,data1b,label1,data2a,data2b,label2,data3a,data3b,label3)
    screenSize = get(0,'ScreenSize');
    figure('Position',[0 screenSize(4) 600 400])
    co = [0 0 1;
        0 1 0;
        1 0 0;
        0 0 0];
    set(gcf,'defaultAxesColorOrder',co)
    set(gcf,'color','w');
    
    if nargin == 4
        scatter(data1a,data1b,'filled','Displayname',label1)
        fileName = strcat(cell2mat(parameter(1)),'-',cell2mat(parameter(2)),'_',label1);
    elseif nargin == 7
        scatter(data1a,data1b,'filled','Displayname',label1)
        hold on
        scatter(data2a,data2b,'filled','Displayname',label2)
        fileName = strcat(cell2mat(parameter(1)),'-',cell2mat(parameter(2)),'_',label1,'_',label2);
    elseif nargin == 10
        scatter(data1a,data1b,'filled','Displayname',label1)
        hold on
        scatter(data2a,data2b,'filled','Displayname',label2)
        scatter(data3a,data3b,'filled','Displayname',label3)
        fileName = strcat(cell2mat(parameter(1)),'-',cell2mat(parameter(2)),'_',label1,'_',label2,'_',label3);
    end
    
    xlabel(parameter(1),'FontSize',14)
    ylabel(parameter(2),'FontSize',14)
    ylim([-1.0 1.0])
    %title(strcat('Observation Position', { ' ' }, parameter(1),'-',parameter(2)),'FontSize',16,'FontWeight', 'normal')
    title(strcat(parameter(1), { ' ' },'vs.',{ ' ' },parameter(2)),'FontSize',16,'FontWeight', 'normal')
    set(gca,'XMinorTick','on','TickDir','out','YMinorTick','on','linewidth',2)
    legend('Location','SouthEastOutside')
    grid on
    ax = gca;
    ax.XAxisLocation = 'origin';
    ax.YAxisLocation = 'origin';
    box on
    
    
    
    print(gcf,'-depsc2', '-loose', strcat(fileName,'.eps'));
    
end
%Plot the histogram, up to 3 data points
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
        
        
        errorbar(xlabels+8*bin_size/16,A_normalized,Aerror,'LineStyle','none','Color','blue','LineWidth',1.5,'HandleVisibility','off')
        
        errorbar(xlabels+8*bin_size/16,B_normalized,Berror,'LineStyle','none','Color','red','LineWidth',1.5,'HandleVisibility','off')
        
        
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