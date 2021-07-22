%Matlab CDF Plotter for MMS
%Andrew Vu 9/19/18
%data=spdfcdfread(filename);
%datainfo=spdfcdfinfo(filename);
clear; close all
warning('off','MATLAB:polyfit:RepeatedPointsOrRescale')
cd '~/Library/Mobile Documents/com~apple~CloudDocs/Research/Analysis/Correlation/'
mms_directory = '/Users/andrewvu/data/mms/';
num_plots = 9;
data_type = 'brst';
plot_gap=1.25;
probe_num = '1';

ssSTDs = 2.00;
minSSDataPoints = 4;

Events = 2:151;



%

eventStructure = struct('eventNumber', zeros(size(Events)),...
    'densitySTD',zeros(size(Events)),...
    'meanDensity',zeros(size(Events)),...
    'relativeMaxDensity',zeros(size(Events)),...
    'Shear',zeros(size(Events)),...
    'vCorr', zeros(size(Events)), ...
    'vxCorr', zeros(size(Events)), ...
    'tempparaCorr', zeros(size(Events)), ...
    'tempperpCorr', zeros(size(Events)), ...
    'tempCorr', zeros(size(Events)));

eventSub = struct('eventNumber', zeros(size(Events)),...
    'densitySTD',zeros(size(Events)),...
    'meanDensity',zeros(size(Events)),...
    'relativeMaxDensity',zeros(size(Events)),...
    'Shear',zeros(size(Events)),...
    'vCorr', zeros(size(Events)), ...
    'vxCorr', zeros(size(Events)), ...
    'tempparaCorr', zeros(size(Events)), ...
    'tempperpCorr', zeros(size(Events)), ...
    'tempCorr', zeros(size(Events)));

eventNoSub = struct('eventNumber', zeros(size(Events)),...
    'densitySTD',zeros(size(Events)),...
    'meanDensity',zeros(size(Events)),...
    'relativeMaxDensity',zeros(size(Events)),...
    'Shear',zeros(size(Events)),...
    'vCorr', zeros(size(Events)), ...
    'vxCorr', zeros(size(Events)), ...
    'tempparaCorr', zeros(size(Events)), ...
    'tempperpCorr', zeros(size(Events)), ...
    'tempCorr', zeros(size(Events)));

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
    
    
    
    [Event_Type, Substructure,threshold_std,~,...
        event_start,event_end,...
        left_OuterEdge, left_InnerEdge,...
        right_InnerEdge, right_OuterEdge,...
        leading_leftmost_date, leading_rightmost_date,...
        trailing_leftmost_date, trailing_rightmost_date,...
        ~, ~,...
        ~, ~] = get_eventTimes(Event_number);
    
    
    
    %     if contains(Event_Type,'HFA')
    %
    %         if Substructure == 0
    %             eventStructure.eventNumber(i) = Event_number;
    %             eventNoSub.eventNumber(i) = Event_number;
    %
    %         elseif Substructure == 1
    %             eventStructure.eventNumber(i) = Event_number;
    %             eventSub.eventNumber(i) = Event_number;
    %         end
    %
    %     else
    %         continue
    %     end
    
    
    
    
    
    formatIn='yyyy-mm-dd HH:MM:SS.FFF';
    tstart = datenum(event_start,formatIn);
    tend = datenum(event_end,formatIn);
    
    
    
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
    
    
    
    %% calculate STD,CV
    %     [coreDensity_STD,coreDensity_CV] = calculate_coreDensitySTD(left_InnerEdge,right_InnerEdge,fpi_i_timedata,fpi_i_ndata);
    %
    %
    %     if contains(Event_Type,'HFA')
    %
    %         if coreDensity_STD <= 1.0
    %             eventStructure.eventNumber(i) = Event_number;
    %             eventNoSub.eventNumber(i) = Event_number;
    %
    %             eventStructure.densitySTD(i) = coreDensity_STD;
    %             eventNoSub.densitySTD(i) = coreDensity_STD;
    %             Substructure = 0;
    %         elseif coreDensity_STD > 1.0
    %             eventStructure.eventNumber(i) = Event_number;
    %             eventSub.eventNumber(i) = Event_number;
    %
    %             eventStructure.densitySTD(i) = coreDensity_STD;
    %             eventSub.densitySTD(i) = coreDensity_STD;
    %             Substructure = 1;
    %         end
    %
    %     else
    %         continue
    %     end
    
    
    
    %% binary for sigmaSigma threshold
    [substructurePresent,~,fpi_coreNdata,~,~,~,~] = calculate_densityStats(fpi_i_timedata,fpi_i_ndata,left_InnerEdge,right_InnerEdge,ssSTDs,minSSDataPoints);
    [coreDensity_Sigma,coreDensity_CV] = calculate_coreDensitySTD(left_InnerEdge,right_InnerEdge,fpi_i_timedata,fpi_i_ndata);
    %
    %
    %
    %     if contains(Event_Type,'HFA')
    %
    %         if substructurePresent == 0
    %             eventStructure.eventNumber(i) = Event_number;
    %             eventNoSub.eventNumber(i) = Event_number;
    %
    %             eventStructure.densitySTD(i) = coreDensity_CV;
    %             eventNoSub.densitySTD(i) = coreDensity_CV;
    %
    %             eventStructure.meanDensity(i) = mean(fpi_coreNdata);
    %             eventStructure.relativeMaxDensity(i) = abs(mean(fpi_coreNdata)-max(fpi_coreNdata))/mean(fpi_coreNdata);
    %             eventNoSub.meanDensity(i) = mean(fpi_coreNdata);
    %             eventNoSub.relativeMaxDensity(i) = abs(mean(fpi_coreNdata)-max(fpi_coreNdata))/mean(fpi_coreNdata);
    %
    %             Substructure = 0;
    %
    %         elseif substructurePresent == 1
    %             eventStructure.eventNumber(i) = Event_number;
    %             eventSub.eventNumber(i) = Event_number;
    %
    %             eventStructure.densitySTD(i) = coreDensity_CV;
    %             eventSub.densitySTD(i) = coreDensity_CV;
    %
    %             eventStructure.meanDensity(i) = mean(fpi_coreNdata);
    %             eventStructure.relativeMaxDensity(i) = abs(mean(fpi_coreNdata)-max(fpi_coreNdata))/mean(fpi_coreNdata);
    %             eventSub.meanDensity(i) = mean(fpi_coreNdata);
    %             eventSub.relativeMaxDensity(i) = abs(mean(fpi_coreNdata)-max(fpi_coreNdata))/mean(fpi_coreNdata);
    %
    %             Substructure = 1;
    %
    %         end
    %
    %     else
    %         continue
    %     end
    %
    
    %%
    
    
    
    
    %Calculations
    
    if threshold_std == 0
        %Current Sheet Normal Calculation Plot, manually.
        [~,B_pre,B_post] = manualCurrentSheet(event_start,event_end,leading_leftmost_date,leading_rightmost_date,...
            trailing_leftmost_date,trailing_rightmost_date,fgm_timedata_srvy,fgm_bdata_srvy);
    else
        %Calculate Current Sheet
        [~,B_pre,B_post] = calculateCurrentSheet(event_start,event_end,fgm_timedata_srvy,fgm_bdata_srvy,threshold_std);
    end
    
    eventStructure.Shear(i) = angle(B_pre,B_post);
    
    %Extra Part for HFA/SHFA
    if eventStructure.Shear(i) > 30
        Substructure = 1;
        eventStructure.eventNumber(i) = Event_number;
        eventSub.eventNumber(i) = Event_number;
        
        eventStructure.densitySTD(i) = coreDensity_CV;
        eventSub.densitySTD(i) = coreDensity_CV;
        
        eventStructure.meanDensity(i) = mean(fpi_coreNdata);
        %         eventStructure.relativeMaxDensity(i) = abs(mean(fpi_coreNdata)-max(fpi_coreNdata))/mean(fpi_coreNdata);
        eventStructure.relativeMaxDensity(i) = max(fpi_coreNdata)/mean(fpi_coreNdata);
        eventSub.meanDensity(i) = mean(fpi_coreNdata);
        %         eventSub.relativeMaxDensity(i) = abs(mean(fpi_coreNdata)-max(fpi_coreNdata))/mean(fpi_coreNdata);
        eventSub.relativeMaxDensity(i) = max(fpi_coreNdata)/mean(fpi_coreNdata);
    else
        Substructure = 0;
        eventStructure.eventNumber(i) = Event_number;
        eventNoSub.eventNumber(i) = Event_number;
        
        eventStructure.densitySTD(i) = coreDensity_CV;
        eventNoSub.densitySTD(i) = coreDensity_CV;
        
        eventStructure.meanDensity(i) = mean(fpi_coreNdata);
        %         eventStructure.relativeMaxDensity(i) = abs(mean(fpi_coreNdata)-max(fpi_coreNdata))/mean(fpi_coreNdata);
        eventStructure.relativeMaxDensity(i) = max(fpi_coreNdata)/mean(fpi_coreNdata);
        eventNoSub.meanDensity(i) = mean(fpi_coreNdata);
        %         eventNoSub.relativeMaxDensity(i) = abs(mean(fpi_coreNdata)-max(fpi_coreNdata))/mean(fpi_coreNdata);
        eventNoSub.relativeMaxDensity(i) = max(fpi_coreNdata)/mean(fpi_coreNdata);
    end
    %End Extra
    
    if Substructure == 0
        eventNoSub.Shear(i) = angle(B_pre,B_post);
    elseif Substructure == 1
        eventSub.Shear(i) = angle(B_pre,B_post);
    end
    
    
    
    %Correlation Coefficients
    [eventStructure.vCorr(i),eventStructure.vxCorr(i),...
        eventStructure.tempparaCorr(i),eventStructure.tempperpCorr(i),...
        eventStructure.tempCorr(i)] = calculate_correlation(left_InnerEdge,right_InnerEdge,...
        fpi_i_timedata,fpi_i_ndata,fpi_i_vdata,fpi_i_tparadata,fpi_i_tperpdata);
    
    
    if Substructure == 0
        eventNoSub.vCorr(i) = eventStructure.vCorr(i);
        eventNoSub.vxCorr(i) = eventStructure.vxCorr(i);
        eventNoSub.tempparaCorr(i) = eventStructure.tempparaCorr(i);
        eventNoSub.tempperpCorr(i) = eventStructure.tempperpCorr(i);
        eventNoSub.tempCorr(i) =  eventStructure.tempCorr(i);
        
    elseif Substructure == 1
        eventSub.vCorr(i) = eventStructure.vCorr(i);
        eventSub.vxCorr(i) = eventStructure.vxCorr(i);
        eventSub.tempparaCorr(i) = eventStructure.tempparaCorr(i);
        eventSub.tempperpCorr(i) = eventStructure.tempperpCorr(i);
        eventSub.tempCorr(i) =  eventStructure.tempCorr(i);
    end
    
    
    %     annotation('textbox',[plot_pos(1)-plot_pos(4), plot_pos(2)-plot_pos(4)*plot_gap*(plot_order-1), plot_pos(3), plot_pos(4)],...
    %         'String',{datestr(tstart,'YYYY mmm dd'),'X-GSE (Re):', 'Y-GSE (Re):', 'Z-GSE (Re):'},...
    %         'VerticalAlignment','Top','Edgecolor','none','FontSize', 14);
    %     %
    %     r_cols = length(mec_r_gsedata);
    %     n_cs = tdnormal(event_start,event_end,fgm_timedata,fgm_bdata,'event');
    %     annotation('textbox',[plot_pos(1)+plot_pos(3)/r_cols+0.05, plot_pos(2)-plot_pos(4)*plot_gap*(plot_order-1), plot_pos(3), plot_pos(4)],...
    %         'String',{'',num2str(round(mec_r_gsedata(1,:)/6371.2*10)/10,'%2g\n')},...
    %         'VerticalAlignment','Top','Edgecolor','none','FontSize', 14);
    %
    %     annotation('textbox',[plot_pos(1)+plot_pos(3)/r_cols+0.05*4 plot_pos(2)-plot_pos(4)*plot_gap*(plot_order-1), plot_pos(3), plot_pos(4)],...
    %         'String',{strcat('n_{cs}: [',num2str(n_cs,'%.4f '),']'),strcat('n_{bs}: [',num2str(shock_normal,'%.4f '),']')},...
    %         'VerticalAlignment','Top','Edgecolor','none','FontSize', 14);
    %
    %     %
    %     % annotation('textbox',[plot_pos(1)+2*(plot_pos(3)/r_cols)+0.025, plot_pos(2)-plot_pos(4)*7.5, plot_pos(3), plot_pos(4)],...
    %     %     'String',{'',num2str(round(mec_r_gsedata(1,:)*10)/10,'%2g\n')},...
    %     %     'VerticalAlignment','Top','Edgecolor','none','FontSize', 14);
    %
    %
    %     %plot_name =  strcat('mms',probe_num,'_ObservationSummary_',date_start(1:19),'_',date_end(1:19),'.pdf');
    %     plot_name =  strcat('mms',probe_num,'_ObservationSummary_','Event_',num2str(Event_number),'.png');
    %     %print(gcf, '-dpdf', '-opengl',plot_name,'-fillpage');
    %     print(gcf,'-dpng','-r300',plot_name);
    %     %movefile(plot_name, '~/Library/Mobile Documents/com~apple~CloudDocs/Research/Analysis')
end
%%
%Get rid of nonzeros because it isn't an HFA
eventStructure.densitySTD(eventStructure.densitySTD==0) = [];
eventStructure.meanDensity(eventStructure.meanDensity==0) = [];
eventStructure.relativeMaxDensity(eventStructure.relativeMaxDensity==0) = [];
eventStructure.Shear(eventStructure.Shear==0) = [];
eventStructure.vCorr(eventStructure.vCorr==0) = [];
eventStructure.vxCorr(eventStructure.vxCorr==0) = [];
eventStructure.tempparaCorr(eventStructure.tempparaCorr==0) = [];
eventStructure.tempperpCorr(eventStructure.tempperpCorr==0) = [];
eventStructure.tempCorr(eventStructure.tempCorr==0) = [];

eventNoSub.densitySTD(eventNoSub.densitySTD==0) = [];
eventNoSub.meanDensity(eventNoSub.meanDensity==0) = [];
eventNoSub.relativeMaxDensity(eventNoSub.relativeMaxDensity==0) = [];
eventNoSub.Shear(eventNoSub.Shear==0) = [];
eventNoSub.vCorr(eventNoSub.vCorr==0) = [];
eventNoSub.vxCorr(eventNoSub.vxCorr==0) = [];
eventNoSub.tempparaCorr(eventNoSub.tempparaCorr==0) = [];
eventNoSub.tempperpCorr(eventNoSub.tempperpCorr==0) = [];
eventNoSub.tempCorr(eventNoSub.tempCorr==0) = [];

eventSub.densitySTD(eventSub.densitySTD==0) = [];
eventSub.meanDensity(eventSub.meanDensity==0) = [];
eventSub.relativeMaxDensity(eventSub.relativeMaxDensity==0) = [];
eventSub.Shear(eventSub.Shear==0) = [];
eventSub.vCorr(eventSub.vCorr==0) = [];
eventSub.vxCorr(eventSub.vxCorr==0) = [];
eventSub.tempparaCorr(eventSub.tempparaCorr==0) = [];
eventSub.tempperpCorr(eventSub.tempperpCorr==0) = [];
eventSub.tempCorr(eventSub.tempCorr==0) = [];

%plot Histograms for all events
plot_histogram({'Core Density CV'},20,eventStructure.densitySTD,'All');
plot_histogram({'V Correlation Coefficient'},20,eventStructure.vCorr,'All');
plot_histogram({'Vx Correlation Coefficient'},20,eventStructure.vxCorr, 'All');
plot_histogram({'Tpara Correlation Coefficient'},20,eventStructure.tempparaCorr, 'All');
plot_histogram({'Tperp Correlation Coefficient'},20,eventStructure.tempperpCorr, 'All');
plot_histogram({'T Correlation Coefficient'},20,eventStructure.tempCorr, 'All');

%plot Histograms for all events
plot_scatter({'Shear Angle';'V Correlation Coefficient'},eventStructure.Shear,eventStructure.vCorr, 'All');
plot_scatter({'Shear Angle';'V_x Correlation Coefficient'},eventStructure.Shear,eventStructure.vxCorr, 'All');
plot_scatter({'Shear Angle';'T_{para} Correlation Coefficient'},eventStructure.Shear,eventStructure.tempparaCorr, 'All');
plot_scatter({'Shear Angle';'T_{perp} Correlation Coefficient'},eventStructure.Shear,eventStructure.tempperpCorr, 'All');
plot_scatter({'Shear Angle';'T Correlation Coefficient'},eventStructure.Shear,eventStructure.tempCorr, 'All');

%plot Histograms for SS and noSS
plot_histogram({'V Correlation Coefficient'},20,eventSub.vCorr,'Sub',eventNoSub.vCorr,'NoSub');
plot_histogram({'Vx Correlation Coefficient'},20,eventSub.vxCorr, 'Sub',eventNoSub.vxCorr,'NoSub');
plot_histogram({'Tpara Correlation Coefficient'},20,eventSub.tempparaCorr, 'Sub',eventNoSub.tempparaCorr,'NoSub');
plot_histogram({'Tperp Correlation Coefficient'},20,eventSub.tempperpCorr, 'Sub',eventNoSub.tempperpCorr,'NoSub');
plot_histogram({'T Correlation Coefficient'},20,eventSub.tempCorr, 'Sub',eventNoSub.tempCorr,'NoSub');

%plot Histograms for SS and noSS
plot_scatter({'Shear Angle';'V Correlation Coefficient'},eventSub.Shear,eventSub.vCorr, 'Sub',eventNoSub.Shear,eventNoSub.vCorr, 'NoSub');
plot_scatter({'Shear Angle';'V_x Correlation Coefficient'},eventSub.Shear,eventSub.vxCorr,'Sub',eventNoSub.Shear,eventNoSub.vxCorr, 'NoSub');
plot_scatter({'Shear Angle';'T_{para} Correlation Coefficient'},eventSub.Shear,eventSub.tempparaCorr,'Sub',eventNoSub.Shear,eventNoSub.tempparaCorr, 'NoSub');
plot_scatter({'Shear Angle';'T_{perp} Correlation Coefficient'},eventSub.Shear,eventSub.tempperpCorr, 'Sub',eventNoSub.Shear,eventNoSub.tempperpCorr, 'NoSub');
plot_scatter({'Shear Angle';'T Correlation Coefficient'},eventSub.Shear,eventSub.tempCorr, 'Sub',eventNoSub.Shear,eventNoSub.tempCorr, 'NoSub');

%plot Histograms for SS and noSS
plot_scatter({'densitySTD';'V Correlation Coefficient'},eventSub.densitySTD,eventSub.vCorr, 'Sub',eventNoSub.densitySTD,eventNoSub.vCorr, 'NoSub');
plot_scatter({'densitySTD';'V_x Correlation Coefficient'},eventSub.densitySTD,eventSub.vxCorr,'Sub',eventNoSub.densitySTD,eventNoSub.vxCorr, 'NoSub');
plot_scatter({'densitySTD';'T_{para} Correlation Coefficient'},eventSub.densitySTD,eventSub.tempparaCorr,'Sub',eventNoSub.densitySTD,eventNoSub.tempparaCorr, 'NoSub');
plot_scatter({'densitySTD';'T_{perp} Correlation Coefficient'},eventSub.densitySTD,eventSub.tempperpCorr, 'Sub',eventNoSub.densitySTD,eventNoSub.tempperpCorr, 'NoSub');
plot_scatter({'densitySTD';'T Correlation Coefficient'},eventSub.densitySTD,eventSub.tempCorr, 'Sub',eventNoSub.densitySTD,eventNoSub.tempCorr, 'NoSub');

%plot Histograms for SS and noSS
plot_histogram({'Mean Density'},20,eventStructure.meanDensity,'All');
plot_histogram({'Mean Density'},20,eventSub.meanDensity,'Sub',eventNoSub.meanDensity,'NoSub');
plot_histogram({'Relative Peak Density'},20,eventSub.relativeMaxDensity,'Sub',eventNoSub.relativeMaxDensity,'NoSub');
% plot_scatter({'densitySTD';'Mean Density'},eventSub.densitySTD,eventSub.meanDensity, 'Sub',eventNoSub.densitySTD,eventNoSub.meanDensity, 'NoSub');






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
