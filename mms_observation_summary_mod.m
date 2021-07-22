%Matlab CDF Plotter for MMS
%Andrew Vu 9/19/18
%data=spdfcdfread(filename);
%datainfo=spdfcdfinfo(filename);
clear
warning('off','MATLAB:polyfit:RepeatedPointsOrRescale')
cd '~/Library/Mobile Documents/com~apple~CloudDocs/Research/Analysis/SummaryPlots/'
mms_directory = '/Users/andrewvu/data/mms/';
num_plots = 9;
data_type = 'brst';
plot_gap=1.25;
probe_num = '1';

ssCount = 0;

ssSTD = 2.00; %1.5 is 87% 68-87-95-99

%1.50 std = 60 events
%2.00 std = 44 events
%2.25 std = 37 events
%2.50 std = 23 events
%3.00 std = 16 events
minSSDataPoints = 4;


for i=154%2:174
    close all
    figure('Position',[1 1 650 850])
    co = [0 0 1;
        0 1 0;
        1 0 0;
        0 0 0];
    set(gcf,'defaultAxesColorOrder',co)
    set(gcf,'color','w');
    
    
    Event_number = i
    
    [~, ~,~,~,...
        date_start,date_end,...
        ~, left_InnerEdge,...
        right_InnerEdge, ~,...
        ~, ~,...
        ~, ~,...
        ~, ~,...
        ~, ~] = get_eventTimes(Event_number);
    
    
    
    
    formatIn='yyyy-mm-dd HH:MM:SS.FFF';
    tstart = datenum(date_start,formatIn);
    tend = datenum(date_end,formatIn);
    
    %November 2019, added data retrieval from .mat file to reduce storage
    eventDataFileName = strcat('MMS1_Data_EventNumber_',num2str(Event_number),'.mat');
    eventDataDirectory = '/Users/andrewvu/Library/Mobile Documents/com~apple~CloudDocs/Research/Event Data/';
    eventDataDirectoryFileName = strcat(eventDataDirectory,eventDataFileName);
    load(eventDataDirectoryFileName)
    % % % % % %
% % % % % %     %%%%%%%%%%%%%%%%%%%%%%%%%%%Position Data%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % %     %load mec data
% % % % % %     [mec_timedata,mec_r_gsedata] = load_mec(date_start,probe_num,'srvy');
% % % % % %     
% % % % % %     % %Find the start and end limits of the event in the data
% % % % % %     % start_index_r = find(mec_timedata > tstart, 1);
% % % % % %     % end_index_r = find(mec_timedata > tend, 1);
% % % % % %     %
% % % % % %     % %Convert the datetime to date String, and then crop to our event timeframe
% % % % % %     % mec_timedata = mec_timedata(start_index_r:end_index_r,1);
% % % % % %     % mec_r_gsedata = mec_r_gsedata(start_index_r:end_index_r,:)/6371.2;
% % % % % %     
% % % % % %     [mec_timedata,mec_r_gsedata,~,~] = crop(mec_timedata,mec_r_gsedata,date_start,date_end);
% % % % % %     
% % % % % %     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % %     
% % % % % %     
% % % % % %     %%%%%%%%%%%%%%%%%%%%%%%%%%Load Data%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % %     
% % % % % %     %load fgm data
% % % % % %     [fgm_timedata, fgm_bdata,~,~]= load_fgm(date_start,date_end,probe_num,data_type);
% % % % % %     
% % % % % %     %Load FPI_e
% % % % % %     [fpi_e_timedata,fpi_e_ndata,fpi_e_vdata,fpi_e_tparadata,fpi_e_tperpdata,...
% % % % % %         fpi_e_edata,fpi_e_espectdata] = load_fpi(date_start,date_end,probe_num,data_type,'e');
% % % % % %     %Load FPI_i
% % % % % %     [fpi_i_timedata,fpi_i_ndata,fpi_i_vdata,fpi_i_tparadata,fpi_i_tperpdata,...
% % % % % %         fpi_i_edata,fpi_i_espectdata,fpi_i_pressdata] = load_fpi(date_start,date_end,probe_num,data_type,'i');
% % % % % %     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%Plot Data%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    plot_order = 1;
    %Plot Magnetic Fields
    plot_fgm_magnetic(date_start,date_end,mms1_fgm_timedata_raw,mms1_fgm_bdata_raw,num_plots,plot_order)
    title('MMS1 Observatory Summary', 'FontSize', 18, 'FontWeight', 'normal')
    plot_pos = get(gca,'Position');
    set(gca,'Position',[plot_pos(1), plot_pos(2)-plot_pos(4)*plot_gap*(plot_order-1), plot_pos(3), plot_pos(4)]);
    plot_order = plot_order+1;
    
    
    
    
    
    %Plot FPI number density
    
    hold on
    plot_fpi_number(date_start,date_end,fpi_i_timedata,(fpi_i_ndata),num_plots,plot_order,'i')
    hold on
    plot_fpi_number(date_start,date_end,fpi_e_timedata,(fpi_e_ndata),num_plots,plot_order,'e')
    
    %Calculate Density Statistics
    [substructurePresent] = plot_densityStats(date_start,date_end,left_InnerEdge,...
        right_InnerEdge,fpi_i_timedata,fpi_i_ndata,num_plots,plot_order,ssSTD,minSSDataPoints);
    if substructurePresent == 1
        ssCount = ssCount + 1
    end
    %Core Density STD
% % % %     [coreDensity_Sigma,coreDensity_CV] = calculate_coreDensitySTD(left_InnerEdge,right_InnerEdge,fpi_i_timedata,fpi_i_ndata);
% % % %     
% % % %     BoxPos = get(gca,'Position');
% % % %     BoxPos(2) = BoxPos(2) + BoxPos(4)/4;
% % % %     annotation('textbox',BoxPos,...
% % % %         'String',{strcat('\sigma=',num2str(coreDensity_Sigma));strcat('cv=',num2str(coreDensity_CV))},...
% % % %         'VerticalAlignment','Top','Edgecolor','none','FontSize', 10);
% % % %     
    
    hold off
    
    set(gca,'Position',[plot_pos(1), plot_pos(2)-plot_pos(4)*plot_gap*(plot_order-1), plot_pos(3), plot_pos(4)]);
    plot_order = plot_order+1;
    
    %Plot FPI velocity data
    plot_fpi_bulkv(date_start,date_end,fpi_i_timedata,fpi_i_vdata,num_plots,plot_order)
    set(gca,'Position',[plot_pos(1), plot_pos(2)-plot_pos(4)*plot_gap*(plot_order-1), plot_pos(3), plot_pos(4)]);
    plot_order = plot_order+1;
    
    %Plot FPI dynamic pressure
%     plot_fpi_dynamic_pressure(date_start,date_end,fpi_i_timedata,fpi_i_ndata,fpi_i_vdata,num_plots,plot_order,'i')
    plot_total_pressure(date_start,date_end,mms1_fgm_timedata_raw,mms1_fgm_bdata_raw,fpi_i_timedata,fpi_i_ndata,fpi_i_vdata,fpi_i_pressdata,num_plots,plot_order,'i')

%     hold on
%     plot_fpi_dynamic_pressure(date_start,date_end,fpi_e_timedata,fpi_e_ndata,fpi_e_vdata,num_plots,plot_order,'e')
%     hold off
    set(gca,'Position',[plot_pos(1), plot_pos(2)-plot_pos(4)*plot_gap*(plot_order-1), plot_pos(3), plot_pos(4)]);
    plot_order = plot_order+1;
    
    %Plot -VxB
%     plot_vcrossb(date_start,date_end,fgm_timedata,fgm_bdata,fpi_i_timedata,fpi_i_vdata,num_plots,plot_order,'n_cs')
%     set(gca,'Position',[plot_pos(1), plot_pos(2)-plot_pos(4)*plot_gap*(plot_order-1), plot_pos(3), plot_pos(4)]);
%     plot_order = plot_order+1;
    %leading edge negative edot n_td means towards the discontinuity and
    %positive edot n_td means towards the discontinuity
    
    shock_normal = calculate_bowshocknormal(date_start,mms1_mec_rdata_raw,fpi_i_timedata,fpi_i_ndata,fpi_i_vdata);
    
    %Plot shock angles
    %plot_shockangle(date_start,date_end,fgm_timedata,fgm_bdata,num_plots,plot_order,shock_normal)
    
    %Plot Entropy
%     plot_entropy(date_start,date_end,fpi_i_timedata,fpi_i_ndata,fpi_i_tparadata,fpi_i_tperpdata,'i',num_plots,plot_order); hold on
%     %plot_entropy(date_start,date_end,fpi_i_timedata,fpi_e_ndata,fpi_e_tparadata,fpi_e_tperpdata,'e',num_plots,plot_order); hold off
%     set(gca,'Position',[plot_pos(1), plot_pos(2)-plot_pos(4)*plot_gap*(plot_order-1), plot_pos(3), plot_pos(4)]);
%     plot_order = plot_order+1;
    
    
    %Plot Temperatures
    plot_fpi_temp(date_start,date_end,fpi_i_timedata,fpi_i_tparadata,fpi_i_tperpdata,num_plots,plot_order,'i')
    hold on
    plot_fpi_temp(date_start,date_end,fpi_e_timedata,fpi_e_tparadata,fpi_e_tperpdata,num_plots,plot_order,'e')
    hold off
    set(gca,'Position',[plot_pos(1), plot_pos(2)-plot_pos(4)*plot_gap*(plot_order-1), plot_pos(3), plot_pos(4)]);
    plot_order = plot_order+1;
    
    %Energy Spectrum for Ions
    plot_fpi_energyspect(date_start,date_end,fpi_i_timedata,fpi_i_edata,fpi_i_espectdata,num_plots,plot_order)
    set(gca,'Position',[plot_pos(1), plot_pos(2)-plot_pos(4)*plot_gap*(plot_order-1), plot_pos(3), plot_pos(4)]);
    plot_order = plot_order+1;
    
    %Energy Spectrum for Electrons
    plot_fpi_energyspect(date_start,date_end,fpi_e_timedata,fpi_e_edata,fpi_e_espectdata,num_plots,plot_order)
    
    set(gca,'Position',[plot_pos(1), plot_pos(2)-plot_pos(4)*plot_gap*(plot_order-1), plot_pos(3), plot_pos(4)]);
    datetick('keeplimits')
    plot_order = plot_order+1;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    annotation('textbox',[plot_pos(1)-plot_pos(4), plot_pos(2)-plot_pos(4)*plot_gap*(plot_order-1), plot_pos(3), plot_pos(4)],...
        'String',{datestr(tstart,'YYYY mmm dd'),'X-GSE (Re):', 'Y-GSE (Re):', 'Z-GSE (Re):'},...
        'VerticalAlignment','Top','Edgecolor','none','FontSize', 14);
    %
    r_cols = length(mms1_mec_rdata_raw);
    n_cs = tdnormal(date_start,date_end,mms1_fgm_timedata_raw,mms1_fgm_bdata_raw,'event');
%     annotation('textbox',[plot_pos(1)+plot_pos(3)/r_cols+0.05, plot_pos(2)-plot_pos(4)*plot_gap*(plot_order-1), plot_pos(3), plot_pos(4)],...
%         'String',{'',num2str(round(mms1_mec_rdata_raw(1,:)/6371.2*10)/10,'%2g\n')},...
%         'VerticalAlignment','Top','Edgecolor','none','FontSize', 14);
%     
%     annotation('textbox',[plot_pos(1)+plot_pos(3)/r_cols+0.05*4 plot_pos(2)-plot_pos(4)*plot_gap*(plot_order-1), plot_pos(3), plot_pos(4)],...
%         'String',{strcat('n_{cs}: [',num2str(n_cs,'%.4f '),']'),strcat('n_{bs}: [',num2str(shock_normal,'%.4f '),']')},...
%         'VerticalAlignment','Top','Edgecolor','none','FontSize', 14);
    
        annotation('textbox',[plot_pos(1)+plot_pos(3)/r_cols+0.05*4, plot_pos(2)-plot_pos(4)*plot_gap*(plot_order-1), plot_pos(3), plot_pos(4)],...
        'String',{'',num2str(round(mms1_mec_rdata_raw(1,:)/6371.2*10)/10,'%2g\n')},...
        'VerticalAlignment','Top','HorizontalAlignment','Left','Edgecolor','none','FontSize', 14);
    
    %
    % annotation('textbox',[plot_pos(1)+2*(plot_pos(3)/r_cols)+0.025, plot_pos(2)-plot_pos(4)*7.5, plot_pos(3), plot_pos(4)],...
    %     'String',{'',num2str(round(mec_r_gsedata(1,:)*10)/10,'%2g\n')},...
    %     'VerticalAlignment','Top','Edgecolor','none','FontSize', 14);
    
    
    %plot_name =  strcat('mms',probe_num,'_ObservationSummary_',date_start(1:19),'_',date_end(1:19),'.pdf');
    plot_name =  strcat('mms',probe_num,'_ObservationSummary_','Event_',num2str(Event_number),'.png');
    %print(gcf, '-dpdf', '-opengl',plot_name,'-fillpage');
    print(gcf,'-dpng','-r300',plot_name);
    %movefile(plot_name, '~/Library/Mobile Documents/com~apple~CloudDocs/Research/Analysis')
end
ssCount