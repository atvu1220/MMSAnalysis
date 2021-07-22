function [shock_normal,n_cs,r_sc] = plot_mms_observation_summary_concise_mod(Event_number,date_start,date_end,...
        mec_timedata,mec_r_gsedata,fgm_timedata, fgm_bdata,...
        fpi_e_timedata,fpi_e_ndata,fpi_e_vdata,fpi_e_pressdata,fpi_e_tparadata,fpi_e_tperpdata,fpi_e_edata,fpi_e_espectdata,...
        fpi_i_timedata,fpi_i_ndata,fpi_i_vdata,fpi_i_pressdata,fpi_i_tparadata,fpi_i_tperpdata,fpi_i_edata,fpi_i_espectdata,...
        n_cs,shear_angle,duration,coreDensity_STD,coreDensity_CV)
    %Matlab CDF Plotter for MMS
    %Andrew Vu 9/19/18
    %data=spdfcdfread(filename);
    %datainfo=spdfcdfinfo(filename);
    figure('Position',[1 1 650 1050])
    co = [0 0 1;
        0 1 0;
        1 0 0;
        0 0 0];
    set(gcf,'defaultAxesColorOrder',co)
    set(gcf,'color','w');
    plot_gap=1.25;
    %cd '~/Library/Mobile Documents/com~apple~CloudDocs/Research/Analysis'
    %figure('PaperPositionMode', 'auto')
    probe_num = '1';
    num_plots = 10;
    data_type = 'brst';
    % date = '2018/03/18';
    
    
    formatIn='yyyy-mm-dd HH:MM:SS.FFF';
    
    tstart = datenum(date_start,formatIn);
    tend = datenum(date_end,formatIn);
    
    if Event_number == 154
        close
        figure('Position',[1 1 650 1350])
        co = [0 0 1;
            0 1 0;
            1 0 0;
            0 0 0];
        set(gcf,'defaultAxesColorOrder',co)
        set(gcf,'color','w');
        num_plots = num_plots+2;
        eventDataFileName = strcat('MMS1_Data_EventNumber_',num2str(Event_number),'.mat');
        %Fucking icloud keeps deleting files every 30 minutes of downloading them
        %eventDataDirectory = '/Users/andrewvu/Library/Mobile Documents/com~apple~CloudDocs/Research/Event Data/';
        eventDataDirectory = '/Users/andrewvu/Library/Mobile Documents/com~apple~CloudDocs/Research/Event Data/';
        eventDataDirectoryFileName = strcat(eventDataDirectory,eventDataFileName);
        load(eventDataDirectoryFileName,'mms1_fgm_timedata_srvy','mms1_fgm_bdata_srvy')
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%Position Data%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %load mec data
    % [mec_timedata,mec_r_gsedata] = load_mec(date_start,probe_num,'srvy');
    
    % %Find the start and end limits of the event in the data
    % start_index_r = find(mec_timedata > tstart, 1);
    % end_index_r = find(mec_timedata > tend, 1);
    %
    % %Convert the datetime to date String, and then crop to our event timeframe
    % mec_timedata = mec_timedata(start_index_r:end_index_r,1);
    % mec_r_gsedata = mec_r_gsedata(start_index_r:end_index_r,:)/6371.2;
    
    [~,mec_r_gsedata,~,~] = crop(mec_timedata,mec_r_gsedata,date_start,date_end);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%Load Data%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    %     %load fgm data
    %     [fgm_timedata, fgm_bdata,~,~]= load_fgm(date_start,probe_num,data_type);
    %
    %     %Load FPI_e
    %     [fpi_e_timedata,fpi_e_ndata,fpi_e_vdata,fpi_e_tparadata,fpi_e_tperpdata,...
    %         fpi_e_edata,fpi_e_espectdata] = load_fpi(date_start,probe_num,data_type,'e');
    %     %Load FPI_i
    %     [fpi_i_timedata,fpi_i_ndata,fpi_i_vdata,fpi_i_tparadata,fpi_i_tperpdata,...
    %         fpi_i_edata,fpi_i_espectdata] = load_fpi(date_start,probe_num,data_type,'i');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%Plot Data%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    plot_order = 1;
    %Plot Magnetic Fields
    plot_fgm_magnetic(date_start,date_end,fgm_timedata,fgm_bdata,num_plots,plot_order)
    if Event_number==154
    else
    title('MMS1 Observatory Summary', 'FontSize', 18, 'FontWeight', 'normal')
    end
    plot_pos = get(gca,'Position');
    set(gca,'Position',[plot_pos(1), plot_pos(2)-plot_pos(4)*plot_gap*(plot_order-1), plot_pos(3), plot_pos(4)]);
    plot_order = plot_order+1;
    
    
    %Plot FPI number density
    plot_fpi_number(date_start,date_end,fpi_i_timedata,fpi_i_ndata,num_plots,plot_order,'i')
    hold on
    plot_fpi_number(date_start,date_end,fpi_e_timedata,fpi_e_ndata,num_plots,plot_order,'e')
    hold off
    set(gca,'Position',[plot_pos(1), plot_pos(2)-plot_pos(4)*plot_gap*(plot_order-1), plot_pos(3), plot_pos(4)]);
    plot_order = plot_order+1;
    
    %Display Sigma and CV
% % % % % %     annotation('textbox',get(gca,'Position'),...
% % % % % %         'String',{strcat('\sigma=',num2str(coreDensity_STD));strcat('cv=',num2str(coreDensity_CV))},...
% % % % % %         'VerticalAlignment','Top','Edgecolor','none','FontSize', 10);
    
    %Plot FPI velocity data
    plot_fpi_bulkv(date_start,date_end,fpi_i_timedata,fpi_i_vdata,num_plots,plot_order)
    set(gca,'Position',[plot_pos(1), plot_pos(2)-plot_pos(4)*plot_gap*(plot_order-1), plot_pos(3), plot_pos(4)]);
    plot_order = plot_order+1;
    
    %Plot FPI velocity data for electrons
    plot_fpi_bulkv(date_start,date_end,fpi_e_timedata,fpi_e_vdata,num_plots,plot_order)
    set(gca,'Position',[plot_pos(1), plot_pos(2)-plot_pos(4)*plot_gap*(plot_order-1), plot_pos(3), plot_pos(4)]);
    plot_order = plot_order+1;
    
    %Plot Current
    plot_fpi_current(date_start,date_end,fpi_e_timedata,fpi_i_timedata,fpi_e_ndata,fpi_e_vdata,fpi_i_vdata,num_plots,plot_order)
    set(gca,'Position',[plot_pos(1), plot_pos(2)-plot_pos(4)*plot_gap*(plot_order-1), plot_pos(3), plot_pos(4)]);
    plot_order = plot_order+1;
    
    %Plot JxB
    plot_JxB(date_start,date_end,fpi_e_timedata,fpi_i_timedata,fpi_e_ndata,fpi_e_vdata,fpi_i_vdata,fgm_timedata,fgm_bdata,num_plots,plot_order)
    set(gca,'Position',[plot_pos(1), plot_pos(2)-plot_pos(4)*plot_gap*(plot_order-1), plot_pos(3), plot_pos(4)]);
    plot_order = plot_order+1;
    

    %Plot FPI dynamic pressure
    plot_total_pressure(date_start,date_end,fgm_timedata,fgm_bdata,fpi_i_timedata,fpi_i_ndata,fpi_i_vdata,fpi_i_pressdata,num_plots,plot_order,'i')
    set(gca,'Position',[plot_pos(1), plot_pos(2)-plot_pos(4)*plot_gap*(plot_order-1), plot_pos(3), plot_pos(4)]);
    plot_order = plot_order+1;
    
    [shock_normal,r_sc] = calculate_bowshocknormal(date_start,mec_r_gsedata,fpi_i_timedata,fpi_i_ndata,fpi_i_vdata);
    
    %Plot Temperatures
    plot_fpi_temp(date_start,date_end,fpi_i_timedata,fpi_i_tparadata,fpi_i_tperpdata,num_plots,plot_order,'i')
    hold on
    plot_fpi_temp(date_start,date_end,fpi_e_timedata,fpi_e_tparadata,fpi_e_tperpdata,num_plots,plot_order,'e')
    hold off
    set(gca,'Position',[plot_pos(1), plot_pos(2)-plot_pos(4)*plot_gap*(plot_order-1), plot_pos(3), plot_pos(4)]);
    plot_order = plot_order+1;
    
    %Energy Spectrum for Ions
    plot_fpi_energyspect(date_start,date_end,fpi_i_timedata,fpi_i_edata,fpi_i_espectdata,num_plots,plot_order)
    ylabel({'Ion';'Energy';'[eV]'},'FontSize', 14)

    set(gca,'Position',[plot_pos(1), plot_pos(2)-plot_pos(4)*plot_gap*(plot_order-1), plot_pos(3), plot_pos(4)]);
    yticks([1,2,3,4])
    plot_order = plot_order+1;
    
    %Energy Spectrum for Electrons
    plot_fpi_energyspect(date_start,date_end,fpi_e_timedata,fpi_e_edata,fpi_e_espectdata,num_plots,plot_order)
    ylabel({'Electron';'Energy';'[eV]'},'FontSize', 14)

    set(gca,'Position',[plot_pos(1), plot_pos(2)-plot_pos(4)*plot_gap*(plot_order-1), plot_pos(3), plot_pos(4)]);
    yticks([1,2,3,4])
    datetick('keeplimits')
    plot_order = plot_order+1;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%For Event 154, Paper Case Study%%%%%
    if Event_number==154
        load('/Users/andrewvu/Library/Mobile Documents/com~apple~CloudDocs/Research/Analyzed Events/Event_Number_154_ForPaper/themisBdataforPaperEvent154.mat','Btimedata','Bdata')
        %MMS Survey Data
        plot_fgm_magnetic('2019-01-07 15:17:00.000','2019-01-07 15:27:00.000',...
            mms1_fgm_timedata_srvy,mms1_fgm_bdata_srvy(:,1:3),num_plots,plot_order)
        ylim([-5.75 5.75])
        ylabel({'MMS1';'B';'[nT]'},'FontSize', 14)
        datetick('keeplimits')
        legend({'B_x', 'B_y', 'B_z','B_t'},'FontSize',10)
        set(gca,'Position',[plot_pos(1), plot_pos(2)-plot_pos(4)*1.01*plot_gap*(plot_order-1), plot_pos(3), plot_pos(4)]);
        plot_order = plot_order+1;
        
        % THEMIS B Survey Data
        plot_fgm_magnetic('2019-01-07 15:07:00.000','2019-01-07 15:17:00.000',...
            Btimedata,Bdata,num_plots,plot_order)
        ylabel({'TH-B';'B';'[nT]'},'FontSize', 14)
        datetick('keeplimits')
        legend({'B_x', 'B_y', 'B_z','B_t'},'FontSize',10)
        set(gca,'Position',[plot_pos(1), plot_pos(2)-plot_pos(4)*1.015*plot_gap*(plot_order-1), plot_pos(3), plot_pos(4)]);
        plot_order = plot_order+1;
        
        colormap(jet)
        %Annotations now have different heights
        annotation('textbox',[plot_pos(1)-1.65*plot_pos(4), plot_pos(2)-plot_pos(4)*0.9875*plot_gap*(plot_order-1), plot_pos(3), plot_pos(4)],...
            'String',{datestr(tstart,'YYYY mmm dd')},...
            'VerticalAlignment','Top','Edgecolor','none','FontSize', 14);

        plot_name =  strcat(num2str(Event_number),'_Event_Number');
        print(gcf,'-dpng','-r900', '-loose', plot_name);
        save(plot_name)

    elseif ((Event_number==7) || (Event_number==48))
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        annotation('textbox',[plot_pos(1)-1.4875*plot_pos(4), plot_pos(2)-plot_pos(4)*0.97*plot_gap*(plot_order-1), plot_pos(3), plot_pos(4)],...
            'String',{datestr(tstart,'YYYY mmm dd')},...
            'VerticalAlignment','Top','Edgecolor','none','FontSize', 14);
        
        plot_name =  strcat(num2str(Event_number),'_Event_Number');
        print(gcf,'-dpng','-r900', '-loose', plot_name);
        save(plot_name)
        
    else
        
        annotation('textbox',[plot_pos(1)-plot_pos(4), plot_pos(2)-plot_pos(4)*plot_gap*(plot_order-1), plot_pos(3), plot_pos(4)],...
            'String',{datestr(tstart,'YYYY mmm dd'),'X-GSE (Re):', 'Y-GSE (Re):', 'Z-GSE (Re):'},...
            'VerticalAlignment','Top','Edgecolor','none','FontSize', 14);
        %
        r_cols = length(mec_r_gsedata);
        %n_cs = tdnormal(date_start,date_end,fgm_timedata,fgm_bdata,'event'); %(-X) GSE for current sheet normal
        annotation('textbox',[plot_pos(1)+plot_pos(3)/r_cols+0.05, plot_pos(2)-plot_pos(4)*plot_gap*(plot_order-1), plot_pos(3), plot_pos(4)],...
            'String',{'',num2str(round(mec_r_gsedata(1,:)/6371.2*10)/10,'%2g\n')},...
            'VerticalAlignment','Top','Edgecolor','none','FontSize', 14);
        
        annotation('textbox',[plot_pos(1)+plot_pos(3)/r_cols+0.05*4 plot_pos(2)-plot_pos(4)*plot_gap*(plot_order-1), plot_pos(3), plot_pos(4)],...
            'String',{strcat('n_{cs}: [',num2str(n_cs,'%.4f '),']'),strcat('n_{bs}: [',num2str(shock_normal,'%.4f '),']'),...
            strcat('Shear Angle (Deg):',num2str(shear_angle,'%2.1f ')),strcat('Duration (s):',num2str(duration,'%2.2f '))},...
            'VerticalAlignment','Top','Edgecolor','none','FontSize', 14);
        
        
        %
        % annotation('textbox',[plot_pos(1)+2*(plot_pos(3)/r_cols)+0.025, plot_pos(2)-plot_pos(4)*7.5, plot_pos(3), plot_pos(4)],...
        %     'String',{'',num2str(round(mec_r_gsedata(1,:)*10)/10,'%2g\n')},...
        %     'VerticalAlignment','Top','Edgecolor','none','FontSize', 14);
        plot_name =  strcat(num2str(Event_number),'_Event_Number');
        %print(gcf, '-dpdf', '-opengl',plot_name,'-fillpage');
        print(gcf,'-dpng','-r300', '-loose', plot_name);
        %movefile(plot_name, '~/Library/Mobile Documents/com~apple~CloudDocs/Research/Analysis')
    end
    

    
end