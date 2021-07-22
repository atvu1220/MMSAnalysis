function [] = plot_total_pressure(date_start,date_end,fgm_timedata,fgm_bdata,fpi_timedata,fpi_ndata,fpi_vdata,fpi_pressure,num_plots,plot_order,specie)
    %plots the plasma parameter with the given start date, end date, datas,
    %number of plots and plot order
    
    formatIn='yyyy-mm-dd HH:MM:SS.FFF';
    tstart = datenum(date_start,formatIn);
    tend = datenum(date_end,formatIn);
    
    
    %Interpolate FGM down to FPI
    [~,fgm_bdata_interp] = interpxyz(fgm_timedata,fgm_bdata,fpi_timedata);
    
    %Find the start and end limits of the event in the data
    start_index = find(fpi_timedata > tstart, 1);
    end_index = find(fpi_timedata > tend, 1);
    
    %crop data
    fpi_timedata = fpi_timedata(start_index:end_index,1);
    fpi_ndata = fpi_ndata(start_index:end_index,:);
    fpi_vdata = fpi_vdata(start_index:end_index,:);
    fpi_pressure = fpi_pressure(:,:,start_index:end_index);
    fgm_bdata = fgm_bdata_interp(start_index:end_index,:);
    
    
    %Calculate dynamic pressure rho*v^2 (in nPa)
    dynamic_pressure = calculate_dynamic_pressure(fpi_ndata,fpi_vdata,specie);
    
    %Calculate magnetic pressure B^2/2mu_o (in nPa)
    magnetic_pressure = calculate_magnetic_pressure(fgm_bdata);
    
    %Calculate thermal pressure px^2 + py^2 + pz^2
    thermal_pressure = calculate_thermal_pressure(fpi_pressure);
    
    %Calculate Total Pressure
    total_pressure = magnetic_pressure + thermal_pressure;% + dynamic_pressure;
    
    subplot(num_plots,1,plot_order);
    
%     plot(fpi_timedata,dynamic_pressure,'LineWidth',1); hold on
%     plot(fpi_timedata,magnetic_pressure,'LineWidth',1);
%     plot(fpi_timedata,thermal_pressure,'LineWidth',1);
%     plot(fpi_timedata,total_pressure,'LineWidth',1);
    
    plot(fpi_timedata,(dynamic_pressure),'LineWidth',1); hold on
    plot(fpi_timedata,(magnetic_pressure),'LineWidth',1);
    plot(fpi_timedata,(thermal_pressure),'LineWidth',1);
    plot(fpi_timedata,(total_pressure),'LineWidth',1);
    
    legend({'Dynamic','Magnetic','Thermal','Total'},'FontSize',10)
    
    legend('boxoff')
    legend('Location','eastoutside')
    colormap('winter');
    datetick
    xlim([fpi_timedata(1) fpi_timedata(end)])
    ylabel({'P';'[nPa]'},'FontSize', 14)
    set(gca, 'XTickLabel', [],'XMinorTick','on','YMinorTick','on','linewidth',1.25)
    
end

