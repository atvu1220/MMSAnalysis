function [] = plot_fpi_dynamic_pressure(date_start,date_end,fpi_timedata,fpi_ndata,fpi_vdata,num_plots,plot_order,specie)
    %plots the plasma parameter with the given start date, end date, datas,
    %number of plots and plot order
    
    formatIn='yyyy-mm-dd HH:MM:SS.FFF';
    tstart = datenum(date_start,formatIn);
    tend = datenum(date_end,formatIn);
    
    
    %Find the start and end limits of the event in the data
    start_index = find(fpi_timedata > tstart, 1);
    end_index = find(fpi_timedata > tend, 1);
    
    %crop data
    fpi_timedata = fpi_timedata(start_index:end_index,1);
    fpi_ndata = fpi_ndata(start_index:end_index,:);
    fpi_vdata = fpi_vdata(start_index:end_index,:);
    
    %Calculate dynamic pressure rho*v^2
    dynamic_pressure = calculate_dynamic_pressure(fpi_ndata,fpi_vdata,specie);
    
    subplot(num_plots,1,plot_order);
    
    
%     semilogy(fpi_timedata,dynamic_pressure)% for difference in i and e
    plot(fpi_timedata,dynamic_pressure,'LineWidth',1)
    if strcmp(specie, 'i')
        legend({'P^i_{ram}'},'FontSize',10)
    else
        legend({'P^i_{ram}','P^e_{ram}'},'FontSize',10)
    end
    %     legend({'N^i','N^e'},'FontSize',14)
    legend('boxoff')
    legend('Location','eastoutside')
    colormap('winter');
    datetick
    xlim([fpi_timedata(1) fpi_timedata(end)])
    ylabel({'P';'[nPa]'},'FontSize', 14)
    set(gca, 'XTickLabel', [],'XMinorTick','on','YMinorTick','on','linewidth',1.25)
    
end

