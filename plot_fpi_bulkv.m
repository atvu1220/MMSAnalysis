function [] = plot_fpi_bulkv(date_start,date_end,fpi_timedata,fpi_vdata,num_plots,plot_order)
    %plots the plasma parameter with the given start date, end date, datas,
    %number of plots and plot order
    
    formatIn='yyyy-mm-dd HH:MM:SS.FFF';
    tstart = datenum(date_start,formatIn);
    tend = datenum(date_end,formatIn);
    
    
    %Find the start and end limits of the event in the data
    start_index = find(fpi_timedata > tstart, 1);
    end_index = find(fpi_timedata > tend, 1);
    
    subplot(num_plots,1,plot_order);
    
    fpi_timedata = fpi_timedata(start_index:end_index,1);
    
    
    
    
    
    
    
    if size(fpi_vdata,2) == 1
        fpi_vdata=(fpi_vdata(start_index:end_index));
        plot(fpi_timedata,fpi_vdata,'LineWidth',1)
        legend({'v'},'FontSize',8)
    else
        %crop data
        fpi_vxdata=(fpi_vdata(start_index:end_index,1));
        fpi_vydata=(fpi_vdata(start_index:end_index,2));
        fpi_vzdata=(fpi_vdata(start_index:end_index,3));
        plot(fpi_timedata,[fpi_vxdata,fpi_vydata,fpi_vzdata],'LineWidth',1)
        legend({'V_x','V_y','V_z'},'FontSize',10)
    end
    
    
    
    legend('boxoff')
    legend('Location','eastoutside')
    colormap('winter');
    datetick
    xlim([fpi_timedata(1) fpi_timedata(end)])
    if (fpi_timedata(2) - fpi_timedata(1) < 1e-6)
        ylabel({'v_{e}';'[km/s]'},'FontSize', 14)
    else
        ylabel({'v_{ion}';'[km/s]'},'FontSize', 14)
    end
    
    set(gca, 'XTickLabel', [],'XMinorTick','on','YMinorTick','on','linewidth',1.25)
    
    
    
    
end

