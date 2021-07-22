function [] = plot_JxB(date_start,date_end,fpi_e_timedata,fpi_i_timedata,fpi_e_ndata,fpi_e_vdata,fpi_i_vdata,fgm_timedata,fgm_bdata,num_plots,plot_order)
    %plots the plasma parameter with the given start date, end date, datas,
    %number of plots and plot order
    
    formatIn='yyyy-mm-dd HH:MM:SS.FFF';
    tstart = datenum(date_start,formatIn);
    tend = datenum(date_end,formatIn);
    
    
    %Find the start and end limits of the event in the data
    start_index = find(fpi_i_timedata > tstart, 1);
    end_index = find(fpi_i_timedata > tend, 1);
    
    subplot(num_plots,1,plot_order);
    
    %Down-Sample B to FPI
    [~,fgm_bdata] = interpxyz(fgm_timedata,fgm_bdata,fpi_i_timedata);
    
    
    [current] = calculate_current(fpi_e_timedata,fpi_i_timedata,fpi_e_ndata,fpi_e_vdata,fpi_i_vdata);
    
    fpi_i_timedata = fpi_i_timedata(start_index:end_index,1);

    %Calculate JxB
    JxB = cross(current,fgm_bdata(:,1:3));
    
    %crop data
    JxB = JxB(start_index:end_index,:);
    plot(fpi_i_timedata,JxB,'LineWidth',1)
    legend({'(JxB)_x','(JxB)_y','(JxB)_z'},'FontSize',10)
    
    legend('boxoff')
    legend('Location','eastoutside')
    colormap('winter');
    datetick
    xlim([fpi_i_timedata(1) fpi_i_timedata(end)])
    ylabel({'JxB';'[nT mA/km^2]'},'FontSize', 12)
    
    
    set(gca, 'XTickLabel', [],'XMinorTick','on','YMinorTick','on','linewidth',1.25)
    
    
    
    
end

