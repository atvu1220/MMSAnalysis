function [] = plot_fgm_magnetic(date_start,date_end,fgm_timedata,fgm_bdata,num_plots,plot_order)
    %plots the plasma parameter with the given start date, end date, datas,
    %number of plots and plot order

    formatIn='yyyy-mm-dd HH:MM:SS.FFF';
    tstart = datenum(date_start,formatIn);
    tend = datenum(date_end,formatIn);
    
    
    %Find the start and end limits of the event in the data
    start_index = find(fgm_timedata >= tstart, 1);
    end_index = find(fgm_timedata >= tend, 1);
    
    %crop data
    fgm_timedata = fgm_timedata(start_index:end_index,1);
    fgm_bdata = fgm_bdata(start_index:end_index,:);
    
    
    subplot(num_plots,1,plot_order);
    

    plot(fgm_timedata,fgm_bdata(:,:),'LineWidth',1)
    if size(fgm_bdata,2) == 1
        legend({'B_t'},'FontSize',8)
    else
        legend({'B_x', 'B_y', 'B_z','B_t'},'FontSize',8)
    end
        
    
    legend('boxoff')
    legend('Location','eastoutside')
    colormap('winter');
    datetick
    xlim([fgm_timedata(1) fgm_timedata(end)])
    ylabel({'B';'[nT]'},'FontSize', 14)
    set(gca, 'XTickLabel', [],'XMinorTick','on','YMinorTick','on','linewidth',1.25)
    
    
end

