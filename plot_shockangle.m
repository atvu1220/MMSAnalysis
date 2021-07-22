function [] = plot_shockangle(date_start,date_end,fgm_timedata,fgm_bdata,num_plots,plot_order,shock_normal)
    %plots the plasma parameter with the given start date, end date, datas,
    %number of plots and plot order

    formatIn='yyyy-mm-dd HH:MM:SS.FFF';
    tstart = datenum(date_start,formatIn);
    tend = datenum(date_end,formatIn);
    
    
    %Find the start and end limits of the event in the data
    start_index = find(fgm_timedata > tstart, 1);
    end_index = find(fgm_timedata > tend, 1);
    
    %crop data
    fgm_timedata = fgm_timedata(start_index:end_index,1);
    fgm_bdata = fgm_bdata(start_index:end_index,:);
    
    
    %arccos(Bu*n/Bu)
    %produces a column of angles
    shock_angle = acosd(abs((sum(fgm_bdata(:,1:3).*shock_normal./vecnorm(fgm_bdata(:,1:3),2,2),2))));
    
    
    
    subplot(num_plots,1,plot_order);
    

    plot(fgm_timedata,shock_angle,'LineWidth',1)
    line([fgm_timedata(1),fgm_timedata(end)],[45,45],'Color','k','LineStyle','--')
    line([fgm_timedata(1),fgm_timedata(end)],[135,135],'Color','k','LineStyle','--')
    %legend({'B_x', 'B_y', 'B_z','B_t'},'FontSize',14)
    %legend('boxoff')
    %legend('Location','eastoutside')
    colormap('winter');
    datetick
    xlim([fgm_timedata(1) fgm_timedata(end)])
    ylabel({'\theta_{Bn}';'[Deg]'},'FontSize', 14)
    set(gca, 'XTickLabel', [],'XMinorTick','on','YMinorTick','on','linewidth',1.25)
    yticks([0 45 90])
    ylim([0 90])
    
    
end

