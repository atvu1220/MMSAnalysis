function [] = plot_fpi_number(date_start,date_end,fpi_timedata,fpi_ndata,num_plots,plot_order,specie)
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
    
    
    subplot(num_plots,1,plot_order);
    
    
    %plot(fpi_timedata,log10(fpi_ndata))
    plot(fpi_timedata,fpi_ndata,'LineWidth',1)
    if strcmp(specie, 'i')
        legend({'N^i'},'FontSize',10)
    elseif strcmp(specie,'se')
        legend({'N^i','Substructure','N^e'},'FontSize',10)
    else
        legend({'N^i','N^e'},'FontSize',10)
    end
    %     legend({'N^i','N^e'},'FontSize',14)
    legend('boxoff')
    legend('Location','eastoutside')
    colormap('winter');
    datetick
    xlim([fpi_timedata(1) fpi_timedata(end)])
    ylabel({'n';'[cm^{-3}]'},'FontSize', 14)
    set(gca, 'XTickLabel', [],'XMinorTick','on','YMinorTick','on','linewidth',1.25)
    
end

