function [] = plot_parameter(date_start,date_end,time,data,num_plots,plot_order,legend_string,yaxis_string)
    %plots the general plasma parameter with the given start date, end date, datas,
    %number of plots and plot order
    
    formatIn='yyyy-mm-dd HH:MM:SS.FFF';
    tstart = datenum(date_start,formatIn);
    tend = datenum(date_end,formatIn);
    
    
    %Find the start and end limits of the event in the data
    start_index = find(time >= tstart, 1);
    end_index = find(time >= tend, 1);
    
    %crop data
    time = time(start_index:end_index,1);
    data = data(start_index:end_index,:);
    
    
    subplot(num_plots,1,plot_order);
    
    
    plot(time,data(:,:),'LineWidth',1)
    if isempty(legend_string)
    else
        legend(legend_string,'FontSize',14)
        legend('boxoff')
        legend('Location','eastoutside')
    end
    colormap('winter');
    datetick
    xlim([time(1) time(end)])
    if isempty(yaxis_string)
    else  
    ylabel(yaxis_string,'FontSize', 14)
    end
    set(gca, 'XTickLabel', [],'XMinorTick','on','YMinorTick','on','linewidth',1.25)
    
    
end

