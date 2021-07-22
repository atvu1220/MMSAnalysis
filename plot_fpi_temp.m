function [] = plot_fpi_temp(date_start,date_end,fpi_timedata,fpi_tparadata,fpi_tperpdata,num_plots,plot_order,specie)
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
    fpi_tparadata=(fpi_tparadata(start_index:end_index));
    fpi_tperpdata=(fpi_tperpdata(start_index:end_index));
    
    
    subplot(num_plots,1,plot_order);
    
    
    plot(fpi_timedata,[fpi_tparadata,fpi_tperpdata])
    if strcmp(specie,'i')
        legend({'T^{i}_{\mid\mid}','T^{i}_{\perp}'},'FontSize',8)
    else
        
        legend({'T^{i}_{\mid\mid}','T^{i}_{\perp}','T^{e}_{\mid\mid}','T^{e}_{\perp}'},'FontSize',8)
    end
    %     legend({'T^{i}_{\mid\mid}','T^{i}_{\perp}','T^{e}_{\mid\mid}','T^{e}_{\perp}'},'FontSize',12)
    legend('boxoff')
    legend('Location','eastoutside')
    colormap('winter');
    datetick
    xlim([fpi_timedata(1) fpi_timedata(end)])
    ylabel({'Temp';'[eV]'},'FontSize', 14)
    
    set(gca, 'YScale','log','XTickLabel', [],'XMinorTick','on','YMinorTick','on','Ytick',[0 10 100 1000],'linewidth',1.25)
    
    
    
end

% l4 = legend(f4,'show');
% l4pos= get(l4, 'Position');
%set(l4, 'Position',[plot_pos(1)+plot_pos(3)+0.01 plot_pos(2)-plot_pos(4)*3.5 l4pos(3) l4pos(4)])
