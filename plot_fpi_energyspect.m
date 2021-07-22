function [] = plot_fpi_energyspect(date_start,date_end,fpi_timedata,fpi_edata,fpi_espectdata,num_plots,plot_order)
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
    fpi_edata=fpi_edata(start_index:end_index,:);
    fpi_espectdata=log10(fpi_espectdata(start_index:end_index,:));
    
    
    subplot(num_plots,1,plot_order);
    
    
    imagesc(fpi_timedata,log10(fpi_edata(1,:)'),fpi_espectdata')
    color_bar = colorbar('Ticks', [3, 4, 5, 6, 7,8],...
        'TickLabels', {'10^3', '10^4', '10^5', '10^6', '10^7','10^8'},'FontSize', 10);
    ylabel(color_bar,{'keV/cm^2 s sr keV'},'FontSize', 8)
    
    %     legend('boxoff')
    %     legend('Location','eastoutside')
    shading interp
    whitejet = [1 1 1; jet];
    colormap(whitejet);
    datetick
    xlim([fpi_timedata(1) fpi_timedata(end)])
    ylabel({'Energy';'[eV]'},'FontSize', 14)
    set(gca,'Ydir','normal', 'XTickLabel', [],'YMinorTick','on','XMinorTick','on','Yticklabels',[10 100 1000 10000],'layer','top','linewidth',1.25)
    yticks([1,2,3,4])
    
    
end

% l4 = legend(f4,'show');
% l4pos= get(l4, 'Position');
%set(l4, 'Position',[plot_pos(1)+plot_pos(3)+0.01 plot_pos(2)-plot_pos(4)*3.5 l4pos(3) l4pos(4)])
