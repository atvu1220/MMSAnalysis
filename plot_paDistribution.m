function [] = plot_paDistribution(event_start,event_end,pa_n,specie,energyChannel,num_plots,plot_order)
    %Function file to calculate Pitch Angle Distributions
    
    [time,pa_bins_values,f_bins,centers] = calculate_paDistribution(event_start,event_end,pa_n,specie,energyChannel);
    
    %Calculte Energy Channel Edges given center value
    d = diff(centers)/2; %differences between bins
    edges = [centers(1)-d(1), centers(1:end-1)+d, centers(end)+d(end)]; %Calculate edges
    energyRange = [edges(energyChannel(1));edges(energyChannel(end)+1)];
    subplot(num_plots,1,plot_order);
    
    pcolor(time,pa_bins_values',log10(f_bins))
    
%     imagesc(time,pa_bins_values',log10(f_bins))
    

    
    shading interp
    whitejet = [1 1 1; jet];
    colormap(whitejet)
    
    xlim([time(1) time(end)])
    %ylabel({'Pitch Angle';sprintf('E %1.f,%1.f,%1.f,%1.f',energyChannel)},'FontSize', 14)
%     ylabel({sprintf('E %1.f,%1.f,%1.f,%1.f',energyChannel)},'FontSize', 10)
%     set(gca, 'XTickLabel', [],'XMinorTick','on','YMinorTick','on','linewidth',1.25)
    datetick('x','keeplimits')
    ylabel({'Pitch Angle';sprintf('%2.f - %2.f eV',energyRange(1),energyRange(end))},'FontSize', 9)
    set(gca,'Ydir','normal', 'XTickLabel', [],'YMinorTick','on','XMinorTick','on','layer','top','linewidth',1.25)
    