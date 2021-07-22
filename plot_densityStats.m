function [substructurePresent] = plot_densityStats(date_start,date_end,left_InnerEdge,...
        right_InnerEdge,fpi_timedata,fpi_ndata,num_plots,plot_order,ssSTDs,minSSDataPoints)
    %Computes the trend line for the core density, then colors in the area of density above the
    %trendline
    
    
    %ssSTDs = 2;
    %minSSDataPoints = 4;
    
    
    
    formatIn='yyyy-mm-dd HH:MM:SS.FFF';
    tstart = datenum(date_start,formatIn);
    tend = datenum(date_end,formatIn);
    
    
    %Find the start and end limits of the event in the data
    start_index = find(fpi_timedata > tstart, 1);
    end_index = find(fpi_timedata > tend, 1);
    
    %crop data
    fpi_timedata = fpi_timedata(start_index:end_index,1);
    fpi_ndata = fpi_ndata(start_index:end_index,:);
   
    
   [substructurePresent,fpi_coreTimedata,fpi_coreNdata,trendDensity,Sigma,substructureStarts,substructureEnds] = calculate_densityStats(fpi_timedata,fpi_ndata,left_InnerEdge,right_InnerEdge,ssSTDs,minSSDataPoints);


    
    substructureIndices = arrayfun(@(s, e) s:e-1, substructureEnds, substructureStarts, 'UniformOutput', false);
    substructureIndices = [substructureIndices{:}];
    [yRange] = ylim; %For base of fill polygon, 0 for linear, and ymin for log plots.
    %substructure = yRange(1).*ones(length(fpi_coreNdata),1);
   substructure = zeros(length(fpi_coreNdata),1);
    substructure(substructureIndices) = 1;
    %substructure(substructureIndex) = 1;
    
    
    subplot(num_plots,1,plot_order); %fill_between_lines = @(X,Y1,Y2,C) fill( [X fliplr(X)],  [Y1 fliplr(Y2)], C );
    
    %plot(fpi_timedata,log10(fpi_ndata))
    
    
    %fill(fpi_coreTimedata,trendDensity+Sigma,'g')
    %area(fpi_coreTimedata,trendDensity+Sigma,'g')
    %fill(fpi_coreTimedata,fpi_coreNdata,'r')
    
    %Fill area with Red of substructure
    
    fillArea = ((fpi_coreNdata).*substructure); %log here
    fillArea(fillArea==0) = yRange(1);
    area(fpi_coreTimedata,fillArea,yRange(1))
    
    
    
    
    
    
    %Plot trendline
    plot(fpi_coreTimedata,(trendDensity),'k','linewidth', 1.5); %log here


    plot(fpi_coreTimedata,(trendDensity+ssSTDs*Sigma),'--k', 'linewidth', 0.75) %log here
    %fill(fpi_coreTimedata,[0;fpi_coreNdata(2:end-1);0],'r')
    %area(fpi_coreTimedata,trendDensity+Sigma,'FaceColor','white','EdgeColor','w')
    %h=get(gca,'Children');
    %set(gca,'Children',[h(4) h(3) h(1) h(2)])
    %legend({'Substructure','','N^i','N^e',},'FontSize',10)
    
    legend({'N^i','N^e','Substructure'},'FontSize',10)
    legend('boxoff')
    legend('Location','eastoutside')
    colormap('winter');
    datetick
    xlim([fpi_timedata(1) fpi_timedata(end)])
    ylabel({'n';'[cm^{-3}]'},'FontSize', 14)
    set(gca, 'XTickLabel', [],'XMinorTick','on','YMinorTick','on','linewidth',1.25)
    
end

