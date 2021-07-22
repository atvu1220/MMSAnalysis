function [n_cs_mean,B_pre,B_post] = manualCurrentSheet(event_start,event_end,leading_leftmost_date, leading_rightmost_date,...
        trailing_leftmost_date,trailing_rightmost_date,fgm_timedata_srvy,fgm_bdata_srvy)
    %Calculates the current sheet normal by mean magnetic field before and
    %after the event. Automatically determines the time interval for mean
    %before and after event by examining how the magnetic field changes.
    %threshold_std > 0 is how much to include based on the total magnetic
    %field averages, smaller is more restrictive.
    %if threshold_std is negative, force the number of minutes before/after
    %event.
    %uses srvy data, thus 16 points per second.
    figure('Position',[0 0 1050 750])
    co = [0 0 1;
        0 1 0;
        1 0 0;
        0 0 0];
    set(gcf,'defaultAxesColorOrder',co)
    set(gcf,'color','w');
    plot_gap=1.25;
    num_plots=4;
    plot_order = 1;
    
    %get indices of event
    [~,~,start_index,end_index] = crop(fgm_timedata_srvy,fgm_bdata_srvy,event_start,event_end);
    [~,~,leading_leftmost_index,leading_rightmost_index] = crop(fgm_timedata_srvy,fgm_bdata_srvy,leading_leftmost_date,leading_rightmost_date);
    [~,~,trailing_leftmost_index,trailing_rightmost_index] = crop(fgm_timedata_srvy,fgm_bdata_srvy,trailing_leftmost_date,trailing_rightmost_date);
    
    B_pre = mean(fgm_bdata_srvy(leading_leftmost_index:leading_rightmost_index,1:3))
    B_post = mean(fgm_bdata_srvy(trailing_leftmost_index:trailing_rightmost_index,1:3))
    
    
    %plot entire range
    formatOut = "yyyy-mm-dd HH:MM:SS.FFF";
    %datestr(mms1_fgm_timedata_srvy(1),formatOut)
    subplot(num_plots,plot_order,1)
    %plot(fgm_timedata_srvy,fgm_bdata_srvy);
    plot_fgm_magnetic(datestr(fgm_timedata_srvy(1),formatOut),datestr(fgm_timedata_srvy(end),formatOut),...
        fgm_timedata_srvy,fgm_bdata_srvy,num_plots,1)
    hold on
    line([(fgm_timedata_srvy(leading_leftmost_index)) (fgm_timedata_srvy(leading_leftmost_index))],...
        get(gca,'YLim'),'Color','k','LineWidth',0.5,'LineStyle','--')
    
    line([(fgm_timedata_srvy(leading_rightmost_index)) (fgm_timedata_srvy(leading_rightmost_index))],...
        get(gca,'YLim'),'Color','k','LineWidth',0.5,'LineStyle','--')
    
    line([(fgm_timedata_srvy(trailing_leftmost_index)) (fgm_timedata_srvy(trailing_leftmost_index))],...
        get(gca,'YLim'),'Color','k','LineWidth',0.5,'LineStyle','--')
    
    line([(fgm_timedata_srvy(trailing_rightmost_index)) (fgm_timedata_srvy(trailing_rightmost_index))],...
        get(gca,'YLim'),'Color','k','LineWidth',0.5,'LineStyle','--')
    
    hold off
    legend({'B_x', 'B_y', 'B_z','B_t'},'FontSize',14)
    
    title('MMS1 FGM Survey Magnetic Fields', 'FontSize', 18, 'FontWeight', 'normal')
    datetick('keeplimits')
    plot_pos = get(gca,'Position');
    set(gca,'Position',[plot_pos(1), plot_pos(2)-plot_pos(4)*plot_gap*(plot_order-1), plot_pos(3), plot_pos(4)]);
    plot_order = plot_order+1;
    
    
    
    %plot pre event
    plot_fgm_magnetic(datestr(fgm_timedata_srvy(leading_leftmost_index),formatOut),datestr(fgm_timedata_srvy(leading_rightmost_index),formatOut),fgm_timedata_srvy,fgm_bdata_srvy,num_plots,plot_order)
    datetick('keeplimits')
    set(gca,'Position',[plot_pos(1), plot_pos(2)-plot_pos(4)*plot_gap*(plot_order-1), plot_pos(3), plot_pos(4)]);
    plot_order = plot_order+1;
    
    
    %plot event range.
    plot_fgm_magnetic(event_start,event_end,fgm_timedata_srvy,fgm_bdata_srvy,num_plots,plot_order)
    hold on
    line([(fgm_timedata_srvy(leading_leftmost_index)) (fgm_timedata_srvy(leading_leftmost_index))],...
        get(gca,'YLim'),'Color','k','LineWidth',0.5,'LineStyle','--')
    
    line([(fgm_timedata_srvy(leading_rightmost_index)) (fgm_timedata_srvy(leading_rightmost_index))],...
        get(gca,'YLim'),'Color','k','LineWidth',0.5,'LineStyle','--')
    
    line([(fgm_timedata_srvy(trailing_leftmost_index)) (fgm_timedata_srvy(trailing_leftmost_index))],...
        get(gca,'YLim'),'Color','k','LineWidth',0.5,'LineStyle','--')
    
    line([(fgm_timedata_srvy(trailing_rightmost_index)) (fgm_timedata_srvy(trailing_rightmost_index))],...
        get(gca,'YLim'),'Color','k','LineWidth',0.5,'LineStyle','--')
    hold off
    legend({'B_x', 'B_y', 'B_z','B_t'},'FontSize',14)
    datetick('keeplimits')
    set(gca,'Position',[plot_pos(1), plot_pos(2)-plot_pos(4)*plot_gap*(plot_order-1), plot_pos(3), plot_pos(4)]);
    plot_order = plot_order+1;
    
    
    %plot post event
    plot_fgm_magnetic(datestr(fgm_timedata_srvy(trailing_leftmost_index),formatOut),datestr(fgm_timedata_srvy(trailing_rightmost_index),formatOut),fgm_timedata_srvy,fgm_bdata_srvy,num_plots,plot_order)
    datetick('keeplimits')
    set(gca,'Position',[plot_pos(1) plot_pos(2)-plot_pos(4)*plot_gap*(plot_order-1), plot_pos(3), plot_pos(4)]);
    plot_order = plot_order+1;
    
    
    %calculate CS normal from means
    n_cs_mean = cross(B_post,B_pre);
    
    n_cs_mean = n_cs_mean/norm(n_cs_mean);
    
    %%%Forcing Positive X
    if n_cs_mean(1) < 0
        n_cs_mean = -n_cs_mean; %force x component to be positive, according to Schwartz et al. 2018.
    end
    
    %n_cs_mean(1) = abs(n_cs_mean(1)); %force x component to be positive, according to Schwartz et al. 2018.
    
    
    
    
    
    
    %calculate CS normal observation summary window, single point means of
    %the initial mean above - initial gap.
    n_cs = cross(mean(fgm_bdata_srvy(end_index-5*16:end_index,1:3)),mean(fgm_bdata_srvy(start_index:start_index+5*16,1:3)));
    
    n_cs = n_cs/norm(n_cs);
    
    %%%Forcing Positive X
    if n_cs(1) < 0
        n_cs = -n_cs; %force x component to be positive, according to Schwartz et al. 2018.
    end
     
    
    
    
%     %% Dsiplay Values
%     annotation('textbox',[plot_pos(1) plot_pos(2)-plot_pos(4)*plot_gap*(plot_order-5), plot_pos(3), plot_pos(4)],...
%         'String',{strcat('n_{cs}: [',num2str(n_cs,'%.4f '),']')},...
%         'VerticalAlignment','Top','Edgecolor','none','FontSize', 14);
%     
%     %     dispaly current sheet using B over the entire interval of 30 minutes.
%     annotation('textbox',[plot_pos(1) plot_pos(2)-plot_pos(4)*plot_gap*(plot_order-3), plot_pos(3), plot_pos(4)],...
%         'String',{strcat('n_{cs}: [',num2str(n_cs_mean,'%.4f '),']')},...
%         'VerticalAlignment','Top','Edgecolor','none','FontSize', 14);
%     
%     %B_pre mean
%     annotation('textbox',[plot_pos(1)+plot_pos(3)/2-plot_pos(4)/2 plot_pos(2)-plot_pos(4)*plot_gap*(plot_order-4), plot_pos(3), plot_pos(4)],...
%         'String',{strcat('B_{pre}: [',num2str(B_pre,'%.4f '),']')},...
%         'VerticalAlignment','Top','Edgecolor','none','FontSize', 11);
%     
%     %B_post mean
%     annotation('textbox',[plot_pos(1)+plot_pos(3)/2-plot_pos(4)/2 plot_pos(2)-plot_pos(4)*plot_gap*(plot_order-2), plot_pos(3), plot_pos(4)],...
%         'String',{strcat('B_{post}: [',num2str(B_post,'%.4f '),']')},...
%         'VerticalAlignment','Top','Edgecolor','none','FontSize', 11);
    
    
    
    %plotting
    orient(gcf,'landscape')
    plot_name =  strcat('2_Current_Sheet_Normal',...
        event_start(1:19),'_','manual','.pdf');
    print(gcf, '-dpdf', '-opengl', plot_name,'-fillpage');
end

