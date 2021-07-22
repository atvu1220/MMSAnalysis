function [n_cs_mean,B_pre,B_post] = calculateCurrentSheet(event_start,event_end,fgm_timedata_srvy,fgm_bdata_srvy,threshold_std)
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
    
    
    
    if threshold_std > 0
        %leading edge mean, upstream of 2 minutes
        %calculate mean
        %calculate variance; (each data points - mean)^2.
        %standard deviation = sqrt(variance);
        %if next point is larger than threshold std, standard deviation, stop.
        interval = 1; %1second interval for means. 0.5 = 30 seconds.
        initial_gap = 3; %number of seconds to start analysis from the summary plot edges, which is start_index and end_index, <5s
        %final_gap =  4; %go back this amount of seconds for the final mean.
        %threshold_std = 3.4;
        
        component = 4; %which bcomponent to analyze. 4=btotal
        leading_leftmost_index = 0;
        leading_rightmost_index = 0;
        
        %Locate beginning index from observation start range to leading shock
        %boundary.
        %start_index_event = start_index+16*initial_gap; %index where event starts.
        
        leading_mean = mean(fgm_bdata_srvy(start_index:start_index+initial_gap*16*interval,4));
        leading_std = std(fgm_bdata_srvy(start_index:start_index+initial_gap*16*interval,4));
        
        for i = initial_gap/interval:30/interval %at most the edge of the event is within 30 seconds of the summary plot
            
            leading_left_index = start_index+16*(i-1)*interval;
            leading_right_index = start_index+16*(i)*interval;
            
            leading_next_mean = mean(fgm_bdata_srvy(leading_left_index:leading_right_index,component));
            
            
            if abs(leading_next_mean-leading_mean) > threshold_std*leading_std %if the diffence between the testing interval mean and the prior mean is greater than the STD, then stop the loop.
                
                %rightmost edge where B is within threshold_std sigma of change to the event
                %edge.
                leading_rightmost_index = leading_left_index-16; %save the index of the rightmost, but minus 1 more interval, to be safe.
                
                break;
                
            elseif abs(leading_next_mean-leading_mean) <= threshold_std*leading_std %if it is less than the STD, then keep going, add it to the interval and calculate new means and STDs
                leading_mean = mean(fgm_bdata_srvy(start_index:leading_right_index,component));
                leading_std = std(fgm_bdata_srvy(start_index:leading_right_index,component));
            end
            
            if i == 30/interval
            leading_rightmost_index = leading_right_index;
            end
        end
        
        
        %Reset means and STD
        leading_mean = mean(fgm_bdata_srvy(leading_rightmost_index-initial_gap*16*interval:leading_rightmost_index,4));
        leading_std = std(fgm_bdata_srvy(leading_rightmost_index-initial_gap*16*interval:leading_rightmost_index,4));
        
        %Revised to restart from rightmost, not start_index
        for i = initial_gap/interval:120/interval %240 half-seconds
            
            leading_left_index = leading_rightmost_index-16*i*interval;
            leading_right_index = leading_rightmost_index-16*(i-1)*interval;
            
            leading_next_mean = mean(fgm_bdata_srvy(leading_left_index:leading_right_index,component));
            
            if abs(leading_next_mean-leading_mean) > threshold_std*leading_std
                
                leading_leftmost_index = leading_right_index+16; %go back one interval, just to be safe.
                B_pre = mean(fgm_bdata_srvy(leading_leftmost_index:leading_rightmost_index,1:3))
                
                break;
                
            elseif abs(leading_next_mean-leading_mean) <= threshold_std*leading_std
                leading_mean = mean(fgm_bdata_srvy(leading_left_index:leading_rightmost_index,component));
                leading_std = std(fgm_bdata_srvy(leading_left_index:leading_rightmost_index,component));
            end
            
            if i==120/interval %if i=120seconds
                leading_leftmost_index = leading_left_index;
                leading_mean = mean(fgm_bdata_srvy(leading_leftmost_index:leading_rightmost_index,component));
                B_pre = mean(fgm_bdata_srvy(leading_leftmost_index:leading_rightmost_index,1:3))
                
            end
        end
        
        
        
        
        
        
        %     for i = 1:120/interval %240 half-seconds
        %
        %         leading_left_index = start_index-16*i*interval;
        %         leading_right_index = start_index-16*(i-1)*interval;
        %
        %         leading_next_mean = mean(fgm_bdata_srvy(leading_left_index:leading_right_index,component));
        %
        %         if abs(leading_next_mean-leading_mean) > threshold_std*leading_std
        %
        %             %leading_mean = mean(fgm_bdata_srvy(leading_right_index:start_index_event+final_gap*16*interval,component));
        %
        %             %leftmost edge where B is within threshold_std sigma of change to the event
        %             %edge.
        %             leading_leftmost_index = leading_right_index; %go back one interval, just to be safe.
        %             B_pre = mean(fgm_bdata_srvy(leading_leftmost_index:leading_rightmost_index,1:3));
        %
        %             break;
        %
        %         elseif abs(leading_next_mean-leading_mean) <= threshold_std*leading_std
        %             leading_mean = mean(fgm_bdata_srvy(leading_left_index:leading_rightmost_index,component));
        %             leading_std = std(fgm_bdata_srvy(leading_left_index:leading_rightmost_index,component));
        %         end
        %
        %         if i==120/interval %if i=120seconds
        %             leading_leftmost_index = leading_left_index;
        %             leading_mean = mean(fgm_bdata_srvy(leading_leftmost_index:leading_rightmost_index,component));
        %             B_pre = mean(fgm_bdata_srvy(leading_leftmost_index:leading_rightmost_index,1:3));
        %
        %         end
        %     end
        %
        
        
        
        
        
        
        
        %trailing edge mean, upstream of up to 2 minutes.
        trailing_rightmost_index = 0; %initialize rightmost index.
        trailing_leftmost_index = 0;
        %     trailing_leftmost_index = end_index-16*initial_gap; %index where event ends.
        
        
        trailing_mean = mean(fgm_bdata_srvy(end_index-initial_gap*16*interval:end_index,component)); %initialize the trailing mean from end of event to the end index of summary plot.
        trailing_std = std(fgm_bdata_srvy(end_index-initial_gap*16*interval:end_index,component)); %initialize the standard deviation of the interval above.
        
        %find the leftmost trailing index
        for i = initial_gap/interval:30/interval %240 half-seconds
            
            trailing_left_index = end_index-16*i*interval;
            trailing_right_index = end_index-16*(i-1)*interval;
            
            trailing_next_mean = mean(fgm_bdata_srvy(trailing_left_index:trailing_right_index,component));
            
            if abs(trailing_next_mean-trailing_mean) > threshold_std*trailing_std
                
                trailing_leftmost_index = trailing_right_index+16; %go back one interval, just to be safe.
                break;
                
            elseif abs(trailing_next_mean-trailing_mean) <= threshold_std*trailing_std
                trailing_mean = mean(fgm_bdata_srvy(trailing_left_index:end_index,component));
                trailing_std = std(fgm_bdata_srvy(trailing_left_index:end_index,component));
            end
            
            if i == 30/interval
            trailing_leftmost_index = trailing_left_index;
            end
        end
        
        %revised for starting at leftmost, not end_index
        trailing_mean = mean(fgm_bdata_srvy(trailing_leftmost_index:trailing_leftmost_index+initial_gap*16*interval,component)); %initialize the trailing mean from end of event to the end index of summary plot.
        trailing_std = std(fgm_bdata_srvy(trailing_leftmost_index:trailing_leftmost_index+initial_gap*16*interval,component)); %initialize the standard deviation of the interval above.
        
        
        
        
        %Trailing mean, after
        for i = initial_gap/interval:120/interval %240 half-seconds
            
            trailing_left_index = trailing_leftmost_index+16*(i-1)*interval; %initialize left index for testing
            trailing_right_index = trailing_leftmost_index+16*i*interval; %initialize right index for testing
            
            trailing_next_mean = mean(fgm_bdata_srvy(trailing_left_index:trailing_right_index,component)); %calculate the mean of B_total between these two left and right indices.
            
            if abs(trailing_next_mean-trailing_mean) > threshold_std*trailing_std %if the diffence between the testing interval mean and the prior mean is greater than the STD, then stop the loop.
                
                %             trailing_mean = mean(fgm_bdata_srvy(trailing_leftmost_index:trailing_left_index-final_gap*16*interval,component)); %calculate mean up to the previous testing interval.
                
                trailing_rightmost_index = trailing_left_index-16; %save the index of the rightmost, but minus 1 more interval, to be safe.
                B_post = mean(fgm_bdata_srvy(trailing_leftmost_index:trailing_rightmost_index,1:3))
                %rightmost edge where B is within threshold_std sigma of change to the event
                %edge.
                
                break;
                
            elseif abs(trailing_next_mean-trailing_mean) <= threshold_std*trailing_std %if it is less than the STD, then keep going, add it to the interval and calculate new means and STDs
                trailing_mean = mean(fgm_bdata_srvy(trailing_leftmost_index:trailing_right_index,component));
                trailing_std = std(fgm_bdata_srvy(trailing_leftmost_index:trailing_right_index,component));
            end
            
            if i ==120/interval %if i=2minutes. we have enough data.
                trailing_rightmost_index = trailing_right_index;
                trailing_mean = mean(fgm_bdata_srvy(trailing_leftmost_index:trailing_rightmost_index,component));
                B_post = mean(fgm_bdata_srvy(trailing_leftmost_index:trailing_right_index,1:3))
                
            end
        end
        
        
        
    else
        %threshold_std is the number of minutes.
        leading_leftmost_index = start_index-16*(-threshold_std)*60;
        if leading_leftmost_index < 0
            leading_leftmost_index = 0;
        end
        
        leading_rightmost_index = start_index+16*5; %5 seconds from the edge of the event.
        
        trailing_leftmost_index = end_index-16*5;
        
        trailing_rightmost_index = end_index+16*(-threshold_std)*60;
        
        if trailing_rightmost_index > length(fgm_bdata_srvy)
            trailing_rightmost_index = length(fgm_bdata_srvy);
        end
        
        B_pre = mean(fgm_bdata_srvy(leading_leftmost_index:leading_rightmost_index,1:3))
        B_post = mean(fgm_bdata_srvy(trailing_leftmost_index:trailing_rightmost_index,1:3))
    end
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    %
    %     %Trailing mean, after
    %     for i = 1:120/interval %240 half-seconds
    %         i
    %         trailing_left_index = end_index+16*(i-1)*interval %initialize left index for testing
    %         trailing_right_index = end_index+16*i*interval %initialize right index for testing
    %
    %         trailing_next_mean = mean(fgm_bdata_srvy(trailing_left_index:trailing_right_index,component)) %calculate the mean of B_total between these two left and right indices.
    %
    %         if abs(trailing_next_mean-trailing_mean) > threshold_std*trailing_std %if the diffence between the testing interval mean and the prior mean is greater than the STD, then stop the loop.
    %
    %             %             trailing_mean = mean(fgm_bdata_srvy(trailing_leftmost_index:trailing_left_index-final_gap*16*interval,component)); %calculate mean up to the previous testing interval.
    %
    %             trailing_rightmost_index = trailing_left_index; %save the index of the rightmost, but minus 1 more interval, to be safe.
    %             B_post = mean(fgm_bdata_srvy(trailing_leftmost_index:trailing_rightmost_index,1:3))
    %             %rightmost edge where B is within threshold_std sigma of change to the event
    %             %edge.
    %
    %             break;
    %
    %         elseif abs(trailing_next_mean-trailing_mean) <= threshold_std*trailing_std %if it is less than the STD, then keep going, add it to the interval and calculate new means and STDs
    %             trailing_mean = mean(fgm_bdata_srvy(trailing_leftmost_index:trailing_right_index,component))
    %             trailing_std = std(fgm_bdata_srvy(trailing_leftmost_index:trailing_right_index,component))
    %         end
    %
    %         if i ==120/interval %if i=2minutes. we have enough data.
    %             trailing_rightmost_index = trailing_right_index;
    %             trailing_mean = mean(fgm_bdata_srvy(trailing_leftmost_index:trailing_rightmost_index,component));
    %             B_post = mean(fgm_bdata_srvy(trailing_leftmost_index:trailing_right_index,1:3));
    %
    %         end
    %     end
    %
    
    
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
    
    title('MMS Current Sheet Magnetic Field Mean Interval ', 'FontSize', 18, 'FontWeight', 'normal')
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
    
    annotation('textbox',[plot_pos(1) plot_pos(2)-plot_pos(4)*plot_gap*(plot_order-5), plot_pos(3), plot_pos(4)],...
        'String',{strcat('n_{cs}: [',num2str(n_cs_mean,'%.4f '),']')},...
        'VerticalAlignment','Top','Edgecolor','none','FontSize', 14);
    
    
    
    %calculate CS normal observation summary window, single point means of
    %the initial mean above - initial gap.
    n_cs = cross(mean(fgm_bdata_srvy(trailing_leftmost_index:end_index,1:3)),mean(fgm_bdata_srvy(start_index:leading_rightmost_index,1:3)));
    
    n_cs = n_cs/norm(n_cs);
    
    %%%Forcing Positive X
    if n_cs(1) < 0
    n_cs = -n_cs; %force x component to be positive, according to Schwartz et al. 2018.
    end
    
    annotation('textbox',[plot_pos(1) plot_pos(2)-plot_pos(4)*plot_gap*(plot_order-3), plot_pos(3), plot_pos(4)],...
        'String',{strcat('n_{cs}: [',num2str(n_cs,'%.4f '),']')},...
        'VerticalAlignment','Top','Edgecolor','none','FontSize', 14);
    
    %B_pre mean
    annotation('textbox',[plot_pos(1)+plot_pos(3)/2-plot_pos(4)/2 plot_pos(2)-plot_pos(4)*plot_gap*(plot_order-4), plot_pos(3), plot_pos(4)],...
        'String',{strcat('B_{pre}: [',num2str(B_pre,'%.4f '),']')},...
        'VerticalAlignment','Top','Edgecolor','none','FontSize', 11);
    
    %B_post mean
    annotation('textbox',[plot_pos(1)+plot_pos(3)/2-plot_pos(4)/2 plot_pos(2)-plot_pos(4)*plot_gap*(plot_order-2), plot_pos(3), plot_pos(4)],...
        'String',{strcat('B_{post}: [',num2str(B_post,'%.4f '),']')},...
        'VerticalAlignment','Top','Edgecolor','none','FontSize', 11);
    %plotting
    orient(gcf,'landscape')
    plot_name =  strcat('2_Current_Sheet_Normal',...
        event_start(1:19),'_',num2str(threshold_std),'.pdf');
    print(gcf, '-dpdf', '-opengl', plot_name,'-fillpage');
end

