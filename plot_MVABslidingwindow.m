function [] = plot_MVABslidingwindow(event_start,event_end,fgm_timedata, fgm_bdata,currentsheet_n,minratio)
    %Andrew Vu 9/25/18
    %cd '~/Library/Mobile Documents/com~apple~CloudDocs/Research/Analysis'
    figure('Position',[1 1 1600 400])
    set(gcf,'color','w');
    co = [0 0 1;
        0 1 0;
        1 0 0;
        0 0 0];
    set(gcf,'defaultAxesColorOrder',co)
    
    probe_num = '1';
    num_plots = 2;
    
    data_type = 'srvy';
    plot_order = 1;
    plot_gap=1;
    window_lower_limit = 4;
    data_points = 16*8; %total on both sides. %16 points per s.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Data Retrieval%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %Find the start and end limits of the event in the data
    [~,~,start_index,end_index] = crop(fgm_timedata,fgm_bdata,event_start,event_end);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%Magnetic Field Data%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    plot_fgm_magnetic(event_start,event_end,fgm_timedata,fgm_bdata,num_plots,plot_order)
    hold on
    
    %Display borders for largest window size on figure
    mid = start_index+round(length(fgm_timedata(start_index:end_index))/2);
    line([(fgm_timedata(mid-data_points/2)) (fgm_timedata(mid-data_points/2))],...
        get(gca,'YLim'),'Color','k','LineWidth',0.5,'LineStyle','--')
    line([(fgm_timedata(mid+data_points/2)) (fgm_timedata(mid+data_points/2))],...
        get(gca,'YLim'),'Color','k','LineWidth',0.5,'LineStyle','--')
    hold off
    
    legend({'B_x', 'B_y', 'B_z','B_{tot}'},'FontSize',14)
    legend('boxoff');
    legend('Location','eastoutside');
    
    %Title
    title_name = strcat('MMS1 MVAB Sliding Window', ',\Deltat_{max}=', num2str(data_points/16),'s');
    title(title_name, 'FontSize', 18, 'FontWeight', 'normal')
    
    %Grab Position of Plot
    plot_pos = get(gca,'Position');
    set(gca,'Position',[plot_pos(1), plot_pos(2)-plot_pos(4)*plot_gap*(plot_order-1), plot_pos(3), plot_pos(4)]);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%Loop over sliding window%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %arbitrary scope set dt_0
    window_range = end_index - start_index+1;
    %srvy has 16 points per second, so 64 points is 4 seconds total
    %data poitns is total, so each side is data_points/2
    window_size = data_points;
    
    nested_lambda=cell(window_range,ceil(window_size/2));
    nested_normal=cell(window_range,ceil(window_size/2));
    nested_time=cell(window_range,ceil(window_size/2));
    nested_lmidmin=zeros(window_range,ceil(window_size/2));
    
    nested_n_min_x=zeros(window_range,ceil(window_size/2));
    nested_n_min_y=zeros(window_range,ceil(window_size/2));
    nested_n_min_z=zeros(window_range,ceil(window_size/2));
    
    t_0_vector = zeros(window_range,1);
    data_point_vector = zeros(1,ceil(window_size/2),1);
    
    for i=1:window_range
        %i is the starting point index
        t_0 = fgm_timedata(i+start_index-1);
        t_0_vector(i) = t_0;
        %window_lower_limit is how many of the first points to skip
        for j=window_lower_limit:ceil(window_size/2)
            data_point_vector(j) = j;
            %j is the number of data points to one side of the center point.
            left_index = (i+start_index)-j;
            right_index = (i+start_index)+j;
            
            bdata_scope = fgm_bdata(left_index:right_index,:);
            
            [~, l, v] = mvab(bdata_scope(:,1:3));
            
            nested_lmidmin(i,j)= l(2)/l(3);
            
            nested_n_min_x(i,j)=v(1,3);
            nested_n_min_y(i,j)=v(2,3);
            nested_n_min_z(i,j)=v(3,3);
            
            nested_lambda{i,j} = l;
            nested_normal{i,j} = v(:,3);
            nested_time{i,j}=t_0;
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Plot 3D image %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    plot_order = plot_order+1;
    subplot(num_plots,1,plot_order);
    
    %Plot one component of the MinVAB normal only where the l2/l3 ratio is
    %greater than 10.
% % % % %             nested_lmidmin(nested_lmidmin < 10 ) = NaN; %Anything below an eigenvalue ratio of less than 10, let it be NaN, and therefore not plotted (colored)
% % % % %             nested_n_min_x_10plus = ~isnan(nested_lmidmin).*nested_n_min_x; %where it is above 10, thus not NaN, let it have color.
% % % % %             nested_n_min_x_10plus(nested_n_min_x_10plus == 0) = NaN; %Anywhere that's zero, let it be white according to colormap
% % % % %             imagesc(fgm_timedata(start_index:end_index,:),2*data_point_vector/16,abs(nested_n_min_x_10plus'))
% % % % % 
% % % % %     c=colorbar;
% % % % % ylabel(c,{'N_x'},'FontSize',14)

    
    nested_lmidmin(nested_lmidmin < minratio ) = NaN;
    nested_n_min_x_10plus = ~isnan(nested_lmidmin).*nested_n_min_x;
    nested_n_min_x_10plus(nested_n_min_x_10plus == 0) = NaN;
    nested_n_min_y_10plus = ~isnan(nested_lmidmin).*nested_n_min_y;
    nested_n_min_y_10plus(nested_n_min_y_10plus == 0) = NaN;
    nested_n_min_z_10plus = ~isnan(nested_lmidmin).*nested_n_min_z;
    nested_n_min_z_10plus(nested_n_min_z_10plus == 0) = NaN;
    nested_angles = zeros(window_range,ceil(window_size/2));
    
    for i=1:window_range
        for j=window_lower_limit:ceil(window_size/2)
            point_normal = [nested_n_min_x_10plus(i,j), nested_n_min_y_10plus(i,j), nested_n_min_z_10plus(i,j)];
            %nested_angles(i,j) = rad2deg(acos(dot(point_normal,n_cs)/(norm(point_normal)*norm(n_cs))));
            if isnan(point_normal) == [1,1,1]
                nested_angles(i,j) = NaN;
            else
                [angle,~] = ambig_angle(point_normal,currentsheet_n);
                nested_angles(i,j) = angle;
            end
         
        end
    end
    
    imagesc(fgm_timedata(start_index:end_index,:),2*data_point_vector/16,abs(nested_angles'))
c=colorbar;
ylabel(c,{'Normal Angle with N_{cs}'},'FontSize',14)


    %plot the angle with the current sheet.

    set(gca,'Ydir','reverse')
    set(gca,'layer','top')
    set(gcf,'color','w');
    box on
    cmap = jet;
    cmap = [1 1 1; cmap];
    colormap(cmap)
    
    xlim([fgm_timedata(start_index) fgm_timedata(end_index)])
    datetick('x','keeplimits')
    set(gca,'TickDir','in');
    xtickangle(-40)
    
    ylabel({'Window Size';'[Seconds]'},'FontSize', 14)
    xlabel(strcat('t_0 position for \lambda >',num2str(minratio)),'FontSize', 14)
    set(gca,'XMinorTick','on','YMinorTick','on','linewidth',1.25);
    set(gca,'Position',[plot_pos(1), plot_pos(2)-plot_pos(4)*plot_gap*(plot_order-1), plot_pos(3), plot_pos(4)]);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    orient(gcf,'landscape')
    plot_name =  strcat('3_SlidingWindow_MVA',...
        event_start(1:19),'_',num2str(data_points),'.pdf');
    print(gcf, '-dpdf', '-opengl', plot_name,'-fillpage');
end

