function [time_ranges] = ...
        plot_Timingslidingwindow(event_start,event_end,...
        mms1_fgm_timedata_raw,mms1_fgm_bdata_raw,...
        mms2_fgm_timedata_raw,mms2_fgm_bdata_raw,...
        mms3_fgm_timedata_raw,mms3_fgm_bdata_raw,...
        mms4_fgm_timedata_raw,mms4_fgm_bdata_raw,...
        mms1_mec_timedata_raw,mms1_mec_rdata_raw,...
        mms2_mec_timedata_raw,mms2_mec_rdata_raw,...
        mms3_mec_timedata_raw,mms3_mec_rdata_raw,...
        mms4_mec_timedata_raw,mms4_mec_rdata_raw,...
        data_points, data_type)
    %Andrew Vu 11/25/18
    %cd '~/Library/Mobile Documents/com~apple~CloudDocs/Research/Analysis'
    tic
    %Minimum angle
    min_angle = 10;
    
    %Number of Plots
    num_plots = 5;
    
    %Vertical Spacing between plots
    plot_gap=1.25;
    
    %Burst or Srvy
    if strcmp(data_type,'brst')
        time_resolution = 128;
    elseif strmp(data_type,'srvy')
        time_resolution = 16;
    end
    
    
    
    %data_points = 256; %total number of brst data points=2s, 1 second on each side,
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Data Interpolation%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Interpolate the Data
    [~,mms1_mec_rdata_interp] = interpxyz(mms1_mec_timedata_raw,mms1_mec_rdata_raw,mms1_fgm_timedata_raw);
    
    [~,mms2_fgm_bdata_interp] = interpxyz(mms2_fgm_timedata_raw,mms2_fgm_bdata_raw,mms1_fgm_timedata_raw);
    [~,mms2_mec_rdata_interp] = interpxyz(mms2_mec_timedata_raw,mms2_mec_rdata_raw,mms1_fgm_timedata_raw);
    
    [~,mms3_fgm_bdata_interp] = interpxyz(mms3_fgm_timedata_raw,mms3_fgm_bdata_raw,mms1_fgm_timedata_raw);
    [~,mms3_mec_rdata_interp] = interpxyz(mms3_mec_timedata_raw,mms3_mec_rdata_raw,mms1_fgm_timedata_raw);
    
    [~,mms4_fgm_bdata_interp] = interpxyz(mms4_fgm_timedata_raw,mms4_fgm_bdata_raw,mms1_fgm_timedata_raw);
    [~,mms4_mec_rdata_interp] = interpxyz(mms4_mec_timedata_raw,mms4_mec_rdata_raw,mms1_fgm_timedata_raw);
    
    %Crop mec data to specific time period
    [~,mms1_mec_rdata,~,~] = crop(mms1_fgm_timedata_raw,mms1_mec_rdata_interp,event_start,event_end);
    [~,mms2_mec_rdata,~,~] = crop(mms1_fgm_timedata_raw,mms2_mec_rdata_interp,event_start,event_end);
    [~,mms3_mec_rdata,~,~] = crop(mms1_fgm_timedata_raw,mms3_mec_rdata_interp,event_start,event_end);
    [~,mms4_mec_rdata,~,~] = crop(mms1_fgm_timedata_raw,mms4_mec_rdata_interp,event_start,event_end);
    
    %Crop fgm data to specific time period
    [mms1_fgm_timedata,mms1_fgm_bdata,~,~] = crop(mms1_fgm_timedata_raw,mms1_fgm_bdata_raw,event_start,event_end);
    [~,mms2_fgm_bdata,~,~] = crop(mms1_fgm_timedata_raw,mms2_fgm_bdata_interp,event_start,event_end);
    [~,mms3_fgm_bdata,~,~] = crop(mms1_fgm_timedata_raw,mms3_fgm_bdata_interp,event_start,event_end);
    [~,mms4_fgm_bdata,~,~] = crop(mms1_fgm_timedata_raw,mms4_fgm_bdata_interp,event_start,event_end);
    %mms_fgm_timedata = [mms1_fgm_timedata,mms2_fgm_timedata,mms3_fgm_timedata,mms4_fgm_timedata];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%
    %%%%%%%%%%%%%%%%%%%%Sliding Window for Timing Method%%%%%%%%%%%%%%%%%%%%
    %Initialize Variables for Loop
    total_data_points = length(mms1_fgm_timedata); %in all of bdata
    
    timing_Bx_n=zeros(total_data_points,3);
    timing_Bx_v=zeros(total_data_points,1);
    timing_Bx_cc=zeros(total_data_points,1);
    timing_Bx_unique=zeros(total_data_points,1);
    
    timing_By_n=zeros(total_data_points,3);
    timing_By_v=zeros(total_data_points,1);
    timing_By_cc=zeros(total_data_points,1);
    timing_By_unique=zeros(total_data_points,1);
    
    timing_Bz_n=zeros(total_data_points,3);
    timing_Bz_v=zeros(total_data_points,1);
    timing_Bz_cc=zeros(total_data_points,1);
    timing_Bz_unique=zeros(total_data_points,1);
    
    %t_0s, loop for t_0s from the start to the end of the event
    parfor t_0=data_points/2+1:total_data_points-data_points/2-1
        
        %Calculate for the n and v for each component at the t_0 with the
        %specified data_points range.
        [n_boundary,v_boundary,ccmean,~, ~, ~, ~, ~, ~, ~,unique_timings] = timing_method(t_0,data_points,...
            mms1_fgm_timedata,mms1_fgm_bdata,...
            mms2_fgm_bdata,...
            mms3_fgm_bdata,...
            mms4_fgm_bdata,...
            mms1_mec_rdata,...
            mms2_mec_rdata,...
            mms3_mec_rdata,...
            mms4_mec_rdata);
        %Store the normals and speeds in arrays.
        
        
        timing_Bx_n(t_0,:) = n_boundary(:,1);
        timing_Bx_v(t_0) = v_boundary(1);
        timing_Bx_cc(t_0) = ccmean(1);
        timing_Bx_unique(t_0) = unique_timings(1);
        
        timing_By_n(t_0,:) = n_boundary(:,2);
        timing_By_v(t_0) = v_boundary(2);
        timing_By_cc(t_0) = ccmean(2);
        timing_By_unique(t_0) = unique_timings(2);
        
        timing_Bz_n(t_0,:) = n_boundary(:,3);
        timing_Bz_v(t_0) = v_boundary(3);
        timing_Bz_cc(t_0) = ccmean(3);
        timing_Bz_unique(t_0) = unique_timings(3);
        
        
        
        
        %only save parameters if the 4SC timings are unique.
        %         if unique_timings(1) == 1
        %             timing_Bx_n(t_0,:) = n_boundary(:,1);
        %             timing_Bx_v(t_0) = v_boundary(1);
        %             timing_Bx_cc(t_0) = ccmean(1);
        %         else
        %             timing_Bx_n(t_0,:) = NaN;
        %             timing_Bx_v(t_0) = NaN;
        %             timing_Bx_cc(t_0) = NaN;
        %         end
        %
        %
        %         if unique_timings(2) == 1
        %             timing_By_n(t_0,:) = n_boundary(:,2);
        %             timing_By_v(t_0) = v_boundary(2);
        %             timing_By_cc(t_0) = ccmean(2);
        %         else
        %             timing_By_n(t_0,:) = NaN;
        %             timing_By_v(t_0) = NaN;
        %             timing_By_cc(t_0) = NaN;
        %         end
        %
        %
        %         if unique_timings(3) == 1
        %             timing_Bz_n(t_0,:) = n_boundary(:,3);
        %             timing_Bz_v(t_0) = v_boundary(3);
        %             timing_Bz_cc(t_0) = ccmean(3);
        %
        %         else
        %             timing_Bz_n(t_0,:) = NaN;
        %             timing_Bz_v(t_0) = NaN;
        %             timing_Bz_cc(t_0) = NaN;
        %         end
        %
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Calculate MVAB and angles
    mvab_eigenvalues = zeros(total_data_points,1);
    mvab_eigenvectors = zeros(total_data_points,3);
    timing_Bx_mva_angle = zeros(total_data_points,1);
    timing_By_mva_angle = zeros(total_data_points,1);
    timing_Bz_mva_angle = zeros(total_data_points,1);
    
    parfor i=1:total_data_points
        
        %Calculate the MVAB for the specified time interval
        if timing_Bx_n(i,:) ~= 0 %don't do calculation if it's zero, outside of event limits.
            
            [~,mvab_l,mvab_n] = mvab(mms1_fgm_bdata(i-data_points/2:i+data_points/2,1:3));
            
            % 1 if ratio is greater than 10.
            mvab_eigenvalues(i) = mvab_l(2)/mvab_l(3) > 10;
            
            %Save Eigenvectors
            mvab_eigenvectors(i,:) = mvab_n(:,3);
            
            %             timing_Bx_mva_angle(i) = mvab_eigenvalues(i)*(ambig_angle(mvab_eigenvectors(i,:),timing_Bx_n(i,:)) < 15);
            %             timing_By_mva_angle(i) = mvab_eigenvalues(i)*(ambig_angle(mvab_eigenvectors(i,:),timing_By_n(i,:)) < 15);
            %             timing_Bz_mva_angle(i) = mvab_eigenvalues(i)*(ambig_angle(mvab_eigenvectors(i,:),timing_Bz_n(i,:)) < 15);
            
            %Calculate angle between timing and MVAB
            timing_Bx_mva_angle(i) = ambig_angle(mvab_eigenvectors(i,:),timing_Bx_n(i,:));
            timing_By_mva_angle(i) = ambig_angle(mvab_eigenvectors(i,:),timing_By_n(i,:));
            timing_Bz_mva_angle(i) = ambig_angle(mvab_eigenvectors(i,:),timing_Bz_n(i,:));
            
            
        end
    end
    %Eigenvalues = 1 if ratio is greater than 10, otherwise NaN.
    mvab_eigenvalues(mvab_eigenvalues == 0) = NaN;
    
    %Minimum angle
    %If angle is greater than 10, equal to NaN;
    timing_Bx_mva_angle(timing_Bx_mva_angle > min_angle) = NaN;
    timing_By_mva_angle(timing_By_mva_angle > min_angle) = NaN;
    timing_Bz_mva_angle(timing_Bz_mva_angle > min_angle) = NaN;
    
    %If 4SC-timings are not unique, = Nan;
    timing_Bx_mva_angle(timing_Bx_unique == 0) = NaN;
    timing_By_mva_angle(timing_By_unique == 0) = NaN;
    timing_Bz_mva_angle(timing_Bz_unique == 0) = NaN;
    
    %Only Non-Nan if less than 10 degrees and eigenvalue ratio is
    %greater than 10.
    timing_Bx_mva_angle = timing_Bx_mva_angle .* mvab_eigenvalues;
    timing_By_mva_angle = timing_By_mva_angle .* mvab_eigenvalues;
    timing_Bz_mva_angle = timing_Bz_mva_angle .* mvab_eigenvalues;
    
    %Plot the number of <10 angles for a given time.
    timing_Ball_mva_angle = zeros(total_data_points,1);
    timing_Ball_mva_angle = ~isnan(timing_Bx_mva_angle) + ~isnan(timing_By_mva_angle) + ~isnan(timing_Bz_mva_angle);
    timing_Ball_mva_angle(timing_Ball_mva_angle == 0 ) = NaN;
    
    
    
    
    
    
    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%ALL Normal Plots%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure('Position',[1 1 1100 700])
    set(gcf,'color','w');
    co = [0 0 1;
        0 1 0;
        1 0 0;
        0 0 0];
    set(gcf,'defaultAxesColorOrder',co)
    plot_order = 1;
    %%%%%%%%%%%%%%%%%%%%%%%%%%Magnetic Field Data%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %plot overall event magnetic field for mms1
    plot_fgm_magnetic(event_start,event_end,mms1_fgm_timedata,mms1_fgm_bdata,num_plots,plot_order)
    
    %Title
    title_name = strcat('MMS Timing Method Sliding Window', ',\Deltat=', num2str(data_points/time_resolution),'s');
    title(title_name, 'FontSize', 18, 'FontWeight', 'normal')
    
    %Grab Position of Plot
    plot_pos = get(gca,'Position');
    %Set Plot Position
    set(gca,'Position',[plot_pos(1), plot_pos(2)-plot_pos(4)*plot_gap*(plot_order-1), plot_pos(3), plot_pos(4)]);
    plot_order = plot_order+1;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%Plot normal Nx for all three B components%%%%%%%%%%%%%%%%%%%%%%%%%%%
    hold on
    plot_parameter(event_start,event_end,mms1_fgm_timedata,[timing_Bx_n(:,1),timing_By_n(:,1),timing_Bz_n(:,1)],num_plots,plot_order,{'B_x','B_y','B_z'},'N_x')
    %     plot_parameter(event_start,event_end,mms1_fgm_timedata,,num_plots,plot_order,'','')
    %     plot_parameter(event_start,event_end,mms1_fgm_timedata,timing_Bx_n(:,3),num_plots,plot_order,,'N')
    hold off
    set(gca,'Position',[plot_pos(1), plot_pos(2)-plot_pos(4)*plot_gap*(plot_order-1), plot_pos(3), plot_pos(4)]);
    plot_order = plot_order+1;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%Plot normal Ny for all three B components%%%%%%%%%%%%%%%%%%%%%%%%%%%
    hold on
    plot_parameter(event_start,event_end,mms1_fgm_timedata,[timing_Bx_n(:,2),timing_By_n(:,2),timing_Bz_n(:,2)],num_plots,plot_order,{'B_x','B_y','B_z'},'N_y')
    %     plot_parameter(event_start,event_end,mms1_fgm_timedata,,num_plots,plot_order,'','')
    %     plot_parameter(event_start,event_end,mms1_fgm_timedata,timing_Bx_n(:,3),num_plots,plot_order,,'N')
    hold off
    set(gca,'Position',[plot_pos(1), plot_pos(2)-plot_pos(4)*plot_gap*(plot_order-1), plot_pos(3), plot_pos(4)]);
    plot_order = plot_order+1;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%Plot normal Nz for all three B components%%%%%%%%%%%%%%%%%%%%%%%%%%%
    hold on
    plot_parameter(event_start,event_end,mms1_fgm_timedata,[timing_Bx_n(:,3),timing_By_n(:,3),timing_Bz_n(:,3)],num_plots,plot_order,{'B_x','B_y','B_z'},'N_z')
    %     plot_parameter(event_start,event_end,mms1_fgm_timedata,,num_plots,plot_order,'','')
    %     plot_parameter(event_start,event_end,mms1_fgm_timedata,timing_Bx_n(:,3),num_plots,plot_order,,'N')
    hold off
    set(gca,'Position',[plot_pos(1), plot_pos(2)-plot_pos(4)*plot_gap*(plot_order-1), plot_pos(3), plot_pos(4)]);
    plot_order = plot_order+1;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%Plot speed for each Component%%%%%%%%%%%%%%%%%%%%%%%%%%%
    plot_parameter(event_start,event_end,mms1_fgm_timedata,[timing_Bx_v,timing_By_v,timing_Bz_v],num_plots,plot_order,{'B_x','B_y','B_z'},{'V';'[km/s]'})
    set(gca,'Position',[plot_pos(1), plot_pos(2)-plot_pos(4)*plot_gap*(plot_order-1), plot_pos(3), plot_pos(4)]);
    plot_order = plot_order+1;
    datetick('keeplimits')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    orient(gcf,'landscape')
    plot_name =  strcat('3_SlidingWindow_Timing_All',...
        event_start(1:19),'_',num2str(data_points),'.pdf');
    print(gcf, '-dpdf', '-opengl', plot_name,'-fillpage');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%Plotting for Bx%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure('Position',[1 1 1100 700])
    set(gcf,'color','w');
    co = [0 0 1;
        0 1 0;
        1 0 0;
        0 0 0];
    set(gcf,'defaultAxesColorOrder',co)
    plot_order = 1;
    %%%%%%%%%%%%%%%%%%%%%%%%%%Magnetic Field Data%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %plot overall event magnetic field for mms1
    plot_fgm_magnetic(event_start,event_end,mms1_fgm_timedata,mms1_fgm_bdata,num_plots,plot_order)
    
    %Title
    title_name = strcat('MMS Timing Method Sliding Window for B_x', ',\Deltat=', num2str(data_points/time_resolution),'s');
    title(title_name, 'FontSize', 18, 'FontWeight', 'normal')
    
    %Grab Position of Plot
    plot_pos = get(gca,'Position');
    %Set Plot Position
    set(gca,'Position',[plot_pos(1), plot_pos(2)-plot_pos(4)*plot_gap*(plot_order-1), plot_pos(3), plot_pos(4)]);
    plot_order = plot_order+1;
    %%%%%%%%%%%%%%%%%%%%%%%Plot Bx for each Spacecraft%%%%%%%%%%%%%%%%%%%%%%%%%%%
    plot_parameter(event_start,event_end,mms1_fgm_timedata,...
        [mms1_fgm_bdata(:,1),mms2_fgm_bdata(:,1),mms3_fgm_bdata(:,1),mms4_fgm_bdata(:,1)],...
        num_plots,plot_order,{'MMS1','MMS2','MMS3','MMS4'},'B_x')
    set(gca,'Position',[plot_pos(1), plot_pos(2)-plot_pos(4)*plot_gap*(plot_order-1), plot_pos(3), plot_pos(4)]);
    plot_order = plot_order+1;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%Plot CCmean for Bx%%%%%%%%%%%%%%%%%%%%%%%%%%%
    plot_parameter(event_start,event_end,mms1_fgm_timedata,timing_Bx_cc,num_plots,plot_order,'','cc_{mean}')
    set(gca,'Position',[plot_pos(1), plot_pos(2)-plot_pos(4)*plot_gap*(plot_order-1), plot_pos(3), plot_pos(4)]);
    plot_order = plot_order+1;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%Plot normal for Bx%%%%%%%%%%%%%%%%%%%%%%%%%%%
    hold on
    plot_parameter(event_start,event_end,mms1_fgm_timedata,[timing_Bx_n(:,1),timing_Bx_n(:,2),timing_Bx_n(:,3)],num_plots,plot_order,{'N_x','N_y','N_z'},'N')
    %     plot_parameter(event_start,event_end,mms1_fgm_timedata,,num_plots,plot_order,'','')
    %     plot_parameter(event_start,event_end,mms1_fgm_timedata,timing_Bx_n(:,3),num_plots,plot_order,,'N')
    hold off
    set(gca,'Position',[plot_pos(1), plot_pos(2)-plot_pos(4)*plot_gap*(plot_order-1), plot_pos(3), plot_pos(4)]);
    plot_order = plot_order+1;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%Plot speed for Bx%%%%%%%%%%%%%%%%%%%%%%%%%%%
    plot_parameter(event_start,event_end,mms1_fgm_timedata,timing_Bx_v(:),num_plots,plot_order,'',{'V';'[km/s]'})
    set(gca,'Position',[plot_pos(1), plot_pos(2)-plot_pos(4)*plot_gap*(plot_order-1), plot_pos(3), plot_pos(4)]);
    plot_order = plot_order+1;
    datetick('keeplimits')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    orient(gcf,'landscape')
    plot_name =  strcat('3_SlidingWindow_Timing_Bx',...
        event_start(1:19),'_',num2str(data_points),'.pdf');
    print(gcf, '-dpdf', '-opengl', plot_name,'-fillpage');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%Plotting for By%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure('Position',[1 1 1100 700])
    set(gcf,'color','w');
    co = [0 0 1;
        0 1 0;
        1 0 0;
        0 0 0];
    set(gcf,'defaultAxesColorOrder',co)
    plot_order=1;
    %%%%%%%%%%%%%%%%%%%%%%%%%%Magnetic Field Data%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %plot overall event magnetic field for mms1
    plot_fgm_magnetic(event_start,event_end,mms1_fgm_timedata,mms1_fgm_bdata,num_plots,plot_order)
    
    %Title
    title_name = strcat('MMS Timing Method Sliding Window for B_y', ',\Deltat=', num2str(data_points/time_resolution),'s');
    title(title_name, 'FontSize', 18, 'FontWeight', 'normal')
    
    %Grab Position of Plot
    plot_pos = get(gca,'Position');
    %Set Plot Position
    set(gca,'Position',[plot_pos(1), plot_pos(2)-plot_pos(4)*plot_gap*(plot_order-1), plot_pos(3), plot_pos(4)]);
    plot_order = plot_order+1;
    %%%%%%%%%%%%%%%%%%%%%%%Plot Bx for each Spacecraft%%%%%%%%%%%%%%%%%%%%%%%%%%%
    plot_parameter(event_start,event_end,mms1_fgm_timedata,...
        [mms1_fgm_bdata(:,2),mms2_fgm_bdata(:,2),mms3_fgm_bdata(:,2),mms4_fgm_bdata(:,2)],...
        num_plots,plot_order,{'MMS1','MMS2','MMS3','MMS4'},'B_y')
    set(gca,'Position',[plot_pos(1), plot_pos(2)-plot_pos(4)*plot_gap*(plot_order-1), plot_pos(3), plot_pos(4)]);
    plot_order = plot_order+1;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%Plot CCmean for By%%%%%%%%%%%%%%%%%%%%%%%%%%%
    plot_parameter(event_start,event_end,mms1_fgm_timedata,timing_By_cc,num_plots,plot_order,'','cc_{mean}')
    set(gca,'Position',[plot_pos(1), plot_pos(2)-plot_pos(4)*plot_gap*(plot_order-1), plot_pos(3), plot_pos(4)]);
    plot_order = plot_order+1;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%Plot normal for By%%%%%%%%%%%%%%%%%%%%%%%%%%%
    hold on
    plot_parameter(event_start,event_end,mms1_fgm_timedata,[timing_By_n(:,1),timing_By_n(:,2),timing_By_n(:,3)],num_plots,plot_order,{'N_x','N_y','N_z'},'N')
    %     plot_parameter(event_start,event_end,mms1_fgm_timedata,,num_plots,plot_order,'','')
    %     plot_parameter(event_start,event_end,mms1_fgm_timedata,timing_Bx_n(:,3),num_plots,plot_order,,'N')
    hold off
    set(gca,'Position',[plot_pos(1), plot_pos(2)-plot_pos(4)*plot_gap*(plot_order-1), plot_pos(3), plot_pos(4)]);
    plot_order = plot_order+1;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%Plot speed for By%%%%%%%%%%%%%%%%%%%%%%%%%%%
    plot_parameter(event_start,event_end,mms1_fgm_timedata,timing_By_v(:),num_plots,plot_order,'',{'V';'[km/s]'})
    set(gca,'Position',[plot_pos(1), plot_pos(2)-plot_pos(4)*plot_gap*(plot_order-1), plot_pos(3), plot_pos(4)]);
    plot_order = plot_order+1;
    datetick('keeplimits')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    orient(gcf,'landscape')
    plot_name =  strcat('3_SlidingWindow_Timing_By',...
        event_start(1:19),'_',num2str(data_points),'.pdf');
    print(gcf, '-dpdf', '-opengl', plot_name,'-fillpage');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%Plotting for Bz%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure('Position',[1 1 1100 700])
    set(gcf,'color','w');
    co = [0 0 1;
        0 1 0;
        1 0 0;
        0 0 0];
    set(gcf,'defaultAxesColorOrder',co)
    plot_order=1;
    %%%%%%%%%%%%%%%%%%%%%%%%%%Magnetic Field Data%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %plot overall event magnetic field for mms1
    plot_fgm_magnetic(event_start,event_end,mms1_fgm_timedata,mms1_fgm_bdata,num_plots,plot_order)
    
    %Title
    title_name = strcat('MMS Timing Method Sliding Window for B_z', ',\Deltat=', num2str(data_points/time_resolution),'s');
    title(title_name, 'FontSize', 18, 'FontWeight', 'normal')
    
    %Grab Position of Plot
    plot_pos = get(gca,'Position');
    %Set Plot Position
    set(gca,'Position',[plot_pos(1), plot_pos(2)-plot_pos(4)*plot_gap*(plot_order-1), plot_pos(3), plot_pos(4)]);
    plot_order = plot_order+1;
    %%%%%%%%%%%%%%%%%%%%%%%Plot Bz for each Spacecraft%%%%%%%%%%%%%%%%%%%%%%%%%%%
    plot_parameter(event_start,event_end,mms1_fgm_timedata,...
        [mms1_fgm_bdata(:,3),mms2_fgm_bdata(:,3),mms3_fgm_bdata(:,3),mms4_fgm_bdata(:,3)],...
        num_plots,plot_order,{'MMS1','MMS2','MMS3','MMS4'},'B_z')
    set(gca,'Position',[plot_pos(1), plot_pos(2)-plot_pos(4)*plot_gap*(plot_order-1), plot_pos(3), plot_pos(4)]);
    plot_order = plot_order+1;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%Plot CCmean for Bz%%%%%%%%%%%%%%%%%%%%%%%%%%%
    plot_parameter(event_start,event_end,mms1_fgm_timedata,timing_Bz_cc,num_plots,plot_order,'','cc_{mean}')
    set(gca,'Position',[plot_pos(1), plot_pos(2)-plot_pos(4)*plot_gap*(plot_order-1), plot_pos(3), plot_pos(4)]);
    plot_order = plot_order+1;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%Plot normal for Bz%%%%%%%%%%%%%%%%%%%%%%%%%%%
    hold on
    plot_parameter(event_start,event_end,mms1_fgm_timedata,[timing_Bz_n(:,1),timing_Bz_n(:,2),timing_Bz_n(:,3)],num_plots,plot_order,{'N_x','N_y','N_z'},'N')
    %     plot_parameter(event_start,event_end,mms1_fgm_timedata,,num_plots,plot_order,'','')
    %     plot_parameter(event_start,event_end,mms1_fgm_timedata,timing_Bx_n(:,3),num_plots,plot_order,,'N')
    hold off
    set(gca,'Position',[plot_pos(1), plot_pos(2)-plot_pos(4)*plot_gap*(plot_order-1), plot_pos(3), plot_pos(4)]);
    plot_order = plot_order+1;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%Plot speed for Bz%%%%%%%%%%%%%%%%%%%%%%%%%%%
    plot_parameter(event_start,event_end,mms1_fgm_timedata,timing_Bz_v(:),num_plots,plot_order,'',{'V';'[km/s]'})
    set(gca,'Position',[plot_pos(1), plot_pos(2)-plot_pos(4)*plot_gap*(plot_order-1), plot_pos(3), plot_pos(4)]);
    plot_order = plot_order+1;
    datetick('keeplimits')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    orient(gcf,'landscape')
    plot_name =  strcat('3_SlidingWindow_Timing_Bz',...
        event_start(1:19),'_',num2str(data_points),'.pdf');
    print(gcf, '-dpdf', '-opengl', plot_name,'-fillpage');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%Plotting for Angles%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure('Position',[1 1 1100 700])
    set(gcf,'color','w');
    co = [0 0 1;
        0 1 0;
        1 0 0;
        0 0 0];
    set(gcf,'defaultAxesColorOrder',co)
    plot_order=1;
    %%%%%%%%%%%%%%%%%%%%%%%%%%Magnetic Field Data%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %plot overall event magnetic field for mms1
    plot_fgm_magnetic(event_start,event_end,mms1_fgm_timedata,mms1_fgm_bdata,num_plots,plot_order)
    datetick('keeplimits')
    %Title
    title_name = strcat('MMS Timing Method & MVAB Angles', ',\Deltat=', num2str(data_points/time_resolution),'s');
    title(title_name, 'FontSize', 18, 'FontWeight', 'normal')
    
    %Grab Position of Plot
    plot_pos = get(gca,'Position');
    %Set Plot Position
    set(gca,'Position',[plot_pos(1), plot_pos(2)-plot_pos(4)*plot_gap*(plot_order-1), plot_pos(3), plot_pos(4)]);
    plot_order = plot_order+1;
    %%%%%%%%%%%%%%%%%%%%%%%Bx angles%%%%%%%%%%%%%%%%%%%%%%%%%%%
    plot_scatter(event_start,event_end,mms1_fgm_timedata,timing_Bx_mva_angle,num_plots,plot_order,'B_x','Angle')
    ylim([0 10])
    set(gca,'Position',[plot_pos(1), plot_pos(2)-plot_pos(4)*plot_gap*(plot_order-1), plot_pos(3), plot_pos(4)]);
    plot_order = plot_order+1;
    datetick('keeplimits')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%By Angles%%%%%%%%%%%%%%%%%%%%%%%%%%%
    plot_scatter(event_start,event_end,mms1_fgm_timedata,timing_By_mva_angle,num_plots,plot_order,'B_y','Angle')
    ylim([0 10])
    set(gca,'Position',[plot_pos(1), plot_pos(2)-plot_pos(4)*plot_gap*(plot_order-1), plot_pos(3), plot_pos(4)]);
    plot_order = plot_order+1;
    datetick('keeplimits')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%Bz angles%%%%%%%%%%%%%%%%%%%%%%%%%%%
    plot_scatter(event_start,event_end,mms1_fgm_timedata,timing_Bz_mva_angle,num_plots,plot_order,'B_z','Angle')
    ylim([0 10])
    set(gca,'Position',[plot_pos(1), plot_pos(2)-plot_pos(4)*plot_gap*(plot_order-1), plot_pos(3), plot_pos(4)]);
    plot_order = plot_order+1;
    datetick('keeplimits')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%All angles less than 10%%%%%%%%%%%%%%%%%%%%%%%%%%%
    plot_scatter(event_start,event_end,mms1_fgm_timedata,timing_Ball_mva_angle,num_plots,plot_order,'','#')
    ylim([0 3])
    set(gca,'Position',[plot_pos(1), plot_pos(2)-plot_pos(4)*plot_gap*(plot_order-1), plot_pos(3), plot_pos(4)]);
    plot_order = plot_order+1;
    datetick('keeplimits')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    orient(gcf,'landscape')
    plot_name =  strcat('3_SlidingWindow_Timing_MVA',...
        event_start(1:19),'_',num2str(data_points),'.pdf');
    print(gcf, '-dpdf', '-opengl', plot_name,'-fillpage');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    formatOut = "yyyy-mm-dd HH:MM:SS.FFF";
    start_datenum = datenum(event_start);
    
    
    
    
    %Display time ranges of angle less than 10
    timing_mva_lessthan10_indices = find(~isnan(timing_Ball_mva_angle));
    time_ranges = cell(length(timing_mva_lessthan10_indices),8,1);
    for i=1:length(timing_mva_lessthan10_indices)
        time_left_index = timing_mva_lessthan10_indices(i) - data_points/2;
        time_right_index = timing_mva_lessthan10_indices(i) + data_points/2;
        time_left_string = datestr(mms1_fgm_timedata(time_left_index),formatOut);
        time_right_string = datestr(mms1_fgm_timedata(time_right_index),formatOut);
        time_ranges{i,1} = {time_left_string};
        time_ranges{i,2} = {time_right_string};
        %Angles
        time_ranges{i,3} = timing_Bx_mva_angle(timing_mva_lessthan10_indices(i));
        time_ranges{i,4} = timing_By_mva_angle(timing_mva_lessthan10_indices(i));
        time_ranges{i,5} = timing_Bz_mva_angle(timing_mva_lessthan10_indices(i));
        %Speeds
        time_ranges{i,6} = timing_Bx_v(timing_mva_lessthan10_indices(i));
        time_ranges{i,7} = timing_By_v(timing_mva_lessthan10_indices(i));
        time_ranges{i,8} = timing_Bz_v(timing_mva_lessthan10_indices(i));
        
        
        
        %For debugging.
                time_ranges{i,12} = mvab_eigenvectors(timing_mva_lessthan10_indices(i),1);
                time_ranges{i,13} = mvab_eigenvectors(timing_mva_lessthan10_indices(i),2);
                time_ranges{i,14} = mvab_eigenvectors(timing_mva_lessthan10_indices(i),3);
        
        
        %Print the normal vector for the smallest angle direction with MVA. Bx By Bz
        B_angles = [timing_Bx_mva_angle(timing_mva_lessthan10_indices(i)), timing_By_mva_angle(timing_mva_lessthan10_indices(i)), ...
            timing_Bz_mva_angle(timing_mva_lessthan10_indices(i))];
        [~,min_B_angles_index]=min(B_angles,[],'omitnan');
        
        switch min_B_angles_index
            case 1
                time_ranges{i,9} = timing_Bx_n(timing_mva_lessthan10_indices(i),1);
                time_ranges{i,10} = timing_Bx_n(timing_mva_lessthan10_indices(i),2);
                time_ranges{i,11} = timing_Bx_n(timing_mva_lessthan10_indices(i),3);
            case 2
                time_ranges{i,9} = timing_By_n(timing_mva_lessthan10_indices(i),1);
                time_ranges{i,10} = timing_By_n(timing_mva_lessthan10_indices(i),2);
                time_ranges{i,11} = timing_By_n(timing_mva_lessthan10_indices(i),3);
            case 3
                time_ranges{i,9} = timing_Bz_n(timing_mva_lessthan10_indices(i),1);
                time_ranges{i,10} = timing_Bz_n(timing_mva_lessthan10_indices(i),2);
                time_ranges{i,11} = timing_Bz_n(timing_mva_lessthan10_indices(i),3);
        end

    end
    %Save to a .dat or .txt file for future reference
    TT = cell2table(time_ranges);
    writetable(TT,strcat('time_ranges',num2str(data_points),'.dat'));
    
    
    %     %Display time ranges for 1 occurence of angle less than 10. out of 3.
    %     timing_mva_lessthan10 = find(timing_Ball_mva_angle == 1);
    %     time_ranges_1 = cell(length(timing_mva_lessthan10),3,1);
    %     for i=1:length(timing_mva_lessthan10)
    %         time_left_index = timing_mva_lessthan10(i)-data_points/2;
    %         time_right_index = timing_mva_lessthan10(i) + data_points/2;
    %         time_left_string = datestr(start_datenum+mms1_fgm_bdata(time_left_index)/86400,formatOut);
    %         time_right_string = datestr(start_datenum+mms1_fgm_bdata(time_right_index)/86400,formatOut);
    %         time_ranges_1{i,1} = {time_left_string};
    %         time_ranges_1{i,2} = {time_right_string};
    %         %time_ranges_1{i,3} = {timing_Ball_mva_angle(timing_mva_lessthan10(i))};
    %     end
    %
    %
    %
    %     %Display time ranges for 2 occurences of angle less than 10. out of 3.
    %     timing_mva_lessthan10 = find(timing_Ball_mva_angle == 2);
    %     time_ranges_2 = cell(length(timing_mva_lessthan10),2,1);
    %     for i=1:length(timing_mva_lessthan10)
    %         time_left_index = timing_mva_lessthan10(i)-data_points/2;
    %         time_right_index = timing_mva_lessthan10(i) + data_points/2;
    %         time_left_string = datestr(start_datenum+mms1_fgm_bdata(time_left_index)/86400,formatOut);
    %         time_right_string = datestr(start_datenum+mms1_fgm_bdata(time_right_index)/86400,formatOut);
    %         time_ranges_2{i,1} = {time_left_string};
    %         time_ranges_2{i,2} = {time_right_string};
    %     end
    %
    %     %Display time ranges for 3 occurences of angle less than 10. out of 3.
    %     timing_mva_lessthan10 = find(timing_Ball_mva_angle == 3);
    %     time_ranges_3 = cell(length(timing_mva_lessthan10),2,1);
    %     for i=1:length(timing_mva_lessthan10)
    %         time_left_index = timing_mva_lessthan10(i)-data_points/2;
    %         time_right_index = timing_mva_lessthan10(i) + data_points/2;
    %         time_left_string = datestr(start_datenum+mms1_fgm_bdata(time_left_index)/86400,formatOut);
    %         time_right_string = datestr(start_datenum+mms1_fgm_bdata(time_right_index)/86400,formatOut);
    %         time_ranges_3{i,1} = {time_left_string};
    %         time_ranges_3{i,2} = {time_right_string};
    %     end
    % formatIn2 = 'yyyy-MM-dd HH:mm:ss.SSS'
    % datetime(event_start,'format',formatIn2)
    % AA = datetime(mms1_fgm_timedata_srvy,'ConvertFrom','datenum')
    toc
end

