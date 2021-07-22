function [] = plot_master(event_analysis,threshold_std,timing_window,sliding_minratio,data_type,...
        Event_number,Event_Type,Substructure,...
        event_start, event_end,...
        left_OuterEdge, left_InnerEdge,...
        right_InnerEdge, right_OuterEdge,...
        leading_leftmost_date, leading_rightmost_date,...
        trailing_leftmost_date,trailing_rightmost_date,...
        leading_start,leading_end,...
        trailing_start, trailing_end,...
        downstream_geometry,upstream_geometry,...
        n)
    %MMS Master Plot
    %Plots summary and boundary analysis entirely on an event
    probe = 1;
    ssSTDs = 3.0000;
    minDataPointsforSS = 1; %normally it is 4.
    ShearCutoff = 15;
    %% Load SW DAta from IDL results
    cd '~/Library/Mobile Documents/com~apple~CloudDocs/Research/Analysis'
    swspeedEvents = importdata('swspeedEvents');
    swVelocityEvents = importdata('swVelocityEvents');
    
    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%Directory Management%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if event_analysis ~= 2
        cd '~/Library/Mobile Documents/com~apple~CloudDocs/Research/Analysis'
        eventFolderName = strcat(event_start(1:10),'/',event_start(12:13),'-',event_start(15:16));
        if ~exist(eventFolderName,'dir')
            mkdir(eventFolderName)
        end
        cd(eventFolderName)
        
    elseif event_analysis == 2
        cd '~/Library/Mobile Documents/com~apple~CloudDocs/Research/Analyzed Events'
        eventFolderName = strcat('Event Number_2021_3/',num2str(Event_number));
        if exist(eventFolderName,'dir')
            rmdir(eventFolderName,'s')
            mkdir(eventFolderName)
        else
            mkdir(eventFolderName)
        end
        cd(eventFolderName)
    end
    %% %Data Retrieval
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Data Retrieval%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if event_analysis == -2 %Just for testing events, don't save data yet.
        %Load MEC Data
        [mms1_mec_timedata_raw, mms1_mec_rdata_raw] = load_mec(event_start,event_end,1,data_type);
        [mms2_mec_timedata_raw, mms2_mec_rdata_raw] = load_mec(event_start,event_end,2,data_type);
        [mms3_mec_timedata_raw, mms3_mec_rdata_raw] = load_mec(event_start,event_end,3,data_type);
        [mms4_mec_timedata_raw, mms4_mec_rdata_raw] = load_mec(event_start,event_end,4,data_type);
        
        %load FGM data
        [mms1_fgm_timedata_raw, mms1_fgm_bdata_raw, ~, ~] = load_fgm(event_start,event_end,1,'brst');
        [mms2_fgm_timedata_raw, mms2_fgm_bdata_raw, ~, ~] = load_fgm(event_start,event_end,2,'brst');
        [mms3_fgm_timedata_raw, mms3_fgm_bdata_raw, ~, ~] = load_fgm(event_start,event_end,3,'brst');
        [mms4_fgm_timedata_raw, mms4_fgm_bdata_raw, ~, ~] = load_fgm(event_start,event_end,4,'brst');
        
        [mms1_fgm_timedata_srvy, mms1_fgm_bdata_srvy] = load_fgm(event_start,event_end,probe,'srvy'); %For Sliding Window
        
        %Load FPI_e
        [fpi_e_timedata,fpi_e_ndata,fpi_e_vdata,fpi_e_tparadata,fpi_e_tperpdata,...
            fpi_e_edata,fpi_e_espectdata,fpi_e_pressdata] = load_fpi(event_start,event_end,probe,'brst','e');
        %Load FPI_i
        [fpi_i_timedata,fpi_i_ndata,fpi_i_vdata,fpi_i_tparadata,fpi_i_tperpdata,...
            fpi_i_edata,fpi_i_espectdata,fpi_i_pressdata] = load_fpi(event_start,event_end,probe,'brst','i');
        
        %Load FPI Dist
        [time,phi_vector,theta_vector,energy_vector,dist] = load_dist(event_start,event_end,probe,'brst','i');
        
        
    else %Once we are past the stage of checking events, any further analysis will save the data
        
        eventDataFileName = strcat('MMS1_Data_EventNumber_',num2str(Event_number),'.mat');
        %Fucking icloud keeps deleting files every 30 minutes of downloading them
        %eventDataDirectory = '/Users/andrewvu/Library/Mobile Documents/com~apple~CloudDocs/Research/Event Data/';
        eventDataDirectory ='/Users/andrewvu/Library/Mobile Documents/com~apple~CloudDocs/Research/Event Data/';
        
        eventDataDirectoryFileName = strcat(eventDataDirectory,eventDataFileName);
        
        if exist(eventDataDirectoryFileName,'file') ~= 2
            
            %Load MEC Data
            [mms1_mec_timedata_raw, mms1_mec_rdata_raw] = load_mec(event_start,event_end,1,data_type);
            [mms2_mec_timedata_raw, mms2_mec_rdata_raw] = load_mec(event_start,event_end,2,data_type);
            [mms3_mec_timedata_raw, mms3_mec_rdata_raw] = load_mec(event_start,event_end,3,data_type);
            [mms4_mec_timedata_raw, mms4_mec_rdata_raw] = load_mec(event_start,event_end,4,data_type);
            
            %load FGM data
            [mms1_fgm_timedata_raw, mms1_fgm_bdata_raw, ~, ~] = load_fgm(event_start,event_end,1,'brst');
            [mms2_fgm_timedata_raw, mms2_fgm_bdata_raw, ~, ~] = load_fgm(event_start,event_end,2,'brst');
            [mms3_fgm_timedata_raw, mms3_fgm_bdata_raw, ~, ~] = load_fgm(event_start,event_end,3,'brst');
            [mms4_fgm_timedata_raw, mms4_fgm_bdata_raw, ~, ~] = load_fgm(event_start,event_end,4,'brst');
            
            [mms1_fgm_timedata_srvy, mms1_fgm_bdata_srvy] = load_fgm(event_start,event_end,probe,'srvy'); %For Sliding Window
            
            %Load FPI_e
            [fpi_e_timedata,fpi_e_ndata,fpi_e_vdata,fpi_e_tparadata,fpi_e_tperpdata,...
                fpi_e_edata,fpi_e_espectdata,fpi_e_pressdata] = load_fpi(event_start,event_end,probe,'brst','e');
            %Load FPI_i
            [fpi_i_timedata,fpi_i_ndata,fpi_i_vdata,fpi_i_tparadata,fpi_i_tperpdata,...
                fpi_i_edata,fpi_i_espectdata,fpi_i_pressdata] = load_fpi(event_start,event_end,probe,'brst','i');
            
            %Load FPI Dist
            [time,phi_vector,theta_vector,energy_vector,dist] = load_dist(event_start,event_end,probe,'brst','i');
            
            
            
            save(eventDataDirectoryFileName,...
                'mms1_mec_timedata_raw',...
                'mms1_mec_rdata_raw',...
                'mms2_mec_timedata_raw',...
                'mms2_mec_rdata_raw',...
                'mms3_mec_timedata_raw',...
                'mms3_mec_rdata_raw',...
                'mms4_mec_timedata_raw',...
                'mms4_mec_rdata_raw',...
                'mms1_fgm_timedata_raw',...
                'mms1_fgm_bdata_raw',...
                'mms2_fgm_timedata_raw',...
                'mms2_fgm_bdata_raw',...
                'mms3_fgm_timedata_raw',...
                'mms3_fgm_bdata_raw',...
                'mms4_fgm_timedata_raw',...
                'mms4_fgm_bdata_raw',...
                'mms1_fgm_timedata_srvy',...
                'mms1_fgm_bdata_srvy',...
                'fpi_e_timedata',...
                'fpi_e_ndata',...
                'fpi_e_vdata',...
                'fpi_e_tparadata',...
                'fpi_e_tperpdata',...
                'fpi_e_edata',...
                'fpi_e_espectdata',...
                'fpi_e_pressdata',...
                'fpi_i_timedata',...
                'fpi_i_ndata',...
                'fpi_i_vdata',...
                'fpi_i_tparadata',...
                'fpi_i_tperpdata',...
                'fpi_i_edata',...
                'fpi_i_espectdata',...
                'fpi_i_pressdata',...
                'time',...
                'phi_vector',...
                'theta_vector',...
                'energy_vector',...
                'dist')
        else
            load(eventDataDirectoryFileName) %#ok<LOAD>
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% %Plots
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Plot Plots%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if threshold_std == 0
        %Current Sheet Normal Calculation Plot, manually.
        [n_cs,B_pre_cs,B_post_cs] = manualCurrentSheet(event_start,event_end,leading_leftmost_date,leading_rightmost_date,...
            trailing_leftmost_date,trailing_rightmost_date,mms1_fgm_timedata_srvy,mms1_fgm_bdata_srvy);
    else
        %Calculate Current Sheet
        [n_cs,B_pre_cs,B_post_cs] = calculateCurrentSheet(event_start,event_end,mms1_fgm_timedata_srvy,mms1_fgm_bdata_srvy,threshold_std);
    end
    
    %Summary Plot
    [bowshock_n,n_cs,r_sc] = plot_mms_observation_summary_mod(event_start,event_end,...
        mms1_mec_timedata_raw,mms1_mec_rdata_raw,mms1_fgm_timedata_raw, mms1_fgm_bdata_raw,...
        fpi_e_timedata,fpi_e_ndata,fpi_e_vdata,fpi_e_tparadata,fpi_e_tperpdata,fpi_e_edata,fpi_e_espectdata,...
        fpi_i_timedata,fpi_i_ndata,fpi_i_vdata,fpi_i_tparadata,fpi_i_tperpdata,fpi_i_edata,fpi_i_espectdata,...
        n_cs)
    
    if event_analysis == -1
        dcm_obj = datacursormode(gcf);
        set(dcm_obj,'UpdateFcn',@getPlotTime)
    end
    
    
    %Either 0, do sliding windows, or 1, do MVA,MDD,Timing, not both., 1.5 is
    %for MVA sliding window only.
    if event_analysis == 0 || event_analysis == 0.5 || event_analysis == 1.5
        %Sliding Window MVA
        plot_MVABslidingwindow(event_start,event_end,mms1_fgm_timedata_srvy, mms1_fgm_bdata_srvy,n_cs,sliding_minratio)
        
        %         figure
        %         subplot(1,4,1)
        %         plot_MVABslidingwindow(event_start,event_end,mms1_fgm_timedata_raw, mms1_fgm_bdata_raw,currentsheet_n,sliding_minratio)
        %         title('MMS1 MVA')
        %         subplot(1,4,2)
        %         plot_MVABslidingwindow(event_start,event_end,mms2_fgm_timedata_raw, mms2_fgm_bdata_raw,currentsheet_n,sliding_minratio)
        %         title('MMS2 MVA')
        %         subplot(1,4,3)
        %         plot_MVABslidingwindow(event_start,event_end,mms3_fgm_timedata_raw, mms3_fgm_bdata_raw,currentsheet_n,sliding_minratio)
        %         title('MMS3 MVA')
        %         subplot(1,4,3)
        %         plot_MVABslidingwindow(event_start,event_end,mms4_fgm_timedata_raw, mms4_fgm_bdata_raw,currentsheet_n,sliding_minratio)
        %         title('MMS4 MVA')
        
        if event_analysis == 1.5
        else
            
            %Sliding Window Timing
            %If timing_window == 0 do 0.5,1,1.5,2,3,4 s intervals.
            %Otherwise do specified interval
            if timing_window == 0
                
                %                 [time_ranges] =  plot_Timingslidingwindow(event_start,event_end,...
                %                     mms1_fgm_timedata_raw,mms1_fgm_bdata_raw,...
                %                     mms2_fgm_timedata_raw,mms2_fgm_bdata_raw,...
                %                     mms3_fgm_timedata_raw,mms3_fgm_bdata_raw,...
                %                     mms4_fgm_timedata_raw,mms4_fgm_bdata_raw,...
                %                     mms1_mec_timedata_raw,mms1_mec_rdata_raw,...
                %                     mms2_mec_timedata_raw,mms2_mec_rdata_raw,...
                %                     mms3_mec_timedata_raw,mms3_mec_rdata_raw,...
                %                     mms4_mec_timedata_raw,mms4_mec_rdata_raw,...
                %                     0.125*256,'brst');
                
                [time_ranges] =  plot_Timingslidingwindow(event_start,event_end,...
                    mms1_fgm_timedata_raw,mms1_fgm_bdata_raw,...
                    mms2_fgm_timedata_raw,mms2_fgm_bdata_raw,...
                    mms3_fgm_timedata_raw,mms3_fgm_bdata_raw,...
                    mms4_fgm_timedata_raw,mms4_fgm_bdata_raw,...
                    mms1_mec_timedata_raw,mms1_mec_rdata_raw,...
                    mms2_mec_timedata_raw,mms2_mec_rdata_raw,...
                    mms3_mec_timedata_raw,mms3_mec_rdata_raw,...
                    mms4_mec_timedata_raw,mms4_mec_rdata_raw,...
                    0.25*256,'brst');
                
                
                [time_ranges] =  plot_Timingslidingwindow(event_start,event_end,...
                    mms1_fgm_timedata_raw,mms1_fgm_bdata_raw,...
                    mms2_fgm_timedata_raw,mms2_fgm_bdata_raw,...
                    mms3_fgm_timedata_raw,mms3_fgm_bdata_raw,...
                    mms4_fgm_timedata_raw,mms4_fgm_bdata_raw,...
                    mms1_mec_timedata_raw,mms1_mec_rdata_raw,...
                    mms2_mec_timedata_raw,mms2_mec_rdata_raw,...
                    mms3_mec_timedata_raw,mms3_mec_rdata_raw,...
                    mms4_mec_timedata_raw,mms4_mec_rdata_raw,...
                    0.375*256,'brst');
                
                %                 [time_ranges] =  plot_Timingslidingwindow(event_start,event_end,...
                %                     mms1_fgm_timedata_raw,mms1_fgm_bdata_raw,...
                %                     mms2_fgm_timedata_raw,mms2_fgm_bdata_raw,...
                %                     mms3_fgm_timedata_raw,mms3_fgm_bdata_raw,...
                %                     mms4_fgm_timedata_raw,mms4_fgm_bdata_raw,...
                %                     mms1_mec_timedata_raw,mms1_mec_rdata_raw,...
                %                     mms2_mec_timedata_raw,mms2_mec_rdata_raw,...
                %                     mms3_mec_timedata_raw,mms3_mec_rdata_raw,...
                %                     mms4_mec_timedata_raw,mms4_mec_rdata_raw,...
                %                     104,'brst');
                
                [time_ranges] =  plot_Timingslidingwindow(event_start,event_end,...
                    mms1_fgm_timedata_raw,mms1_fgm_bdata_raw,...
                    mms2_fgm_timedata_raw,mms2_fgm_bdata_raw,...
                    mms3_fgm_timedata_raw,mms3_fgm_bdata_raw,...
                    mms4_fgm_timedata_raw,mms4_fgm_bdata_raw,...
                    mms1_mec_timedata_raw,mms1_mec_rdata_raw,...
                    mms2_mec_timedata_raw,mms2_mec_rdata_raw,...
                    mms3_mec_timedata_raw,mms3_mec_rdata_raw,...
                    mms4_mec_timedata_raw,mms4_mec_rdata_raw,...
                    0.5*256,'brst');
                
                [time_ranges] =  plot_Timingslidingwindow(event_start,event_end,...
                    mms1_fgm_timedata_raw,mms1_fgm_bdata_raw,...
                    mms2_fgm_timedata_raw,mms2_fgm_bdata_raw,...
                    mms3_fgm_timedata_raw,mms3_fgm_bdata_raw,...
                    mms4_fgm_timedata_raw,mms4_fgm_bdata_raw,...
                    mms1_mec_timedata_raw,mms1_mec_rdata_raw,...
                    mms2_mec_timedata_raw,mms2_mec_rdata_raw,...
                    mms3_mec_timedata_raw,mms3_mec_rdata_raw,...
                    mms4_mec_timedata_raw,mms4_mec_rdata_raw,...
                    0.75*256,'brst');
                
                [time_ranges] =  plot_Timingslidingwindow(event_start,event_end,...
                    mms1_fgm_timedata_raw,mms1_fgm_bdata_raw,...
                    mms2_fgm_timedata_raw,mms2_fgm_bdata_raw,...
                    mms3_fgm_timedata_raw,mms3_fgm_bdata_raw,...
                    mms4_fgm_timedata_raw,mms4_fgm_bdata_raw,...
                    mms1_mec_timedata_raw,mms1_mec_rdata_raw,...
                    mms2_mec_timedata_raw,mms2_mec_rdata_raw,...
                    mms3_mec_timedata_raw,mms3_mec_rdata_raw,...
                    mms4_mec_timedata_raw,mms4_mec_rdata_raw,...
                    1*256,'brst');
                
                [time_ranges] =  plot_Timingslidingwindow(event_start,event_end,...
                    mms1_fgm_timedata_raw,mms1_fgm_bdata_raw,...
                    mms2_fgm_timedata_raw,mms2_fgm_bdata_raw,...
                    mms3_fgm_timedata_raw,mms3_fgm_bdata_raw,...
                    mms4_fgm_timedata_raw,mms4_fgm_bdata_raw,...
                    mms1_mec_timedata_raw,mms1_mec_rdata_raw,...
                    mms2_mec_timedata_raw,mms2_mec_rdata_raw,...
                    mms3_mec_timedata_raw,mms3_mec_rdata_raw,...
                    mms4_mec_timedata_raw,mms4_mec_rdata_raw,...
                    1.25*256,'brst');
                
                [time_ranges] =  plot_Timingslidingwindow(event_start,event_end,...
                    mms1_fgm_timedata_raw,mms1_fgm_bdata_raw,...
                    mms2_fgm_timedata_raw,mms2_fgm_bdata_raw,...
                    mms3_fgm_timedata_raw,mms3_fgm_bdata_raw,...
                    mms4_fgm_timedata_raw,mms4_fgm_bdata_raw,...
                    mms1_mec_timedata_raw,mms1_mec_rdata_raw,...
                    mms2_mec_timedata_raw,mms2_mec_rdata_raw,...
                    mms3_mec_timedata_raw,mms3_mec_rdata_raw,...
                    mms4_mec_timedata_raw,mms4_mec_rdata_raw,...
                    1.5*256,'brst');
                
                [time_ranges] =  plot_Timingslidingwindow(event_start,event_end,...
                    mms1_fgm_timedata_raw,mms1_fgm_bdata_raw,...
                    mms2_fgm_timedata_raw,mms2_fgm_bdata_raw,...
                    mms3_fgm_timedata_raw,mms3_fgm_bdata_raw,...
                    mms4_fgm_timedata_raw,mms4_fgm_bdata_raw,...
                    mms1_mec_timedata_raw,mms1_mec_rdata_raw,...
                    mms2_mec_timedata_raw,mms2_mec_rdata_raw,...
                    mms3_mec_timedata_raw,mms3_mec_rdata_raw,...
                    mms4_mec_timedata_raw,mms4_mec_rdata_raw,...
                    2*256,'brst');
                
                [time_ranges] =  plot_Timingslidingwindow(event_start,event_end,...
                    mms1_fgm_timedata_raw,mms1_fgm_bdata_raw,...
                    mms2_fgm_timedata_raw,mms2_fgm_bdata_raw,...
                    mms3_fgm_timedata_raw,mms3_fgm_bdata_raw,...
                    mms4_fgm_timedata_raw,mms4_fgm_bdata_raw,...
                    mms1_mec_timedata_raw,mms1_mec_rdata_raw,...
                    mms2_mec_timedata_raw,mms2_mec_rdata_raw,...
                    mms3_mec_timedata_raw,mms3_mec_rdata_raw,...
                    mms4_mec_timedata_raw,mms4_mec_rdata_raw,...
                    2.5*256,'brst');
                
                [time_ranges] =  plot_Timingslidingwindow(event_start,event_end,...
                    mms1_fgm_timedata_raw,mms1_fgm_bdata_raw,...
                    mms2_fgm_timedata_raw,mms2_fgm_bdata_raw,...
                    mms3_fgm_timedata_raw,mms3_fgm_bdata_raw,...
                    mms4_fgm_timedata_raw,mms4_fgm_bdata_raw,...
                    mms1_mec_timedata_raw,mms1_mec_rdata_raw,...
                    mms2_mec_timedata_raw,mms2_mec_rdata_raw,...
                    mms3_mec_timedata_raw,mms3_mec_rdata_raw,...
                    mms4_mec_timedata_raw,mms4_mec_rdata_raw,...
                    3*256,'brst');
                
                [time_ranges] =  plot_Timingslidingwindow(event_start,event_end,...
                    mms1_fgm_timedata_raw,mms1_fgm_bdata_raw,...
                    mms2_fgm_timedata_raw,mms2_fgm_bdata_raw,...
                    mms3_fgm_timedata_raw,mms3_fgm_bdata_raw,...
                    mms4_fgm_timedata_raw,mms4_fgm_bdata_raw,...
                    mms1_mec_timedata_raw,mms1_mec_rdata_raw,...
                    mms2_mec_timedata_raw,mms2_mec_rdata_raw,...
                    mms3_mec_timedata_raw,mms3_mec_rdata_raw,...
                    mms4_mec_timedata_raw,mms4_mec_rdata_raw,...
                    4*256,'brst');
                
                
                %
                %             else
                %
                %                 [time_ranges] =  plot_Timingslidingwindow(event_start,event_end,...
                %                     mms1_fgm_timedata_raw,mms1_fgm_bdata_raw,...
                %                     mms2_fgm_timedata_raw,mms2_fgm_bdata_raw,...
                %                     mms3_fgm_timedata_raw,mms3_fgm_bdata_raw,...
                %                     mms4_fgm_timedata_raw,mms4_fgm_bdata_raw,...
                %                     mms1_mec_timedata_raw,mms1_mec_rdata_raw,...
                %                     mms2_mec_timedata_raw,mms2_mec_rdata_raw,...
                %                     mms3_mec_timedata_raw,mms3_mec_rdata_raw,...
                %                     mms4_mec_timedata_raw,mms4_mec_rdata_raw,...
                %                     timing_window,'brst');
            end
            
        end
    end
    
    
    if event_analysis == 1 || event_analysis == 0.5 || event_analysis == 1.5 || event_analysis == 2
        plot_MVABslidingwindow(event_start,event_end,mms1_fgm_timedata_srvy, mms1_fgm_bdata_srvy,n_cs,sliding_minratio);
        
        %Inner Leading Boundary
        [MVA_n1,timing_n1,timing_v1,MVAB_timing_angle1] = plot_complete_boundary_analysis(event_start,event_end,leading_start,leading_end,left_InnerEdge,right_InnerEdge,...
            mms1_fgm_timedata_raw,mms1_fgm_bdata_raw,...
            mms2_fgm_timedata_raw,mms2_fgm_bdata_raw,...
            mms3_fgm_timedata_raw,mms3_fgm_bdata_raw,...
            mms4_fgm_timedata_raw,mms4_fgm_bdata_raw,...
            mms1_mec_timedata_raw,mms1_mec_rdata_raw,...
            mms2_mec_timedata_raw,mms2_mec_rdata_raw,...
            mms3_mec_timedata_raw,mms3_mec_rdata_raw,...
            mms4_mec_timedata_raw,mms4_mec_rdata_raw)
        
        %Inner Trailing Boundary
        [MVA_n2,timing_n2,timing_v2,MVAB_timing_angle2] = plot_complete_boundary_analysis(event_start,event_end,trailing_start,trailing_end,left_InnerEdge,right_InnerEdge,...
            mms1_fgm_timedata_raw,mms1_fgm_bdata_raw,...
            mms2_fgm_timedata_raw,mms2_fgm_bdata_raw,...
            mms3_fgm_timedata_raw,mms3_fgm_bdata_raw,...
            mms4_fgm_timedata_raw,mms4_fgm_bdata_raw,...
            mms1_mec_timedata_raw,mms1_mec_rdata_raw,...
            mms2_mec_timedata_raw,mms2_mec_rdata_raw,...
            mms3_mec_timedata_raw,mms3_mec_rdata_raw,...
            mms4_mec_timedata_raw,mms4_mec_rdata_raw)
        
        
        
    end
    %% %Event Parameters
    if event_analysis == 1.5 || event_analysis == 2 || event_analysis == 3 || event_analysis == 4
        if event_analysis == 4
            timing_n1=[0,0,0];
            timing_n2=[0,0,0];
            timing_v1=0;
            timing_v2=0;
            MVA_n1=[0,0,0];
            MVA_n2=[0,0,0];
            MVAB_timing_angle1=0;
            MVAB_timing_angle2=0;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%Event Parameters%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if Event_number==154
            n_cs = [0.28, -0.36, -0.89]
        end %This is because if we use the MMS TD normal, the length is 0.334....mostly Z, we use the TH-B normal instead
        totalDuration = findDuration(left_OuterEdge,right_OuterEdge);
        coreDuration = findDuration(left_InnerEdge,right_InnerEdge);
        
        %Position Magnitudes
        Rxy = sqrt(r_sc(1)^2 + r_sc(2)^2);
        Ryz = sqrt(r_sc(2)^2 + r_sc(3)^2);
        Rmag = sqrt(r_sc(1)^2 + r_sc(2)^2 + r_sc(3)^2);
        %currentsheet_n=currentsheet_n.*[-1,1,1]; %-X GSE for CS normal,
        %currentsheet_n(1)=abs(currentsheet_n(1)); %+X GSE for CS normal
        
        %%% Current sheet is aligned to +X direction.,according to Schwartz et al. 2018
        %n_cs * v_Sw < 0, since v_sw is always negative, n_cs has to always point +X
        %Find average values before and after event
        [B_pre,B_post] = calculate_prepostAverages(mms1_fgm_timedata_raw,mms1_fgm_bdata_raw(:,1:3),event_start,event_end,'B',5)
        %         [B_pre,B_post] = pre_post(mms1_fgm_timedata_raw,mms1_fgm_bdata_raw,event_start,event_end);
        B_pre_norm = norm(B_pre);
        B_post_norm = norm(B_post);
        
        [V_pre,V_post] = calculate_prepostAverages(fpi_i_timedata,fpi_i_vdata,event_start,event_end,'i',5)
        V_pre_norm = norm(V_pre);
        V_post_norm = norm(V_post);
        [n_pre,n_post] = calculate_prepostAverages(fpi_i_timedata,fpi_i_ndata,event_start,event_end,'i',5);
        E_pre = -cross(V_pre,B_pre(1:3));
        E_post= -cross(V_post,B_post(1:3));
        
        
        %Load OMNI Data averages
        [~,Omni_Bx,Omni_By,Omni_Bz,Omni_B,Omni_Vx,Omni_Vy,Omni_Vz,Omni_V,Omni_n,~,~,Omni_beta,Omni_M,Omni_MGSM] = load_omni(event_start,event_start,n);
        
        %Fractional difference between reflected and solar wind speeds
        %         reflectedDownSpeed = downDensityRatio*(V_pre_norm - downstreamSpeed) + downstreamSpeed
        %         reflectedUpSpeed = upDensityRatio*(V_post_norm - upstreamSpeed) + upstreamSpeed
        %
        %         fractionalDownRatio = (downstreamSpeed - reflectedDownSpeed)/downstreamSpeed
        %         fractionalUpRatio = (upstreamSpeed - reflectedUpSpeed)/upstreamSpeed
        
        %[downstreamSpeed,upstreamSpeed,~,~] = calculate_updownSW(event_start,event_end,Event_number)
        downstreamSpeed = swspeedEvents.data(Event_number-1,2);
        upstreamSpeed = swspeedEvents.data(Event_number-1,3);
        downstreamDensity = swspeedEvents.data(Event_number-1,4);
        upstreamDensity = swspeedEvents.data(Event_number-1,5);
        downVx = swVelocityEvents.data(Event_number-1,2);
        downVy = swVelocityEvents.data(Event_number-1,3);
        downVz = swVelocityEvents.data(Event_number-1,4);
        upVx = swVelocityEvents.data(Event_number-1,5);
        upVy = swVelocityEvents.data(Event_number-1,6);
        upVz = swVelocityEvents.data(Event_number-1,7);
        
        V_down = [downVx,downVy,downVz];
        V_down_norm = norm(V_down);
        
        V_up= [upVx,upVy,upVz];
        V_up_norm = norm(V_up);
        
        
        
        
        
        
        %Cone Angle
        %         ConeAngle_pre = asind( B_pre(1)/norm(B_pre))
        %         ConeAngle_post = asind( B_post(1)/norm(B_post))
        
        ConeAngle_pre = angle( B_pre_cs,[1;0;0])
        %ConeAngle_pre = acosd( norm(B_pre)/[1,0,0])
        ConeAngle_post = angle( B_post_cs,[1;0;0])
        %ConeAngle_post = acosd( norm(B_post)/[1,0,0])
        
        ConeAngleV_pre = angle(V_pre,B_pre_cs)
        ConeAngleV_post = angle(V_post,B_post_cs)
        
        %Clock Angle
        ClockAngle_pre = atand(B_pre_cs(2)/B_pre_cs(3));
        ClockAngle_post = atand(B_post_cs(2)/B_post_cs(3));
        
        %         DiscontinuityRatio1 = abs(diff([norm(B_pre),norm(B_post)]))/max([norm(B_pre),norm(B_post)])
        %         DiscontinuityRatio2 = abs(max([dot(B_pre,currentsheet_n),dot(B_post,currentsheet_n)]))/max([norm(B_pre),norm(B_post)])
        
        %Mach Number
        mass_ion = 1.6726219e-27; %kg
        mu_o = 1.25663706e-6;
        
        VA_pre = (norm(B_pre)*1e-9/(mass_ion*n_pre*1e6*mu_o)^(1/2))/1e3
        M_pre = norm(V_down)/VA_pre
        %                 M_pre = downstreamSpeed/VA_pre
        
        VA_post = (norm(B_post)*1e-9/(mass_ion*n_post*1e6*mu_o)^(1/2))/1e3
        M_post = norm(V_up)/VA_post
        %                 M_post = upstreamSpeed/VA_post
        
        %Current Sheet Speeds, Use V_sw before Event
        V_pre_ncs  = dot(V_pre,n_cs);
        V_post_ncs = dot(V_post,n_cs);
        V_ncs      = V_pre_ncs
        V_ncs      = V_post_ncs
        
        
        
        
        %magnetic shear angle, if less than 30, probably a SHFA, no
        %discontinuity.
        Shear_angle = angle(B_pre_cs,B_post_cs)
        
        %Calculate closest distance, and distance to current sheet connection, and shock normal
        %         [shock_normal,HFAtoBS_Distance,bs_pos,closest_shock_normal,ClosestDistance,closest_point,Intersection] = calculate_MerkaBS_CS(...
        %             mms1_fgm_timedata_raw,...
        %             mms1_mec_timedata_raw,mms1_mec_rdata_raw,...
        %             mms2_mec_timedata_raw,mms2_mec_rdata_raw,...
        %             mms3_mec_timedata_raw,mms3_mec_rdata_raw,...
        %             mms4_mec_timedata_raw,mms4_mec_rdata_raw,...
        %             event_start,event_end,...
        %             left_InnerEdge,right_InnerEdge,...
        %             leading_start,leading_end,...
        %             trailing_start,trailing_end,...
        %             n_cs,V_pre,M_pre)
        
        %Calculate closest distance, and distance to current sheet connection, and shock normal
        [shock_normal,HFAtoBS_Distance,bs_pos,closest_shock_normal,ClosestDistance,closest_point] = calculate_MerkaBS_CS2(...
            mms1_fgm_timedata_raw,...
            mms1_mec_timedata_raw,mms1_mec_rdata_raw,...
            mms2_mec_timedata_raw,mms2_mec_rdata_raw,...
            mms3_mec_timedata_raw,mms3_mec_rdata_raw,...
            mms4_mec_timedata_raw,mms4_mec_rdata_raw,...
            left_InnerEdge,right_InnerEdge,...
            n_cs,downstreamDensity,V_down,M_pre)
        
        
        %CAlculate MAhx number along bow shock normal
        downstreamSWspeedAlongBSNormal = abs(dot(shock_normal,V_down));
        upstreamSWspeedAlongBSNormal = abs(dot(shock_normal,V_up));
        downstreamMBS_pre = downstreamSWspeedAlongBSNormal./VA_pre;
        upstreamMBS_pre = upstreamSWspeedAlongBSNormal./VA_post;
        
        
        %Magnetosonic Mach number
        downstreamMGSMBS_pre =  downstreamSWspeedAlongBSNormal./V_down_norm .* Omni_MGSM;
        upstreamMGSMBS_post =  upstreamSWspeedAlongBSNormal./V_up_norm .* Omni_MGSM;
        


        %Omni data
        Omni_Ma = Omni_M;
        Omni_Ms = Omni_MGSM;
        Omni_MaBS = Omni_M.*closest_shock_normal(1);
        Omni_MsBS = Omni_MGSM.*closest_shock_normal(1);
        Omni_thetaBn = angle(closest_shock_normal,[Omni_Bx,Omni_By,Omni_Bz]); 
        Omni_coneAngle = acosd(Omni_Bx/Omni_B);
        
        
        
        Intersection=1;
        %Calculate event edges interscetion with BS and with each other, as well as the shock normal
        %for each edge and the event height from the intersection point of the two planes
        %         [~,leading_shock_normal,...
        %             ~,trailing_shock_normal,...
        %             leading_DistanceFromIntersectiontoBS,trailing_DistanceFromIntersectiontoBS,...
        %             event_width,event_height] = calculate_MerkaBS_EventEdges(...
        %             mms1_fgm_timedata_raw,...
        %             mms1_mec_timedata_raw,mms1_mec_rdata_raw,...
        %             mms2_mec_timedata_raw,mms2_mec_rdata_raw,...
        %             mms3_mec_timedata_raw,mms3_mec_rdata_raw,...
        %             mms4_mec_timedata_raw,mms4_mec_rdata_raw,...
        %             event_start,event_end,...
        %             left_InnerEdge,right_InnerEdge,...
        %             leading_start,leading_end,...
        %             trailing_start,trailing_end,...
        %             timing_n1,timing_n2,...
        %             V_pre,M_pre,timing_v1)
        
        [~,leading_shock_normal,...
            ~,trailing_shock_normal,...
            leading_DistanceFromIntersectiontoBS,trailing_DistanceFromIntersectiontoBS,...
            event_width,event_height] = calculate_MerkaBS_EventEdges2(...
            mms1_fgm_timedata_raw,...
            mms1_mec_timedata_raw,mms1_mec_rdata_raw,...
            mms2_mec_timedata_raw,mms2_mec_rdata_raw,...
            mms3_mec_timedata_raw,mms3_mec_rdata_raw,...
            mms4_mec_timedata_raw,mms4_mec_rdata_raw,...
            event_start,event_end,...
            left_InnerEdge,right_InnerEdge,...
            leading_start,leading_end,...
            trailing_start,trailing_end,...
            timing_n1,timing_n2,...
            downstreamDensity,V_down,M_pre,timing_v1)
        
        
        
        edgesAngle = 180-angle(timing_n1,timing_n2);
        
        %event edges normals are pointing inward, fix back to +x
        if timing_n1(1) < 0
            timing_n1 = -timing_n1; %force x component to be positive, according to Schwartz et al. 2018.
            timing_v1 = -timing_v1;
        end
        
        if timing_n2(1) < 0
            timing_n2 = -timing_n2; %force x component to be positive, according to Schwartz et al. 2018.
            timing_v2 = -timing_v2;
        end
        
        %shock angles outside of each trailing each at the bs connection position
        shockAngle_leading_edge = acosd(dot(B_pre_cs,leading_shock_normal)/norm(B_pre_cs))
        shockAngle_trailing_edge = acosd(dot(B_post_cs,trailing_shock_normal)/norm(B_post_cs))
        
        %Caclulate shock angles via intersection of B vector at BS surface starting from MMS
        %position
        [downstream_shocknormal,downstream_shockAngle] = calculate_BBSShockAngle(...
            mms1_mec_timedata_raw,mms1_mec_rdata_raw,...
            mms1_fgm_timedata_raw,...
            left_OuterEdge,B_pre_cs,...
            downstreamDensity,V_down,M_pre)
        
        [upstream_shocknormal,upstream_shockAngle] = calculate_BBSShockAngle(...
            mms1_mec_timedata_raw,mms1_mec_rdata_raw,...
            mms1_fgm_timedata_raw,...
            right_OuterEdge,B_post_cs,...
            upstreamDensity,V_up,M_post)
        
        %Event area, in the plane of spacecraft trajectory
        heron_s = (event_width + leading_DistanceFromIntersectiontoBS + trailing_DistanceFromIntersectiontoBS)/2;
        event_area = sqrt(heron_s * (heron_s-event_width) * (heron_s-leading_DistanceFromIntersectiontoBS) * (heron_s-trailing_DistanceFromIntersectiontoBS))
        check_event_edges_angle = acosd(( leading_DistanceFromIntersectiontoBS^2 + trailing_DistanceFromIntersectiontoBS^2 - event_width^2)/(2*leading_DistanceFromIntersectiontoBS*trailing_DistanceFromIntersectiontoBS))
        
        
        %Shock Angles before and after event
        %If HFA, use CS-BS connection shock normal, if SHFA or FB, use Closest Shock Normal
        if Shear_angle >= ShearCutoff
            Theta_B_pre_n   = acosd(dot(B_pre_cs,shock_normal)/norm(B_pre_cs))
            Theta_B_post_n  = acosd(dot(B_post_cs,shock_normal)/norm(B_post_cs))
            Theta_BpreBpost = acosd(dot(B_pre_cs,B_post_cs)/(norm(B_pre_cs)*norm(B_post_cs)));
        else
            Theta_B_pre_n   = acosd(dot(B_pre_cs,closest_shock_normal)/norm(B_pre_cs))
            Theta_B_post_n  = acosd(dot(B_post_cs,closest_shock_normal)/norm(B_post_cs))
            Theta_BpreBpost = acosd(dot(B_pre_cs,B_post_cs)/(norm(B_pre_cs)*norm(B_post_cs)));
        end
        %         Theta_B_pre_n   = acosd(dot(B_pre_cs,shock_normal)/norm(B_pre_cs))
        %         Theta_B_post_n  = acosd(dot(B_post_cs,shock_normal)/norm(B_post_cs))
        %         Theta_BpreBpost = acosd(dot(B_pre_cs,B_post_cs)/(norm(B_pre_cs)*norm(B_post_cs)));
        
        if strcmp(downstream_geometry,'ll')
            downstream_geometry = 'Quasi-Parallel';
        elseif strcmp(downstream_geometry,'pp')
            downstream_geometry = 'Quasi-Perpendicular';
        end
        if strcmp(upstream_geometry,'ll')
            upstream_geometry = 'Quasi-Parallel';
        elseif strcmp(upstream_geometry,'pp')
            upstream_geometry = 'Quasi-Perpendicular';
        end
        
        
        %Electric field Direction to CS
        E_pre_dot_cs_n  = dot(E_pre,n_cs) %positive is towards
        E_post_dot_cs_n = dot(E_post,n_cs) %positive is away
        %E_post_dot_cs_n = (-1)*dot(E_post,currentsheet_n) %positive is towards
        
        
        %Electric field Direction to CS
        E_pre_cs_n_angle  = angle(E_pre,n_cs) %Less than 90 is away
        E_post_cs_n_angle = angle(E_post,n_cs) %Greater than 90 is towards
        
        %Angle between CS and BS
        if Shear_angle >= ShearCutoff
            Theta_csbs = acosd(dot(shock_normal,n_cs))
        else
            Theta_csbs = acosd(dot(closest_shock_normal,n_cs))
        end
        %         %Just use same method for SHFAs as HFAs.
        %         Theta_csbs = acosd(dot(shock_normal,n_cs))
        
        %Transversal Speed
        if Shear_angle >= ShearCutoff
            V_tr = V_ncs/((sind(Theta_csbs))^2)*(n_cs-shock_normal*cosd(Theta_csbs))
        else
            V_tr = V_ncs/((sind(Theta_csbs))^2)*(n_cs-closest_shock_normal*cosd(Theta_csbs))
        end
        %         V_tr = V_ncs/((sind(Theta_csbs))^2)*(n_cs-shock_normal*cosd(Theta_csbs))
        mag_V_tr = norm(V_tr)
        
        
        
        %angle with solar wind velocity
        Theta_csswpre=acosd(dot(V_pre,n_cs)/norm(V_pre));
        Theta_csswpost=acosd(dot(V_post,n_cs)/norm(V_post));
        %
        %         if Shear_angle >= ShearCutoff
        %             Theta_bsswpre=acosd(dot(V_pre,shock_normal)/norm(V_pre));
        %             Theta_bsswpost=acosd(dot(V_post,shock_normal)/norm(V_post));
        %         else
        %             Theta_bsswpre=acosd(dot(V_pre,closest_shock_normal)/norm(V_pre));
        %             Theta_bsswpost=acosd(dot(V_post,closest_shock_normal)/norm(V_post));
        %         end
        
        
        %Pre Vg
        if Shear_angle >= ShearCutoff
            V_g_pre = 2 * abs(dot(V_pre,shock_normal) * sind(Theta_B_pre_n));
        else
            V_g_pre = 2 * abs(dot(V_pre,closest_shock_normal) * sind(Theta_B_pre_n));
        end
        %         V_g_pre = 2 * abs(dot(V_pre,shock_normal) * sind(Theta_B_pre_n));
        V_tr_pre_over_V_g = abs(mag_V_tr/V_g_pre)
        
        %V_tr_pre_over_V_g_2 = abs(cosd(Theta_csswpre))/(2*cosd(Theta_bsswpre)*sind(Theta_B_pre_n)*sind(Theta_csbs))
        
        %Post Vg
        if Shear_angle >= ShearCutoff
            V_g_post = 2*abs(dot(V_post,shock_normal) * sind(Theta_B_post_n));
        else
            V_g_post = 2*abs(dot(V_post,closest_shock_normal) * sind(Theta_B_post_n));
        end
        %         V_g_post = 2*abs(dot(V_post,shock_normal) * sind(Theta_B_post_n));
        V_tr_post_over_V_g = abs(mag_V_tr/V_g_post)
        
        %V_tr_post_over_V_g_2 = abs(cosd(Theta_csswpost))/(2*cosd(Theta_bsswpost)*sind(Theta_B_post_n)*sind(Theta_csbs))
        
        
        %size
        Size = mag_V_tr*coreDuration/6371.2
        %         Size = mag_V_tr*totalDuration/6371.2
        
        
        
        
        Vn_1 = timing_n1.*timing_v1
        Vn_2 = timing_n2.*timing_v2
        
        
        
        
        %Event Boundary Shock Angle
        Leading_shockAngle = angle(B_pre(1:3),timing_n1)
        Trailing_shockAngle = angle(B_post(1:3),timing_n2)
        
        
        Delta_V =(Vn_2-Vn_1)
        %HFA_exp_V = abs(dot(Delta_V,N_cs))
        %Delta_V=(392-325);
        HFA_exp_V=(dot(Delta_V,n_cs))
        
        Event_age = Size/abs(HFA_exp_V)*6371.2
        Distance_traveled = Event_age*mag_V_tr/6371.2
        
        %Expansion in SW frame 2019October
        V_solarwindFrame = (V_pre+V_post)./2;
        V1sw = timing_v1 - dot(V_solarwindFrame,timing_n1);
        V2sw = timing_v2 - dot(V_solarwindFrame,timing_n2);
        Delta_Vsw = V2sw*timing_n2-V1sw*timing_n1
        Delta_Vsw_mag = norm(Delta_Vsw)
        HFA_exp_Vsw = (dot(Delta_Vsw,n_cs))
        
        %Use Solar Wind Speed
        %         V_solarwindFrame = ([Omni_Vx,Omni_Vy,Omni_Vz])
        V1 = timing_v1 - V_solarwindFrame*timing_n1'
        V2 = timing_v2 - V_solarwindFrame*timing_n2'
        %V2 = (V2.*timing_n2)*timing_n1' %Projected onto coordinates of V1
        BoundariesExpansion = V2+V1 %positive is contracting, negative is expanding
        BoundariesExpansion_Minus = V2-V1 %THIS IS THE MAIN BOUNDARY EXPANSION SPEED
        
        
        V1b = (timing_v1.*timing_n1 - V_solarwindFrame)
        V2b = (timing_v2.*timing_n2 - V_solarwindFrame)
        V1b_mag = dot(V1b,timing_n1)
        V2b_mag = dot(V2b,timing_n2)
        %V2 = (V2.*timing_n2)*timing_n1' %Projected onto coordinates of V1
        BoundariesExpansionb = V2b_mag+V1b_mag %positive is contracting, negative is expanding
        BoundariesExpansionb_Minus = V2b_mag-V1b_mag
        
        
        %This gives the least number of ocntracting events.
        %V_solarwindFramec = [-(downstreamSpeed+upstreamSpeed)./2, 0 ,0 ]; %Since changing to IDL sw speeds, we got more negative contraction
        %V_solarwindFramec = (V_pre+V_post)./2
        V_solarwindFramec = (V_down+V_up)./2
        
        [~,indexofMaxV] =  max([V_down_norm,V_up_norm]);%,downstreamSpeed,upstreamSpeed]);
        

            if indexofMaxV == 1
                V_solarwindFramec = V_down
            elseif indexofMaxV == 2
                V_solarwindFramec = V_up
%             elseif indexofMaxV == 3
%                  V_solarwindFramec = [downstreamSpeed,0,0]
%             elseif indexofMaxV == 4
%                  V_solarwindFramec = [upstreamSpeed,0,0]
            end
        
        V1c = (timing_v1.*timing_n1 - V_solarwindFramec)
        V2c = (timing_v2.*timing_n2 - V_solarwindFramec)
        V1c_mag = dot(V1c,-timing_n1)
        V2c_mag = dot(V2c,timing_n2)
        %V2 = (V2.*timing_n2)*timing_n1' %Projected onto coordinates of V1
        BoundariesExpansionc = +V2c_mag+V1c_mag %positive is contracting, negative is expanding
        BoundariesExpansionc_Minus = V2c_mag-V1c_mag
        
        
        %Boundaries Expansion August30,2020
        n_leading = timing_n1;
        v_leading = timing_v1;
        leading_velocity = v_leading.*n_leading;
        if Shear_angle > ShearCutoff
            
            ncs_leading_angle = angle(-n_cs,leading_velocity)
            if (ncs_leading_angle < 90) && (n_leading(1) > 0 ) %If Less than 90 degrees, make nleading negative (same as n_cs)
                
                n_leading = - n_leading;
                v_leading  = - v_leading;
                
            elseif (ncs_leading_angle > 90) && (n_leading(1) < 0 ) %If Greater than 90 degrees, making nleading positive (opposite of n_cs)
                n_leading = - n_leading;
                v_leading   = - v_leading;
                
            end
            if n_leading(1).*v_leading > 0
                v_leading = -(v_leading);
            end
        end
        
        %Expansion Speed Calculation
        V1_2020 = (v_leading - dot(V_down,n_leading,2));
        V2_2020 = (-abs(timing_v2) - dot(V_up,timing_n2,2));
        BoundariesExpansion_2020 = V1_2020 + V2_2020;
        
        
        
        Size2 = abs(BoundariesExpansionc*coreDuration/6371.2)
        %Size along n_cs, 2/20/2020
        Age2 = Size/dot(V2c-V1c,n_cs)*6371.2
        Age2p = Size/dot(V2c+V1c,n_cs)*6371.2
        
        theta_n1_sunearth = angle([1,0,0],timing_n1)
        theta_n2_sunearth = angle([1,0,0],timing_n2)
        
        %Timing velocity in Vn1 frame
        Vt1 = timing_v1
        Vt2 = (Vn_2)*timing_n1'
        BoundariesExpansion2 = (timing_v2.*timing_n2)*timing_n1'-timing_v1
        
        
        %Core Density STD
        [coreDensity_Sigma,coreDensity_CV] = calculate_coreDensitySTD(left_InnerEdge,right_InnerEdge,fpi_i_timedata,fpi_i_ndata)
        
        [tempCorr] = calculate_coreTempCorr(left_InnerEdge,right_InnerEdge,fpi_i_timedata,fpi_i_tparadata,fpi_i_tperpdata)
        %         [core_BmagCorr,core_vxCorr,core_vmagCorr,core_tparaCorr,core_tperpCorr] = calculate_coreCorr(left_InnerEdge,right_InnerEdge,mms1_fgm_timedata_raw,mms1_fgm_bdata_raw,fpi_i_timedata,fpi_i_ndata,fpi_i_vdata,fpi_i_tparadata,fpi_i_tperpdata)
        
        % %         %size
        % %         Size = mag_V_tr*coreDuration/6371.2
        % %
        % %         Vn_1 = timing_n1.*timing_v1
        % %         Vn_2 = (timing_n2.*timing_v2).*timing_n1
        % %
        % %         Delta_V =(Vn_2-Vn_1)
        % %         %HFA_exp_V = abs(dot(Delta_V,N_cs))
        % %         %Delta_V=(392-325);
        % %         HFA_exp_V=(dot(Delta_V,currentsheet_n))
        % %
        % %         %Xiao et al. 2015
        % %         V_solarwindFrame = (V_pre+V_post)./2;
        % %         V1 = timing_v1 - V_solarwindFrame*timing_n1';
        % %         V2 = timing_v2 - V_solarwindFrame*timing_n2';
        % %         V2 = (V2.*timing_n2)*timing_n1'; %Projected onto coordinates of V1
        % %         theta_n1_sunearth = angle([1,0,0],timing_n1);
        % %         theta_n2_sunearth = angle([1,0,0],timing_n2);
        % %         BoundariesExpansion = V2-V1; %positive is contracting, negative is expanding
        % %         BoundariesExpansion2 = (timing_v2.*timing_n2)*timing_n1'-timing_v1;
        
                %Substructure Analysis
        [Substructure,~,~,~,~,~,~,...
            NCoreratio,CoreDeflectionAngle,CoreVratio,CoreDPratio,CoreTratio,nVcorr,nBcorr,...
            maxNSWratio,SWMaxDeflectionAngle,maxVSWratio,maxDPSWratio,maxTSWratio,...
            maxBcoreratio,maxBSWratio,...
            SSSize,SSSize_in_electronScales,SS_speed,SS_duration,durationRatio,core_ion_gyroradius,SSstart_indexFraction,SSSizeBulkFlow,SSSizeBulkVflow_ionGyro,...
            CoreSWDeflectionAngle,maxTecoreratio,maxTeSWRatio,SS_i_nTCorr,SS_e_nTCorr,SSSize_in_IonScales,maxNsigma,ionGyroradius,ionInertiallength,...
            SSSizeBulkVflow_electronGyro,electronGyroradius,electronInertiallength] = calculate_densityStats_2020(...
            fpi_i_timedata,fpi_i_ndata,fpi_i_vdata,fpi_i_tparadata,fpi_i_tperpdata,...
            mms1_fgm_timedata_srvy,mms1_fgm_bdata_srvy,fpi_e_timedata,fpi_e_ndata,fpi_e_tperpdata,fpi_e_tparadata,...
            event_start,event_end,left_OuterEdge,right_OuterEdge,left_InnerEdge,right_InnerEdge,...
            ssSTDs,minDataPointsforSS,...
            mms1_fgm_timedata_raw,mms1_fgm_bdata_raw,mms1_mec_timedata_raw,mms1_mec_rdata_raw,...
            mms2_fgm_timedata_raw,mms2_fgm_bdata_raw,mms2_mec_timedata_raw,mms2_mec_rdata_raw,...
            mms3_fgm_timedata_raw,mms3_fgm_bdata_raw,mms3_mec_timedata_raw,mms3_mec_rdata_raw,...
            mms4_fgm_timedata_raw,mms4_fgm_bdata_raw,mms4_mec_timedata_raw,mms4_mec_rdata_raw,...
            V_down,V_up,downstreamDensity,upstreamDensity)
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% %Write to Excel Spreadsheet; Event Parameters Database
        cd '~/Library/Mobile Documents/com~apple~CloudDocs/Research/Analysis'
        databaseFilename = 'EventParameterDatabase_All_Events_2021_3.xlsx';
        %         databaseFilename = 'EventParameterDatabase.xlsx';
        
        % if Substructure == 0
        %    sub = 'no';
        % else
        %     sub = 'yes';
        % end
        
        
        %         %Check if magnetic shear is less than 30, then SHFA
        %         if strcmp(Event_Type,'HFA') || strcmp(Event_Type,'SHFA')
        %             if Shear_angle < ShearCutoff
        %                 Event_Type = 'SHFA';
        %             else
        %                 Event_Type = 'HFA';
        %             end
        %         elseif strcmp(Event_Type,'Foreshock Bubble')
        %         end
        %
        
        %Check if magnetic shear is less than 30, then SHFA
        if strcmp(Event_Type,'HFA') || strcmp(Event_Type,'SHFA') || strcmp(Event_Type,'Foreshock Bubble')
            if Shear_angle < ShearCutoff
                Event_Type = 'SHFA';
            else
%                 Event_Type = 'HFA';
            end
        end
        
        
        %Criteria 1
        %If density STD is below 1.0, then there is no substructure
        %         if coreDensity_Sigma <= 1.0
        %             Substructure = 0;
        %         else
        %             Substructure = 1;
        %         end
        
        %Criteria 2
        %If density CV is below 0.25, then there is no substructure
        
        %         if coreDensity_CV <= 0.25
        %             Substructure = 0;
        %         else
        %             Substructure = 1;
        %         end
        
        

        
        
        
        if Substructure == 1
            SizeIonLength = Size*6371.2/(core_ion_gyroradius)
            Size2IonLength = Size2*6371.2/(core_ion_gyroradius)
            CoreSSBeginIonLength1 = SSstart_indexFraction*SizeIonLength
            CoreSSBeginIonLength2 = SSstart_indexFraction*Size2IonLength
            SSSize_durationRatioSize = durationRatio*Size
        else
            SizeIonLength = {[]}
            Size2IonLength = {[]}
            CoreSSBeginIonLength1 = {[]}
            CoreSSBeginIonLength2 = {[]}
            SSSize_durationRatioSize = {[]}
        end
        
        %Check for substructures
        %         [dynamic_pressure] = calculate_dynamic_pressure(fpi_i_ndata,fpi_i_vdata,'i');
        %         [~,dynamic_pressure,~,~] = crop(fpi_i_timedata,dynamic_pressure,left_InnerEdge,right_InnerEdge);
        %         mad(dynamic_pressure)
        %         mad(dynamic_pressure,1)
        
        %
        %     [~,n,~,~] = crop(fpi_i_timedata,fpi_i_ndata,left_InnerEdge,right_InnerEdge);
        %     [~,vx,~,~] = crop(fpi_i_timedata,fpi_i_vdata(:,1),left_InnerEdge,right_InnerEdge);
        %     [~,vmag,~,~] = crop(fpi_i_timedata,fpi_i_vdata(:,1:3),left_InnerEdge,right_InnerEdge);
        %
        %     vmag = (vmag(:,1).^2 + vmag(:,2).^2 + vmag(:,3).^2).^(1/2);
        %
        %     n_vx_corrcoeff = corrcoef(n,vx)
        %     n_absvx_corrcoeff = corrcoef(n,abs(vx))
        %     n_vmag_corcoeff = corrcoef(n,vmag)
        
        
        
        
        T = table(Event_number,convertCharsToStrings(Event_Type),Substructure,convertCharsToStrings(event_start),convertCharsToStrings(event_end),r_sc,n_cs, bowshock_n,Theta_csbs,Shear_angle,...
            convertCharsToStrings(leading_start),convertCharsToStrings(leading_end),MVA_n1,timing_n1,timing_v1,MVAB_timing_angle1,...
            convertCharsToStrings(trailing_start),convertCharsToStrings(trailing_end),MVA_n2,timing_n2,timing_v2,MVAB_timing_angle2,...
            B_pre,B_post,Theta_BpreBpost,V_pre,V_post,n_pre,n_post,Theta_B_pre_n,Theta_B_post_n,ConeAngle_pre,ConeAngle_post,...
            VA_pre,M_pre,VA_post,M_post,...
            E_pre_cs_n_angle,E_post_cs_n_angle,E_pre_dot_cs_n,E_post_dot_cs_n,...
            V_ncs,V_tr,mag_V_tr,...
            V_g_pre,V_tr_pre_over_V_g,V_g_post,V_tr_post_over_V_g,...
            Vn_1,Vn_2, Delta_V, Size, HFA_exp_V,...
            totalDuration,coreDuration,Event_age,Distance_traveled,...
            V1,V2,BoundariesExpansion,Vt1,Vt2,BoundariesExpansion2,theta_n1_sunearth,theta_n2_sunearth,...
            Leading_shockAngle,Trailing_shockAngle,edgesAngle,...
            coreDensity_Sigma,tempCorr,coreDensity_CV,...
            upstreamSpeed,downstreamSpeed,...
            HFAtoBS_Distance,ClosestDistance,shock_normal,closest_shock_normal,bs_pos,closest_point,...
            shockAngle_leading_edge,shockAngle_trailing_edge,...
            leading_DistanceFromIntersectiontoBS,trailing_DistanceFromIntersectiontoBS,...
            event_width,event_height,event_area,check_event_edges_angle,...
            V_pre_norm,V_post_norm,B_pre_norm,B_post_norm,...
            Rxy,Ryz,Rmag,...
            B_pre_cs,B_post_cs,norm(B_pre_cs),norm(B_post_cs),...
            downstream_shockAngle,upstream_shockAngle,...
            Omni_B,Omni_V,Omni_n,Omni_beta,Omni_M,...
            Size2, convertCharsToStrings(downstream_geometry),convertCharsToStrings(upstream_geometry),...
            V1sw,V2sw,Delta_Vsw,Delta_Vsw_mag,HFA_exp_Vsw,BoundariesExpansion_Minus,Intersection,...
            V1b,V2b,V1b_mag,V2b_mag,BoundariesExpansionb,BoundariesExpansionb_Minus,...
            V_solarwindFramec,V1c,V2c,V1c_mag,V2c_mag,BoundariesExpansionc,BoundariesExpansionc_Minus,...
            NCoreratio,CoreDeflectionAngle,CoreVratio,CoreDPratio,CoreTratio,nVcorr,...
            maxNSWratio,SWMaxDeflectionAngle,maxVSWratio,maxDPSWratio,maxTSWratio,...
            maxBcoreratio,maxBSWratio,nBcorr,SSSize,SSSize_in_electronScales,SS_speed,SS_duration,durationRatio,...
            SizeIonLength,Size2IonLength,SSstart_indexFraction,CoreSSBeginIonLength1,CoreSSBeginIonLength2,Age2,...
            SSSizeBulkFlow,SSSizeBulkVflow_ionGyro,SSSize_durationRatioSize,Age2p,downstreamDensity,upstreamDensity,...
            CoreSWDeflectionAngle,maxTecoreratio,maxTeSWRatio,SS_i_nTCorr,SS_e_nTCorr,SSSize_in_IonScales,maxNsigma,ionGyroradius,ionInertiallength,...
            SSSizeBulkVflow_electronGyro,electronGyroradius,electronInertiallength,V_down,V_up,...
            BoundariesExpansion_2020,V1_2020,V2_2020,...
            Omni_M,downstreamSWspeedAlongBSNormal,upstreamSWspeedAlongBSNormal,downstreamMBS_pre,upstreamMBS_pre,Omni_MGSM,downstreamMGSMBS_pre,upstreamMGSMBS_post,...
            Omni_Ma,Omni_Ms,Omni_MaBS,Omni_MsBS,Omni_thetaBn,Omni_coneAngle)
        
        

        %downDensityRatio,upDensityRatio,reflectedDownSpeed,reflectedUpSpeed,fractionalDownRatio,fractionalUpRatio,...
        %n_vx_corrcoeff(1,2),n_absvx_corrcoeff(1,2),n_vmag_corcoeff(1,2),...
        
        if Event_number == 0
            %Don't write to Excel, this is a practice run.
        else
            %Write to Excel on the next line, which corresponds to Event number.
            EventRow = strcat('A',num2str(Event_number));
            writetable(T,databaseFilename,'Sheet',1,'WriteVariableNames',false,'Range',EventRow) %Master sheet of all events
            
            %Separate events into Type,  substructure/no substructure tabs
            %only contain SHFA and HFAs. Tabs 4 and 5 and 6 are HFA, SHFA,
            %FB, respectively
            pause(5)
            
            switch Event_Type
                case 'HFA'
                    writetable(T,databaseFilename,'Sheet',4,'WriteVariableNames',false,'Range',EventRow)
                    %Substructure yes/no sheet division
                    if Substructure == 1 %yes substructure
                        writetable(T,databaseFilename,'Sheet',2,'WriteVariableNames',false,'Range',EventRow)
                        writetable(T,databaseFilename,'Sheet',5,'WriteVariableNames',false,'Range',EventRow)
                    elseif Substructure == 0 %no substructure
                        writetable(T,databaseFilename,'Sheet',3,'WriteVariableNames',false,'Range',EventRow)
                        writetable(T,databaseFilename,'Sheet',6,'WriteVariableNames',false,'Range',EventRow)
                    end
                case 'SHFA'
                    writetable(T,databaseFilename,'Sheet',7,'WriteVariableNames',false,'Range',EventRow)
                    if Substructure == 1 %yes substructure
                        writetable(T,databaseFilename,'Sheet',2,'WriteVariableNames',false,'Range',EventRow)
                        writetable(T,databaseFilename,'Sheet',8,'WriteVariableNames',false,'Range',EventRow)
                    elseif Substructure == 0 %no substructure
                        writetable(T,databaseFilename,'Sheet',3,'WriteVariableNames',false,'Range',EventRow)
                        writetable(T,databaseFilename,'Sheet',9,'WriteVariableNames',false,'Range',EventRow)
                    end
                case 'Foreshock Bubble'
                    writetable(T,databaseFilename,'Sheet',10,'WriteVariableNames',false,'Range',EventRow)
                    if Substructure == 1 %yes substructure
                        writetable(T,databaseFilename,'Sheet',2,'WriteVariableNames',false,'Range',EventRow)
                        writetable(T,databaseFilename,'Sheet',11,'WriteVariableNames',false,'Range',EventRow)
                    elseif Substructure == 0 %no substructure
                        writetable(T,databaseFilename,'Sheet',3,'WriteVariableNames',false,'Range',EventRow)
                        writetable(T,databaseFilename,'Sheet',12,'WriteVariableNames',false,'Range',EventRow)
                    end
                    
            end
            pause(5)
            %Plot concise FPI summary for quick view, all events in one directory.
            cd '~/Library/Mobile Documents/com~apple~CloudDocs/Research/Analyzed Events/All Events_2021_3'
            %             cd '~/Library/Mobile Documents/com~apple~CloudDocs/Research/Analyzed Events/All Events'
            [~,~,~] = plot_mms_observation_summary_concise_mod(Event_number,event_start,event_end,...
                mms1_mec_timedata_raw,mms1_mec_rdata_raw,mms1_fgm_timedata_raw, mms1_fgm_bdata_raw,...
                fpi_e_timedata,fpi_e_ndata,fpi_e_vdata,fpi_e_pressdata,fpi_e_tparadata,fpi_e_tperpdata,fpi_e_edata,fpi_e_espectdata,...
                fpi_i_timedata,fpi_i_ndata,fpi_i_vdata,fpi_i_pressdata,fpi_i_tparadata,fpi_i_tperpdata,fpi_i_edata,fpi_i_espectdata,...
                n_cs,Shear_angle,totalDuration,coreDensity_Sigma,coreDensity_CV);
        end
        
        
        
        
    end
