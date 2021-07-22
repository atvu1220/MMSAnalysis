function [MVA_n,timing_n,timing_v,MVA_timing_angle] = plot_complete_boundary_analysis(event_start,event_end,date_start,date_end,left_InnerEdge,right_InnerEdge,...
        mms1_fgm_timedata_raw,mms1_fgm_bdata_raw,...
        mms2_fgm_timedata_raw,mms2_fgm_bdata_raw,...
        mms3_fgm_timedata_raw,mms3_fgm_bdata_raw,...
        mms4_fgm_timedata_raw,mms4_fgm_bdata_raw,...
        mms1_mec_timedata_raw,mms1_mec_rdata_raw,...
        mms2_mec_timedata_raw,mms2_mec_rdata_raw,...
        mms3_mec_timedata_raw,mms3_mec_rdata_raw,...
        mms4_mec_timedata_raw,mms4_mec_rdata_raw)
    %Andrew Vu 9/25/18
    %version 2 has fixed the problem of the time interval of the "centered" mva
    %and in addition, calculates the total G matrix and averages it over the
    %interval -- rather than the average of the normal vectors.
    %data=spdfcdfread(filename);
    %datainfo=spdfcdfinfo(filename);
    %clear
    %%
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%Load FGM/Mec Data%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %addpath ~/Library/'Mobile Documents'/com~apple~CloudDocs/Research
    %addpath ~/data/
    %cd '~/Library/Mobile Documents/com~apple~CloudDocs/Research/Analysis'
    figure('Position',[0 0 800 700])
    movegui('northwest')
    mms_directory = '/Users/andrewvu/data/mms/';
    probe_num = '1';
    num_plots = 8;
    data_type = 'brst';
    
    plot_pdf = 'yes';
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%Event Date Input%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Enter 0 for placeholders for dates.
    %Overall Event Start and End Times
    
    
    %Specific Portion of the Event for Analysis
    % event_start = '2018-03-01 01:03:45.000';
    % event_end = '2018-03-01 01:04:20.000';
    %
    % date_start = '2018-03-01 01:03:54.600';
    % date_end = '2018-03-01 01:03:55.100';
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%Adjustable Parameters%%%%%%%%%%%%%%%%%%%%%%%%%
    avg_i = 0; %for removing fluctuations from calibrationb errors in G
    fluctuate_mdd = 0;%1 for force all vectors to stay within their same sign
    fluctuate_mva = 0;
    MVAB_window_size = 0; %this is the number of data points on either side of t_0; enter 0 for default window size of the time interval
    
    %Plot the minimum, intermediate, or maximum vectors from eigenvalue problem
    %1 is max, 2 is mid, 3 is min
    mddmm = 1;
    mvamm = 3;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%Center +/- sides Method%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % %center of discontinuity
    % center_time = '2018-03-12 07:36:15.000';
    % %dt is window size, on each side.
    % dt = 1;
    % %Create left and Right Initial Bounds
    % formatIn='yyyy-MM-dd HH:mm:ss.SSS';
    % center_datetime = datetime(center_time,'InputFormat',formatIn);
    % date_start = (center_datetime - seconds(dt));
    % date_end = (center_datetime + seconds(dt));
    %
    % formatIn='yyyy-mm-dd HH:MM:SS.FFF';
    % date_start = datestr(date_start,formatIn);
    % date_end = datestr(date_end,formatIn);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %Set the Strings for which eigenvector to use
    if mddmm==1
        mddl = '{max}';
    elseif mddmm==2
        mddl = '{mid}';
    elseif mddmm==3
        mddl = '{min}';
    end
    
    if mvamm==1
        mval = '{max}';
    elseif mvamm ==2
        mval = '{mid}';
    elseif mvamm==3
        mval = '{min}';
    end
    
    formatIn='yyyy-mm-dd HH:MM:SS.FFF';
    tstart = datenum(date_start,formatIn);
    tend = datenum(date_end,formatIn);
    %%
    % %%%%%%%%%%%%%%%%%%%%%%%%%Load FGM/Mec Data%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % [mms1_fgm_timedata_raw, mms1_fgm_bdata_raw, ~, ~] = load_fgm(event_start,1,data_type);
    % [mms2_fgm_timedata_raw, mms2_fgm_bdata_raw, ~, ~] = load_fgm(event_start,2,data_type);
    % [mms3_fgm_timedata_raw, mms3_fgm_bdata_raw, ~, ~] = load_fgm(event_start,3,data_type);
    % [mms4_fgm_timedata_raw, mms4_fgm_bdata_raw, ~, ~] = load_fgm(event_start,4,data_type);
    %
    % %Mec data has fewer data points than FGM Epheremis
    % [mms1_mec_timedata_raw, mms1_mec_rdata_raw] = load_mec(event_start,1,'brst');
    % [mms2_mec_timedata_raw, mms2_mec_rdata_raw] = load_mec(event_start,2,'brst');
    % [mms3_mec_timedata_raw, mms3_mec_rdata_raw] = load_mec(event_start,3,'brst');
    % [mms4_mec_timedata_raw, mms4_mec_rdata_raw] = load_mec(event_start,4,'brst');
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%Interpolate position and bdata with bdata from mms1 %%%%%%%%%%%%%
    %Interpolate the Data
    
    [~,mms1_mec_rdata_interp] = interpxyz(mms1_mec_timedata_raw,mms1_mec_rdata_raw,mms1_fgm_timedata_raw);
    
    [~,mms2_fgm_bdata_interp] = interpxyz(mms2_fgm_timedata_raw,mms2_fgm_bdata_raw,mms1_fgm_timedata_raw);
    [~,mms2_mec_rdata_interp] = interpxyz(mms2_mec_timedata_raw,mms2_mec_rdata_raw,mms1_fgm_timedata_raw);
    
    [~,mms3_fgm_bdata_interp] = interpxyz(mms3_fgm_timedata_raw,mms3_fgm_bdata_raw,mms1_fgm_timedata_raw);
    [~,mms3_mec_rdata_interp] = interpxyz(mms3_mec_timedata_raw,mms3_mec_rdata_raw,mms1_fgm_timedata_raw);
    
    [~,mms4_fgm_bdata_interp] = interpxyz(mms4_fgm_timedata_raw,mms4_fgm_bdata_raw,mms1_fgm_timedata_raw);
    [~,mms4_mec_rdata_interp] = interpxyz(mms4_mec_timedata_raw,mms4_mec_rdata_raw,mms1_fgm_timedata_raw);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %%%%%%%%%%%%%%%%%Data Cropping and Matrix Formation%%%%%%%%%%%%%%%%%%%%%%%%
    %Find the start and end limits of the event in the data, index of, in
    %datenum format
    
    %These indices are for plotting the entire event
    event_start_index = find(mms1_fgm_timedata_raw >= datenum(event_start,formatIn), 1);
    event_end_index = find(mms1_fgm_timedata_raw >= datenum(event_end,formatIn), 1);
    
    %These indices are for the analysis in the specified region
    start_index = find(mms1_fgm_timedata_raw >= tstart, 1);
    end_index = find(mms1_fgm_timedata_raw >= tend, 1);
    
    % %Crop mec data to specific time period
    [~,mms1_mec_rdata,~,~] = crop(mms1_fgm_timedata_raw,mms1_mec_rdata_interp,date_start,date_end);
    [~,mms2_mec_rdata,~,~] = crop(mms1_fgm_timedata_raw,mms2_mec_rdata_interp,date_start,date_end);
    [~,mms3_mec_rdata,~,~] = crop(mms1_fgm_timedata_raw,mms3_mec_rdata_interp,date_start,date_end);
    [~,mms4_mec_rdata,~,~] = crop(mms1_fgm_timedata_raw,mms4_mec_rdata_interp,date_start,date_end);
    
    %Crop fgm data to specific time period
    [mms1_fgm_timedata,mms1_fgm_bdata,~,~] = crop(mms1_fgm_timedata_raw,mms1_fgm_bdata_raw,date_start,date_end);
    [mms2_fgm_timedata,mms2_fgm_bdata,~,~] = crop(mms1_fgm_timedata_raw,mms2_fgm_bdata_interp,date_start,date_end);
    [mms3_fgm_timedata,mms3_fgm_bdata,~,~] = crop(mms1_fgm_timedata_raw,mms3_fgm_bdata_interp,date_start,date_end);
    [mms4_fgm_timedata,mms4_fgm_bdata,~,~] = crop(mms1_fgm_timedata_raw,mms4_fgm_bdata_interp,date_start,date_end);
    mms_fgm_timedata = [mms1_fgm_timedata,mms2_fgm_timedata,mms3_fgm_timedata,mms4_fgm_timedata];
    
    
    
    %group components of B and r together for faster loops.
    bx=[mms1_fgm_bdata(:,1) mms2_fgm_bdata(:,1) mms3_fgm_bdata(:,1) mms4_fgm_bdata(:,1)];
    by=[mms1_fgm_bdata(:,2) mms2_fgm_bdata(:,2) mms3_fgm_bdata(:,2) mms4_fgm_bdata(:,2)];
    bz=[mms1_fgm_bdata(:,3) mms2_fgm_bdata(:,3) mms3_fgm_bdata(:,3) mms4_fgm_bdata(:,3)];
    
    % bx=[smooth(mms1_fgm_bdata(:,1)) smooth(mms2_fgm_bdata(:,1)) smooth(mms3_fgm_bdata(:,1)) smooth(mms4_fgm_bdata(:,1))];
    % by=[smooth(mms1_fgm_bdata(:,2)) smooth(mms2_fgm_bdata(:,2)) smooth(mms3_fgm_bdata(:,2)) smooth(mms4_fgm_bdata(:,2))];
    % bz=[smooth(mms1_fgm_bdata(:,3)) smooth(mms2_fgm_bdata(:,3)) smooth(mms3_fgm_bdata(:,3)) smooth(mms4_fgm_bdata(:,3))];
    
    rx=[mms1_mec_rdata(:,1) mms2_mec_rdata(:,1) mms3_mec_rdata(:,1) mms4_mec_rdata(:,1)];
    ry=[mms1_mec_rdata(:,2) mms2_mec_rdata(:,2) mms3_mec_rdata(:,2) mms4_mec_rdata(:,2)];
    rz=[mms1_mec_rdata(:,3) mms2_mec_rdata(:,3) mms3_mec_rdata(:,3) mms4_mec_rdata(:,3)];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%MDD Calculations%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Calculate the Center of the Tetrahedron
    rx_0 = mean(rx,2);
    ry_0 = mean(ry,2);
    rz_0 = mean(rz,2);
    
    %Calculate the relative distance of each probe from the center
    rx = rx-rx_0;
    ry = ry-ry_0;
    rz = rz-rz_0;
    
    %Calculate components Volumetric Tensor
    Rxx = (1/4)*sum(rx.^2,2);
    Rxy = (1/4)*sum(rx.*ry,2);
    Rxz = (1/4)*sum(rx.*rz,2);
    Ryy = (1/4)*sum(ry.^2,2);
    Ryz = (1/4)*sum(ry.*rz,2);
    Rzz = (1/4)*sum(rz.^2,2);
    
    %initialize volumetric tensor matrix
    R=zeros(3,3,length(rx));
    R_inv = zeros(3,3,length(rx));
    
    
    %Loop for creating Volumetric Tensor
    for i=1:length(rx)
        R(1,1,i) = Rxx(i);
        R(1,2,i) = Rxy(i);
        R(1,3,i) = Rxz(i);
        R(2,1,i) = Rxy(i);
        R(2,2,i) = Ryy(i);
        R(2,3,i) = Ryz(i);
        R(3,1,i) = Rxz(i);
        R(3,2,i) = Ryz(i);
        R(3,3,i) = Rzz(i);
        
        %Inverse Matrix
        R_inv(:,:,i) = inv(R(:,:,i));
    end
    
    %For combination of probes
    probes = [1,2,3,4];
    Summ_probes = nchoosek(probes,2);
    
    
    % function [grad_B] = gradB(time,Summ_probes,bx,by,bz,rx,ry,rz,R_inv)
    
    %Initialize Matrices and Vectors for G and its eigenvalues/eigenvectors
    grad_B = zeros(3,3,length(rx));
    eig_lambda = zeros(3,3,length(rx));
    eig_vectors = zeros(3,3,length(rx));
    
    
    %Find the grad_B matrix first, and also its eigenvectors and eigenvalues
    for i=1:length(rx)
        %(Grad B)*(Grad_B)T is symmetric
        grad_B(:,:,i) = gradB(i,Summ_probes,probes,bx,by,bz,rx,ry,rz,R_inv);
        L = grad_B(:,:,i)*grad_B(:,:,i)';
        [eig_v, eig_l] = eig(L);
        
        %Order the vectors and eigenvalues from largest to smallest, for all
        %times
        [eig_l,order] = sort(diag(eig_l),'descend');
        eig_v = eig_v(:,order);
        
        %Store each set of eigenvectors and eigenvalues in the 3D array, for
        %each time.
        eig_lambda(:,:,i) = diag(eig_l);
        eig_vectors(:,:,i) = eig_v;
    end
    %For Calibration Error Elimination, Denton et al. 2010
    %Find the largest value for the largest eigenvalue, corresponds to largest G.
    Gmax_index = find(squeeze(eig_lambda(1,1,:)) == max(squeeze(eig_lambda(1,1,:))));
    
    %Calculate the components of the average G in the time interval enclosing
    %the largest G matrix, or largest eigenvalue
    Gxx = mean(squeeze(grad_B(1,1,Gmax_index-avg_i:Gmax_index+avg_i)));
    Gyx = mean(squeeze(grad_B(1,2,Gmax_index-avg_i:Gmax_index+avg_i)));
    Gzx = mean(squeeze(grad_B(1,3,Gmax_index-avg_i:Gmax_index+avg_i)));
    
    Gxy = mean(squeeze(grad_B(2,1,Gmax_index-avg_i:Gmax_index+avg_i)));
    Gyy = mean(squeeze(grad_B(2,2,Gmax_index-avg_i:Gmax_index+avg_i)));
    Gzy = mean(squeeze(grad_B(2,3,Gmax_index-avg_i:Gmax_index+avg_i)));
    
    Gxz = mean(squeeze(grad_B(3,1,Gmax_index-avg_i:Gmax_index+avg_i)));
    Gyz = mean(squeeze(grad_B(3,2,Gmax_index-avg_i:Gmax_index+avg_i)));
    Gzz = mean(squeeze(grad_B(3,3,Gmax_index-avg_i:Gmax_index+avg_i)));
    
    %Average grad_B matrix in the time interval of largest G, or largest eigenvalue
    G_0 = abs([Gxx, Gyx, Gzx;...
        Gxy, Gyy, Gzy;...
        Gxz, Gyz, Gzz]);
    
    %Initialize matrices for MDD error calculation
    div_B = zeros(length(rx),1);
    curl_B = zeros(length(rx),3);
    
    
    %Calculate the Eigenvalues and Eigenvectors of L = GG', where we have
    %subtracted the average value of G in the largest G interval
    for i=1:length(rx)
        %Form G*GT and calculate its eigenvalues and eigenvectors
        grad_B(:,:,i) = grad_B(:,:,i) - G_0;
        L = (grad_B(:,:,i))*(grad_B(:,:,i))';
        [eig_v, eig_l] = eig(L);
        
        %Order the eigenvectors and eigenvalues from largest to smallest
        [eig_l,order] = sort(diag(eig_l),'descend');
        eig_v = eig_v(:,order);
        
        %Store the eigenvectors and eigenvalues for each time
        eig_lambda(:,:,i) = diag(eig_l);
        eig_vectors(:,:,i) = eig_v;
        
        %Calculate the MDD error for each time, divB and curlB
        div_B(i,1) = trace(grad_B(:,:,i));
        curl_B(i,1) = grad_B(2,3,i) - grad_B(3,2,i);
        curl_B(i,2) = grad_B(3,1,i) - grad_B(1,3,i);
        curl_B(i,3) = grad_B(1,2,i) - grad_B(2,1,i);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%
    %%%%%%%%%%%%%%%%%%%%%%%Sliding Window for MVAB normals%%%%%%%%%%%%%%%%%%%%%
    t_0_range = end_index - start_index+1;
    
    %If there is no inputted window size, the sliding window MVAB analysis will default to the same window size
    %specified for the MDD analysis
    if MVAB_window_size == 0
        MVAB_window_size = floor(end_index - start_index);
    end
    
    %initialize loop matrices
    MVAB_lambda_all=zeros(3,t_0_range);
    MVAB_normal_all=zeros(3,t_0_range);
    MVAB_time=zeros(1,t_0_range);
    
    %Sliding Window MVAB for every point same as MDD analysis
    for i=1:t_0_range
        t_0 = mms1_fgm_timedata_raw(i+start_index-1); %center point, loop through all of these between the start/end indices
        left_index = i+start_index-MVAB_window_size/2; %determine the left boundary from the center point, half of window size.
        right_index = i+start_index+MVAB_window_size/2; %determine the right boundary from the center point, half of window size.
        
        bdata_scope = mms1_fgm_bdata_raw(left_index:right_index,1:3); %Grab the Bdata and splice it to just within the left and right boundaries.
        
        [~, l, v] = mvab(bdata_scope(:,1:3)); %Calculate the MVAB from the spliced time interval
        
        MVAB_lambda_all(:,i) = l; %Store the eigenvalues
        MVAB_normal_all(:,i) = v(:,mvamm);%Store only eigenvector chosen at the beginning of the script
        
        MVAB_time(1,i)=t_0;%Store the center point(time) used for MVAB analysis
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Unfluctuate%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %If needed and specified in the beginning, flip the equal and opposite
    %signs of the eigenvectors. Calculate the mean, and determine the overal
    %sign of the eigenvectors
    
    if fluctuate_mdd == 1
        mdd_nx = Unfluctuate(squeeze(eig_vectors(1,mddmm,:)));
        mdd_ny = Unfluctuate(squeeze(eig_vectors(2,mddmm,:)));
        mdd_nz = Unfluctuate(squeeze(eig_vectors(3,mddmm,:)));
    elseif fluctuate_mdd == 0
        mdd_nx = (squeeze(eig_vectors(1,mddmm,:)));
        mdd_ny = (squeeze(eig_vectors(2,mddmm,:)));
        mdd_nz = (squeeze(eig_vectors(3,mddmm,:)));
        
    end
    
    if fluctuate_mva == 1
        MVAB_normal_all(1,:) = Unfluctuate(MVAB_normal_all(1,:));
        MVAB_normal_all(2,:) = Unfluctuate(MVAB_normal_all(2,:));
        MVAB_normal_all(3,:) = Unfluctuate(MVAB_normal_all(3,:));
    elseif fluctuate_mva == 0
        MVAB_normal_all(1,:) = (MVAB_normal_all(1,:));
        MVAB_normal_all(2,:) = (MVAB_normal_all(2,:));
        MVAB_normal_all(3,:) = (MVAB_normal_all(3,:));
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%
    %%%%%%%%%%%%%%%%%%%%%%%Average Angle Calculations%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %%%%%%%%%%%%%%%%%Mean MDD and Center MVA Angle Calculation %%%%%%%%%%%%%%%%
    
    
    %Calulation for the mean G_matrix of all time.
    AA = mean(squeeze(grad_B(1,1,:)));
    BB = mean(squeeze(grad_B(1,2,:)));
    CC = mean(squeeze(grad_B(1,3,:)));
    DD = mean(squeeze(grad_B(2,1,:)));
    EE = mean(squeeze(grad_B(2,2,:)));
    FF = mean(squeeze(grad_B(2,3,:)));
    GG = mean(squeeze(grad_B(3,1,:)));
    HH = mean(squeeze(grad_B(3,2,:)));
    II = mean(squeeze(grad_B(3,3,:)));
    
    %Form the mean G_matrix of all time
    KK = [AA, BB, CC; DD, EE, FF; GG, HH, II];
    %Calculate the eigenvalues and eigenvectors for the mean GGT
    [VV,DD] = eig(KK*KK');
    
    %Order the eigenvectors and eigenvalues from largest to smallest
    [~,order] = sort(diag(DD),'descend');
    VV = VV(:,order);
    
    %Store the specified eigenvector
    Mean_grad_B_MDD_normal = VV(:,mddmm);
    
    formatOut = "yyyy-mm-dd HH:MM:SS.FFF";
    datestr(mms1_fgm_timedata(1),formatOut)
    datestr(mms1_fgm_timedata(end),formatOut)
    
    %Calculate the MVAB for the specified time interval
    [~,mvab_l,mvab_v] = mvab(mms1_fgm_bdata(:,1:3));
    
    
    %Store the eigenvalues and eigenvectors
    MVA_normal = ([mvab_v(1,mvamm), mvab_v(2,mvamm), mvab_v(3,mvamm)]);
    MVA_eigenvalues = (mvab_l)' %#ok<NOPTS>
    
    MVA_min_normal = ([mvab_v(1,3), mvab_v(2,3), mvab_v(3,3)]) %Store the Minimum Eigenvector direction, specifically.
    
    %Because MVA is ambiguous for one direction, determine the best normal
    %direction to compare with MDD. The best normal direction is the one with
    %the smallest angle with MDD.
    
    % % % % % %Calculate all 4 possibilities
    % % % % % MVA_ambig(1)=angle(MVA_normal.*[-1,1,1],Mean_grad_B_MDD_normal);
    % % % % % MVA_ambig(2)=angle(MVA_normal.*[1,-1,1],Mean_grad_B_MDD_normal);
    % % % % % MVA_ambig(3)=angle(MVA_normal.*[1,1,-1],Mean_grad_B_MDD_normal);
    % % % % % MVA_ambig(4)=angle(MVA_normal.*[1,1,1],Mean_grad_B_MDD_normal);
    % % % % %
    % % % % % %Find the normal with the smallest angle with MDD, and replace the edited normal with the old
    % % % % % %normal.
    % % % % % switch find(MVA_ambig==min(MVA_ambig))
    % % % % %     case 1
    % % % % %         MVA_normal = MVA_normal.*[-1,1,1];
    % % % % %     case 2
    % % % % %         MVA_normal = MVA_normal.*[1,-1,1];
    % % % % %     case 3
    % % % % %         MVA_normal = MVA_normal.*[1,1,-1];
    % % % % %     case 4
    % % % % %         MVA_normal = MVA_normal; %#ok<ASGSL>
    % % % % % end
    % % % % %
    % % % % % %The smallest angle between MDD and MVA, possible
    % % % % % mean_MDD_dot_MVA_center = min(MVA_ambig);
    mean_MDD_dot_MVA_center = ambig_angle(MVA_normal,Mean_grad_B_MDD_normal);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %%%%%%%%%%%%%%%%%%%%Calculation of Angle at every point%%%%%%%%%%%%%%%%%%%%
    %Calculations for Point by Point Analysis of MDD and MVA
    %Calculate average of each component for the MDDB Normal
    mean_x = mean(mdd_nx,'omitnan');
    mean_y = mean(mdd_ny,'omitnan');
    mean_z = mean(mdd_nz,'omitnan');
    % mean_x = mean(mdd_nx,'omitnan');
    % mean_y = mean(mdd_ny,'omitnan');
    % mean_z = mean(mdd_nz,'omitnan');
    MDD_mean_normal = ([mean_x, mean_y, mean_z])./norm([mean_x, mean_y, mean_z]);
    
    %Set up the point by points of MDD and MVA
    MDDB_normal_all = [mdd_nx, mdd_ny, mdd_nz]; %MDDB normal at every point
    MVAB_normal_all = MVAB_normal_all'; %MVAB normals at every point
    
    %Initialize the matrices
    MVA_dot_MDD_angle_all = zeros(length(rx),1);
    MVA_dot_MDD_center_all = zeros(length(rx),1);
    MVA_ambig_dot_center= zeros(length(rx),4);
    
    for i=1:length(rx)
        MVA_dot_MDD_angle_all(i,1) = ambig_angle(MVAB_normal_all(i,:),MDDB_normal_all(i,:));
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%
    %%Custom Plotting
    % plots = [B_all';'B_portion';'MDD_Normal';'MDD_Lambda';'MVA_Normal';'MVA_Lambda';'Angle'];
    %
    % for i=1:length(plots)
    %     switch plot(i)
    %         case 'B_all'
    %         case 'B_portion'
    %         case 'MDD_Normal'
    %         case 'MDD_Lambda'
    %         case 'MVA_Normal'
    %         case 'MVA_Lambda'
    %         case 'Angle'
    %     end
    %
    %
    % end
    
    
    %%
    %%%%%%%%%%%%%%%%%%%%%%%Plot Magnetic Field Data%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    f1 = subplot(num_plots,1,1);
    set(gcf,'color','w');
    co = [0 0 1;
        0 1 0;
        1 0 0;
        0 0 0];
    set(gcf,'defaultAxesColorOrder',co)
    
    plot(mms1_fgm_timedata_raw(event_start_index:event_end_index,1),mms1_fgm_bdata_raw(event_start_index:event_end_index,1:4))
    
    line([(mms1_fgm_timedata_raw(start_index)) (mms1_fgm_timedata_raw(start_index))],...
        get(gca,'YLim'),'Color','k','LineWidth',0.5,'LineStyle','--')
    line([(mms1_fgm_timedata_raw(end_index)) (mms1_fgm_timedata_raw(end_index))],...
        get(gca,'YLim'),'Color','k','LineWidth',0.5,'LineStyle','--')
    
    box on
    datetick
    set(gca,'layer','top')
    legend({'B_x', 'B_y', 'B_z','B_{tot}'},'FontSize',10)
    legend('boxoff');
    legend('Location','eastoutside');
    ylabel({'B';'[nT]'},'FontSize', 14)
    xlim([mms1_fgm_timedata_raw(event_start_index), mms1_fgm_timedata_raw(event_end_index)])
    set(f1,'XMinorTick','on','YMinorTick','on')
    
    title_name = 'MMS MDDB&MVAB Normals';
    title(title_name, 'FontSize', 18, 'FontWeight', 'normal')
    
    
    plot_pos = get(f1,'Position');
    set(f1,'Position',[plot_pos(1), plot_pos(2), plot_pos(3), plot_pos(4)]);
    
    l1= findobj(gcf, 'Type', 'Legend');
    l1pos=get(l1, 'Position');
    set(l1,'Position',[plot_pos(1)+plot_pos(3)-2, plot_pos(2)+plot_pos(4)*0.5-l1pos(4)*0.5, l1pos(3), l1pos(4)]);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%Plot Portion of B%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    f2 = subplot(num_plots,1,2);
    
    plot(mms1_fgm_timedata_raw(start_index:end_index),mms1_fgm_bdata_raw(start_index:end_index,1:3))
    box on
    
    set(gca,'layer','top')
    legend({'B_x', 'B_y', 'B_z'},'FontSize',10)
    legend('boxoff');
    legend('Location','eastoutside');
    datetick('x','SS.FFF','keepticks')
    
    ylabel({'B';'[nT]'},'FontSize', 14)
    datetick
    xlim([mms1_fgm_timedata(1) mms1_fgm_timedata(end)])
    
    set(f2,'XMinorTick','on','YMinorTick','on','XTickLabel',[])
    
    set(f2,'Position',[plot_pos(1), plot_pos(2)-1*plot_pos(4)-0.05, plot_pos(3), plot_pos(4)]);
    
    l1= findobj(gcf, 'Type', 'Legend');
    l1pos=get(l1(2), 'Position');
    set(l1(2),'Position',[plot_pos(1)+plot_pos(3)+0.010, plot_pos(2)+1*plot_pos(4)-0.1, l1pos(3), l1pos(4)]);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Plot MDD Normals%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    f3 = subplot(num_plots,1,3);
    
    %Minimum
    plot(mms1_fgm_timedata,mdd_nx)
    hold on
    plot(mms1_fgm_timedata,mdd_ny)
    plot(mms1_fgm_timedata,mdd_nz)
    hold off
    
    box on
    datetick
    set(gca,'layer','top')
    legend({'N_x', 'N_y', 'N_z'},'FontSize',10)
    legend('boxoff');
    legend('Location','eastoutside');
    
    ylabel(strcat('MDD_',mddl),'FontSize', 14)
    
    
    xlim([mms1_fgm_timedata(1) mms1_fgm_timedata(end)])
    
    set(f3,'XMinorTick','on','YMinorTick','on','XTickLabel',[])
    
    set(f3,'Position',[plot_pos(1), plot_pos(2)-2*plot_pos(4)-0.075, plot_pos(3), plot_pos(4)*1]);
    
    l1= findobj(gcf, 'Type', 'Legend');
    l1pos=get(l1(3), 'Position');
    %set(l1(3),'Position',[plot_pos(1)+plot_pos(3)+0.075, plot_pos(2)+2*plot_pos(4)+0.15, l1pos(3), l1pos(4)]);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Plot MDD Eigenvalues%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    f4 = subplot(num_plots,1,4);
    
    plot(mms1_fgm_timedata,sqrt(squeeze(eig_lambda(1,1,:))./squeeze(eig_lambda(2,2,:))))
    hold on
    plot(mms1_fgm_timedata,sqrt(squeeze(eig_lambda(2,2,:))/squeeze(eig_lambda(3,3,:))),'g')
    %plot(mms1_fgm_btimedata,sqrt(squeeze(eig_lambda(3,3,:))/squeeze(eig_lambda(1,1,:))),'r')
    hold off
    
    datetick
    set(gca,'layer','top')
    legend({'\lambda_{max}/\lambda_{mid}', '\lambda_{mid}/\lambda_{min}', '\lambda_{min}/\lambda_{max}'},'FontSize',10)
    legend('boxoff');
    legend('Location','eastoutside');
    
    ylabel(strcat('\lambda_{MDD}'),'FontSize', 14)
    
    xlim([mms1_fgm_timedata(1) mms1_fgm_timedata(end)])
    
    set(f4,'XMinorTick','on','YMinorTick','on','XTickLabel',[])
    
    set(f4,'Position',[plot_pos(1), plot_pos(2)-3*plot_pos(4)-0.1, plot_pos(3), plot_pos(4)*1]);
    
    % l1= findobj(gcf, 'Type', 'Legend');
    % l1pos=get(l1(4), 'Position');
    % set(l1(4),'Position',[plot_pos(1)+plot_pos(3)+0.1, plot_pos(2)+3*plot_pos(4), l1pos(3), l1pos(4)]);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Plot MDD Error%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    f5 = subplot(num_plots,1,5);
    
    plot(mms1_fgm_timedata,abs(div_B./(sqrt(curl_B(:,1).^2+curl_B(:,2).^2+curl_B(:,3).^2))))
    hold on
    plot(mms1_fgm_timedata,abs(div_B./(norm(G_0))))
    line(xlim(), [0.4,0.4], 'LineWidth', 0.25, 'Color', 'k','LineStyle','--');
    hold off
    
    datetick
    set(gca,'layer','top')
    legend({'\nabla\cdotB/\nablaxB', '\nabla\cdotB/max(\nablaB)'},'FontSize',10)
    legend('boxoff');
    legend('Location','eastoutside');
    %xlim([datetime(tstart,'ConvertFrom','datenum') datetime(tend,'ConvertFrom','datenum')]);
    
    ylabel(strcat('MDD_{Error}'),'FontSize', 14)
    
    xlim([mms1_fgm_timedata(1) mms1_fgm_timedata(end)])
    
    set(f5,'XMinorTick','on','YMinorTick','on','XTickLabel',[])
    
    set(f5,'Position',[plot_pos(1), plot_pos(2)-4*plot_pos(4)-0.125, plot_pos(3), plot_pos(4)*1]);
    
    % l1= findobj(gcf, 'Type', 'Legend');
    % l1pos=get(l1(4), 'Position');
    % set(l1(4),'Position',[plot_pos(1)+plot_pos(3)+0.1, plot_pos(2)+3*plot_pos(4), l1pos(3), l1pos(4)]);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Plot MVAB Normals%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    f6 = subplot(num_plots,1,6);
    
    plot(mms1_fgm_timedata,MVAB_normal_all')
    
    datetick
    set(gca,'layer','top')
    legend({'N_x', 'N_y', 'N_z'},'FontSize',10)
    legend('boxoff');
    legend('Location','eastoutside');
    
    ylabel(strcat('MVAB_',mval),'FontSize', 14)
    
    xlim([mms1_fgm_timedata(1) mms1_fgm_timedata(end)])
    
    set(f6,'XMinorTick','on','YMinorTick','on','XTickLabel',[])
    
    set(f6,'Position',[plot_pos(1), plot_pos(2)-5*plot_pos(4)-0.15, plot_pos(3), plot_pos(4)*1]);
    
    % l1= findobj(gcf, 'Type', 'Legend');
    % l1pos=get(l1(4), 'Position');
    % set(l1(4),'Position',[plot_pos(1)+plot_pos(3)+0.1, plot_pos(2)+3*plot_pos(4), l1pos(3), l1pos(4)]);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Plot MVAB Lambdas%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    f7 = subplot(num_plots,1,7);
    MVAB_lambda_12 = MVAB_lambda_all(1,:)./MVAB_lambda_all(2,:);
    % MVAB_lambda_12(MVAB_lambda_12 > 25) = NaN;
    
    MVAB_lambda_23 = MVAB_lambda_all(2,:)./MVAB_lambda_all(3,:);
    % MVAB_lambda_23(MVAB_lambda_23 > 25) = NaN;
    
    plot(mms1_fgm_timedata,MVAB_lambda_12)
    hold on
    plot(mms1_fgm_timedata,MVAB_lambda_23)
    hold off
    
    datetick
    set(gca,'layer','top')
    legend({'\lambda_{max}/\lambda_{min}','\lambda_{mid}/\lambda_{min}'},'FontSize',10)
    legend('boxoff');
    legend('Location','eastoutside');
    %xlim([datetime(tstart,'ConvertFrom','datenum') datetime(tend,'ConvertFrom','datenum')]);
    
    ylabel(strcat('\lambda_{MVAB}'),'FontSize', 14)
    % ylim([0 25])
    xlim([mms1_fgm_timedata(1) mms1_fgm_timedata(end)])
    
    set(f7,'XMinorTick','on','YMinorTick','on','XTickLabel',[])
    
    set(f7,'Position',[plot_pos(1), plot_pos(2)-6*plot_pos(4)-0.175, plot_pos(3), plot_pos(4)*1]);
    
    % l1= findobj(gcf, 'Type', 'Legend');
    % l1pos=get(l1(4), 'Position');
    % set(l1(4),'Position',[plot_pos(1)+plot_pos(3)+0.1, plot_pos(2)+3*plot_pos(4), l1pos(3), l1pos(4)]);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Plot Angles%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    f8 = subplot(num_plots,1,8);
    
    plot(mms1_fgm_timedata,MVA_dot_MDD_angle_all)
    box on
    datetick('x','SS.FFF','keeplimits')
    set(gca,'layer','top')
    legend({'Sliding MVA'},'FontSize',10)
    legend('boxoff');
    legend('Location','eastoutside');
    
    ylabel({'Angle'; '[Deg]'},'FontSize', 14)
    
    xlim([mms1_fgm_timedata(1) mms1_fgm_timedata(end)])
    set(f8,'XMinorTick','on','YMinorTick','on')
    
    set(f8,'Position',[plot_pos(1), plot_pos(2)-7*plot_pos(4)-0.2, plot_pos(3), plot_pos(4)*1]);
    
    
    % l1= findobj(gcf, 'Type', 'Legend');
    % l1pos=get(l1(5), 'Position');
    % set(l1(8),'Position',[plot_pos(1)+plot_pos(3)+0.1, plot_pos(2)+4*plot_pos(4), l1pos(3), l1pos(4)]);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%Text on Plot%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % annotation('textbox',[plot_pos(1), plot_pos(2)-8*plot_pos(4)-.18, plot_pos(3), plot_pos(4)*1],...
    %     'String',{date_start(1:10)},...
    %     'VerticalAlignment','Bottom','Edgecolor','none','FontSize', 14);
    
    
    annotation('textbox',[plot_pos(1), plot_pos(2)-8*plot_pos(4)-.18, plot_pos(3), plot_pos(4)*1],...
        'String',strcat({date_start},' - ',{date_end(11:end)}),...
        'VerticalAlignment','Bottom','Edgecolor','none','FontSize', 14);
    
    % annotation('textbox',[plot_pos(1), plot_pos(2)-8*plot_pos(4)-.20, plot_pos(3), plot_pos(4)*1],...
    %     'String',{date_end},...
    %     'VerticalAlignment','Bottom','Edgecolor','none','FontSize', 14);
    
    % annotation('textbox',[plot_pos(1), plot_pos(2)-8*plot_pos(4)-0.180, plot_pos(3), plot_pos(4)*1],...
    %     'String',strcat('(Median MDD)\cdot(Centered MVA) Angle= ', ' ', num2str(median(MVA_dot_MDD_angle_all),'%2.2f'),char(176)),...
    %     'VerticalAlignment','Bottom','Edgecolor','none','FontSize', 14);
    
    annotation('textbox',[plot_pos(1), plot_pos(2)-8*plot_pos(4)-0.20, plot_pos(3), plot_pos(4)*1],...
        'String',strcat('(Average MDD)\cdot(Center MVA) Angle= ', ' ', num2str(mean_MDD_dot_MVA_center,'%2.2f'),char(176)),...
        'VerticalAlignment','Bottom','Edgecolor','none','FontSize', 14);
    
    
    % annotation('textbox',[plot_pos(1), plot_pos(2)-8*plot_pos(4)-0.20, plot_pos(3), plot_pos(4)*1],...
    %     'String',strcat('(Mean MDD)\cdot(Centered MVA) Angle= ', ' ', num2str(MDD_dot_MVA_center,'%2.2f'),char(176)),...
    %     'VerticalAlignment','Bottom','Edgecolor','none','FontSize', 14);
    
    annotation('textbox',[plot_pos(1)+0.55*plot_pos(3), plot_pos(2)-8*plot_pos(4)-0.180, plot_pos(3), plot_pos(4)*1],...
        'String',strcat('Average MDD Normal = ', ' ', num2str(Mean_grad_B_MDD_normal',' [%0.4f, %0.4f, %0.4f]')),...
        'VerticalAlignment','Bottom','Edgecolor','none','FontSize', 14);
    
    % annotation('textbox',[plot_pos(1)+0.55*plot_pos(3), plot_pos(2)-8*plot_pos(4)-0.20, plot_pos(3), plot_pos(4)*1],...
    %     'String',strcat('Mean MDD Normal = ', ' ', num2str(MDD_mean_normal,' [%0.4f, %0.4f, %0.4f]')),...
    %     'VerticalAlignment','Bottom','Edgecolor','none','FontSize', 14);
    
    annotation('textbox',[plot_pos(1)+0.55*plot_pos(3), plot_pos(2)-8*plot_pos(4)-0.20, plot_pos(3), plot_pos(4)*1],...
        'String',strcat('Center MVA Normal = ', ' ', num2str(MVA_min_normal,' [%0.4f, %0.4f, %0.4f]')),...
        'VerticalAlignment','Bottom','Edgecolor','none','FontSize', 14);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if strcmp(plot_pdf,'yes')
        orient(gcf,'landscape')
        
        %plot_name =  strcat('mms','_mdd_',mddl,'_mva_',mval,'_angles2',...
        %    date_start(1:23),'_avg_i',num2str(avg_i),'_fluctmdd',num2str(fluctuate_mdd),'_fluctmva',num2str(fluctuate_mva),'_MVAB_window',num2str(MVAB_window_size),'.pdf');
        plot_name =  strcat('4_BoundaryAnalysis_MVA-MDD_',date_start(1:23),'.pdf');
        print(gcf, '-dpdf',plot_name,'-fillpage');
        
        %movefile(plot_name, '~/Library/Mobile Documents/com~apple~CloudDocs/Research/Analysis')
    end
    
    
    
    figure('Position',[800 700 850 700])
    set(gcf,'color','w');
    co = [0 0 1;
        0 1 0;
        1 0 0;
        0 0 0];
    set(gcf,'defaultAxesColorOrder',co)
    plot_order = 1;
    plot_gap = 1.4;
    num_plots = 6;
    data_type='brst';
    
    %Dates for Event Boundaries
    
    
    
    %%
    %Timing Method
    %Calculate midpoint of data, then number of data points
    all_data = length(mms1_fgm_bdata);
    half_data = floor(all_data/2);
    quarter_data = floor(half_data/2);
    
    %Center index of bdata
    t_0 = 1 + half_data;
    [n_boundary,v_boundary,~,cc12, cc13, cc14, cc12max, cc13max, cc14max, time_boundary,unique_timings] = ...
        timing_method(t_0,all_data,...
        mms1_fgm_timedata,mms1_fgm_bdata,...
        mms2_fgm_bdata,...
        mms3_fgm_bdata,...
        mms4_fgm_bdata,...
        mms1_mec_rdata,...
        mms2_mec_rdata,...
        mms3_mec_rdata,...
        mms4_mec_rdata);
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%Correlation of Data%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %Separate the data into quadrants,
%     %we really only test the correlation in the 2 & 3rd quadriles of mms1 bdata
%     all_data = length(mms1_fgm_timedata); %
%     half_data = floor(all_data/2); %index for half data
%     quarter_data = floor(half_data/2); %index for quarter data
%     
%     %Initialize correlation vectors, xyz, for each spacecraft pairings
%     %We only correlate half the data's worth of points
%     cc12 = zeros(half_data,3);
%     cc13 = zeros(half_data,3);
%     cc14 = zeros(half_data,3);
%     
%    
%     %loop over i=0 to data_points/2 for bxyz for mms234
%     for i=1:half_data %loop for each subsequent data point.
%         %timing window is 1/4 - 3/4 of time range of mms1_btime
%         for j=1:3 %loop for each xyz
%             %Correlation for this point,i, for this component,j, for MMS1 and MMS2
%             temp_matrix = corrcoef(mms1_fgm_bdata(quarter_data:all_data-quarter_data,j),...
%                 mms2_fgm_bdata(i:all_data-half_data+i,j));
%             cc12(i,j) =  temp_matrix(1,2);
%             
%             %Correlation for this point,i, for this component,j, for MMS1 and MMS3
%             temp_matrix = corrcoef(mms1_fgm_bdata(quarter_data:all_data-quarter_data,j),...
%                 mms3_fgm_bdata(i:all_data-half_data+i,j));
%             cc13(i,j) = temp_matrix(1,2);
%             
%             %Correlation for this point,i, for this component,j, for MMS1 and MMS4
%             temp_matrix = corrcoef(mms1_fgm_bdata(quarter_data:all_data-quarter_data,j),...
%                 mms4_fgm_bdata(i:all_data-half_data+i,j));
%             cc14(i,j) = temp_matrix(1,2);
%             
%         end
%     end
%     
%     %find the best index for the starting point of the time range for the most
%     %correlation, highest correlation coefficient "most same" after sliding.
%     %for each component and each spacecraft, find the largest correlation
%     %coefficient, closer to 1 is better. we can only choose the normal
%     %components from one direction, x y or z, so should choose the direction that has the
%     %mean highest correlation coefficient for all spacecraft pairs.
%     
%     %each i is a component x y z
%     for i=1:3
%         [cc12max(i),cc12max_index(i)] = max(cc12(:,i));
%         [cc13max(i),cc13max_index(i)] = max(cc13(:,i));
%         [cc14max(i),cc14max_index(i)] = max(cc14(:,i));
%     end

%     
%     %the index is the starting point and the ending point is i+half_data, this
%     %max index has the largest correlation coefficient, thus we should use this
%     %for the timing method
%     
%     
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %%
%     %%%%%%%%%%%%%%%%%%%Calculation through Timing Method%%%%%%%%%%%%%%%%%%%%%%%
%     %Initialize matrices for all three components' parameters time_boundary =
%     %zeros(4,3);
%     n_boundary = zeros(3,3);
%     v_boundary = zeros(1,3);
%     
%     
%     for i=1:3
%         %subscript is the starting point for btime and analysis, chosen because
%         %of its highest correlation coefficient.
%         index1=quarter_data;
%         index2=cc12max_index(i);
%         index3=cc13max_index(i);
%         index4=cc14max_index(i);
%         
%         %time of crossing of the discontinity by each spacecraft
%         time1=mms1_fgm_timedata(index1);
%         time2=mms1_fgm_timedata(index2);
%         time3=mms1_fgm_timedata(index3);
%         time4=mms1_fgm_timedata(index4);
%         
%         
%         %calculation of normal and speed
%         %Calculate the time difference between spacecrafts
%         time_vector = 86400*[time1-time2;
%             time1-time3;
%             time1-time4;
%             time2-time3;
%             time2-time4;
%             time3-time4]'; %86400 is from datenum to seconds, datenum is in days.
%         
%         %Distance vector of each spacecraft
%         r1=mms1_mec_rdata(index1,1:3);
%         r2=mms2_mec_rdata(index2,1:3);
%         r3=mms3_mec_rdata(index3,1:3);
%         r4=mms4_mec_rdata(index4,1:3);
%         
%         %Calculate the distance difference between spacecrafts
%         r_matrix = [r1-r2;
%             r1-r3;
%             r1-r4;
%             r2-r3;
%             r2-r4;
%             r3-r4];
%         
%         %Calculate volumetric matrix
%         %rearrange positions for easier calculation of R_alphabeta
%         rx = [r1(1),r2(1),r3(1),r4(1)];
%         ry = [r1(2),r2(2),r3(2),r4(2)];
%         rz = [r1(3),r2(3),r3(3),r4(3)];
%         
%         %Calculate the Center of the Tetrahedron
%         rx_0 = mean(rx,2);
%         ry_0 = mean(ry,2);
%         rz_0 = mean(rz,2);
%         
%         %Calculate the relative distance of each probe from the center
%         rx = rx-rx_0;
%         ry = ry-ry_0;
%         rz = rz-rz_0;
%         
%         %Calculate components Volumetric Tensor
%         Rxx = (1/4)*sum(rx.^2,2);
%         Rxy = (1/4)*sum(rx.*ry,2);
%         Rxz = (1/4)*sum(rx.*rz,2);
%         Ryy = (1/4)*sum(ry.^2,2);
%         Ryz = (1/4)*sum(ry.*rz,2);
%         Rzz = (1/4)*sum(rz.^2,2);
%         
%         %initialize volumetric tensor matrix
%         R=zeros(3,3);
%         %Set values of Volumetric tensor
%         R(1,1) = Rxx;
%         R(1,2) = Rxy;
%         R(1,3) = Rxz;
%         R(2,1) = Rxy;
%         R(2,2) = Ryy;
%         R(2,3) = Ryz;
%         R(3,1) = Rxz;
%         R(3,2) = Ryz;
%         R(3,3) = Rzz;
%         
%         
%         m_l= (1/4^2)*(time_vector*r_matrix)/R; %Section 12.1.2 Analysis Methods of Multispacecrafts
%         mag_m_l = norm(m_l);
%         
%         vv = 1/mag_m_l; %speed
%         nn = m_l*vv; %normal
%         
%         %Save variables, columns are for each component
%         time_boundary(:,i) = [datetime(time1,'ConvertFrom','datenum','format','yyyy-MM-dd/HH:mm:ss.SSSS')
%             datetime(time2,'ConvertFrom','datenum','format','yyyy-MM-dd/HH:mm:ss.SSSS')
%             datetime(time3,'ConvertFrom','datenum','format','yyyy-MM-dd/HH:mm:ss.SSSS')
%             datetime(time4,'ConvertFrom','datenum','format','yyyy-MM-dd/HH:mm:ss.SSSS')];
%         
%         n_boundary(:,i) = nn;
%         v_boundary(i) = vv;
%     end
%     
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%
    %%%%%%%%%%%%%%%%%%%%%%%Print Normal And Speed%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %From one CC plot, find the three highest cc for each spacecraft pairing,
    %the closer to 1 these numbers are, the better.
    %the component that has all three numbers closest to 1 is the best for the
    %timing method and whose normal and speeds we should use.
    
    %find the best component to use based on mean correlation coefficients for
    %each spacecraft pairing
    [cc_max_component_value,cc_max_component_index] = ...
        max([mean([cc12max(1),cc13max(1),cc14max(1)]);...
        mean([cc12max(2),cc13max(2),cc14max(2)]);...
        mean([cc12max(3),cc13max(3),cc14max(3)])]);
    
    %display normal and speeds
    MVAB_timing_angle_all = [ambig_angle(MVA_min_normal,n_boundary(:,1)),...
        ambig_angle(MVA_min_normal,n_boundary(:,2)),...
        ambig_angle(MVA_min_normal,n_boundary(:,3))];
    MDD_timing_angle_all = [angle((Mean_grad_B_MDD_normal),(n_boundary(:,1))),...
        angle((Mean_grad_B_MDD_normal),(n_boundary(:,2))),...
        angle((Mean_grad_B_MDD_normal),(n_boundary(:,3)))];
    
    %Use the minimal angle that has 4 unique spacecraft crossings.
    unique_MVAB_timing_angle_all = MVAB_timing_angle_all.*unique_timings;
    unique_MVAB_timing_angle_all(unique_MVAB_timing_angle_all == 0) = NaN;
    
    
    [~,MVAB_timing_angle_index] = min(unique_MVAB_timing_angle_all);
    timing_n = n_boundary(:,MVAB_timing_angle_index)'
    timing_v = v_boundary(MVAB_timing_angle_index)'
    % print_timing_n = sprintf(strcat(num2str(timing_n,'[%.4f %.4f %.4f]')))
    %
    % print_timing_v = sprintf(strcat(num2str(timing_v,'%3.1f'),'km/s'))
    
    [MVA_timing_angle,MVA_n]= ambig_angle(MVA_min_normal,timing_n);
    
    
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    %%Update 2019 August for angle between edges, we choose inward facing and let speed determine
    %%direction
   
    %Find center position of MMS within core of event
    left_InnerEdge_index = find(mms1_fgm_timedata_raw >= datenum(left_InnerEdge,formatIn), 1);
    right_InnerEdge_index = find(mms1_fgm_timedata_raw >= datenum(right_InnerEdge,formatIn), 1);
    
    CenterofCore_TetrahedronPosition = mean([...
        mms1_mec_rdata_interp(left_InnerEdge_index:right_InnerEdge_index,:);...
        mms2_mec_rdata_interp(left_InnerEdge_index:right_InnerEdge_index,:);...
        mms3_mec_rdata_interp(left_InnerEdge_index:right_InnerEdge_index,:);...
        mms4_mec_rdata_interp(left_InnerEdge_index:right_InnerEdge_index,:)...
        ]);
    
    CenterofBoundary_TetrahedronPosition = mean([...
        mms1_mec_rdata;...
        mms2_mec_rdata;...
        mms3_mec_rdata;...
        mms4_mec_rdata...
        ]);
    
    %Inward pointing vector from boundary position
    inwardVector = CenterofCore_TetrahedronPosition - CenterofBoundary_TetrahedronPosition;
    
    %choose the normal that has the smallest angle with the inward pointing vector
    [~,indexforInwardNormal] = min ( [angle(inwardVector,timing_n) angle(inwardVector,-timing_n) ] );
    
    if indexforInwardNormal == 2
        timing_n = -timing_n;
        timing_v = -timing_v;
    end
    

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Plotting%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Plot the total magnetic field for each spacecraft
    plot_fgm_magnetic(date_start,date_end,mms1_fgm_timedata,mms1_fgm_bdata(:,1:1),num_plots,plot_order)
    hold on
    plot_fgm_magnetic(date_start,date_end,mms1_fgm_timedata,mms2_fgm_bdata(:,1:1),num_plots,plot_order)
    plot_fgm_magnetic(date_start,date_end,mms1_fgm_timedata,mms3_fgm_bdata(:,1:1),num_plots,plot_order)
    plot_fgm_magnetic(date_start,date_end,mms1_fgm_timedata,mms4_fgm_bdata(:,1:1),num_plots,plot_order)
    hold off
    ylabel({'B_x';'[nT]'},'FontSize', 14)
    legend('off')
    
    title('MMS1 Timing Method: Correlation Coefficients', 'FontSize', 18, 'FontWeight', 'normal')
    legend({'1', '2', '3','4'},'FontSize',14)
    legend('boxoff')
    legend('Location','eastoutside')
    datetick
    xlim([mms1_fgm_timedata(1) mms1_fgm_timedata(end)])
    set(gca,'XMinorTick','on','YMinorTick','on','linewidth',1.25)
    plot_pos = get(gca,'Position');
    set(gca,'Position',[plot_pos(1), plot_pos(2)-plot_pos(4)*plot_gap*(plot_order-1), plot_pos(3), plot_pos(4)]);
    
    plot_order = plot_order+1;
    
    time_delay = 0.0078*linspace(-quarter_data,quarter_data,half_data);%7.8ms for brst data
    
    %Plot the time delay for x
    subplot(num_plots,1,plot_order)
    plot(time_delay,[cc12(:,1),cc13(:,1),cc14(:,1)])
    legend({'1-2', '1-3', '1-4'},'FontSize',14)
    legend('boxoff')
    legend('Location','eastoutside')
    set(gca,'XTicklabel',[],'XMinorTick','on','YMinorTick','on','linewidth',1.25)
    ylabel({'cc_x'},'FontSize', 14)
    %xlabel({'Time Delay [s]'},'FontSize', 14)
    set(gca,'Position',[plot_pos(1), plot_pos(2)-plot_pos(4)*plot_gap*(plot_order-1), plot_pos(3), plot_pos(4)]);
    
    plot_order = plot_order+1;
    
    %Plot the time delay for y
    subplot(num_plots,1,plot_order)
    plot(time_delay,[cc12(:,2),cc13(:,2),cc14(:,2)])
    legend({'1-2', '1-3', '1-4'},'FontSize',14)
    legend('boxoff')
    legend('Location','eastoutside')
    set(gca,'XTicklabel',[],'XMinorTick','on','YMinorTick','on','linewidth',1.25)
    ylabel({'cc_y'},'FontSize', 14)
    %xlabel({'Time Delay [s]'},'FontSize', 14)
    set(gca,'Position',[plot_pos(1), plot_pos(2)-plot_pos(4)*plot_gap*(plot_order-1)+(((plot_pos(4)*(plot_gap-1.1)*1))), plot_pos(3), plot_pos(4)]);
    
    plot_order = plot_order+1;
    
    %Plot the time delay for z
    subplot(num_plots,1,plot_order)
    plot(time_delay,[cc12(:,3),cc13(:,3),cc14(:,3)])
    legend({'1-2', '1-3', '1-4'},'FontSize',14)
    legend('boxoff')
    legend('Location','eastoutside')
    set(gca,'XMinorTick','on','YMinorTick','on','linewidth',1.25)
    ylabel({'cc_z'},'FontSize', 14)
    xlabel({'Time Delay [s]'},'FontSize', 12)
    set(gca,'Position',[plot_pos(1), plot_pos(2)-plot_pos(4)*plot_gap*(plot_order-1)+(((plot_pos(4)*(plot_gap-1.1)*2))), plot_pos(3), plot_pos(4)]);
    
    plot_order = plot_order+1;
    
    %Plot the max correlation coefficient for each component
    subplot(num_plots,1,plot_order)
    plot([cc12max(1),cc13max(1),cc14max(1)])
    hold on
    plot([cc12max(2),cc13max(2),cc14max(2)])
    plot([cc12max(3),cc13max(3),cc14max(3)])
    hold off
    legend({'x', 'y', 'z'},'FontSize',14)
    legend('boxoff')
    legend('Location','eastoutside')
    colormap('winter');
    set(gca,'XMinorTick','off','YMinorTick','on','linewidth',1.25)
    xticks([1,2,3])
    xticklabels({'1-2','1-3','1-4'})
    xlim([0.9 3.1])
    ylabel({'cc_{max}'},'FontSize', 14)
    xlabel({'Spacecraft Pair'},'FontSize', 12)
    set(gca,'Position',[plot_pos(1), plot_pos(2)-plot_pos(4)*plot_gap*(plot_order-1)+(((plot_pos(4)*(plot_gap-1.25)*3))), plot_pos(3), plot_pos(4)]);
    
    plot_order = plot_order+1;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%Annotations%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %Date
    date_range = {date_start;date_end};
    annotation('textbox',[0, plot_pos(2)-plot_pos(4)*plot_gap*(plot_order-1.25), plot_pos(3), plot_pos(4)],...
        'String',date_range,...
        'VerticalAlignment','Top','Edgecolor','none','FontSize', 12);
    
    %Angles
    annotation('textbox',[0, plot_pos(2)-plot_pos(4)*plot_gap*(plot_order-1.25)-0.05, plot_pos(3), plot_pos(4)],...
        'String',strcat('MVA Angle:', num2str(MVAB_timing_angle_all')),...
        'VerticalAlignment','Top','Edgecolor','none','FontSize', 12);
    annotation('textbox',[0, plot_pos(2)-plot_pos(4)*plot_gap*(plot_order-1.25)-0.125, plot_pos(3), plot_pos(4)],...
        'String',strcat('MDD Angle:', num2str(MDD_timing_angle_all')),...
        'VerticalAlignment','Top','Edgecolor','none','FontSize', 12);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%Anotations for X%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Mean cc
    annotation('textbox',[plot_pos(1)+plot_pos(4), plot_pos(2)-plot_pos(4)*plot_gap*(plot_order-1.25), plot_pos(3), plot_pos(4)],...
        'String',strcat('cc_{x}mean: ', num2str(mean([cc12max(1),cc13max(1),cc14max(1)]))),...
        'VerticalAlignment','Top','Edgecolor','none','FontSize', 12);
    %Timings for each spacecraft
    annotation('textbox',[plot_pos(1)+plot_pos(4), plot_pos(2)-plot_pos(4)*plot_gap*(plot_order-1.25)-0.025, plot_pos(3), plot_pos(4)],...
        'String',strcat('Timings:', datestr(time_boundary(:,1),'HHMMSS.FFF')),...
        'VerticalAlignment','Top','Edgecolor','none','FontSize', 12);
    %Normal
    annotation('textbox',[plot_pos(1)+plot_pos(4), plot_pos(2)-plot_pos(4)*plot_gap*(plot_order-1.25)-0.12, plot_pos(3), plot_pos(4)],...
        'String',strcat('Normal: ',' ', num2str(n_boundary(:,1)','[%0.4f, %0.4f, %0.4f]')),...
        'VerticalAlignment','Top','Edgecolor','none','FontSize', 12);
    %speed
    annotation('textbox',[plot_pos(1)+plot_pos(4), plot_pos(2)-plot_pos(4)*plot_gap*(plot_order-1.25)-0.14, plot_pos(3), plot_pos(4)],...
        'String',strcat('Speed: ', num2str(v_boundary(1),'%3.2f'),'km/s'),...
        'VerticalAlignment','Top','Edgecolor','none','FontSize', 12);
    
    % annotation('textbox',[plot_pos(1)+plot_pos(4), plot_pos(2)-plot_pos(4)*plot_gap*(plot_order-1.25)-0.14, plot_pos(3), plot_pos(4)],...
    %     'String',strcat('Speed: ', num2str(v_boundary(1),'%3.2f'),'km/s'),...
    %     'VerticalAlignment','Top','Edgecolor','none','FontSize', 12);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%Anotations for Y%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Mean cc
    annotation('textbox',[plot_pos(1)+3.5*plot_pos(4), plot_pos(2)-plot_pos(4)*plot_gap*(plot_order-1.25), plot_pos(3), plot_pos(4)],...
        'String',strcat('cc_{y}mean: ', num2str(mean([cc12max(2),cc13max(2),cc14max(2)]))),...
        'VerticalAlignment','Top','Edgecolor','none','FontSize', 12);
    %Timings for each spacecraft
    annotation('textbox',[plot_pos(1)+3.5*plot_pos(4), plot_pos(2)-plot_pos(4)*plot_gap*(plot_order-1.25)-0.025, plot_pos(3), plot_pos(4)],...
        'String',strcat('Timings:', datestr(time_boundary(:,2),'HHMMSS.FFF')),...
        'VerticalAlignment','Top','Edgecolor','none','FontSize', 12);
    %Normal
    annotation('textbox',[plot_pos(1)+3.5*plot_pos(4), plot_pos(2)-plot_pos(4)*plot_gap*(plot_order-1.25)-0.12, plot_pos(3), plot_pos(4)],...
        'String',strcat('Normal: ',' ', num2str(n_boundary(:,2)','[%0.4f, %0.4f, %0.4f]')),...
        'VerticalAlignment','Top','Edgecolor','none','FontSize', 12);
    %speed
    annotation('textbox',[plot_pos(1)+3.5*plot_pos(4), plot_pos(2)-plot_pos(4)*plot_gap*(plot_order-1.25)-0.14, plot_pos(3), plot_pos(4)],...
        'String',strcat('Speed: ', num2str(v_boundary(2),'%3.2f'),'km/s'),...
        'VerticalAlignment','Top','Edgecolor','none','FontSize', 12);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%Anotations for z%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Mean cc
    annotation('textbox',[plot_pos(1)+6*plot_pos(4), plot_pos(2)-plot_pos(4)*plot_gap*(plot_order-1.25), plot_pos(3), plot_pos(4)],...
        'String',strcat('cc_{z}mean: ', num2str(mean([cc12max(3),cc13max(3),cc14max(3)]))),...
        'VerticalAlignment','Top','Edgecolor','none','FontSize', 12);
    %Timings for each spacecraft
    annotation('textbox',[plot_pos(1)+6*plot_pos(4), plot_pos(2)-plot_pos(4)*plot_gap*(plot_order-1.25)-0.025, plot_pos(3), plot_pos(4)],...
        'String',strcat('Timings:', datestr(time_boundary(:,3),'HHMMSS.FFF')),...
        'VerticalAlignment','Top','Edgecolor','none','FontSize', 12);
    %Normal
    annotation('textbox',[plot_pos(1)+6*plot_pos(4), plot_pos(2)-plot_pos(4)*plot_gap*(plot_order-1.25)-0.12, plot_pos(3), plot_pos(4)],...
        'String',strcat('Normal: ',' ', num2str(n_boundary(:,3)','[%0.4f, %0.4f, %0.4f]')),...
        'VerticalAlignment','Top','Edgecolor','none','FontSize', 12);
    %speed
    annotation('textbox',[plot_pos(1)+6*plot_pos(4), plot_pos(2)-plot_pos(4)*plot_gap*(plot_order-1.25)-0.14, plot_pos(3), plot_pos(4)],...
        'String',strcat('Speed: ', num2str(v_boundary(3),'%3.2f'),'km/s'),...
        'VerticalAlignment','Top','Edgecolor','none','FontSize', 12);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
    
    %Save Plots
    if strcmp(plot_pdf,'yes')
        orient(gcf,'landscape')
        plot_name =  strcat('4_BoundaryAnalysis_Timing_',date_start(1:23),'.pdf');
        print(gcf, '-dpdf',plot_name,'-fillpage');
        %movefile(plot_name, '~/Library/Mobile Documents/com~apple~CloudDocs/Research/Analysis')
    end
end

