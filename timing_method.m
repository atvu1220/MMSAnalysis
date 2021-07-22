function [n_boundary,v_boundary,ccmean,cc12, cc13, cc14, cc12max, cc13max, cc14max, time_boundary,unique_timings] = timing_method(t_0,data_points,...
        mms1_fgm_timedata,mms1_fgm_bdata,...
        mms2_fgm_bdata,...
        mms3_fgm_bdata,...
        mms4_fgm_bdata,...
        mms1_mec_rdata,...
        mms2_mec_rdata,...
        mms3_mec_rdata,...
        mms4_mec_rdata) %Data should be interpolated, quadratic or spline with mms1
    %Timing Method Calculation for a specific t_0, center time point.
    %Create quadratures of the data from mms1_btime
    half_data = floor(data_points/2); %Half of the data points 
    quarter_data = floor(half_data/2); %Half of Half of the data points
    
    
    %Initialize correlation vectors, xyz, for each spacecraft pairings
        %We only correlate half the data's worth of points
        cc12 = zeros(half_data,3);
        cc13 = zeros(half_data,3);
        cc14 = zeros(half_data,3);
        
        
        %loop over i=0 to data_points/2,half_data, for bxyz for mms234
        for i=1:half_data %loop for each subsequent data point.
            %timing window is 1/4 - 3/4 of time range of mms1_bdata, and
            %the range to compare with other sc is 0 - 1/2 of time range of
            %mms2,3,4_bdata
            for j=1:3 %loop for each xyz
                %Correlation for this point,i, for this component,j, for MMS1 and MMS2
                temp_matrix = corrcoef(mms1_fgm_bdata(t_0-quarter_data:t_0+quarter_data,j),...
                    mms2_fgm_bdata(i+t_0-half_data:t_0+i,j));
                cc12(i,j) =  temp_matrix(1,2);
                
                %Correlation for this point,i, for this component,j, for MMS1 and MMS3
                temp_matrix = corrcoef(mms1_fgm_bdata(t_0-quarter_data:t_0+quarter_data,j),...
                    mms3_fgm_bdata(i+t_0-half_data:t_0+i,j));
                cc13(i,j) = temp_matrix(1,2);
                
                %Correlation for this point,i, for this component,j, for MMS1 and MMS4
                temp_matrix = corrcoef(mms1_fgm_bdata(t_0-quarter_data:t_0+quarter_data,j),...
                    mms4_fgm_bdata(i+t_0-half_data:t_0+i,j));
                cc14(i,j) = temp_matrix(1,2);
                
            end
        end
        %i is many of data points from t_0 is the best correlation.
        
        %find the best index for the starting point of the time range for the most
        %correlation, highest correlation coefficient "most same" after sliding.
        %for each component and each spacecraft, find the largest correlation
        %coefficient, closer to 1 is better. we can only choose the normal
        %components from one direction, x y or z, so should choose the direction that has the
        %mean highest correlation coefficient for all spacecraft pairs.
        
        %each j is a component x y z
        for j=1:3
            [cc12max(j),cc12max_index(j)] = max(cc12(:,j));
            [cc13max(j),cc13max_index(j)] = max(cc13(:,j));
            [cc14max(j),cc14max_index(j)] = max(cc14(:,j));
        end
        
        for j=1:3 %store max mean correlation coefficients.
        ccmean(j) = mean([cc12max(j),cc13max(j),cc14max(j)]);
        end
        %the index is the starting point and the ending point is i+half_data, this
        %max index has the largest correlation coefficient, thus we should use this
        %for the timing method
        
        %%
        %%%%%%%%%%%%%%%%%%%Calculation through Timing Method%%%%%%%%%%%%%%%%%%%%%%%
        %Initialize matrices for all three components' parameters
        %zeros(4,3);
        n_boundary = zeros(3,3);
        v_boundary = zeros(1,3);
        %time_boundary = zeros(4,3);
        unique_timings = zeros(1,3);
        
        
        for j=1:3
            %subscript is the starting point for btime and analysis, chosen because
            %of its highest correlation coefficient.
            index1=quarter_data;
            index2=cc12max_index(j);
            index3=cc13max_index(j);
            index4=cc14max_index(j);
            
            %time of crossing of the discontinity by each spacecraft
            time1=mms1_fgm_timedata(t_0+index1);
            time2=mms1_fgm_timedata(t_0+index2);
            time3=mms1_fgm_timedata(t_0+index3);
            time4=mms1_fgm_timedata(t_0+index4);
            
            
            %calculation of normal and speed
            %Calculate the time difference between spacecrafts
            time_vector = 86400*[time1-time2;
                time1-time3;
                time1-time4;
                time2-time3;
                time2-time4;
                time3-time4]'; %86400 is from datenum to seconds, datenum is in days.
            
            %Distance vector of each spacecraft
            r1=mms1_mec_rdata(t_0+index1,1:3);
            r2=mms2_mec_rdata(t_0+index2,1:3);
            r3=mms3_mec_rdata(t_0+index3,1:3);
            r4=mms4_mec_rdata(t_0+index4,1:3);
            
            %Calculate the distance difference between spacecrafts
            r_matrix = [r1-r2;
                r1-r3;
                r1-r4;
                r2-r3;
                r2-r4;
                r3-r4];
            
            %Calculate volumetric matrix
            %rearrange positions for easier calculation of R_alphabeta
            rx = [r1(1),r2(1),r3(1),r4(1)];
            ry = [r1(2),r2(2),r3(2),r4(2)];
            rz = [r1(3),r2(3),r3(3),r4(3)];
            
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
            R=zeros(3,3);
            %Set values of Volumetric tensor
            R(1,1) = Rxx;
            R(1,2) = Rxy;
            R(1,3) = Rxz;
            R(2,1) = Rxy;
            R(2,2) = Ryy;
            R(2,3) = Ryz;
            R(3,1) = Rxz;
            R(3,2) = Ryz;
            R(3,3) = Rzz;
            
            
            m_l= (1/4^2)*(time_vector*r_matrix)/R; %Section 12.1.2 Analysis Methods of Multispacecrafts, EQ. 12.13
            mag_m_l = norm(m_l);
            
            v = 1/mag_m_l; %speed
            n = m_l*v; %normal
            
            %Save variables, columns are for each component
            time_boundary(:,j) = [datetime(time1,'ConvertFrom','datenum','format','yyyy-MM-dd/HH:mm:ss.SSSS')
                datetime(time2,'ConvertFrom','datenum','format','yyyy-MM-dd/HH:mm:ss.SSSS')
                datetime(time3,'ConvertFrom','datenum','format','yyyy-MM-dd/HH:mm:ss.SSSS')
                datetime(time4,'ConvertFrom','datenum','format','yyyy-MM-dd/HH:mm:ss.SSSS')];
            
            %saves a matrix of 0 for not unique, and 1 for 4 unique timings
            %for each data point. [1,1,1] means Bx By and Bz are all
            %unique)
            unique_timings(j) = length(unique([time1,time2,time3,time4])) == 4;
            
        
            %for each t_0
            n_boundary(:,j) = n;
            v_boundary(j) = v;
        end
    
    
    
    
end

