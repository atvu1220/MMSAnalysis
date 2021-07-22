function [Timedata,Bxdata,Bydata,Bzdata,Bmagdata,Vxdata,Vydata,Vzdata,Speeddata,ndata,Tdata,Pdata,Betadata,Mdata,MGSMdata] = load_omni(date_start,date_end,n)
    %Load Solar Wind Magnetic Field Data OMNI Data that has been time shifted to nose of BS for all files within directory.
    %date_start and date_end are in 'YYYY-MM-DD' or 'YYYY-MM-DD HH:MM:SS'
    %n is only used when we are retrieving OMNI values for a particular event. n is the number of
    %additional data points before and after the event times to use for averaging.
    if nargin == 2
        n=1; %Default number of data points
    end
    %% Indices for Parameters (5min)
    %     EpochData_index = 1;
    %     Timeshift_index = 11;
    %
    %     Bx_index = 15;
    %     By_index = 16;
    %     Bz_index = 17;
    %
    %
    %     Speed_index = 22;
    %     Vx_index = 23;
    %     Vy_index = 24;
    %     Vz_index = 25;
    %
    %     ionDensity_index = 26;
    %     Temp_index = 27;
    %     Pressure_index = 28;
    %     E_index = 29;
    %     Beta_index = 30;
    %     MachNumber_index = 31;
    
    %% Indices for Parameters (1min)
    %1minute averaged indices
    EpochData_index = 1;
    Timeshift_index = 11;
    
    Bx_index = 16;
    By_index = 17;
    Bz_index = 18;
    
    
    Speed_index = 23;
    Vx_index = 24;
    Vy_index = 25;
    Vz_index = 26;
    
    ionDensity_index = 27;
    Temp_index = 28;
    Pressure_index = 29;
    E_index = 30;
    Beta_index = 31;
    MachNumber_index = 32;
    MGSMachNumber_index = 33;
    
    %% Retrieve Filenames and Calculate the number of months
    data_directory = '~/data/OMNI/'; %Change this to your file directory of OMNI CDF files. All the files should be in the same directory level
    %Search for files in directory
    directory_structure = dir(data_directory);
    
    %Find the exact filename for all files in directory
    all_filenames = {directory_structure(:,1).name};
    
    %Search for files that only contain omni
    data_filenames = all_filenames(contains(all_filenames,'omni'));
    
    
    %Use only the files whose filenames are within the time range given by input arguments
    %date_start and date_end
    index_start = find(contains(data_filenames,strcat(date_start(1:4),date_start(6:7))));
    index_end = find(contains(data_filenames,strcat(date_end(1:4),date_end(6:7))));
    data_filenames = data_filenames(index_start:index_end);
    
    months_of_data = size(data_filenames,2); %Omni data is given in monthly files
    
    clear all_filenames
    clear directory_structure
    
    %% Load First Month of Data
    %Load the first month of Data
    currentDataFilename = strcat(data_directory,char(data_filenames(1)));
    data_info = spdfcdfinfo(currentDataFilename); %Load the CDF file
    data_variables = data_info.Variables; %Get the structure for all variable names
    
    
    EpochData_label = data_variables(EpochData_index,1); %Epoch
    
    BxData_label = data_variables(Bx_index,1); %Bx GSE
    ByData_label = data_variables(By_index,1); %By GSE
    BzData_label = data_variables(Bz_index,1); %Bz GSE
    
    SpeedData_label = data_variables(Speed_index,1); %V GSE
    VxData_label = data_variables(Vx_index,1); %Vx GSE
    VyData_label = data_variables(Vy_index,1); %Vy GSE
    VzData_label = data_variables(Vz_index,1); %Vz GSE
    
    ionDensityData_label = data_variables(ionDensity_index,1); %Ion density
    TempData_label = data_variables(Temp_index,1); %Temp
    PressureData_label = data_variables(Pressure_index,1); %Pressure
    BetaData_label = data_variables(Beta_index,1); %Plasma Beta
    MachNumberData_label = data_variables(MachNumber_index,1); %MachNumber
    MGSMachNumberData_label = data_variables(MGSMachNumber_index,1);
    
    %Load the data with the correct index from the structure.
    Timedata = spdfcdfread(currentDataFilename, 'Variables', EpochData_label, 'ConvertEpochToDatenum', true);
    
    Bxdata = spdfcdfread(currentDataFilename,'Variables', BxData_label);
    Bydata = spdfcdfread(currentDataFilename,'Variables', ByData_label);
    Bzdata = spdfcdfread(currentDataFilename,'Variables', BzData_label);
    
    Speeddata = spdfcdfread(currentDataFilename,'Variables', SpeedData_label);
    Vxdata = spdfcdfread(currentDataFilename,'Variables', VxData_label);
    Vydata = spdfcdfread(currentDataFilename,'Variables', VyData_label);
    Vzdata = spdfcdfread(currentDataFilename,'Variables', VzData_label);
    
    ndata = spdfcdfread(currentDataFilename,'Variables', ionDensityData_label);
    Tdata = spdfcdfread(currentDataFilename,'Variables', TempData_label);
    Pdata = spdfcdfread(currentDataFilename,'Variables', PressureData_label);
    Betadata = spdfcdfread(currentDataFilename,'Variables', BetaData_label);
    Mdata = spdfcdfread(currentDataFilename,'Variables', MachNumberData_label);
    MGSMdata = spdfcdfread(currentDataFilename,'Variables', MGSMachNumberData_label);
    
    %Change for the data type to allow for future calculations
    Timedata = double(Timedata);
    
    Bxdata = double(Bxdata);
    Bydata = double(Bydata);
    Bzdata = double(Bzdata);
    
    Speeddata = double(Speeddata);
    Vxdata = double(Vxdata);
    Vydata = double(Vydata);
    Vzdata = double(Vzdata);
    
    ndata = double(ndata);
    Tdata = double(Tdata);
    Pdata = double(Pdata);
    Betadata = double(Betadata);
    Mdata = double(Mdata);
    MGSMdata = double(MGSMdata);
    
    %% Single date Data Retrieval for Events or Monthly Data Retrieval for Normalization
    %After loading the first file is completed, either continue to load all the data from other months within time interval if
    %the time between date_start and date_end is greater than 1 month OR
    %if we are just looking for the OMNI values for a particular event (couple of minutes), then we
    %just pick out the time intervals of the event within the month's worth of data.
    
    
    if months_of_data==1   
        %% Single Event Data
        %For each single event, find the mean parameters within n data points of event time
        %Crop Data
        formatIn='yyyy-mm-dd HH:MM:SS.FFF';
        
        %Convert the time from the data file from Epoch to our date format
        timeStart = datenum(date_start,formatIn);
        timeEnd = datenum(date_end,formatIn);
        
        %Find the start and end indices of the data for the event times
        start_index = find(Timedata <= timeStart, 1,'last');
        end_index = find(Timedata >= timeEnd, 1);
        
        %Crop the month's worth of Time data to just the times near the Event Times, then convert
        %from DateTime to DateString
        Timedata = Timedata(start_index-n:end_index+n,1);
        Timedata = datestr(Timedata,'yyyy-mm-dd HH:MM:SS');
        
        
        %Also crop the values for other parameters
        Bxdata = (Bxdata(start_index-n:end_index+n,:));
        Bydata = (Bydata(start_index-n:end_index+n,:));
        Bzdata = (Bzdata(start_index-n:end_index+n,:));
        Bmagdata = vecnorm([Bxdata Bydata Bzdata],2,2);
        
        Vxdata = (Vxdata(start_index-n:end_index+n,:));
        Vydata = (Vydata(start_index-n:end_index+n,:));
        Vzdata = (Vzdata(start_index-n:end_index+n,:));
        Speeddata = (Speeddata(start_index-n:end_index+n,:));
        
        ndata = (ndata(start_index-n:end_index+n,:));
        Tdata = (Tdata(start_index-n:end_index+n,:));
        Pdata = (Pdata(start_index-n:end_index+n,:));
        Betadata = (Betadata(start_index-n:end_index+n,:));
        Mdata = (Mdata(start_index-n:end_index+n,:));
        MGSMdata = (MGSMdata(start_index-n:end_index+n,:));
        
        
        
        %Get rid of Extreme data values
        OmniDataMatrix = [Bxdata,Bydata,Bzdata,Bmagdata,Vxdata,Vydata,Vzdata,Speeddata,ndata,Tdata,Pdata,Betadata,Mdata,MGSMdata];
        
        for i=1:size(OmniDataMatrix,1)
            if ((OmniDataMatrix(i,1) >= 1e3) || (OmniDataMatrix(i,2) >= 1e3) || (OmniDataMatrix(i,3)  >= 1e3) || (OmniDataMatrix(i,4)  >= 1e4) || ...
                    (OmniDataMatrix(i,5)  >= 1e4) || (OmniDataMatrix(i,6)  >= 1e4) || (OmniDataMatrix(i,7)  >= 1e4) || (OmniDataMatrix(i,8)  >= 1e4) || ...
                    (OmniDataMatrix(i,9)  >= 900) || (OmniDataMatrix(i,10)  >= 9999990) || (OmniDataMatrix(i,11)  >= 99.000) || (OmniDataMatrix(i,12)  >= 999) || (OmniDataMatrix(i,13)  >= 999))
                OmniDataMatrix(i,:) = NaN;
                Timedata(i,:) = NaN;
            end
        end
        
        %Get the indices where NaN, then clear the entire row
        Timedata(any(isnan(OmniDataMatrix),2), :) = [];
        OmniDataMatrix(any(isnan(OmniDataMatrix),2), :) = [];
        
        
        %Divide up the matrix into 1 variable per vector.
        Bxdata=OmniDataMatrix(:,1);
        Bydata=OmniDataMatrix(:,2);
        Bzdata=OmniDataMatrix(:,3);
        Bmagdata=OmniDataMatrix(:,4);
        Vxdata=OmniDataMatrix(:,5);
        Vydata=OmniDataMatrix(:,6);
        Vzdata=OmniDataMatrix(:,7);
        Speeddata=OmniDataMatrix(:,8);
        ndata=OmniDataMatrix(:,9);
        Tdata=OmniDataMatrix(:,10);
        Pdata=OmniDataMatrix(:,11);
        Betadata=OmniDataMatrix(:,12);
        Mdata=OmniDataMatrix(:,13);
        MGSMdata=OmniDataMatrix(:,14);
        
        %Take the mean of all data points that are n data points from the event start and end times.
        Bxdata = mean(Bxdata);
        Bydata = mean(Bydata);
        Bzdata = mean(Bzdata);
        Bmagdata = mean(Bmagdata);
        
        Vxdata = mean(Vxdata);
        Vydata = mean(Vydata);
        Vzdata = mean(Vzdata);
        Speeddata = mean(Speeddata);
        
        ndata = mean(ndata);
        Tdata = mean(Tdata);
        Pdata = mean(Pdata);
        Betadata = mean(Betadata);
        Mdata = mean(Mdata);
        MGSMdata = mean(MGSMdata);
            
    else
        %% Multiple Months for Normalization
        %Otherwise, we take multiple days for normalization
        %Loop through all months
        for i=2:months_of_data
            currentDataFilename = strcat(data_directory,char(data_filenames(i)));
            currentDataFilename(27:34)
            Timedata_additional = spdfcdfread(currentDataFilename, 'Variables', EpochData_label, 'ConvertEpochToDatenum', true);
            
            Bxdata_additional = spdfcdfread(currentDataFilename,'Variables', BxData_label);
            Bydata_additional = spdfcdfread(currentDataFilename,'Variables', ByData_label);
            Bzdata_additional = spdfcdfread(currentDataFilename,'Variables', BzData_label);
            
            Speeddata_additional = spdfcdfread(currentDataFilename,'Variables', SpeedData_label);
            Vxdata_additional = spdfcdfread(currentDataFilename,'Variables', VxData_label);
            Vydata_additional = spdfcdfread(currentDataFilename,'Variables', VyData_label);
            Vzdata_additional = spdfcdfread(currentDataFilename,'Variables', VzData_label);
            
            ndata_additional = spdfcdfread(currentDataFilename,'Variables', ionDensityData_label);
            Tdata_additional = spdfcdfread(currentDataFilename,'Variables', TempData_label);
            Pdata_additional = spdfcdfread(currentDataFilename,'Variables', PressureData_label);
            Betadata_additional = spdfcdfread(currentDataFilename,'Variables', BetaData_label);
            Mdata_additional = spdfcdfread(currentDataFilename,'Variables', MachNumberData_label);
            MGSMdata_additional = spdfcdfread(currentDataFilename,'Variables', MGSMachNumberData_label);

            
            
            %Convert to Double
            Timedata_additional = double(Timedata_additional);
            
            Bxdata_additional = double(Bxdata_additional);
            Bydata_additional = double(Bydata_additional);
            Bzdata_additional = double(Bzdata_additional);
            
            Speeddata_additional = double(Speeddata_additional);
            Vxdata_additional = double(Vxdata_additional);
            Vydata_additional = double(Vydata_additional);
            Vzdata_additional = double(Vzdata_additional);
            
            ndata_additional = double(ndata_additional);
            Tdata_additional = double(Tdata_additional);
            Pdata_additional = double(Pdata_additional);
            Betadata_additional = double(Betadata_additional);
            Mdata_additional = double(Mdata_additional);
            MGSMdata_additional = double(MGSMdata_additional);
            
            %append this file's data to the main data array
            Timedata = [Timedata;Timedata_additional];
            
            Bxdata = [Bxdata;Bxdata_additional];
            Bydata = [Bydata;Bydata_additional];
            Bzdata = [Bzdata;Bzdata_additional];
            
            Vxdata = [Vxdata;Vxdata_additional];
            Vydata = [Vydata;Vydata_additional];
            Vzdata = [Vzdata;Vzdata_additional];
            Speeddata = [Speeddata;Speeddata_additional];
            
            ndata = [ndata;ndata_additional];
            Tdata = [Tdata;Tdata_additional];
            Pdata = [Pdata;Pdata_additional];
            Betadata = [Betadata;Betadata_additional];
            Mdata = [Mdata;Mdata_additional];
            MGSMdata = [MGSMdata;MGSMdata_additional];
            
            clear Timedata_additional;
            clear Bxdata_additional;
            clear Bydata_additional;
            clear Bzdata_additional;
            clear Vxdata_additional;
            clear Vydata_additional;
            clear Vzdata_additional;
            clear Speeddata_additional;
            clear ndata_additional;
            clear Tdata_additional;
            clear Pdata_additional;
            clear Betadata_additional;
            clear Mdata_additional;
            clear MGSMdata_additional;
            
        end
        Bmagdata = vecnorm([Bxdata Bydata Bzdata],2,2);
        Timedata = datestr(Timedata,'yyyy-mm-dd HH:MM:SS');
        
        %Get rid of Extreme data values
        OmniDataMatrix = [Bxdata,Bydata,Bzdata,Bmagdata,Vxdata,Vydata,Vzdata,Speeddata,ndata,Tdata,Pdata,Betadata,Mdata,MGSMdata];
        
        for i=1:length(OmniDataMatrix)
            if ((OmniDataMatrix(i,1) >= 1e3) || (OmniDataMatrix(i,2) >= 1e3) || (OmniDataMatrix(i,3)  >= 1e3) || (OmniDataMatrix(i,4)  >= 1e4) || ...
                    (OmniDataMatrix(i,5)  >= 1e4) || (OmniDataMatrix(i,6)  >= 1e4) || (OmniDataMatrix(i,7)  >= 1e4) || (OmniDataMatrix(i,8)  >= 1e4) || ...
                    (OmniDataMatrix(i,9)  >= 900) || (OmniDataMatrix(i,10)  >= 9999990) || (OmniDataMatrix(i,11)  >= 99.000) || (OmniDataMatrix(i,12)  >= 999) || (OmniDataMatrix(i,13)  >= 999))
                OmniDataMatrix(i,:) = NaN;
                Timedata(i,:) = NaN;
            end
        end
        
        %Get the indices where NaN, then clear the entire row
        Timedata(any(isnan(OmniDataMatrix),2), :) = [];
        OmniDataMatrix(any(isnan(OmniDataMatrix),2), :) = [];
        
        
        Bxdata=OmniDataMatrix(:,1);
        Bydata=OmniDataMatrix(:,2);
        Bzdata=OmniDataMatrix(:,3);
        Bmagdata=OmniDataMatrix(:,4);
        Vxdata=OmniDataMatrix(:,5);
        Vydata=OmniDataMatrix(:,6);
        Vzdata=OmniDataMatrix(:,7);
        Speeddata=OmniDataMatrix(:,8);
        ndata=OmniDataMatrix(:,9);
        Tdata=OmniDataMatrix(:,10);
        Pdata=OmniDataMatrix(:,11);
        Betadata=OmniDataMatrix(:,12);
        Mdata=OmniDataMatrix(:,13);
        MGSMdata=OmniDataMatrix(:,13);
        
        %For saving the data to a MATLAB File.
                OMNI_Bmag = Bmagdata;
                OMNI_Vmag = Speeddata;
                OMNI_n = ndata;
                OMNI_T = Tdata;
                OMNI_P = Pdata;
                OMNI_Beta = Betadata;
                OMNI_Mach = Mdata;
                OMNI_MGSMach = MGSMdata;
                cd  '/Users/andrewvu/Library/Mobile Documents/com~apple~CloudDocs/Research/MATLAB Analysis'
                save('OMNI_2017_2019',...
                    'OMNI_Bmag','OMNI_Vmag','OMNI_n','OMNI_T','OMNI_P','OMNI_Beta','OMNI_Mach','OMNI_MGSMach')
        
    end
end