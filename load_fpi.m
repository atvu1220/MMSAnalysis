function [tdata,ndata,vdata,tparadata,tperpdata,energydata,espectdata,pressdata] = load_fpi(date_start,date_end,probe,data_type,specie)
    %Splice date details
    year = date_start(1:4);
    month = date_start(6:7);
    day = date_start(9:10);
    hour = date_start(12:13);
    minute = date_start(15:16);
    
    %Create Directory String
    data_directory = '~/data/mms/';
    if specie == 'e'
        probe_instrument_datatype_directory = strcat('mms',num2str(probe),'/fpi/',data_type,'/l2/des-moms/');
    elseif specie == 'i'
        probe_instrument_datatype_directory = strcat('mms',num2str(probe),'/fpi/',data_type,'/l2/dis-moms/');
        
    end
    
    %Directory differences for brst and srvy
    if data_type == 'brst'
        date_directory = strcat(year,'/',month,'/',day,'/');
        buffer = 128*3; %3s buffer for brst
    else
        date_directory = strcat(year,'/',month,'/');
        buffer = 16*60*5; %5m buffer for srvy.
    end
    %Complete Directory String
    total_directory = strcat(data_directory,probe_instrument_datatype_directory,date_directory);
    
    %Search for files in directory
    directory_structure = dir(total_directory);
    
    %Find the exact filename if there are multiple files in directory
    all_filenames = {directory_structure(:,1).name};
    
    
    
    
    
    number_of_filenames = size(all_filenames,2); %subtract 2 for .. and .,
    
    if number_of_filenames == 2
        disp('No data for entered event times')
        quit
    elseif number_of_filenames > 2
        %Find the first mms data file
        filter = contains(all_filenames, strcat(year,month,day));
        file_index = find(filter==1,1);
        
        %First dataset
        file_name = directory_structure(file_index,1).name; % . .. first-file second-file;; 3 is the first file.
        file_directory = directory_structure(file_index,1).folder;
        
        %Do a Preliminary info on the dataset for the data set labels
        file_directory_name = strcat(file_directory,'/',file_name);
        data_info = spdfcdfinfo(file_directory_name);
        data_variables = data_info.Variables;
        
        %Load data labels for each specie
        if specie == 'e'
            ndata_label = data_variables(23,1);
            vdata_label = data_variables(29,1);
            
            pressdata_label = data_variables(33,1);
            
            tparadata_label = data_variables(45,1);
            tperpdata_label = data_variables(46,1);
            
            edata_label = data_variables(43,1);%gauges
            espectdata_label = data_variables(22,1);
        elseif specie == 'i'
            ndata_label = data_variables(19,1);
            vdata_label = data_variables(25,1);
            
            pressdata_label = data_variables(29,1);
            
            tparadata_label = data_variables(40,1);
            tperpdata_label = data_variables(41,1);
            
            edata_label = data_variables(38,1);%gauges
            espectdata_label = data_variables(16,1);
        end
        
        %Load timedata as datenum
        tdata = spdfcdfread(file_directory_name, 'Variables', 'Epoch', 'ConvertEpochToDatenum', true);
        
        %load FPI data
        ndata = spdfcdfread(file_directory_name,'Variables', ndata_label);
        vdata = spdfcdfread(file_directory_name,'Variables', vdata_label);
        
        pressdata = spdfcdfread(file_directory_name,'Variables', pressdata_label);
        
        tparadata = spdfcdfread(file_directory_name,'Variables', tparadata_label);
        tperpdata = spdfcdfread(file_directory_name,'Variables', tperpdata_label);
        
        energydata = spdfcdfread(file_directory_name,'Variables', edata_label);
        espectdata = spdfcdfread(file_directory_name,'Variables', espectdata_label);
        
        %Change from single to double datatype for math
        tdata = double(tdata);
        ndata = double(ndata);
        
        vdata = double(vdata);
        pressdata = double(pressdata);
        tparadata = double(tparadata);
        tperpdata = double(tperpdata);
        
        energydata = double(energydata);
        espectdata = double(espectdata);
        
    end
    
    if number_of_filenames >= 4
        for i=4:number_of_filenames
            %Create additional filename-directory strings
            file_name = directory_structure(i,1).name; % . .. first-file second-file;; 3 is the first file.
            file_directory = directory_structure(i,1).folder;
            
            %Do a Preliminary info on the dataset for the data set labels
            file_directory_name = strcat(file_directory,'/',file_name);
            data_info = spdfcdfinfo(file_directory_name);
            data_variables = data_info.Variables;
            
            
            
            %Load data labels for each specie
            if specie == 'e'
                ndata_label = data_variables(23,1);
                vdata_label = data_variables(29,1);
                
                pressdata_label = data_variables(33,1);
                
                tparadata_label = data_variables(45,1);
                tperpdata_label = data_variables(46,1);
                
                edata_label = data_variables(43,1);%gauges
                espectdata_label = data_variables(22,1);
            elseif specie == 'i'
                ndata_label = data_variables(19,1);
                vdata_label = data_variables(25,1);
                
                pressdata_label = data_variables(29,1);
                
                tparadata_label = data_variables(40,1);
                tperpdata_label = data_variables(41,1);
                
                edata_label = data_variables(38,1);%gauges
                espectdata_label = data_variables(16,1);
            end
            
            %Load timedata as datenum
            tdata_additional = spdfcdfread(file_directory_name, 'Variables', 'Epoch', 'ConvertEpochToDatenum', true);
            
            %load FPI data
            ndata_additional = spdfcdfread(file_directory_name,'Variables', ndata_label);
            vdata_additional = spdfcdfread(file_directory_name,'Variables', vdata_label);
            
            pressdata_additional = spdfcdfread(file_directory_name,'Variables', pressdata_label);
            
            tparadata_additional = spdfcdfread(file_directory_name,'Variables', tparadata_label);
            tperpdata_additional = spdfcdfread(file_directory_name,'Variables', tperpdata_label);
            
            energydata_additional = spdfcdfread(file_directory_name,'Variables', edata_label);
            espectdata_additional = spdfcdfread(file_directory_name,'Variables', espectdata_label);
            
            %convert to double type for more precision
            tdata_additional = double(tdata_additional);
            
            ndata_additional = double(ndata_additional);
            vdata_additional = double(vdata_additional);
            
            pressdata_additional = double(pressdata_additional);
            
            tparadata_additional = double(tparadata_additional);
            tperpdata_additional = double(tperpdata_additional);
            
            energydata_additional = double(energydata_additional);
            espectdata_additional = double(espectdata_additional);
            
            
            %append this file's data to the main data array
            tdata = [tdata; tdata_additional];
            
            ndata = [ndata; ndata_additional];
            vdata = [vdata; vdata_additional];
            
            %pressdata = [pressdata; pressdata_additional];
            pressdata = cat(3,pressdata,pressdata_additional);
            
            tparadata = [tparadata; tparadata_additional];
            tperpdata = [tperpdata; tperpdata_additional];
            
            energydata = [energydata; energydata_additional];
            espectdata = [espectdata; espectdata_additional];
            
            clear tdata_additional;
            clear ndata_additional;
            clear vdata_additional;
            clear pressdata_additional;
            clear tparadata_additional;
            clear tperpdata_additional;
            clear energydata_additional;
            clear espectdata_additional;
            
            
            
        end
        
    end
    
    
    %Splicing to reduce the large file size
    %Convert event times to datenum
    formatIn='yyyy-mm-dd HH:MM:SS.FFF';
    tstart = datenum(date_start,formatIn);
    tend = datenum(date_end,formatIn);
    
    %Splice data
    % %Find the start and end limits of the event in the data, then add an
    % extra 64 data points = 1/2second for interpolation methods.
    start_index = find(tdata >= tstart, 1)-buffer;
    end_index = find(tdata >= tend, 1)+buffer;
    
    %Restrict it to the data if indices are outside of bounds
    if start_index < 0
        start_index = 1;
    end
    if end_index > length(tdata)
        end_index = length(tdata);
    end
    
    %Convert the datetime to date String, and then crop to our event timeframe
    tdata = tdata(start_index:end_index,1);
    
    ndata = ndata(start_index:end_index,:);
    vdata = vdata(start_index:end_index,:);
    
    pressdata = pressdata(:,:,start_index:end_index);
    %pressdata = reshape(pressdata,length(pressdata),3,3);
    
    tparadata = tparadata(start_index:end_index,:);
    tperpdata = tperpdata(start_index:end_index,:);
    
    energydata = energydata(start_index:end_index,:);
    espectdata = espectdata(start_index:end_index,:);
    
end

%     filter = contains(all_filenames, strcat(year,month,day));
%
%     %narrow it down from day,hour,minute for the exact file
%     if sum(filter) == 1
%         file_index = find(filter==1);
%     elseif sum(filter) > 1
%         filter = contains(all_filenames, strcat(year,month,day,hour));
%         if sum(filter) == 1
%             file_index = find(filter == 1);
%         elseif sum(filter) > 1
%             %             filter = contains(all_filenames, strcat(year,month,day,hour,minute));
%             %             if sum(filter) == 1
%             %                 file_index = find(filter ==1);
%
%             for i=1:length(filter)
%                 if filter(i) == 1
%                     file_name = directory_structure(i,1).name;
%                     file_directory = directory_structure(i,1).folder;
%                     file_directory_name = strcat(file_directory,'/',file_name);
%                     temp_timedata = spdfcdfread(file_directory_name, 'Variables', 'Epoch', 'ConvertEpochToDatenum', true);
%                     if temp_timedata(1) < datenum(date_start) && temp_timedata(end) > datenum(date_start)
%                         file_index = i;
%                         break;
%                         %%%to be continued with 2 filenames combined.
%                     end
%                 else
%                     continue;
%
%                 end
%
%
%             end
%
%
%
%
%
%             %             end
%         end
%     end

%
%     %Now that we have the right file, create the directory listing
%     file_name = directory_structure(file_index,1).name;
%     file_directory = directory_structure(file_index,1).folder;
%
%     %Do a Preliminary info on the dataset for the data set labels
%     file_directory_name = strcat(file_directory,'/',file_name);
%     data_info = spdfcdfinfo(file_directory_name);
%     data_variables = data_info.Variables;
%
%     %B and R labels, from spdfcdfinfo
%     %     tdata_label = data_variables(1,1);
%
%
%
%     if specie == 'e'
%         ndata_label = data_variables(23,1);
%         vdata_label = data_variables(29,1);
%
%         tparadata_label = data_variables(45,1);
%         tperpdata_label = data_variables(46,1);
%
%         edata_label = data_variables(43,1);%gauges
%         espectdata_label = data_variables(22,1);
%     elseif specie == 'i'
%         ndata_label = data_variables(19,1);
%         vdata_label = data_variables(25,1);
%
%         tparadata_label = data_variables(40,1);
%         tperpdata_label = data_variables(41,1);
%
%         edata_label = data_variables(38,1);%gauges
%         espectdata_label = data_variables(16,1);
%     end
%
%
%
%
%
%     %Load timedata as datenum
%     tdata = spdfcdfread(file_directory_name, 'Variables', 'Epoch', 'ConvertEpochToDatenum', true);
%
%     %load FGM data
%     ndata = spdfcdfread(file_directory_name,'Variables', ndata_label);
%     vdata = spdfcdfread(file_directory_name,'Variables', vdata_label);
%
%     tparadata = spdfcdfread(file_directory_name,'Variables', tparadata_label);
%     tperpdata = spdfcdfread(file_directory_name,'Variables', tperpdata_label);
%
%     energydata = spdfcdfread(file_directory_name,'Variables', edata_label);
%     espectdata = spdfcdfread(file_directory_name,'Variables', espectdata_label);
%
%     %Change from single to double datatype for math
%     tdata = double(tdata);
%     ndata = double(ndata);
%
%
%     vdata = double(vdata);
%     tparadata = double(tparadata);
%     tperpdata = double(tperpdata);
%
%     energydata = double(energydata);
%     espectdata = double(espectdata);
%

%1Epoch
%2gse brst
%3gsm brst
%7Epochr
%8r_rse
%9r_gsm