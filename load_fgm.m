function [Btimedata,Bdata,rtimedata,rdata] = load_fgm(date_start,date_end,probe,data_type,Btype)
    
    if isduration(date_start) == 1
        formatOut='yyyy-mm-dd HH:MM:SS.FFF';
        date_start = datestr(date_start,formatOut);
        date_end = datestr(date_end,formatOut);
    end
    
    %Splice date details
    year = date_start(1:4);
    month = date_start(6:7);
    day = date_start(9:10);
    hour = date_start(12:13);
    minute = date_start(15:16);
    
    if nargin == 5
        switch Btype
            case 'GSE'
                Btype_index = 2;
            case 'GSM'
                Btype_index = 3;
            case 'DMPA'
                Btype_index = 4;
            case 'BCS'
                Btype_index = 5;
        end
        
    else
        Btype_index = 2; %GSE
    end
    
    
    
    %Create Directory String
    data_directory = '~/data/mms/';
    probe_instrument_datatype_directory = strcat('mms',num2str(probe),'/fgm/',data_type,'/l2/');
    
    
    
    %Directory differences for brst and srvy
    if data_type == 'brst'
        date_directory = strcat(year,'/',month,'/',day,'/');
        buffer = 128*3; %3s buffer for brst
    else
        date_directory = strcat(year,'/',month,'/');
        buffer = 2*16*60*15; %15 minute buffer for srvy
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
        
        %B and R labels, from spdfcdfinfo
        Bdata_label = data_variables(Btype_index,1);
        rdata_label = data_variables(8,1);
        
        
        %Load timedata as datenum
        Btimedata = spdfcdfread(file_directory_name, 'Variables', 'Epoch', 'ConvertEpochToDatenum', true);
        rtimedata = spdfcdfread(file_directory_name, 'Variables', 'Epoch_state', 'ConvertEpochToDatenum', true);
        %Themis use,  Btimedata = spdftt2000todatenum(Btimedata)
      
        %load FGM data
        Bdata = spdfcdfread(file_directory_name,'Variables', Bdata_label);
        rdata = spdfcdfread(file_directory_name,'Variables', rdata_label);
        
        %Change from single to double datatype for math
        Btimedata = double(Btimedata);
        rtimedata = double(rtimedata);
        Bdata = double(Bdata);
        rdata = double(rdata);
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
            
            %B and R labels, from spdfcdfinfo
            Bdata_label = data_variables(Btype_index,1);
            rdata_label = data_variables(8,1);
            
            
            %Load timedata as datenum
            Btimedata_additional = spdfcdfread(file_directory_name, 'Variables', 'Epoch', 'ConvertEpochToDatenum', true);
            rtimedata_additional = spdfcdfread(file_directory_name, 'Variables', 'Epoch_state', 'ConvertEpochToDatenum', true);
            
            %load FGM data
            Bdata_additional = spdfcdfread(file_directory_name,'Variables', Bdata_label);
            rdata_additional = spdfcdfread(file_directory_name,'Variables', rdata_label);
            
            %convert to double type for more precision
            Btimedata_additional = double(Btimedata_additional);
            rtimedata_additional = double(rtimedata_additional);
            
            Bdata_additional = double(Bdata_additional);
            rdata_additional = double(rdata_additional);
            
            
            %append this file's data to the main data array
            Btimedata = [Btimedata; Btimedata_additional];
            rtimedata = [rtimedata; rtimedata_additional];
            
            Bdata = [Bdata; Bdata_additional];
            rdata = [rdata; rdata_additional];
            
            clear Btimedata_additional;
            clear rtimedata_additional
            clear Bdata_additional;
            clear rdata_additional;
        end
        
    end
    
    [Btimedata,sortedIndices] = sort(Btimedata,'ascend');
%    rtimedata = rtimedata(sortedIndices);
    Bdata = [Bdata(sortedIndices,1),Bdata(sortedIndices,2),Bdata(sortedIndices,3),Bdata(sortedIndices,4)];

    %Splicing to reduce the large file size
    %Convert event times to datenum
    formatIn='yyyy-mm-dd HH:MM:SS.FFF';
    tstart = datenum(date_start,formatIn);
    tend = datenum(date_end,formatIn);
    
    %Splice Bdata
    % %Find the start and end limits of the event in the data, then add an
    % extra 64 data points = 1/2second for interpolation methods.
    start_index = find(Btimedata >= tstart, 1)-buffer;
    end_index = find(Btimedata >= tend, 1)+buffer;
    
    %Restrict it to the data if indices are outside of bounds
    if start_index < 0
        start_index = 1;
    end
    if end_index > length(Btimedata)
        end_index = length(Btimedata);
    end
    
    
    %Convert the datetime to date String, and then crop to our event timeframe
    Btimedata = Btimedata(start_index:end_index,1);
    Bdata = Bdata(start_index:end_index,:);
    
    %Splice rdata
    % %Find the start and end limits of the event in the data
    start_index = find(rtimedata >= tstart, 1);
    end_index = find(rtimedata >= tend, 1);
    
    %Restrict it to the data if indices are outside of bounds
    if start_index < 0
        start_index = 1;
    end
    if end_index > length(rtimedata)
        end_index = length(rtimedata);
    end
    
    %Convert the datetime to date String, and then crop to our event timeframe
    rtimedata = rtimedata(start_index:end_index,1);
    rdata = rdata(start_index:end_index,:);
    
    
    
    
    %if there are more data files, retrieve them and combine them into one big
    %array
    
    % filter = contains(all_filenames, strcat(year,month,day));
    
    
    
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
    %     Bdata_label = data_variables(2,1);
    %     rdata_label = data_variables(8,1);
    %
    %
    %     %Load timedata as datenum
    %     Btimedata = spdfcdfread(file_directory_name, 'Variables', 'Epoch', 'ConvertEpochToDatenum', true);
    %     rtimedata = spdfcdfread(file_directory_name, 'Variables', 'Epoch_state', 'ConvertEpochToDatenum', true);
    %
    %     %load FGM data
    %     Bdata = spdfcdfread(file_directory_name,'Variables', Bdata_label);
    %     rdata = spdfcdfread(file_directory_name,'Variables', rdata_label);
    %
    %     %Change from single to double datatype for math
    %     Btimedata = double(Btimedata);
    %     rtimedata = double(rtimedata);
    %     Bdata = double(Bdata);
    %     rdata = double(rdata);
    
    
    
    
    
    
end

%1Epoch
%2gse brst
%3gsm brst
%7Epochr
%8r_rse
%9r_gsm