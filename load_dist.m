function [tdata,phidata,thetadata,energydata,distdata,start_index,end_index] = load_dist(date_start,date_end,probe,data_type,specie)
    
    
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
    
    %Create Directory String
    data_directory = '~/data/mms/';
    if specie == 'e'
        probe_instrument_datatype_directory = strcat('mms',num2str(probe),'/fpi/',data_type,'/l2/des-dist/');
    elseif specie == 'i'
        probe_instrument_datatype_directory = strcat('mms',num2str(probe),'/fpi/',data_type,'/l2/dis-dist/');
        
    end
    
    %Directory differences for brst and srvy
    if strcmp(data_type,'brst')
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
        i=3;
        while i > 2
            file_name = directory_structure(i,1).name; % . .. first-file second-file;; 3 is the first file.
            file_directory = directory_structure(i,1).folder;
            file_directory_name = strcat(file_directory,'/',file_name);
            tdata = spdfcdfread(file_directory_name, 'Variables', 'Epoch', 'ConvertEpochToDatenum', true);
            
            
            %Splicing to reduce the large file size
            %Convert event times to datenum
            formatIn='yyyy-mm-dd HH:MM:SS.FFF';
            tstart = datenum(date_start,formatIn);
            tend = datenum(date_end,formatIn);
            
            
            %Find start and end indices.
            start_index = find(tdata >= tstart, 1);
            end_index = find(tdata >= tend, 1);
            
            if isempty(start_index) || isempty(end_index)
                i=i+1;
            
            elseif start_index > 0 || end_index < length(tdata)
%                 file_index = i;
                break;
            elseif i==size(filter,2)
                disp('No data in all of folder matches this date')
                break;
            end
            
        end
        
%         %First dataset
%         file_name = directory_structure(file_index,1).name; % . .. first-file second-file;; 3 is the first file.
%         file_directory = directory_structure(file_index,1).folder;
%         
%         %Do a Preliminary info on the dataset for the data set labels
%         file_directory_name = strcat(file_directory,'/',file_name);

        data_info = spdfcdfinfo(file_directory_name);
        data_variables = data_info.Variables;
        
        %Load data labels for each specie
        if specie == 'e'
            phidata_label = data_variables(9,1);
            thetadata_label = data_variables(17,1);
            energydata_label = data_variables(19,1);
            distdata_label = data_variables(11,1);
        elseif specie == 'i'
            phidata_label = data_variables(9,1);
            thetadata_label = data_variables(17,1);
            energydata_label = data_variables(19,1);
            distdata_label = data_variables(11,1);
        end
        
        %Load timedata as datenum
%         tdata = spdfcdfread(file_directory_name, 'Variables', 'Epoch', 'ConvertEpochToDatenum', true);
        
        %load FPI data
        phidata = spdfcdfread(file_directory_name,'Variables', phidata_label);
        thetadata = spdfcdfread(file_directory_name,'Variables', thetadata_label);
        energydata = spdfcdfread(file_directory_name,'Variables', energydata_label);
        distdata = spdfcdfread(file_directory_name,'Variables', distdata_label);
        
        
        %Change from single to double datatype for math
        phidata = double(phidata);
        thetadata = double(thetadata);
        
        energydata = double(energydata);
        distdata = double(distdata);
        
    end
    %
    %     if number_of_filenames >= 4
    %         for i=4:number_of_filenames
    %             %Create additional filename-directory strings
    %             file_name = directory_structure(i,1).name; % . .. first-file second-file;; 3 is the first file.
    %             file_directory = directory_structure(i,1).folder;
    %
    %             %Do a Preliminary info on the dataset for the data set labels
    %             file_directory_name = strcat(file_directory,'/',file_name);
    %             data_info = spdfcdfinfo(file_directory_name);
    %             data_variables = data_info.Variables;
    %
    %
    %
    %         %Load data labels for each specie
    %         if specie == 'e'
    %             phidata_label = data_variables(9,1);
    %             thetadata_label = data_variables(17,1);
    %             energydata_label = data_variables(19,1);
    %             distdata_label = data_variables(11,1);
    %         elseif specie == 'i'
    %             phidata_label = data_variables(9,1);
    %             thetadata_label = data_variables(17,1);
    %             energydata_label = data_variables(19,1);
    %             distdata_label = data_variables(11,1);
    %         end
    %
    %         %Load timedata as datenum
    %         tdata_additional = spdfcdfread(file_directory_name, 'Variables', 'Epoch', 'ConvertEpochToDatenum', true);
    %
    %         %load FPI data
    %         phidata_additional = spdfcdfread(file_directory_name,'Variables', phidata_label);
    %         thetadata_additional = spdfcdfread(file_directory_name,'Variables', thetadata_label);
    %         energydata_additional = spdfcdfread(file_directory_name,'Variables', energydata_label);
    %         distdata_additional = spdfcdfread(file_directory_name,'Variables', distdata_label);
    %
    %         %convert to double type for more precision
    %         phidata_additional = double(phidata_additional);
    %         thetadata_additional = double(thetadata_additional);
    %         energydata_additional = double(energydata_additional);
    %         distdata_additional = double(distdata_additional);
    %
    %         %append this file's data to the main data array
    %         tdata =[tdata;tdata_additional];
    %         phidata = [phidata; phidata_additional];
    %         thetadata = [thetadata; thetadata_additional];
    %         energydata = [energydata; energydata_additional];
    %         distdata = [distdata; distdata_additional];
    %
    %             clear tdata_additional
    %             clear phidata_additional;
    %             clear thetadata_additional;
    %             clear energydata_additional;
    %             clear distdata_additional;
    %
    %
    %
    %         end
    
% end


%Convert the datetime to date String, and then crop to our event timeframe
tdata = tdata(start_index:end_index);
phidata = phidata(start_index:end_index,:);
energydata = energydata(start_index:end_index,:);
distdata = distdata(:,:,:,start_index:end_index);

%Change Angles
phidata = mod(phidata + 180,360); %ensure between [0 and 360)
thetadata = thetadata - 90; %range from -90 to 90.



end
