%Given an event number, load the data, plot the density, choose ambient densities, choose duration of each boundary
clear
close all
tic

%Event Numbers
Typical7=  [2;3;5;6;7;8;9;11;12;14;15;16;17;18;20;21;22;24;25;28;29;30;31;32;34;35;36;37;38;39;41;42;43;47;48;50;51;52;54;55;56;57;58;60;61;65;66;67;68;76;77;78;79;80;81;82;83;84;85;87;93;94;97;98;99;100;104;105;106;111;116;117;121;122;123;124;127;130;132;133;135;136;138;142;143;152;153;154;157;158;159;160;161;162;164;165;166;167;168;169;170;...
    179;183;184;185;187;188;189;191;199;200;201;202;203;204;205]';

%Create matrix for all values
%Event Number, downstream Ambient value, downstream value, downstream compression ratio,
%upstream ambient value, upstream value, upstream compression ratio,
%comparison ratio


Event=[];
if ~isempty(Event)
    Typical7 = Event;
end
BoundaryCompressionRatios = zeros(length(Typical7),12);
for i=1:length(Typical7)
    Event_number = Typical7(i);
    %Load the data for the event
    
    [Event_Type, Substructure,threshold_std,data_type,...
        event_start,event_end,...
        left_OuterEdge, left_InnerEdge,...
        right_InnerEdge, right_OuterEdge,...
        leading_leftmost_date, leading_rightmost_date,...
        trailing_leftmost_date, trailing_rightmost_date,...
        leading_start, leading_end,...
        trailing_start, trailing_end,...
        downstream_geometry,upstream_geometry] = get_eventTimes(Event_number);
    
    
    eventDataFileName = strcat('MMS1_Data_EventNumber_',num2str(Event_number),'.mat');
    eventDataDirectory = '/Users/andrewvu/Library/Mobile Documents/com~apple~CloudDocs/Research/Event Data/';
    
    eventDataDirectoryFileName = strcat(eventDataDirectory,eventDataFileName);
    load(eventDataDirectoryFileName) %#ok<LOAD>
    
    
    %Plot the Density
    if length(Typical7) ==1
        figure('Position',[1 1 1050 650])
        co = [0 0 1;
            0 1 0;
            1 0 0;
            0 0 0];
        set(gcf,'defaultAxesColorOrder',co)
        set(gcf,'color','w');
        plot_fpi_number(event_start,event_end,fpi_e_timedata,fpi_e_ndata,1,1,'e');datetick('keeplimits')
    end
    %Check if there's two boundaries, if not skip.
    
    %Choose ambient density
    [preMean,postMean] = calculate_ambientMean(fpi_e_timedata,fpi_e_ndata,event_start,left_OuterEdge,left_InnerEdge,right_InnerEdge,right_OuterEdge,event_end);
    
    
    %Choose downstream outer edge, then inner edge
    downstreamDuration = datenum(split(between(datetime(left_OuterEdge),datetime(left_InnerEdge)),'time'))*60*60*24;
    downstreamMax = calculate_max(fpi_e_timedata,fpi_e_ndata,left_OuterEdge,left_InnerEdge);
    downstreamCompressionRatio = downstreamMax/preMean;
    downstreamValue = downstreamCompressionRatio*downstreamDuration;
    
    
    
    %Choose upstream outer edge, then inner edge
    upstreamDuration = datenum(split(between(datetime(right_InnerEdge),datetime(right_OuterEdge)),'time'))*60*60*24;
    upstreamMax = calculate_max(fpi_e_timedata,fpi_e_ndata,right_InnerEdge,right_OuterEdge);
    upstreamCompressionRatio = upstreamMax/postMean;
    upstreamValue = upstreamCompressionRatio*upstreamDuration;
    
    
    %comparisonValue between boundaries
    comparisonRatio = max([downstreamValue/upstreamValue, upstreamValue/downstreamValue]);
    %SCalculate and save data to matrix
    
    BoundaryCompressionRatios(i,:) = [Event_number,...
        preMean,downstreamDuration,downstreamMax,downstreamCompressionRatio,downstreamValue,...
        postMean,upstreamDuration,upstreamMax,upstreamCompressionRatio,upstreamValue,...
        comparisonRatio];
    
end
BoundaryCompressionRatios(2:6)
BoundaryCompressionRatios(7:11)
BoundaryCompressionRatios(12)
% hist(BoundaryCompressionRatios(:,end))
toc
function [maxinBoundary] = calculate_max(fpi_timedata,fpi_ndata,left,right)
%Calculate max or min.
formatIn='yyyy-mm-dd HH:MM:SS.FFF';

%Core Time limits
Start = datenum(left,formatIn);
End = datenum(right,formatIn);

Start_index = find(fpi_timedata > Start, 1);
End_index = find(fpi_timedata > End, 1);

fpi_boundaryNdata = fpi_ndata(Start_index:End_index,1);

maxinBoundary = max(fpi_boundaryNdata);
end

function [preMean,postMean] = calculate_ambientMean(timedata,data,event_start,left_OuterEdge,left_InnerEdge,right_InnerEdge,right_OuterEdge,event_end)
%Calculate max or min.
formatIn='yyyy-mm-dd HH:MM:SS.FFF';

% Time limits

leading_start = datenum(left_OuterEdge,formatIn);
leading_start_index = find(timedata > leading_start, 1);
trailing_end = datenum(right_OuterEdge,formatIn);
trailing_end_index = find(timedata > trailing_end, 1);

eventStart = datenum(event_start,formatIn);
eventStart_index = find(timedata > eventStart, 1);
eventEnd = datenum(event_end,formatIn);
eventEnd_index = find(timedata > eventEnd, 1);

preMean = mean(data(eventStart_index:leading_start_index,:));
postMean = mean(data(trailing_end_index:eventEnd_index,:));

end