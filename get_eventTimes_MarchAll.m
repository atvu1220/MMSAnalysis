function [Event_Type, Substructure,threshold_std,data_type,...
          event_start,event_end,...
          leading_leftmost_date, leading_rightmost_date,...
          trailing_leftmost_date, trailing_rightmost_date,...
          leading_start, leading_end,...
          trailing_start, trailing_end] = get_eventTimes_MarchAll(Event_number)
    %%
    
    
    %timing window %estimate it based on the time length of the gradient.
    threshold_std = 0;
    leading_leftmost_date   = '0';
    leading_rightmost_date  = '0';
    trailing_leftmost_date  = '0';
    trailing_rightmost_date = '0';
    
    data_type = 'srvy';
    switch Event_number
        %% % %Event 2
        case 2 %compression of wave activity after event, averages therefore would be inaccurate
            Event_Type = 'SHFA';
            Substructure = 0;
            threshold_std = 2.125;
            
            event_start = '2018-03-01 01:03:44.000';
            event_end = '2018-03-01 01:04:19.000';
            
            leading_start = '2018-03-01 01:03:54.600';
            leading_end = '2018-03-01 01:03:55.100';
            
            trailing_start = '2018-03-01 01:04:10.500';
            trailing_end = '2018-03-01 01:04:11.250';
            
            
            %% % %Event 3
        case 3 %Use Mec Srvy
            Event_Type = 'HFA';
            Substructure = 0;
            threshold_std = 2.75;
            
            event_start = '2015-12-28 05:26:48.000';
            event_end = '2015-12-28 05:27:36.000';
            
            leading_start = '2015-12-28 05:26:57.750';
            leading_end = '2015-12-28 05:26:59.000';
            
            trailing_start = '2015-12-28 05:27:25.750';
            trailing_end = '2015-12-28 05:27:27.500';
            
            %% % %Event 4
        case 4
            Event_Type = 'HFA';
            Substructure = 1;
            threshold_std = 0;
            
            event_start    = '2017-11-19 11:31:15.000';
            event_end      = '2017-11-19 11:32:40.000';
            
            leading_leftmost_date   = '2017-11-19 11:29:40.000';
            leading_rightmost_date  = '2017-11-19 11:31:30.000';
            trailing_leftmost_date  = '2017-11-19 11:32:20.000';
            trailing_rightmost_date = '2017-11-19 11:33:20.000';
            
            %             leading_start  = '2017-11-19 11:31:35.082';
            %             leading_end    = '2017-11-19 11:31:35.582';
            leading_start  = '2017-11-19 11:31:39.777';
            leading_end    = '2017-11-19 11:31:40.277';
            
            
            trailing_start = '2017-11-19 11:32:08.833';
            trailing_end   = '2017-11-19 11:32:09.833';
            %% % %Event 5
        case 5
            Event_Type = 'HFA';
            Substructure = 0; %srvy mec
            threshold_std = 2.75;
            
            event_start = '2017-11-10 17:25:45.000';
            event_end = '2017-11-10 17:27:57.000';
            
            leading_start = '2017-11-10 17:26:17.500';
            leading_end = '2017-11-10 17:26:19.000';
            
            trailing_start = '2017-11-10 17:27:38.750';
            trailing_end = '2017-11-10 17:27:40.000';
            %% % %Event 6
        case 6
            Event_Type = 'HFA';
            Substructure = 1;
            threshold_std=3;
            event_start      = '2017-12-18 11:59:58.000';
            event_end        = '2017-12-18 12:01:20.000';
            
            leading_start    = '2017-12-18 12:00:21.500';
            leading_end      = '2017-12-18 12:00:23.500';
            
            trailing_start   = '2017-12-18 12:01:00.375';
            trailing_end     = '2017-12-18 12:01:01.875';
            %% % %Event 7
        case 7
            Event_Type = 'HFA';
            Substructure = 1;
            threshold_std = 2.5;
            event_start    = '2018-01-09 08:34:27.000';
            event_end      = '2018-01-09 08:35:02.000';
            
            leading_start  = '2018-01-09 08:34:34.750';
            leading_end    = '2018-01-09 08:34:35.250';
            
            trailing_start = '2018-01-09 08:34:49.500';
            trailing_end   = '2018-01-09 08:34:51.000';
            %% % %Event 8
        case 8
            Event_Type = 'HFA';
            Substructure = 0;
            threshold_std = 2.75;
            
            event_start    = '2018-01-12 01:50:21.000';
            event_end      = '2018-01-12 01:52:30.000';
            
            leading_start  = '2018-01-12 01:50:37.250';
            leading_end    = '2018-01-12 01:50:38.000';
            
            trailing_start = '2018-01-12 01:51:47.375';
            trailing_end   = '2018-01-12 01:51:48.250';
            %% % %Event 9
        case 9 %lots of wave activity
            Event_Type = 'SHFA';
            Substructure = 1;
            threshold_std = 2.75;
            event_start    = '2018-01-14 23:06:50.000';
            event_end      = '2018-01-14 23:08:05.000';
            
            leading_start  = '2018-01-14 23:07:02.000';
            leading_end    = '2018-01-14 23:07:02.500';
            
            trailing_start = '2018-01-14 23:07:48.250';
            trailing_end   = '2018-01-14 23:07:51.000';
            %% % %Event 10
        case 10 %mec srvy is broken, use brst.
            Event_Type = 'SHFA';
            data_type = 'brst';
            Substructure = 1;  %waves before event/after
            threshold_std = 2.75;
            event_start    = '2018-01-26 08:23:05.000';
            event_end      = '2018-01-26 08:23:35.000';
            
            leading_start  = '2018-01-26 08:23:15.500';
            leading_end    = '2018-01-26 08:23:17.250';
            
            trailing_start = '2018-01-26 08:23:25.250';
            trailing_end   = '2018-01-26 08:23:26.250';
            %% % %Event 11
        case 11
            Event_Type = 'HFA';
            Substructure = 0;
            threshold_std = 5;
            
            event_start = '2017-12-29 18:10:55.000';
            event_end = '2017-12-29 18:12:25.000';
            
            leading_start  = '2017-12-29 18:11:16.750';
            leading_end    = '2017-12-29 18:11:18.500';
            
            trailing_start = '2017-12-29 18:11:52.750';
            trailing_end   = '2017-12-29 18:11:54.750';
            %% % %Event 12
        case 12
            Event_Type = 'HFA';
            Substructure = 1;
            threshold_std = 2.75;
            
            event_start = '2018-01-29 18:34:30.000';
            event_end = '2018-01-29 18:35:30.000';
            
            leading_start  = '2018-01-29 18:34:57.000';
            leading_end    = '2018-01-29 18:34:57.750';
            
            trailing_start = '2018-01-29 18:35:09.000';
            trailing_end   = '2018-01-29 18:35:12.000';
            %% % %Event 13
        case 13
            Event_Type = 'SHFA';
            Substructure = 1;
            threshold_std = 2.75;
            
            event_start = '2018-03-18 22:00:04.000';
            event_end = '2018-03-18 22:00:52.000';
            
            leading_start  = '2018-03-18 22:00:26.000';
            leading_end    = '2018-03-18 22:00:27.000';
            
            trailing_start = '2018-03-18 22:00:32.500';
            trailing_end   = '2018-03-18 22:00:33.500';
            %% % %Event 14
        case 14
            Event_Type = 'HFA';
            Substructure = 1;
            threshold_std = 2.75;
            
            event_start = '2018-03-26 10:00:40.000';
            event_end = '2018-03-26 10:01:20.000';
            
            leading_start  = '2018-03-26 10:00:57.850';
            leading_end    = '2018-03-26 10:01:01.850';
            
            trailing_start = '2018-03-26 10:01:07.257';
            trailing_end   = '2018-03-26 10:01:09.257';
            %% % %Event 15
        case 15
            Event_Type = 'HFA';
            Substructure = 0;
            threshold_std = 3.25;
            
            event_start    = '2017-12-18 12:55:30.000';
            event_end      = '2017-12-18 12:57:15.000';
            
            leading_start  = '2017-12-18 12:55:58.753';
            leading_end    = '2017-12-18 12:56:01.753';
            
            trailing_start = '2017-12-18 12:56:40.636';
            trailing_end   = '2017-12-18 12:56:42.636';
            %% % %Event 16
        case 16
            Event_Type = 'HFA';
            Substructure = 0;
            threshold_std = 3.75;
            
            event_start = '2018-04-01 01:10:34.000';
            event_end = '2018-04-01 01:11:22.000';
            
            leading_start  = '2018-04-01 01:10:50.500';
            leading_end    = '2018-04-01 01:10:54.500';
            
            trailing_start = '2018-04-01 01:10:58.000';
            trailing_end   = '2018-04-01 01:11:00.000';
            %% % %Event 17
        case 17
            Event_Type = 'HFA';
            Substructure = 1;
            threshold_std = 2.25;
            
            event_start    = '2018-04-27 19:47:30.000';
            event_end      = '2018-04-27 19:48:22.000';
            
            leading_start  = '2018-04-27 19:47:45.500';
            leading_end    = '2018-04-27 19:47:47.500';
            
            trailing_start = '2018-04-27 19:47:58.750';
            trailing_end   = '2018-04-27 19:47:59.750';
            %% % %Event 18
        case 18
            Event_Type = 'HFA';
            Substructure = 0;
            threshold_std = 2.5;
            
            event_start    = '2017-12-29 17:20:35.000';
            event_end      = '2017-12-29 17:22:10.000';
            
            leading_start  = '2017-12-29 17:21:01.570';
            leading_end    = '2017-12-29 17:21:02.570';
            
            trailing_start = '2017-12-29 17:21:34.461';
            trailing_end   = '2017-12-29 17:21:34.961';
            %% % %Event 19 %foreshock bubble
        case 19
            Event_Type = 'Foreshock Bubble';
            Substructure = 0;
            threshold_std = 2.75;
            
            event_start    = '2018-04-03 13:57:25.000';
            event_end      = '2018-04-03 13:58:35.000';
            
            leading_start  = '2018-04-03 13:57:32.095';
            leading_end    = '2018-04-03 13:57:34.095';
            
            trailing_start = '2018-04-03 13:58:03.228';
            trailing_end   = '2018-04-03 13:58:05.228';
            %% % %Event 20
        case 20
            Event_Type = 'HFA';
            Substructure = 0;
            threshold_std = 0;
            
            event_start    = '2017-12-17 17:52:30.000';
            event_end      = '2017-12-17 17:53:35.000';
            
            leading_leftmost_date   = '2017-12-17 17:59:00.000';
            leading_rightmost_date  = '2017-12-17 18:04:00.000';
            trailing_leftmost_date  = '2017-12-17 17:38:00.000';
            trailing_rightmost_date = '2017-12-17 17:43:00.000';
            
            leading_start  = '2017-12-17 17:52:49.066';
            leading_end    = '2017-12-17 17:52:50.066';
            
            trailing_start = '2017-12-17 17:53:17.278';
            trailing_end   = '2017-12-17 17:53:19.278';
            %% % %Event 21
        case 21
            Event_Type = 'HFA';
            Substructure = 1;
            threshold_std = 2.00;
            
            event_start    = '2018-01-31 23:16:35.000';
            event_end      = '2018-01-31 23:17:25.000';
            
            leading_start  = '2018-01-31 23:16:48.670';
            leading_end    = '2018-01-31 23:16:50.171';
            
            trailing_start = '2018-01-31 23:17:12.499';
            trailing_end   = '2018-01-31 23:17:14.499';
            %% % %Event 22
        case 22
            Event_Type = 'Foreshock Bubble';
            Substructure = 1;
            threshold_std = 0;
            
            %             event_start    = '2017-12-29 19:11:21.000';
            event_start    = '2017-12-29 19:11:00.000';
            event_end      = '2017-12-29 19:11:50.000';
            
            leading_leftmost_date   = '2017-12-29 19:11:00.000';
            leading_rightmost_date  = '2017-12-29 19:11:05.000';
            trailing_leftmost_date  = '2017-12-29 19:13:00.000';
            trailing_rightmost_date = '2017-12-29 19:15:00.000';
            
            %             leading_start  = '2017-12-29 19:11:24.880';
            %             leading_end    = '2017-12-29 19:11:25.880';
            
            leading_start  = '2017-12-29 19:11:12.576';
            leading_end    = '2017-12-29 19:11:13.076';
            
            trailing_start = '2017-12-29 19:11:31.974';
            trailing_end   = '2017-12-29 19:11:33.974';
            %% % %Event 23 %SHFA
        case 23
            Event_Type = 'SHFA';
            Substructure = 0;
            threshold_std = 2.75;
            
            event_start    = '2017-11-16 12:13:00.000';
            event_end      = '2017-11-16 12:14:00.000';
            
            leading_start  = '2017-11-16 12:13:22.457';
            leading_end    = '2017-11-16 12:13:23.457';
            
            
            trailing_start = '2017-11-16 12:13:37.090';
            trailing_end   = '2017-11-16 12:13:37.590';
            %% % %Event 24
        case 24
            Event_Type = 'HFA';
            Substructure = 1;
            threshold_std = 3.25;
            
            event_start    = '2015-12-28 04:27:04.000';
            event_end      = '2015-12-28 04:28:34.000';
            
            leading_start  = '2015-12-28 04:27:26.374';
            leading_end    = '2015-12-28 04:27:26.874';
            
            trailing_start = '2015-12-28 04:27:58.413';
            trailing_end   = '2015-12-28 04:27:59.413';
            %% % %Event 25
        case 25
            Event_Type = 'SHFA';
            Substructure = 0;
            threshold_std = 0;
            
            event_start    = '2018-12-30 22:26:43.000';
            event_end      = '2018-12-30 22:28:00.000';
            
            leading_leftmost_date   = '2018-12-30 22:12:00.000';
            leading_rightmost_date  = '2018-12-30 22:26:25.000';
            trailing_leftmost_date  = '2018-12-30 22:28:00.000';
            trailing_rightmost_date = '2018-12-30 22:42:00.000';
            
            leading_start  = '2018-12-30 22:26:54.445';
            leading_end    = '2018-12-30 22:26:58.446';
            
            trailing_start = '2018-12-30 22:27:35.094';
            trailing_end   = '2018-12-30 22:27:39.095';
            
            %% % %Event 26
        case 26
            Event_Type = 'HFA';
            Substructure = 0; %slight substrcutures that divides interior into two regions
            threshold_std = 1.50;
            
            event_start    = '2017-11-19 07:43:20.000';
            event_end      = '2017-11-19 07:44:08.000';
            
            leading_start  = '2017-11-19 07:43:27.580';
            leading_end    = '2017-11-19 07:43:28.580';
            
            trailing_start = '2017-11-19 07:43:58.041';
            trailing_end   = '2017-11-19 07:43:59.041';
            %% % %Event 27
        case 27
            Event_Type = 'HFA';
            Substructure = 0;
            threshold_std = 2.25;
            
            event_start    = '2017-10-21 21:18:55.000';
            event_end      = '2017-10-21 21:20:00.000';
            
            leading_start  = '2017-10-21 21:19:28.038';
            leading_end    = '2017-10-21 21:19:29.038';
            
            trailing_start = '2017-10-21 21:19:32.882';
            trailing_end   = '2017-10-21 21:19:34.882';
            %% % %Event 28
        case 28
            Event_Type = 'HFA';
            Substructure = 0;
            threshold_std = 0;
            
            event_start    = '2018-01-31 23:17:18.000';
            event_end      = '2018-01-31 23:19:00.000';
            
            %             leading_leftmost_date  = '2018-01-31 23:09:30.000';
            %             leading_rightmost_date = '2018-01-31 23:11:50.000';
            
            leading_leftmost_date  = '2018-01-31 23:17:20.000';
            leading_rightmost_date = '2018-01-31 23:17:26.000';
            trailing_leftmost_date   = '2018-01-31 23:19:03.000';
            trailing_rightmost_date  = '2018-01-31 23:21:03.000';
            
            
            leading_start  = '2018-01-31 23:17:30.273';
            leading_end    = '2018-01-31 23:17:34.273';
            
            trailing_start = '2018-01-31 23:18:11.765';
            trailing_end   = '2018-01-31 23:18:13.765';
            %% % %Event 29
        case 29
            Event_Type = 'HFA';
            Substructure = 1;
            threshold_std = 2.75;
            
            event_start    = '2017-11-20 08:06:30.000';
            event_end      = '2017-11-20 08:09:00.000';
            
            leading_start  = '2017-11-20 08:07:06.555';
            leading_end    = '2017-11-20 08:07:07.555';
            
            trailing_start = '2017-11-20 08:08:14.931';
            trailing_end   = '2017-11-20 08:08:15.931';
            %% % %Event 30
        case 30
            Event_Type = 'Foreshock Bubble';
            Substructure = 0;
            threshold_std = 3.5;
            
            event_start    = '2017-12-01 14:11:20.000';
            event_end      = '2017-12-01 14:12:20.000';
            
            leading_start  = '2017-12-01 14:11:34.275';
            leading_end    = '2017-12-01 14:11:38.275';
            
            trailing_start = '2017-12-01 14:12:01.643';
            trailing_end   = '2017-12-01 14:12:04.643';
            %% % %Event 31
        case 31
            Event_Type = 'Foreshock Bubble';
            Substructure = 0;
            threshold_std = 2.25;
            
            event_start    = '2017-12-01 14:12:40.000';
            event_end      = '2017-12-01 14:13:40.000';
            
            leading_start  = '2017-12-01 14:12:55.440';
            leading_end    = '2017-12-01 14:12:55.940';
            
            trailing_start = '2017-12-01 14:13:14.355';
            trailing_end   = '2017-12-01 14:13:15.355';
            %% % %Event 32
        case 32
            Event_Type = 'HFA';
            Substructure = 0;
            threshold_std = 3.75;
            
            event_start    = '2017-12-10 00:42:10.000';
            event_end      = '2017-12-10 00:44:50.000';
            
            leading_start  = '2017-12-10 00:42:24.729';
            leading_end    = '2017-12-10 00:42:28.729';
            
            trailing_start = '2017-12-10 00:44:11.441';
            trailing_end   = '2017-12-10 00:44:13.441';
            %% % %Event 33
        case 33
            Event_Type = 'Foreshock Bubble';
            Substructure = 1;
            threshold_std = 2.00;
            
            event_start    = '2017-12-15 14:22:05.000';
            event_end      = '2017-12-15 14:24:00.000';
            
            leading_start  = '2017-12-15 14:22:38.947';
            leading_end    = '2017-12-15 14:22:41.947';
            
            trailing_start = '2017-12-15 14:23:07.338';
            trailing_end   = '2017-12-15 14:23:09.338';
            %% % %Event 34
        case 34
            Event_Type = 'Foreshock Bubble';
            Substructure = 0;
            threshold_std = 2.75;
            
            event_start    = '2017-12-18 12:35:30.000';
            event_end      = '2017-12-18 12:36:30.000';
            
            leading_start  = '2017-12-18 12:35:54.962';
            leading_end    = '2017-12-18 12:35:55.962';
            
            trailing_start = '2017-12-18 12:36:12.478';
            trailing_end   = '2017-12-18 12:36:14.478';
            %% % %Event 35
        case 35
            Event_Type = 'HFA';
            Substructure = 0;
            threshold_std = 1.5;
            
            event_start    = '2017-12-18 14:10:20.000';
            event_end      = '2017-12-18 14:11:30.000';
            
            leading_start  = '2017-12-18 14:10:57.318';
            leading_end    = '2017-12-18 14:11:00.318';
            
            trailing_start = '2017-12-18 14:11:07.443';
            trailing_end   = '2017-12-18 14:11:09.443';
            %% % %Event 36
        case 36
            Event_Type = 'Foreshock Bubble';
            Substructure = 0;
            threshold_std = 0;
            
            event_start    = '2018-02-04 10:36:50.000';
            event_end      = '2018-02-04 10:38:35.000';
            
            leading_leftmost_date  = '2018-02-04 10:31:00.000';
            leading_rightmost_date = '2018-02-04 10:37:00.000';
            trailing_leftmost_date   = '2018-02-04 10:38:30.000';
            trailing_rightmost_date  = '2018-02-04 10:50:00.000';
            
            %2018-02-04 10:38:04.504	2018-02-04 10:38:06.004
            leading_start  = '2018-02-04 10:37:47.433';
            leading_end    = '2018-02-04 10:37:48.933';
            
            trailing_start = '2018-02-04 10:38:05.824';
            trailing_end   = '2018-02-04 10:38:07.824';
            %% % %Event 37
        case 37
            Event_Type = 'HFA';
            Substructure = 1;
            threshold_std = 0;
            
            event_start    = '2018-03-04 16:30:06.000';
            event_end      = '2018-03-04 16:30:46.000';
            
            leading_leftmost_date  = '2018-03-04 16:30:05.000';
            leading_rightmost_date = '2018-03-04 16:30:13.500';
            trailing_leftmost_date   = '2018-03-04 16:30:40.000';
            trailing_rightmost_date  = '2018-03-04 16:31:12.000';
            
            leading_start  = '2018-03-04 16:30:13.345';
            leading_end    = '2018-03-04 16:30:17.345';
            
            trailing_start = '2018-03-04 16:30:28.588';
            trailing_end   = '2018-03-04 16:30:32.588';
            %% % %Event 38
        case 38
            Event_Type = 'HFA';
            Substructure = 1;
            threshold_std = 0;
            
            event_start    = '2018-04-12 06:35:00.000';
            event_end      = '2018-04-12 06:36:00.000';
            
            leading_leftmost_date  = '2018-04-12 06:33:04.000';
            leading_rightmost_date = '2018-04-12 06:35:04.000';
            trailing_leftmost_date   = '2018-04-12 06:35:57.500';
            trailing_rightmost_date  = '2018-04-12 06:36:00.000';
            
            leading_start  = '2018-04-12 06:35:18.453';
            leading_end    = '2018-04-12 06:35:21.454';
            
            trailing_start = '2018-04-12 06:35:44.196';
            trailing_end   = '2018-04-12 06:35:45.696';
            %% % %Event 39
        case 39
            Event_Type = 'Foreshock Bubble';
            Substructure = 0;
            threshold_std = 0;
            
            event_start    = '2017-01-14 06:49:15.000';
            event_end      = '2017-01-14 06:50:05.000';
            
            leading_leftmost_date  = '2017-01-14 06:34:20.000';
            leading_rightmost_date = '2017-01-14 06:40:20.000';
            trailing_leftmost_date   = '2017-01-14 06:50:01.000';
            trailing_rightmost_date  = '2017-01-14 06:51:55.000';
            
            leading_start  = '2017-01-14 06:49:23.917';
            leading_end    = '2017-01-14 06:49:25.917';
            
            trailing_start = '2017-01-14 06:49:43.902';
            trailing_end   = '2017-01-14 06:49:47.902';
            %             trailing_start = '2017-01-14 06:49:44.011';
            %             trailing_end   = '2017-01-14 06:49:48.011';
            %% % %Event 40
        case 40
            Event_Type = 'SHFA';
            Substructure = 0;
            %             threshold_std = 2.85;
            threshold_std = -2;
            
            event_start    = '2017-12-18 14:28:40.000';
            event_end      = '2017-12-18 14:30:35.000';
            
            leading_start  = '2017-12-18 14:29:24.583';
            leading_end    = '2017-12-18 14:29:28.584';
            
            trailing_start = '2017-12-18 14:29:51.435';
            trailing_end   = '2017-12-18 14:29:54.435';
            %% % %Event 41
        case 41
            Event_Type = 'SHFA';
            Substructure = 1;
            threshold_std = 0;
            
            event_start    = '2018-02-17 19:45:35.000';
            event_end      = '2018-02-17 19:46:50.000';
            
            leading_leftmost_date   = '2018-02-17 19:45:35.000';
            leading_rightmost_date  = '2018-02-17 19:45:43.000';
            trailing_leftmost_date  = '2018-02-17 19:46:45.000';
            trailing_rightmost_date = '2018-02-17 19:47:30.000';
            
            leading_start  = '2018-02-17 19:46:07.594';
            leading_end    = '2018-02-17 19:46:10.594';
            
            trailing_start = '2018-02-17 19:46:22.852';
            trailing_end   = '2018-02-17 19:46:26.852';
            
            %% % %Event 42
        case 42
            Event_Type = 'HFA';
            Substructure = 0;
            threshold_std = 2.75;
            % threshold_std = -2;
            
            event_start    = '2017-12-29 17:29:45.000';
            event_end      = '2017-12-29 17:30:35.000';
            
            leading_start  = '2017-12-29 17:30:03.874';
            leading_end    = '2017-12-29 17:30:05.874';
            
            trailing_start = '2017-12-29 17:30:15.570';
            trailing_end   = '2017-12-29 17:30:17.570';
            %% % %Event 43 %SHFA
        case 43
            Event_Type = 'SHFA';
            Substructure = 0;
            threshold_std = 2.25;
            % threshold_std = -2;
            
            event_start    = '2017-12-29 19:22:40.000';
            event_end      = '2017-12-29 19:23:15.000';
            
            leading_start  = '2017-12-29 19:22:51.828';
            leading_end    = '2017-12-29 19:22:55.828';
            
            trailing_start = '2017-12-29 19:23:07.625';
            trailing_end   = '2017-12-29 19:23:08.625';
            %% % %Event 44, connected 44-45
        case 44
            Event_Type = 'HFA';
            Substructure = 0;
            threshold_std = 4;
            %threshold_std = -2;
            
            event_start    = '2018-01-09 09:18:45.000';
            event_end      = '2018-01-09 09:20:00.000';
            
            leading_start  = '2018-01-09 09:18:52.976';
            leading_end    = '2018-01-09 09:18:55.976';
            
            trailing_start = '2018-01-09 09:19:12.593';
            trailing_end   = '2018-01-09 09:19:16.593';
            %% % %Event 45, connected 44-45
        case 45
            Event_Type = 'HFA';
            Substructure = 1;
            threshold_std = 4;
            % threshold_std = -2;
            
            event_start    = '2018-01-09 09:18:45.000';
            event_end      = '2018-01-09 09:20:00.000';
            
            leading_start  = '2018-01-09 09:19:18.953';
            leading_end    = '2018-01-09 09:19:22.953';
            
            trailing_start = '2018-01-09 09:19:43.469';
            trailing_end   = '2018-01-09 09:19:45.469';
            %% % %Event 46
        case 46
            Event_Type = 'HFA';
            Substructure = 1;
            threshold_std = 2.75;
            
            event_start    = '2018-03-29 03:25:45.000';
            event_end      = '2018-03-29 03:26:35.000';
            
            leading_leftmost_date   = '2018-03-29 03:17:40.000';
            leading_rightmost_date  = '2018-03-29 03:24:40.000';
            trailing_leftmost_date  = '2018-03-29 03:26:40.000';
            trailing_rightmost_date = '2018-03-29 03:28:20.000';
            
            leading_start  = '2018-03-29 03:25:53.481';
            leading_end    = '2018-03-29 03:25:54.481';
            
            trailing_start = '2018-03-29 03:26:23.013';
            trailing_end   = '2018-03-29 03:26:25.013';
            %% % %Event 47
        case 47
            Event_Type = 'Foreshock Bubble';
            Substructure = 0;
            %threshold_std = 5;
            threshold_std = -2.65;
            
            event_start    = '2018-01-23 09:11:50.000';
            event_end      = '2018-01-23 09:12:50.000';
            
            leading_start  = '2018-01-23 09:12:05.643';
            leading_end    = '2018-01-23 09:12:06.643';
            
            trailing_start = '2018-01-23 09:12:22.557';
            trailing_end   = '2018-01-23 09:12:26.557';
            %% % %Event 48
        case 48
            Event_Type = 'HFA';
            Substructure = 1;
            threshold_std = 1;
            
            event_start    = '2018-02-20 13:42:00.000';
            event_end      = '2018-02-20 13:46:30.000';
            
            leading_leftmost_date   = '2018-02-20 13:30:00.000';
            leading_rightmost_date  = '2018-02-20 13:42:00.000';
            trailing_leftmost_date  = '2018-02-20 13:46:00.000';
            trailing_rightmost_date = '2018-02-20 13:52:55.000';
            
            leading_start  = '2018-02-20 13:42:56.566';
            leading_end    = '2018-02-20 13:42:57.066';
            
            %             trailing_start = '2018-02-20 13:43:37.215';
            %             trailing_end   = '2018-02-20 13:43:39.215';
            
            trailing_start = '2018-02-20 13:45:45.904';
            trailing_end   = '2018-02-20 13:45:48.904';
            %% % %Event 49, connected with Event 50
        case 49
            Event_Type = 'Foreshock Bubble';
            Substructure = 0;
            threshold_std = 0;
            
            event_start    = '2018-01-29 03:35:45.000';
            event_end      = '2018-01-29 03:36:45.000';
            
            leading_leftmost_date   = '2018-01-29 03:20:45.000';
            leading_rightmost_date  = '2018-01-29 03:35:50.000';
            trailing_leftmost_date  = '2018-01-29 03:36:25.000';
            trailing_rightmost_date = '2018-01-29 03:36:35.000';
            
            leading_start  = '2018-01-29 03:35:51.472';
            leading_end    = '2018-01-29 03:35:54.472';
            
            trailing_start = '2018-01-29 03:36:15.027';
            trailing_end   = '2018-01-29 03:36:17.027';
            %% % %Event 50, connected with Event 49
        case 50
            Event_Type = 'HFA';
            Substructure = 1;
            %threshold_std = 5;
            threshold_std = 0;
            
            %Overall Events
            event_start    = '2018-03-01 01:05:40.000';
            event_end      = '2018-03-01 01:06:20.000';
            
            leading_leftmost_date   = '2018-03-01 01:05:25.000';
            leading_rightmost_date  = '2018-03-01 01:05:45.000';
            trailing_leftmost_date  = '2018-03-01 01:06:18.000';
            trailing_rightmost_date = '2018-03-01 01:06:26.500';
            
            leading_start  = '2018-03-01 01:05:55.875';
            leading_end    = '2018-03-01 01:05:57.875';
            
            trailing_start = '2018-03-01 01:06:09.922';
            trailing_end   = '2018-03-01 01:06:10.922';
            
            
            %% % %Event 51, connected with Event 50,49,48,2
        case 51
            Event_Type = 'Foreshock Bubble';
            Substructure = 1;
            %threshold_std = 5;
            threshold_std = 2.75;
            
            event_start    = '2018-03-01 01:06:15.000';
            event_end      = '2018-03-01 01:07:36.000';
            
            leading_start  = '2018-03-01 01:06:44.360';
            leading_end    = '2018-03-01 01:06:48.360';
            
            trailing_start = '2018-03-01 01:07:01.517';
            trailing_end   = '2018-03-01 01:07:03.517';
            %% % %Event 52
        case 52
            Event_Type = 'HFA';
            Substructure = 1;
            threshold_std = 0;
            
            event_start    = '2017-12-18 12:04:50.000';
            event_end      = '2017-12-18 12:06:35.000';
            
            leading_leftmost_date   = '2017-12-18 12:04:25.000';
            leading_rightmost_date  = '2017-12-18 12:04:45.000';
            trailing_leftmost_date  = '2017-12-18 12:06:16.000';
            trailing_rightmost_date = '2017-12-18 12:06:29.000';
            
            leading_start  = '2017-12-18 12:05:24.420';
            leading_end    = '2017-12-18 12:05:25.420';
            
            %             leading_start  = '2017-12-18 12:05:24.100';
            %             leading_end    = '2017-12-18 12:05:25.600';
            
            %             trailing_start = '2017-12-18 12:06:01.827';
            %             trailing_end   = '2017-12-18 12:06:03.827';
            
            trailing_start = '2017-12-18 12:06:02.194';
            trailing_end   = '2017-12-18 12:06:03.694';
            
            %% % %Event 53 %Foreshock Bubble
        case 53
            Event_Type = 'Foreshock Bubble';
            Substructure = 0;
            %threshold_std = 5;
            threshold_std = -0.5;
            
            event_start    = '2018-04-03 13:57:10.000';
            event_end      = '2018-04-03 13:58:40.000';
            
            leading_start  = '2018-04-03 13:57:32.259';
            leading_end    = '2018-04-03 13:57:34.259';
            
            trailing_start = '2018-04-03 13:58:01.447';
            trailing_end   = '2018-04-03 13:58:04.447';
            %% % %Event 54
        case 54
            Event_Type = 'Foreshock Bubble';
            Substructure = 1;
            threshold_std = 3.75;
            
            event_start    = '2018-03-30 08:40:45.000';
            event_end      = '2018-03-30 08:42:45.000';
            
            leading_start  = '2018-03-30 08:41:00.583';
            leading_end    = '2018-03-30 08:41:04.583';
            
            %             trailing_start = '2018-03-30 08:42:18.459';
            %             trailing_end   = '2018-03-30 08:42:21.460';
            
            trailing_start = '2018-03-30 08:42:25.569';
            trailing_end   = '2018-03-30 08:42:27.569';
            %% % %Event 55
        case 55
            Event_Type = 'SHFA';
            Substructure = 1;
            threshold_std = 4.75;
            
            
            event_start    = '2018-03-18 00:10:10.000';
            event_end      = '2018-03-18 00:11:05.000';
            
            leading_start  = '2018-03-18 00:10:38.597';
            leading_end    = '2018-03-18 00:10:39.597';
            
            trailing_start = '2018-03-18 00:10:47.699';
            trailing_end   = '2018-03-18 00:10:49.699';
            %% % %Event 56 %SHFA
        case 56
            Event_Type = 'SHFA';
            Substructure = 1;
            threshold_std = 2.75;
            
            event_start    = '2018-03-01 00:31:40.000';
            event_end      = '2018-03-01 00:32:25.000';
            
            leading_start  = '2018-03-01 00:31:51.590';
            leading_end    = '2018-03-01 00:31:54.590';
            
            trailing_start = '2018-03-01 00:32:09.020';
            trailing_end   = '2018-03-01 00:32:10.020';
            %% % %Event 57
        case 57
            Event_Type = 'SHFA';
            Substructure = 0;
            threshold_std = 1.75;
            
            event_start    = '2018-03-01 00:29:50.000';
            event_end      = '2018-03-01 00:31:10.000';
            
            leading_start  = '2018-03-01 00:30:23.776';
            leading_end    = '2018-03-01 00:30:27.776';
            
            trailing_start = '2018-03-01 00:30:52.550';
            trailing_end   = '2018-03-01 00:30:55.550';
            %% % %Event 58
        case 58
            Event_Type = 'HFA';
            Substructure = 1;
            threshold_std = 0;
            
            event_start    = '2018-02-21 11:34:30.000';
            event_end      = '2018-02-21 11:37:20.000';
            
            leading_leftmost_date   = '2018-02-21 11:25:00.000';
            leading_rightmost_date  = '2018-02-21 11:30:00.000';
            trailing_leftmost_date  = '2018-02-21 11:37:15.000';
            trailing_rightmost_date = '2018-02-21 11:40:00.000';
            
            leading_start  = '2018-02-21 11:35:03.498';
            leading_end    = '2018-02-21 11:35:04.498';
            
            trailing_start = '2018-02-21 11:36:21.054';
            trailing_end   = '2018-02-21 11:36:25.054';
            %% % %Event 59
        case 59
            Event_Type = 'HFA';
            Substructure = 0;
            threshold_std = 0;
            
            event_start    = '2017-12-15 14:24:50.000';
            event_end      = '2017-12-15 14:26:10.000';
            
            leading_leftmost_date   = '2017-12-15 14:23:30.000';
            leading_rightmost_date  = '2017-12-15 14:24:20.000';
            trailing_leftmost_date  = '2017-12-15 14:26:40.000';
            trailing_rightmost_date = '2017-12-15 14:28:40.000';
            
            %             leading_start  = '2017-12-15 14:25:21.754';
            %             leading_end    = '2017-12-15 14:25:24.754';
            
            leading_start  = '2017-12-15 14:25:24.286';
            leading_end    = '2017-12-15 14:25:26.286';
            
            trailing_start = '2017-12-15 14:25:30.559';
            trailing_end   = '2017-12-15 14:25:31.559';
            %% % %Event 60
        case 60
            Event_Type = 'HFA';
            Substructure = 1;
            threshold_std = 0;
            
            event_start    = '2017-12-15 14:24:50.000';
            event_end      = '2017-12-15 14:26:10.000';
            
            leading_leftmost_date   = '2017-12-15 14:23:30.000';
            leading_rightmost_date  = '2017-12-15 14:24:20.000';
            trailing_leftmost_date  = '2017-12-15 14:26:40.000';
            trailing_rightmost_date = '2017-12-15 14:28:40.000';
            
            leading_start  = '2017-12-15 14:25:31.403';
            leading_end    = '2017-12-15 14:25:31.903';
            
            trailing_start = '2017-12-15 14:25:36.809';
            trailing_end   = '2017-12-15 14:25:37.809';
            %% % %Event 61
        case 61
            Event_Type = 'HFA';
            Substructure = 1;
            threshold_std = 0;
            
            event_start    = '2017-12-15 14:24:50.000';
            event_end      = '2017-12-15 14:26:10.000';
            
            leading_leftmost_date   = '2017-12-15 14:23:30.000';
            leading_rightmost_date  = '2017-12-15 14:24:20.000';
            trailing_leftmost_date  = '2017-12-15 14:26:40.000';
            trailing_rightmost_date = '2017-12-15 14:28:40.000';
            
            leading_start  = '2017-12-15 14:25:37.216';
            leading_end    = '2017-12-15 14:25:38.216';
            
            trailing_start = '2017-12-15 14:25:44.223';
            trailing_end   = '2017-12-15 14:25:47.224';
            %% % %Event 62
        case 62
            Event_Type = 'HFA';
            Substructure = 0;
            threshold_std = 0;
            
            event_start    = '2018-01-07 05:50:30.000';
            event_end      = '2018-01-07 05:53:15.000';
            
            leading_leftmost_date   = '2018-01-07 05:40:00.000';
            leading_rightmost_date  = '2018-01-07 05:45:00.000';
            trailing_leftmost_date  = '2018-01-07 05:52:50.000';
            trailing_rightmost_date = '2018-01-07 05:53:14.000';
            
            leading_start  = '2018-01-07 05:51:10.652';
            leading_end    = '2018-01-07 05:51:13.653';
            
            trailing_start = '2018-01-07 05:51:33.309';
            trailing_end   = '2018-01-07 05:51:35.309';
            %% % %Event 63
        case 63
            Event_Type = 'HFA';
            Substructure = 0;
            threshold_std = 0;
            
            event_start    = '2018-01-07 05:50:30.000';
            event_end      = '2018-01-07 05:53:15.000';
            
            leading_leftmost_date   = '2018-01-07 05:40:00.000';
            leading_rightmost_date  = '2018-01-07 05:45:00.000';
            trailing_leftmost_date  = '2018-01-07 05:52:50.000';
            trailing_rightmost_date = '2018-01-07 05:53:14.000';
            
            leading_start  = '2018-01-07 05:52:08.216';
            leading_end    = '2018-01-07 05:52:11.216';
            
            trailing_start = '2018-01-07 05:52:33.036';
            trailing_end   = '2018-01-07 05:52:35.036';
            %% % %Event 64
        case 64
            Event_Type = 'Foreshock Bubble';
            Substructure = 0;
            threshold_std = 0;
            
            event_start    = '2018-01-07 05:53:00.000';
            event_end      = '2018-01-07 05:54:10.000';
            
            leading_leftmost_date   = '2018-01-07 05:52:50.000';
            leading_rightmost_date  = '2018-01-07 05:53:14.000';
            trailing_leftmost_date  = '2018-01-07 05:53:50.000';
            %             trailing_rightmost_date = '2018-01-07 06:09:10.000';
            trailing_rightmost_date = '2018-01-07 05:54:06.000';
            
            leading_start  = '2018-01-07 05:53:13.826';
            leading_end    = '2018-01-07 05:53:15.826';
            
            trailing_start = '2018-01-07 05:53:23.365';
            trailing_end   = '2018-01-07 05:53:26.365';
            %% % %Event 65
        case 65
            Event_Type = 'Foreshock Bubble';
            Substructure = 0;
            threshold_std = 2.75;
            
            event_start    = '2018-01-15 16:54:30.000';
            event_end      = '2018-01-15 16:55:10.000';
            
            leading_leftmost_date   = '2018-01-15 16:54:30.000';
            leading_rightmost_date  = '2018-01-15 16:54:30.000';
            trailing_leftmost_date  = '2018-01-15 16:55:10.000';
            trailing_rightmost_date = '2018-01-15 16:55:10.000';
            
            leading_start  = '2018-01-15 16:54:35.726';
            leading_end    = '2018-01-15 16:54:39.726';
            
            trailing_start = '2018-01-15 16:54:54.546';
            trailing_end   = '2018-01-15 16:54:55.546';
            %% % %Event 66
        case 66
            Event_Type = 'Foreshock Bubble';
            Substructure = 0;
            threshold_std = 0;
            
            event_start    = '2018-01-12 02:59:15.000';
            event_end      = '2018-01-12 03:00:10.000';
            
            leading_leftmost_date   = '2018-01-12 02:45:30.000';
            leading_rightmost_date  = '2018-01-12 02:59:20.000';
            trailing_leftmost_date  = '2018-01-12 03:00:10.000';
            trailing_rightmost_date = '2018-01-12 03:07:00.000';
            
            leading_start  = '2018-01-12 02:59:23.250';
            leading_end    = '2018-01-12 02:59:34.250';
            
            trailing_start = '2018-01-12 02:59:51.027';
            trailing_end   = '2018-01-12 02:59:54.027';
            %% % %Event 67
        case 67
            Event_Type = 'HFA';
            Substructure = 1;
            threshold_std = 0;
            
            event_start    = '2017-12-18 12:08:17.000';
            event_end      = '2017-12-18 12:10:55.000';
            
            leading_leftmost_date   = '2017-12-18 12:08:28.000';
            leading_rightmost_date  = '2017-12-18 12:08:32.000';
            trailing_leftmost_date  = '2017-12-18 12:12:30.000';
            trailing_rightmost_date = '2017-12-18 12:21:00.000';
            
            leading_start  = '2017-12-18 12:08:58.127';
            leading_end    = '2017-12-18 12:09:00.127';
            
            trailing_start = '2017-12-18 12:09:50.065';
            trailing_end   = '2017-12-18 12:09:53.065';
            %% % %Event 68
        case 68
            Event_Type = 'HFA';
            Substructure = 1;
            threshold_std = 0;
            
            event_start    = '2018-10-25 02:10:00.000';
            event_end      = '2018-10-25 02:10:40.000';
            
            leading_leftmost_date   = '2018-10-25 02:02:00.000';
            leading_rightmost_date  = '2018-10-25 02:09:00.000';
            %             trailing_leftmost_date  = '2018-10-25 02:11:40.000';
            %             trailing_rightmost_date = '2018-10-25 02:20:00.000';
            trailing_leftmost_date  = '2018-10-25 02:11:30.000';
            trailing_rightmost_date = '2018-10-25 02:12:10.000';
            
            leading_start  = '2018-10-25 02:10:08.860';
            leading_end    = '2018-10-25 02:10:10.360';
            
            trailing_start = '2018-10-25 02:10:26.712';
            trailing_end   = '2018-10-25 02:10:28.212';
            %% % %Event 69 ['','']
        case 69
            Event_Type = 'HFA';
            Substructure = 1;
            threshold_std = 0;
            
            event_start    = '2018-11-11 22:03:40.000';
            event_end      = '2018-11-11 22:05:20.000';
            
            leading_leftmost_date   = '2018-11-11 21:58:40.000';
            leading_rightmost_date  = '2018-11-11 22:02:30.000';
            trailing_leftmost_date  = '2018-11-11 22:05:20.000';
            trailing_rightmost_date = '2018-11-11 22:06:05.000';
            
            leading_start  = '2018-11-11 22:04:17.450';
            leading_end    = '2018-11-11 22:04:19.450';
            
            trailing_start = '2018-11-11 22:04:48.817';
            trailing_end   = '2018-11-11 22:04:52.817';
            %% % %Event 70 ['','']
        case 70
            Event_Type = 'HFA';
            Substructure = 1;
            threshold_std = 0;
            
            event_start    = '2018-11-14 15:12:20.000';
            event_end      = '2018-11-14 15:13:40.000';
            
            leading_leftmost_date   = '2018-11-14 15:10:55.000';
            leading_rightmost_date  = '2018-11-14 15:11:50.000';
            trailing_leftmost_date  = '2018-11-14 15:14:40.000';
            trailing_rightmost_date = '2018-11-14 15:17:20.000';
            
            leading_start  = '2018-11-14 15:12:46.208';
            leading_end    = '2018-11-14 15:12:49.208';
            
            trailing_start = '2018-11-14 15:13:15.068';
            trailing_end   = '2018-11-14 15:13:17.068';
            %% % %Event 71
        case 71
            Event_Type = 'SHFA';
            Substructure = 0;
            threshold_std = 0;
            
            event_start    = '2018-11-15 23:01:00.000';
            event_end      = '2018-11-15 23:02:30.000';
            
            %             leading_leftmost_date   = '2018-11-15 22:47:00.000';
            %             leading_rightmost_date  = '2018-11-15 22:53:00.000';
            %             trailing_leftmost_date  = '2018-11-15 23:07:30.000';
            %             trailing_rightmost_date = '2018-11-15 23:17:00.000';
            
            leading_leftmost_date   = '2018-11-15 22:58:10.000';
            leading_rightmost_date  = '2018-11-15 22:59:20.000';
            trailing_leftmost_date  = '2018-11-15 23:05:40.000';
            trailing_rightmost_date = '2018-11-15 23:06:30.000';
            
            leading_start  = '2018-11-15 23:01:42.615';
            leading_end    = '2018-11-15 23:01:46.615';
            
            %             trailing_start = '2018-11-15 23:01:54.263';
            %             trailing_end   = '2018-11-15 23:01:55.263';
            trailing_start = '2018-11-15 23:01:50.467';
            trailing_end   = '2018-11-15 23:01:51.967';
            %% % %Event 72 ['','']
        case 72
            Event_Type = 'Foreshock Bubble';
            Substructure = 0;
            threshold_std = 0;
            
            event_start    = '2018-11-17 09:48:40.000';
            event_end      = '2018-11-17 09:49:05.000';
            
            leading_leftmost_date   = '2018-11-17 09:48:25.000';
            leading_rightmost_date  = '2018-11-17 09:48:42.000';
            trailing_leftmost_date  = '2018-11-17 09:57:30.000';
            trailing_rightmost_date = '2018-11-17 10:04:00.000';
            
            %             leading_start  = '2018-11-17 09:48:47.411';
            %             leading_end    = '2018-11-17 09:48:50.411';
            leading_start  = '2018-11-17 09:48:48.161';
            leading_end    = '2018-11-17 09:48:50.161';
            
            trailing_start = '2018-11-17 09:48:51.340';
            trailing_end   = '2018-11-17 09:48:54.340';
            %% % %Event 73 ['','']
        case 73
            Event_Type = 'HFA';
            Substructure = 0;
            threshold_std = 0;
            
            event_start    = '2018-11-24 09:53:45.000';
            event_end      = '2018-11-24 09:55:20.000';
            
            leading_leftmost_date   = '2018-11-24 09:53:40.000';
            leading_rightmost_date  = '2018-11-24 09:53:55.000';
            trailing_leftmost_date  = '2018-11-24 10:07:30.000';
            trailing_rightmost_date = '2018-11-24 10:10:00.000';
            
            leading_start  = '2018-11-24 09:54:15.457';
            leading_end    = '2018-11-24 09:54:19.457';
            
            trailing_start = '2018-11-24 09:54:53.293';
            trailing_end   = '2018-11-24 09:54:57.293';
            %% % %Event 74  ['','']
        case 74
            Event_Type = 'Foreshock Bubble';
            Substructure = 0;
            threshold_std = 0;
            
            event_start    = '2018-11-24 10:03:30.000';
            event_end      = '2018-11-24 10:04:45.000';
            
            leading_leftmost_date   = '2018-11-24 10:02:00.000';
            leading_rightmost_date  = '2018-11-24 10:02:10.000';
            trailing_leftmost_date  = '2018-11-24 10:04:40.000';
            trailing_rightmost_date = '2018-11-24 10:05:40.000';
            
            leading_start  = '2018-11-24 10:03:54.629';
            leading_end    = '2018-11-24 10:03:55.629';
            
            trailing_start = '2018-11-24 10:04:13.707';
            trailing_end   = '2018-11-24 10:04:16.707';
            %% % %Event 75 ['','']
        case 75
            Event_Type = 'HFA';
            Substructure = 1;
            threshold_std = 0;
            
            event_start    = '2018-12-11 06:52:03.000';
            event_end      = '2018-12-11 06:53:35.000';
            
            leading_leftmost_date   = '2018-12-11 06:51:00.000';
            leading_rightmost_date  = '2018-12-11 06:51:30.000';
            trailing_leftmost_date  = '2018-12-11 06:53:40.000';
            trailing_rightmost_date = '2018-12-11 06:57:00.000';
            
            leading_start  = '2018-12-11 06:52:18.399';
            leading_end    = '2018-12-11 06:52:22.399';
            
            trailing_start = '2018-12-11 06:52:59.618';
            trailing_end   = '2018-12-11 06:53:02.618';
            %% % %Event 76 ['','']
        case 76
            Event_Type = 'HFA';
            Substructure = 1;
            threshold_std = 0;
            
            event_start    = '2018-12-12 05:27:35.000';
            event_end      = '2018-12-12 05:28:40.000';
            
            leading_leftmost_date   = '2018-12-12 05:27:10.000';
            leading_rightmost_date  = '2018-12-12 05:27:36.000';
            trailing_leftmost_date  = '2018-12-12 05:28:38.000';
            trailing_rightmost_date = '2018-12-12 05:28:56.000';
            
            leading_start  = '2018-12-12 05:27:36.829';
            leading_end    = '2018-12-12 05:27:40.829';
            
            trailing_start = '2018-12-12 05:28:24.126';
            trailing_end   = '2018-12-12 05:28:27.126';
            
            %% % %Event 77 ['','']
        case 77
            Event_Type = 'Foreshock Bubble';
            Substructure = 0;
            threshold_std = 0;
            
            event_start    = '2018-12-14 03:31:40.000';
            event_end      = '2018-12-14 03:32:40.000';
            
            leading_leftmost_date   = '2018-12-14 03:31:15.000';
            leading_rightmost_date  = '2018-12-14 03:32:00.000';
            trailing_leftmost_date  = '2018-12-14 03:32:36.000';
            trailing_rightmost_date = '2018-12-14 03:44:00.000';
            
            leading_start  = '2018-12-14 03:32:01.024';
            leading_end    = '2018-12-14 03:32:02.524';
            
            trailing_start = '2018-12-14 03:32:20.001';
            trailing_end   = '2018-12-14 03:32:24.001';
            %% % %Event 78 ['','']
        case 78
            Event_Type = 'Foreshock Bubble';
            Substructure = 0;
            threshold_std = 0;
            
            event_start    = '2018-12-14 03:44:25.000';
            event_end      = '2018-12-14 03:45:10.000';
            
            leading_leftmost_date   = '2018-12-14 03:35:00.000';
            leading_rightmost_date  = '2018-12-14 03:44:10.000';
            trailing_leftmost_date  = '2018-12-14 03:45:30.000';
            trailing_rightmost_date = '2018-12-14 03:50:50.000';
            
            leading_start  = '2018-12-14 03:44:44.129';
            leading_end    = '2018-12-14 03:44:45.629';
            
            trailing_start = '2018-12-14 03:44:57.895';
            trailing_end   = '2018-12-14 03:44:59.895';
            %% % %Event 79 ['','']
        case 79
            Event_Type = 'SHFA';
            Substructure = 0;
            threshold_std = 0;
            
            event_start    = '2018-12-14 04:01:40.000';
            event_end      = '2018-12-14 04:02:50.000';
            
            leading_leftmost_date   = '2018-12-14 04:00:10.000';
            leading_rightmost_date  = '2018-12-14 04:01:20.000';
            trailing_leftmost_date  = '2018-12-14 04:15:10.000';
            trailing_rightmost_date = '2018-12-14 04:17:10.000';
            
            leading_start  = '2018-12-14 04:01:59.237';
            leading_end    = '2018-12-14 04:02:03.237';
            
            trailing_start = '2018-12-14 04:02:28.699';
            trailing_end   = '2018-12-14 04:02:31.699';
            %% % %Event 80 ['','']
        case 80
            Event_Type = 'Foreshock Bubble';
            Substructure = 0;
            threshold_std = 0;
            
            event_start    = '2018-12-14 04:20:23.000';
            event_end      = '2018-12-14 04:22:00.000';
            
            leading_leftmost_date   = '2018-12-14 04:19:10.000';
            leading_rightmost_date  = '2018-12-14 04:20:05.000';
            trailing_leftmost_date  = '2018-12-14 04:21:54.000';
            trailing_rightmost_date = '2018-12-14 04:22:02.000';
            
            leading_start  = '2018-12-14 04:21:01.277';
            leading_end    = '2018-12-14 04:21:05.277';
            
            trailing_start = '2018-12-14 04:21:30.621';
            trailing_end   = '2018-12-14 04:21:33.621';
            %% % %Event 81 ['','']
        case 81
            Event_Type = 'SHFA';
            Substructure = 0;
            threshold_std = 0;
            
            event_start    = '2018-12-10 04:10:05.000';
            event_end      = '2018-12-10 04:11:05.000';
            
            leading_leftmost_date   = '2018-12-10 04:04:00.000';
            leading_rightmost_date  = '2018-12-10 04:10:05.000';
            trailing_leftmost_date  = '2018-12-10 04:10:55.000';
            trailing_rightmost_date = '2018-12-10 04:18:00.000';
            
            leading_start  = '2018-12-10 04:10:23.072';
            leading_end    = '2018-12-10 04:10:24.572';
            
            trailing_start = '2018-12-10 04:10:37.447';
            trailing_end   = '2018-12-10 04:10:40.447';
            
            %% % %Event 82
        case 82
            Event_Type = 'SHFA';
            Substructure = 0;
            threshold_std = 0;
            
            event_start    = '2018-12-10 04:18:40.000';
            event_end      = '2018-12-10 04:20:30.000';
            
            leading_leftmost_date   = '2018-12-10 04:18:40.000';
            leading_rightmost_date  = '2018-12-10 04:18:49.500';
            trailing_leftmost_date  = '2018-12-10 04:20:57.000';
            trailing_rightmost_date = '2018-12-10 04:21:08.500';
            
            leading_start  = '2018-12-10 04:18:51.532';
            leading_end    = '2018-12-10 04:18:55.532';
            
            trailing_start = '2018-12-10 04:19:40.408';
            trailing_end   = '2018-12-10 04:19:42.408';
            %% % %Event 83 ['','']
        case 83
            Event_Type = 'HFA';
            Substructure = 0;
            threshold_std = 0;
            
            event_start    = '2018-12-10 04:40:55.000';
            event_end      = '2018-12-10 04:42:00.000';
            
            leading_leftmost_date   = '2018-12-10 04:40:57.000';
            leading_rightmost_date  = '2018-12-10 04:40:59.000';
            trailing_leftmost_date  = '2018-12-10 04:41:57.000';
            trailing_rightmost_date = '2018-12-10 04:42:07.000';
            
            leading_start  = '2018-12-10 04:41:08.208';
            leading_end    = '2018-12-10 04:41:11.208';
            
            trailing_start = '2018-12-10 04:41:37.169';
            trailing_end   = '2018-12-10 04:41:38.169';
            %% % %Event 84 ['','']
        case 84
            Event_Type = 'HFA';
            Substructure = 1;
            threshold_std = 0;
            
            event_start    = '2018-12-10 04:52:40.000';
            event_end      = '2018-12-10 04:53:50.000';
            
            leading_leftmost_date   = '2018-12-10 04:48:30.000';
            leading_rightmost_date  = '2018-12-10 04:49:06.000';
            trailing_leftmost_date  = '2018-12-10 04:56:42.000';
            trailing_rightmost_date = '2018-12-10 04:58:27.000';
            
            %             leading_start  = '2018-12-10 04:53:03.561';
            %             leading_end    = '2018-12-10 04:53:05.061';
            %
            %             leading_start  = '2018-12-10 04:53:03.585';
            %             leading_end    = '2018-12-10 04:53:05.085';
            
            leading_start  = '2018-12-10 04:53:05.397';
            leading_end    = '2018-12-10 04:53:06.897';
            
            trailing_start = '2018-12-10 04:53:30.031';
            trailing_end   = '2018-12-10 04:53:33.031';
            %% % %Event 85 ['','']
        case 85
            Event_Type = 'HFA';
            Substructure = 1;
            threshold_std = 0;
            
            event_start    = '2018-12-10 05:21:40.000';
            event_end      = '2018-12-10 05:22:40.000';
            
            leading_leftmost_date   = '2018-12-10 05:13:20.000';
            leading_rightmost_date  = '2018-12-10 05:15:30.000';
            trailing_leftmost_date  = '2018-12-10 05:22:30.000';
            trailing_rightmost_date = '2018-12-10 05:22:32.500';
            
            leading_start  = '2018-12-10 05:21:49.320';
            leading_end    = '2018-12-10 05:21:50.820';
            
            trailing_start = '2018-12-10 05:22:28.360';
            trailing_end   = '2018-12-10 05:22:30.860';
            
            %% % %Event 86 ['','']
        case 86
            Event_Type = 'HFA';
            Substructure = 0;
            threshold_std = 0;
            
            event_start    = '2018-12-10 05:23:20.000';
            event_end      = '2018-12-10 05:25:10.000';
            
            leading_leftmost_date   = '2018-12-10 05:23:11.000';
            leading_rightmost_date  = '2018-12-10 05:23:22.500';
            trailing_leftmost_date  = '2018-12-10 05:25:00.000';
            trailing_rightmost_date = '2018-12-10 05:25:40.000';
            
            %             leading_start  = '2018-12-10 05:23:28.962';
            %             leading_end    = '2018-12-10 05:23:30.962';
            
            leading_start  = '2018-12-10 05:23:30.087';
            leading_end    = '2018-12-10 05:23:32.587';
            
            trailing_start = '2018-12-10 05:24:41.409';
            trailing_end   = '2018-12-10 05:24:43.409';
            
            %% % %Event 87 ['','']
        case 87
            Event_Type = 'HFA';
            Substructure = 0;
            threshold_std = 0;
            
            event_start    = '2018-12-10 05:25:30.000';
            event_end      = '2018-12-10 05:27:35.000';
            
            leading_leftmost_date   = '2018-12-10 05:25:03.500';
            leading_rightmost_date  = '2018-12-10 05:25:43.000';
            trailing_leftmost_date  = '2018-12-10 05:27:27.000';
            trailing_rightmost_date = '2018-12-10 05:27:35.500';
            
            leading_start  = '2018-12-10 05:25:44.370';
            leading_end    = '2018-12-10 05:25:45.870';
            
            trailing_start = '2018-12-10 05:26:17.222';
            trailing_end   = '2018-12-10 05:26:18.723';
            
            %% % %Event 88 ['','']
        case 88
            Event_Type = 'HFA';
            Substructure = 1;
            threshold_std = 0;
            
            event_start    = '2018-12-10 05:48:25.000';
            event_end      = '2018-12-10 05:49:16.000';
            
            leading_leftmost_date   = '2018-12-10 05:43:00.000';
            leading_rightmost_date  = '2018-12-10 05:48:20.000';
            trailing_leftmost_date  = '2018-12-10 05:49:11.500';
            trailing_rightmost_date = '2018-12-10 05:49:14.000';
            
            leading_start  = '2018-12-10 05:48:46.906';
            leading_end    = '2018-12-10 05:48:49.906';
            
            trailing_start = '2018-12-10 05:49:04.593';
            trailing_end   = '2018-12-10 05:49:07.594';
            
            %% % %Event 89 ['','']
        case 89
            Event_Type = 'SHFA';
            Substructure = 0;
            threshold_std = 0;
            
            event_start    = '2018-12-10 05:49:11.000';
            event_end      = '2018-12-10 05:50:05.000';
            
            leading_leftmost_date   = '2018-12-10 05:49:11.500';
            leading_rightmost_date  = '2018-12-10 05:49:14.000';
            trailing_leftmost_date  = '2018-12-10 05:51:45.000';
            trailing_rightmost_date = '2018-12-10 05:53:30.000';
            
            %             leading_start  = '2018-12-10 05:49:24.078';
            %             leading_end    = '2018-12-10 05:49:24.578';
            %
            %             leading_start  = '2018-12-10 05:49:23.492';
            %             leading_end    = '2018-12-10 05:49:23.992';
            
            leading_start  = '2018-12-10 05:49:23.711';
            leading_end    = '2018-12-10 05:49:25.211';
            
            %             trailing_start = '2018-12-10 05:49:45.703';
            %             trailing_end   = '2018-12-10 05:49:48.703';
            
            trailing_start = '2018-12-10 05:49:45.836';
            trailing_end   = '2018-12-10 05:49:48.336';
            
            %% % %Event 90 ['','']
        case 90
            Event_Type = 'HFA';
            Substructure = 1;
            threshold_std = 0;
            
            event_start    = '2018-12-10 06:27:45.000';
            event_end      = '2018-12-10 06:29:15.000';
            
            leading_leftmost_date   = '2018-12-10 06:27:51.000';
            leading_rightmost_date  = '2018-12-10 06:27:53.000';
            
            %             trailing_leftmost_date  = '2018-12-10 06:32:18.000';
            %             trailing_rightmost_date = '2018-12-10 06:34:50.000';
            
            trailing_leftmost_date  = '2018-12-10 06:28:56.000';
            trailing_rightmost_date = '2018-12-10 06:28:57.500';
            
            leading_start  = '2018-12-10 06:28:02.025';
            leading_end    = '2018-12-10 06:28:04.525';
            
            trailing_start = '2018-12-10 06:28:39.588';
            trailing_end   = '2018-12-10 06:28:42.588';
            
            %% % %Event 91  ['','']
        case 91
            Event_Type = 'HFA';
            Substructure = 0;
            threshold_std = 0;
            
            event_start    = '2018-12-10 06:47:55.000';
            event_end      = '2018-12-10 06:48:35.000';
            
            leading_leftmost_date   = '2018-12-10 06:37:40.000';
            leading_rightmost_date  = '2018-12-10 06:44:30.000';
            trailing_leftmost_date  = '2018-12-10 06:48:34.000';
            trailing_rightmost_date = '2018-12-10 06:48:51.000';
            
            leading_start  = '2018-12-10 06:48:00.308';
            leading_end    = '2018-12-10 06:48:02.808';
            
            trailing_start = '2018-12-10 06:48:23.558';
            trailing_end   = '2018-12-10 06:48:26.558';
            %% % %Event 92  ['','']
        case 92
            Event_Type = 'HFA';
            Substructure = 1;
            threshold_std = 0;
            
            event_start    = '2018-12-10 06:50:35.000';
            event_end      = '2018-12-10 06:51:55.000';
            
            leading_leftmost_date   = '2018-12-10 06:50:38.500';
            leading_rightmost_date  = '2018-12-10 06:50:45.000';
            trailing_leftmost_date  = '2018-12-10 06:52:07.500';
            trailing_rightmost_date = '2018-12-10 06:52:15.000';
            
            leading_start  = '2018-12-10 06:50:45.216';
            leading_end    = '2018-12-10 06:50:47.716';
            
            trailing_start = '2018-12-10 06:51:20.084';
            trailing_end   = '2018-12-10 06:51:22.584';
            %% % %Event 93  ['','']
        case 93
            Event_Type = 'HFA';
            Substructure = 0;
            threshold_std = 0;
            
            event_start    = '2018-12-10 07:39:10.000';
            event_end      = '2018-12-10 07:40:25.000';
            
            leading_leftmost_date   = '2018-12-10 07:39:07.000';
            leading_rightmost_date  = '2018-12-10 07:39:12.000';
            trailing_leftmost_date  = '2018-12-10 07:40:17.000';
            trailing_rightmost_date = '2018-12-10 07:40:22.000';
            
            %             leading_start  = '2018-12-10 07:39:33.297';
            %             leading_end    = '2018-12-10 07:39:35.797';
            
            leading_start  = '2018-12-10 07:39:33.961';
            leading_end    = '2018-12-10 07:39:35.961';
            
            trailing_start = '2018-12-10 07:40:11.353';
            trailing_end   = '2018-12-10 07:40:13.353';
            
            %% % %Event 94  ['','']
        case 94
            Event_Type = 'SHFA';
            Substructure = 0;
            threshold_std = 0;
            
            event_start    = '2018-12-10 07:40:20.000';
            event_end      = '2018-12-10 07:41:15.000';
            
            leading_leftmost_date   = '2018-12-10 07:40:25.000';
            leading_rightmost_date  = '2018-12-10 07:40:29.000';
            trailing_leftmost_date  = '2018-12-10 07:41:40.000';
            trailing_rightmost_date = '2018-12-10 07:42:08.000';
            
            leading_start  = '2018-12-10 07:40:38.579';
            leading_end    = '2018-12-10 07:40:40.580';
            
            trailing_start = '2018-12-10 07:41:02.213';
            trailing_end   = '2018-12-10 07:41:04.713';
            
            %% % %Event 95  ['','']
        case 95
            Event_Type = 'HFA';
            Substructure = 1;
            threshold_std = 0;
            
            event_start    = '2018-12-10 08:30:03.000';
            event_end      = '2018-12-10 08:31:52.500';
            
            leading_leftmost_date   = '2018-12-10 08:30:05.500';
            leading_rightmost_date  = '2018-12-10 08:30:10.000';
            trailing_leftmost_date  = '2018-12-10 08:31:50.000';
            trailing_rightmost_date = '2018-12-10 08:31:54.000';
            
            leading_start  = '2018-12-10 08:30:14.708';
            leading_end    = '2018-12-10 08:30:16.208';
            
            trailing_start = '2018-12-10 08:30:29.216';
            trailing_end   = '2018-12-10 08:30:30.716';
            %% % %Event 96
            
            
            
            
            
            
            
        case 96
            Event_Type = 'Foreshock Bubble';
            Substructure = 0;
            threshold_std = 0;
            
            event_start    = '2018-12-30 22:26:43.000';
            event_end      = '2018-12-30 22:28:00.000';
            
            leading_leftmost_date   = '2018-12-30 22:12:00.000';
            leading_rightmost_date  = '2018-12-30 22:26:25.000';
            trailing_leftmost_date  = '2018-12-30 22:28:00.000';
            trailing_rightmost_date = '2018-12-30 22:42:00.000';
            
            leading_start  = '2018-12-30 22:26:54.445';
            leading_end    = '2018-12-30 22:26:58.446';
            
            trailing_start = '2018-12-30 22:27:35.094';
            trailing_end   = '2018-12-30 22:27:39.095';
            
            %% % %Event 97
        case 97
            Event_Type = 'Foreshock Bubble';
            Substructure = 0;
            threshold_std = 0;
            
            event_start    = '2018-12-30 22:26:43.000';
            event_end      = '2018-12-30 22:28:00.000';
            
            leading_leftmost_date   = '2018-12-30 22:12:00.000';
            leading_rightmost_date  = '2018-12-30 22:26:25.000';
            trailing_leftmost_date  = '2018-12-30 22:28:00.000';
            trailing_rightmost_date = '2018-12-30 22:42:00.000';
            
            leading_start  = '2018-12-30 22:26:54.445';
            leading_end    = '2018-12-30 22:26:58.446';
            
            trailing_start = '2018-12-30 22:27:35.094';
            trailing_end   = '2018-12-30 22:27:39.095';
    end
end
%%

% Event_Type = 'Foreshock Bubble'; %Age is too long, expansion speed is very slow
%             Substructure = 1;
%             threshold_std = 0;
%
%             event_start    = '2016-12-08 09:59:50.000';
%             event_end      = '2016-12-08 10:01:05.000';
%
%             leading_leftmost_date   = '2016-12-08 09:59:53.500';
%             leading_rightmost_date  = '2016-12-08 10:00:00.000';
%             trailing_leftmost_date  = '2016-12-08 10:01:02.000';
%             trailing_rightmost_date = '2016-12-08 10:01:07.000';
%
% %             leading_start  = '2016-12-08 10:00:03.785';
% %             leading_end    = '2016-12-08 10:00:04.785';
% %
% %             trailing_start = '2016-12-08 10:00:40.059';
% %             trailing_end   = '2016-12-08 10:00:41.559';
%
%
%             leading_start  = '2016-12-08 10:00:23.223';
%             leading_end    = '2016-12-08 10:00:25.223';
%
%             trailing_start = '2016-12-08 10:00:40.677';
%             trailing_end   = '2016-12-08 10:00:41.677';

%             %% % %Event 78
%         case 78;
%             Event_Type = 'Foreshock Bubble';
%             Substructure = 0;
%             threshold_std = 0;
%
%             event_start    = '2017-01-14 06:49:15.000';
%             event_end      = '2017-01-14 06:50:05.000';
%
%             leading_leftmost_date  = '2017-01-14 06:34:20.000';
%             leading_rightmost_date = '2017-01-14 06:40:20.000';
%             trailing_leftmost_date   = '2017-01-14 06:50:01.000';
%             trailing_rightmost_date  = '2017-01-14 06:51:55.000';
%
%             leading_start  = '2017-01-14 06:49:23.917';
%             leading_end    = '2017-01-14 06:49:25.917';
%
%             trailing_start = '2017-01-14 06:49:43.941';
%             trailing_end   = '2017-01-14 06:49:47.941';
