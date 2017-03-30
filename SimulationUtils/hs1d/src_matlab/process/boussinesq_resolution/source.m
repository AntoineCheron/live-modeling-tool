classdef source
% class storing the forcing due to rainfall (or recharge following the problem)
    properties(Access=public)
        type                % type of source term: steady, periodical, data_based, random, etc.
        time                % time_properties object storing the time properties of the simulation
        recharge_chronicle  % Ntx1 array storing recharge amounts per unit/time
        ETP_chronicle       % Ntx1 array storing ETP amounts per unit/time
        recharge_rate       % Recharge rate if steady or Recharge Amplitude if periodical or square
        period              % period in case of sinusoidal or square recharge
        recharge_mean       % averaged recharge over all the chronicle
    end
    
    methods(Access=public)
        % constructor
        function obj=source(type,period)
            % period in days
            if(nargin<1) type='steady'; end
            if(strcmp(type,'periodical') && nargin<2) period=10; end
            obj.type=type;
            switch type
                case 'steady'
                    obj.period=Inf;
                case 'random'
                    obj.period=0;
                case 'data_based'
                    obj.period=-1;
                case 'periodical'
                    obj.period=period;
                    obj.period=time_unit.time_to_seconds(obj.period,'day');
                case 'square'
                    obj.period=period;
                    obj.period=time_unit.time_to_seconds(obj.period,'day');
                case 'unit_hydrograph'
                    obj.period=-2;
                otherwise
                    fprintf(strcat('The type requested, i.e. ',type,{' '},'does not correspond to a correct recharge type'));
            end
            obj.ETP_chronicle=nan;
        end
        
        function obj=set_recharge_chronicle(obj,t,recharge_rate)
            obj.time=t;
%             obj.time=time_unit.time_to_seconds(obj.time,unit);
            obj.recharge_rate=recharge_rate*1e-3/(3600*24);
            switch obj.period
                case Inf
                    obj.recharge_chronicle=obj.recharge_rate*ones(obj.time.get_size);
                case 0
                    obj.recharge_chronicle=obj.set_random_recharge;
                case -1
                    obj.recharge_chronicle=nan;
                    fprintf('Data has not been uploaded because they are not supplied in arguments with this function. \n Call function upload_recharge_chronicle \n');
                case -2
                    size_time=obj.time.get_size;
                    tmax_total=obj.time.get_tmax;
                    tmax=3600*24;
                    tmax=round(tmax*size_time(2)/tmax_total);
                    obj.recharge_chronicle=0*ones(size_time);
                    obj.recharge_chronicle(1:tmax)=obj.recharge_rate;
                otherwise
                    [t_chronicle,~]=obj.time.get_properties;
                    if(strcmp(obj.type,'periodical'))
                        obj.recharge_chronicle=obj.recharge_rate*(1+cos(2*pi*t_chronicle/obj.period));
                    elseif(strcmp(obj.type,'square'))
                        Int_=floor(t_chronicle/obj.period);
                        Odd_rest=mod(Int_,2);
                        obj.recharge_chronicle=obj.recharge_rate*Odd_rest;
                        % not to have stiff forcing the rain decrease from Recharge_rate to 0 in 60 s
                        rem_stiff=rem(t_chronicle,2*obj.period);
                        Bool1=rem_stiff<60;
                        obj.recharge_chronicle(Bool1)=obj.recharge_rate*(1-(rem_stiff(Bool1))/60); 
                        % not to have stiff forcing either when the rain starts so goes from 0 to Recharge_rate in 60s
                        rem_stiff=rem(t_chronicle+obj.period,2*obj.period);
                        Bool1=rem_stiff<60;
                        obj.recharge_chronicle(Bool1)=obj.recharge_rate*((rem_stiff(Bool1))/60); 
                    end
            end
            obj.recharge_mean=nanmean(obj.recharge_chronicle);
        end
        
        function Recharge=compute_recharge_rate(obj,t)
            switch obj.period
                case Inf
                    Recharge=obj.recharge_rate*ones(size(t));
                case 0
                    [t_chronicle,~]=obj.time.get_properties;
                    Recharge=interpn(t_chronicle,obj.recharge_chronicle,t,'linear');
                case -1
                    [t_chronicle,~]=obj.time.get_properties;
%                     Recharge=interpn(t_chronicle,obj.recharge_chronicle,t,'linear');
                    if(length(t)>1 && length(t)<length(t_chronicle) && sum(t_chronicle(1:length(t))-t)==0)
                        fprintf('WARNING: in compute_recharge_rate method time_results has not the same size as t_chronicle \n');
                        t_chronicle=t;
                    end
                    if(length(t)==length(t_chronicle) || length(t)==1)
                        time_pos=find(t_chronicle-t<0,1,'last');
                        time_neg=find(t-t_chronicle<0,1,'first');
                        time_null=find(t-t_chronicle==0);
                        if(~isempty(time_null))
                            Recharge=obj.recharge_chronicle(time_null);
                        else
                            t1=t_chronicle(time_pos);
                            t2=t_chronicle(time_neg);
                            Recharge1=obj.recharge_chronicle(time_pos); Recharge2=obj.recharge_chronicle(time_neg);
                            Recharge=(Recharge1*(t2-t)+Recharge2*(t-t1))/(t2-t1);
                        end
                    else
                        Recharge=interpn(t_chronicle,obj.recharge_chronicle,t,'linear');
                    end
                case -2
                    t2=obj.time.get_properties;
                    pos_tmax=find(obj.recharge_chronicle>0,1,'last');
                    tmax=t2(pos_tmax);
                    Bool1=t<tmax;
                    Bool2=t-tmax<60;
                    Bool3=Bool2 & ~Bool1;
                    Recharge=obj.recharge_rate.*Bool1+obj.recharge_rate*(1-(t-tmax)/60).*Bool3;
                otherwise
                    if(strcmp(obj.type,'periodical'))
                        Recharge=obj.recharge_rate*(1+cos(2*pi*t/obj.period));
                    elseif(strcmp(obj.type,'square'))
                        Int_=floor(t/obj.period);
                        Odd_rest=mod(Int_,2);
                        Recharge=obj.recharge_rate*Odd_rest;
                        % not to have stiff forcing the rain decrease from Recharge_rate to 0 in 60 s
                        rem_stiff=rem(t,2*obj.period);
                        Bool1=rem_stiff<60; %Bool2=rem_stiff>0; Bool3= Bool1 & Bool2;
                        Recharge(Bool1)=obj.recharge_rate*(1-(rem_stiff(Bool1))/60); 
                        % not to have stiff forcing either when the rain starts so goes from 0 to Recharge_rate in 60s
                        rem_stiff=rem(t+obj.period,2*obj.period);
                        Bool1=rem_stiff<60;
                        Recharge(Bool1)=obj.recharge_rate*((rem_stiff(Bool1))/60);    
                    end
            end
        end
        
        function ETP=compute_ETP_rate(obj,t)
            switch obj.period
                case -1
                    if(isnan(obj.ETP_chronicle))
                        ETP=nan;
                    else
                        [t_chronicle,~]=obj.time.get_properties;
                        %                     Recharge=interpn(t_chronicle,obj.recharge_chronicle,t,'linear');
                        time_pos=find(t_chronicle-t<0,1,'last');
                        time_neg=find(t-t_chronicle<0,1,'first');
                        time_null=find(t-t_chronicle==0);
                        if(~isempty(time_null))
                            ETP=obj.ETP_chronicle(time_null);
                        else
                            t1=t_chronicle(time_pos);
                            t2=t_chronicle(time_neg);
                            ETP1=obj.ETP_chronicle(time_pos); ETP2=obj.ETP_chronicle(time_neg);
                            ETP=(ETP1*(t2-t)+ETP2*(t-t1))/(t2-t1);
                        end
                    end
                otherwise
                    ETP=nan;
            end
        end
        
        function recharge_chronicle=set_random_recharge(obj,t)
            recharge_chronicle=0;
        end
        
        function [t,obj]=upload_recharge_chronicle(obj,Pluvio_directory,ratio_P_R)
            if(nargin<3) ratio_P_R=1; end
            [Date,Pluvio,ETP]=obj.read_data(Pluvio_directory);
            [t,obj]=obj.set_recharge_chronicle_data_based(Date,ratio_P_R,Pluvio,'mm/h',ETP);
        end
        
        function [t,obj]=set_recharge_chronicle_data_based(obj,Date,ratio_P_R,Pluvio,unity,ETP)
            % by default assume pluvio in mm/h 
            if(nargin<5) unity='mm/h'; end
            if(nargin<6) ETP=nan; end
            tmin=Date(1); tmax=Date(end); Nt=length(Date); time_unity_type='day';
            t=time_properties(tmin,tmax,Nt,time_unity_type);
            if(strcmp(unity,'mm/h'))
                obj.recharge_chronicle=ratio_P_R*Pluvio/(1000*3600);
            elseif(strcmp(unity,'mm/d'))
                obj.recharge_chronicle=ratio_P_R*Pluvio/(1000*3600*24);
            elseif(strcmp(unity,'m/s'))
                obj.recharge_chronicle=ratio_P_R*Pluvio;
            else
                fprintf('Problem with Pluvio units to set recharge_chronicle property in source object');
            end
            obj.ETP_chronicle=ETP/(1000*3600);
            obj.time=t;
            obj.recharge_mean=nanmean(obj.recharge_chronicle);
        end
        
%         function [Date_hourly_refined,Pluvio_hourly_refined,ETP_hourly_refined]=read_data(obj,Pluvio_directory)
        function [Date_hourly,Pluvio_hourly,ETP_hourly]=read_data(obj,Pluvio_directory)
%         function [Date_daily_refined,Pluvio_daily_refined,ETP_daily_refined]=read_data(obj,Pluvio_directory)
%         function [Date_daily,Pluvio_daily,ETP_daily]=read_data(obj,Pluvio_directory)
            fid = fopen(Pluvio_directory);

            if fid>0
                % note how we skip the header lines and use the delimiter
                data = textscan(fid,'%s %s %s %s','Delimiter',';');
                
                % close the file
                fclose(fid);
                Pluvio_ETP=nan(length(data{1}),4);
                % date Pluvio Q
                for i=1:length(data{1})
                    if(strcmp(data{2}{i},''))
                        Pluvio_ETP(i,1)=datenum(data{1}{i},'dd/mmm/yyyy');
                    else
                        data{5}{i}=[data{1}{i},' ',data{2}{i}];
                        Pluvio_ETP(i,1)=datenum(data{5}{i},'dd/mmm/yyyy HH:MM:SS');
                    end
                    Pluvio_ETP(i,end)=datenum(data{1}{i},'dd/mmm/yyyy');
                    if(strcmp(data{3}{i},''))
                        Pluvio_ETP(i,2)=nan;
                    else
                        Pluvio_ETP(i,2)=str2num(data{3}{i});
                    end
                    if(strcmp(data{4}{i},''))
                        Pluvio_ETP(i,3)=nan;
                    else
                        Pluvio_ETP(i,3)=str2num(data{4}{i});
                    end
                end
            end
            
            [Y,M,D,H] = datevec (Pluvio_ETP(:,1));
            [c0,~,c1] = unique([Y,M,D,H],'rows');
            Pluvio_hourly=accumarray(c1,Pluvio_ETP(:,2),[],@nansum);
            ETP_hourly=accumarray(c1,Pluvio_ETP(:,3),[],@sum);
            Date_hourly=datenum(c0(:,1),c0(:,2),c0(:,3),c0(:,4),0,0);
            Date_hourly=Date_hourly-Date_hourly(1);
            Date_hourly=Date_hourly';
            Pluvio_hourly=Pluvio_hourly';
            ETP_hourly=ETP_hourly';
            Date_hourly_refined=Date_hourly(1):(Date_hourly(end)-Date_hourly(1))/(10*length(Date_hourly)-1):Date_hourly(end);
            Pluvio_hourly_refined=interpn(Date_hourly,Pluvio_hourly,Date_hourly_refined);
            ETP_hourly_refined=interpn(Date_hourly,ETP_hourly,Date_hourly_refined);
            % to take into account ETP
            Pluvio_hourly_refined=Pluvio_hourly_refined;
            
            [c0,~,c1] = unique([Y,M,D],'rows');
            Pluvio_daily=accumarray(c1,Pluvio_ETP(:,2),[],@nansum);
            ETP_daily=accumarray(c1,Pluvio_ETP(:,3),[],@sum);
            Date_daily=datenum(c0(:,1),c0(:,2),c0(:,3));
            Date_daily=Date_daily-Date_daily(1);
            Date_daily=Date_daily';
            Pluvio_daily=Pluvio_daily';
            ETP_daily=ETP_daily';
            % convert Pluvio_daily in mm/h
            Pluvio_daily=Pluvio_daily/24;
            ETP_daily=ETP_daily/24;
            Date_daily_refined=Date_daily(1):(Date_daily(end)-Date_daily(1)+1)/(240*(length(Date_daily))-1):(Date_daily(end)+1);
            Pluvio_daily_refined=interpn(Date_daily,Pluvio_daily,floor(Date_daily_refined));
            ETP_daily_refined=interpn(Date_daily,ETP_daily,floor(Date_daily_refined));
            Pluvio_daily_refined=Pluvio_daily_refined(1:end-1);
            ETP_daily_refined=ETP_daily_refined(1:end-1);
            Date_daily_refined=Date_daily_refined(1:end-1);
            %
            %     [c0,~,c1] = unique([Y],'rows');
            %     Pluvio_yearly=accumarray(c1,Pluvio_PF(:,2),[],@nansum);
            %     Date_yearly=datenum(c0(:,1));
            ETP_nan_values=sum(isnan(Pluvio_ETP(:,3)));
            if(ETP_nan_values==length(Pluvio_ETP(:,3)))
                ETP_daily_refined=nan;
                ETP_hourly_refined=nan;
                ETP_hourly=nan;
                ETP_daily=nan;
            end
            
        end
        
        function save_recharge_chronicle(obj,file_output,strfolder)
            % save one recharge chronicle
            folder_create(file_output);
            filename=strcat(file_output,'\hydrologic.input');
            [t,unit]=obj.time.get_properties;
            M=nan(length(obj.recharge_chronicle),2); M(:,1)=t'; M(:,2)=(obj.recharge_chronicle)';
            fid = fopen(filename, 'w');
            a=strfind(file_output,'\'); 
            recharge_type=file_output(a(end)+1:end); 
            string_char=sprintf([recharge_type,' \n']);
            string_char2=sprintf(strfolder);
            fprintf(fid, string_char2);
            fprintf(fid, string_char);
            string_char=sprintf('t(s)\trecharge_chronicle(m.s-1)\n');
            fprintf(fid, string_char);
            fclose(fid);
            dlmwrite(filename,M, '-append', 'precision', '%E','delimiter','\t');
            save([file_output,'\source.mat'],'obj');
            obj.plot_save_recharge_chronicle(file_output,t);
            close all;
        end
                
        function plot_save_recharge_chronicle(obj,filename,t)
            figure; hold on;
            plot(t/(24*3600),obj.recharge_chronicle*1000*24*3600,'o-');
            title('recharge chronicle');
            xlabel('time [days]');
            ylabel('infiltration [mm/d]');
            savefig([filename,'\recharge_chronicle.fig']); 
            print([filename,'\recharge_chronicle.png'],'-dpng'); 
        end
    end
    methods(Static)
        function obj_set=generate_hydrologic_forcing(folder_output)
            % generate synthetic hydrologic forcing
            time_unity_type='day';          % unity of time properties
            tmin=0;                         % minimum time for the simulation
            tmax=35;                        % maximum time for the simulation
            Nt=tmax*24;                         % Number of time steps where the solution will be exported (so here 1/day)
            % 1/ steady source terms
            recharge_rate=2.5*0.3;              % in mm/d mean recharge rate
            t=time_properties(tmin,tmax,Nt,time_unity_type);
            recharge=source('steady');
            recharge=recharge.set_recharge_chronicle(t,recharge_rate);
            obj_set{1}=recharge;
            strfolder{1}='Synthetic infiltration data(steady case) \n';
            strfolder2{1}='Synthetic1';
            % 2/ square signal
            recharge_rate=5*0.3;              % in mm/d mean recharge rate
            period=5;                       % period in days
            t=time_properties(tmin,tmax,Nt,time_unity_type);
            recharge=source('square',period);
            recharge=recharge.set_recharge_chronicle(t,recharge_rate);
            obj_set{2}=recharge;
            strfolder{2}='Synthetic infiltration data(square case) \n';
            strfolder2{2}='Synthetic2';
            % 3/ real rainfall 1
            recharge=source('data_based');
            Pluvio_directory=which('dataQPluvioPF3.txt');
            ratio_P_R=0.3;
            [~,recharge]=recharge.upload_recharge_chronicle(Pluvio_directory,ratio_P_R);
            obj_set{3}=recharge;
            strfolder{3}='Real infiltration data(case 1) \n';
            strfolder2{3}='Real1';
            % 4/ real rainfall 2
            recharge=source('data_based');
            Pluvio_directory=which('dataQPluvioPF12.txt');
            [~,recharge]=recharge.upload_recharge_chronicle(Pluvio_directory,ratio_P_R);
            obj_set{4}=recharge;
            strfolder{4}='Real infiltration data(case 2) \n';
            strfolder2{4}='Real2';
            % 5/ real rainfall 3 (1 year)
            recharge=source('data_based');
            Pluvio_directory=which('dataQPluvioPF13.txt');
            [~,recharge]=recharge.upload_recharge_chronicle(Pluvio_directory,ratio_P_R);
            obj_set{5}=recharge;
            strfolder{5}='Real infiltration data(case 3) \n';
            strfolder2{5}='Real3';
            
            
            % save recharge chronicle
            for i=1:length(obj_set)
                folder_output_specific=[folder_output,'\',strfolder2{i}];
                obj_set{i}.save_recharge_chronicle(folder_output_specific,strfolder{i});
            end
        end
        
        function [t,R]=read_recharge_chronicle(source_file_path)
            fid=fopen(source_file_path);
            if(fid>0)
                data_recharge=dlmread(source_file_path,'\t',3,0);
                t=data_recharge(:,1);
                R=data_recharge(:,2);
            else
                fprintf('Cannot open source file where lies recharge chronicle \n');
                t=nan;
                R=nan;
            end
        end
        
        function obj=load_source_class_from_txt_file(source_file_path)
            [t,R]=source.read_recharge_chronicle(source_file_path);
            t=t/(24*3600);
            obj=source('data_based',-1);
            [~,obj]=obj.set_recharge_chronicle_data_based(t,1,R,'m/s',nan);
        end
    end
end