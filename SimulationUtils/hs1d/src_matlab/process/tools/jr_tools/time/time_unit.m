classdef time_unit
    % time units
    properties(Access=protected)
        % unit of time
        unit
    end
    
    methods(Access=public)
        % constructor
        function obj=time_unit(unit)
            if(nargin<1); obj.unit='sec'; else obj.unit=unit; end
        end
    end

    methods(Static)
        
        function t=time_to_years(t,time_unit)
            %transforms travel times in time_unit to years
            days_in_year = 365;
            hours_in_year = 24 * days_in_year;
            mins_in_year = 60 * hours_in_year;
            secs_in_year = 60 * mins_in_year;
            
            if (strcmp (time_unit, 'sec'))
                t(:) = t(:)/ secs_in_year;
            elseif (strcmp (time_unit, 'min'))
                t(:) = t(:)/ mins_in_year;
            elseif (strcmp (time_unit, 'hour'))
                t(:) = t(:)/ hours_in_year;
            elseif (strcmp (time_unit, 'day'))
                t(:) = t(:)/ days_in_year;
            elseif (strcmp (time_unit, 'year'))
                t(:) = t(:)/ 1 ;
            else
                fprintf('UNIT UNKNOWN %s\n', time_unit); 
            end
        end
        
        function t=time_to_seconds(t,time_unit)
            %transforms travel times in time_unit to seconds
            secs_in_mins = 60;
            secs_in_hours = 60 * secs_in_mins;
            secs_in_days = 24 * secs_in_hours;
            secs_in_year = 365 * secs_in_days;
            
            if (strcmp (time_unit, 'year'))
                t(:) = t(:) * secs_in_year;
            elseif (strcmp (time_unit, 'day'))
                t(:) = t(:) * secs_in_days;
            elseif (strcmp (time_unit, 'hour'))
                t(:) = t(:) * secs_in_hours;
            elseif (strcmp (time_unit, 'min'))
                t(:) = t(:) * secs_in_mins;
            elseif (strcmp (time_unit, 'sec'))
                t(:) = t(:) * 1 ;
            else
                fprintf('UNIT UNKNOWN %s\n', unit); 
            end
        end
        
        function t=time_to_days(t,time_unit)
            %transforms travel times in time_unit to days
            days_in_year = 365;
            hours_in_day = 24;
            mins_in_day = 60 * hours_in_day;
            secs_in_day = 60 * mins_in_day;
            
            if (strcmp (time_unit, 'sec'))
                t(:) = t(:)/ secs_in_day;
            elseif (strcmp (time_unit, 'min'))
                t(:) = t(:)/ mins_in_day;
            elseif (strcmp (time_unit, 'hour'))
                t(:) = t(:)/ hours_in_day;
            elseif (strcmp (time_unit, 'day'))
                t(:) = t(:) / 1;
            elseif (strcmp (time_unit, 'year'))
                t(:) = t(:)* days_in_year ;
            else
                fprintf('UNIT UNKNOWN %s\n', time_unit); 
            end
        end
    end


end