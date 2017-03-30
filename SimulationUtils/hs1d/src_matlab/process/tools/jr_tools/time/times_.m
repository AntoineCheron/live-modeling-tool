classdef times_ < time_unit
    % times
    properties(Access=private)
        t
    end
    
    methods(Access=public)
        
        % constructor
        function obj=times_(t,unit)
            if(nargin<1); obj.t=[]; obj.unit='sec'; else obj.t=t; obj.unit=unit; end
        end
        
        % Time in years (whatever the current unit)
        function t=time_in_years(obj)
            t=time_unit.time_to_years(obj.t,obj.unit); 
        end
        
        % Time in seconds (whatever the current unit)
        function t=time_in_seconds(obj)
            t=time_unit.time_to_seconds(obj.t,obj.unit); 
        end
        
        % Time in the unit given in parameter
        function obj=time_in_unit(obj,unit_new)
            if(strcmp(unit_new,'sec'))
                obj.t=obj.time_in_seconds;
                obj.unit='sec'; 
            elseif(strcmp(unit_new,'year'))
                obj.t=obj.time_in_years;
                obj.unit='year'; 
            else
                fprintf('unsupported time option %s\n', obj1.unit);
            end
        end
        
        % Gets time differential (with the same dimension)
        function dy=dt(obj)
            dy=diff(obj.t); 
            %#JR dy(end+1)=dy(end); 
        end
        
        % Function "+"
        function r = plus(obj1,obj2)
            if(strcmp(obj1.unit,obj2.unit)==0)
                obj2=obj2.time_in_unit(obj1.unit);
            end
            r = times_(obj1.t+obj2.t,obj1.unit);
        end
        
        % Function "-"
        function r = minus(obj1,obj2)
            if(strcmp(obj1.unit,obj2.unit)==0)
                obj2=obj2.time_in_unit(obj1.unit);
            end
            r = times_(obj1.t-obj2.t,obj1.unit);
        end
        
    end
    
    methods(Static)
        % Test function 
        function test
            t1=times_(10,'year'); 
            t2=times_(5,'year'); 
            u=t1+t2; 
            u2=t1-t2;
            t3=times_(10,'year'); 
            t4=times_(time_unit.time_to_seconds(5,'year'),'sec'); 
            v=t3+t4; 
            v2=t3-t4;
            t5=times_(10,'year');
            t6=times_(time_unit.time_to_seconds(5,'year'),'sec'); 
            w=t5+t6; 
            w2=t5-t6; 
            w=times_(w.time_in_years,'year'); 
            w2=times_(w2.time_in_years,'year'); 
            if(u.t==v.t); fprintf('success\n'); else fprintf('failure\n'); end 
            if(u.t==w.t); fprintf('success\n'); else fprintf('failure\n'); end 
            if(u2.t==v2.t); fprintf('success\n'); else fprintf('failure\n'); end 
            if(u2.t==w2.t); fprintf('success\n'); else fprintf('failure\n'); end 
        end
    end
    
end