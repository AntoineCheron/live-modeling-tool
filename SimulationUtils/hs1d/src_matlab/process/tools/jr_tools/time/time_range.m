classdef time_range
    % Time Range and discretization
    properties(Access=private)
        % minimum time
        tmin
        % maximum time
        tmax
        % number of time 
        nt
    end
    
    methods(Access=public)
        % constructor
        function obj=time_range(tmin,tmax,nt)
            % Default parameter values
            if(nargin<1)
                tmin=0.001; 
                tmax=100; 
                nt=200;
            end
            
            % Sets time units to seconds
            obj.tmin=tmin; 
            obj.tmax=tmax;
            obj.nt=nt; 
        end
        
        % Gets time
        function y=t(obj)
            % Decompositions of the last century in years
            y=obj.tmin:(obj.tmax-obj.tmin)/obj.nt:obj.tmax; % y=y'; 
        end
        
        % Gets time differential (with the same dimension)
        function dy=dt(obj)
            dy=diff(obj.t); 
            %#JR dy(end+1)=dy(end); 
        end
    end

 


end