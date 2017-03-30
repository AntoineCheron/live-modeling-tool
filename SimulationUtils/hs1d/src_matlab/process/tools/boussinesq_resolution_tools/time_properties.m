classdef time_properties
    properties(Access=private)
        tmin
        tmax
        Nt
        t
        unit
    end
    
    methods(Access=public)
        function obj=time_properties(tmin,tmax,Nt,unit)
            tmin=time_unit.time_to_seconds(tmin,unit);
            tmax=time_unit.time_to_seconds(tmax,unit);
            
            obj.tmin=tmin;
            obj.tmax=tmax;
            obj.Nt=Nt;
            obj.t=obj.tmin:(obj.tmax-obj.tmin)/(Nt-1):obj.tmax;
            obj.unit=unit;
        end
    
        function tmin=get_tmin(obj)
            tmin=obj.tmin;
        end
        
        function tmax=get_tmax(obj)
            tmax=obj.tmax;
        end
        
        function unit=get_unit(obj)
            unit=obj.unit;
        end
        
        function Nt=get_Nt(obj)
            Nt=obj.Nt;
        end
        
        function size_=get_size(obj)
            size_=size(obj.t);
        end
        
        function [t,unit]=get_properties(obj)
            t=obj.t;
            unit=obj.unit;
        end
    end
end