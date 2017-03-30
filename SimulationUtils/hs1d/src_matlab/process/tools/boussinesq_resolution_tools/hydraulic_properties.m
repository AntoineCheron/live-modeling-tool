classdef hydraulic_properties
    properties(Access=public)
%         i
        f
        k
    end
    
    methods(Access=public)
        function obj=hydraulic_properties(f,k)
%             obj.i=i;
            obj.f=f;
            obj.k=k;
        end
        
        function [f,k]=get_hydraulic_properties(obj)
%             i=obj.i;
            f=obj.f;
            k=obj.k;
        end
    end
end