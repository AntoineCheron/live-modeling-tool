classdef ode_options
    properties(Access=private)
        maxstep 
        RelTol
        AbsTol
    end
    
    methods(Access=public)
        % constructor
        function obj=ode_options
        end
        
        function obj=instantiate(obj,maxstep,RelTol,AbsTol)
            obj.maxstep=maxstep;
            obj.RelTol=RelTol;
            obj.AbsTol=AbsTol;
        end
        
        function [maxstep,RelTol,AbsTol]=get_parameters(obj)
            maxstep=obj.maxstep;
            RelTol=obj.RelTol;
            AbsTol=obj.AbsTol;
        end
    end
end