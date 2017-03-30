classdef initial_conditions
    properties(Access=private)
        percentage_loaded
        S
        Q
        QS
    end
    
    methods(Access=public)
        function obj=initial_conditions(Smax,percentage_loaded)
            obj.percentage_loaded=percentage_loaded;
            obj=obj.set_initial(Smax);
        end
    end
    
    methods(Access=private)
        function obj=set_initial(obj,Smax)
            % initial conditions at time t=0
            obj.Sin=obj.percentage_loaded*Smax;
            obj.Qin=obj.compute_Q_from_S(obj.Sin,0);
            obj.QSin=zeros(length(obj.x),1);
        end
    end
end