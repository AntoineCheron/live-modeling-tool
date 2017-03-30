classdef parameter_grid
    % Grid of paraemeters 
    
    properties(Access=private)
        min_    % Minimal parameter value
        max_    % Maximal parameter value
        n       % Discretization of parameters
        lin_log % Scale lin or log or customed
    end
    
    properties (Constant)
        % Scale type
        lin=0;
        log=1;
        custom=-1;
    end
    
    methods(Access=public)
        
        %% Constructor
        function obj=parameter_grid(min_,max_,n,lin_log)
            obj.min_=min_; 
            obj.max_=max_; 
            obj.n=n;
            obj.lin_log=lin_log; 
        end
        
        function lin_log=get_lin_log(obj)
            lin_log=obj.lin_log;
        end
        
        %% Discretization 
        function val=discretize_interval(obj)
            % Discretize interval [min_,max] in n regular spaces
            %   lin_log='log': does it on a logarithmic basis
            
            if(obj.lin_log==parameter_grid.lin)
                val=obj.min_:(obj.max_-obj.min_)/obj.n:obj.max_;
            elseif(obj.lin_log==parameter_grid.custom)
                val=[obj.min_,obj.max_];
            else
                obj.min_=log10(obj.min_); obj.max_=log10(obj.max_);
                val=obj.min_:(obj.max_-obj.min_)/obj.n:obj.max_;
                val=10.^(val);
            end
            
        end
        
    end
    
    
end