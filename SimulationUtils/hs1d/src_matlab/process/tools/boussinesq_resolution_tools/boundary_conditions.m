classdef boundary_conditions
    properties
        type    % 2 value string array defining the type (S, Q or free) of the boundary conditions for xmin et xmax 
        value   % 2 value array giving the value of the boundary condition (Slim or Qlim)
    end
    
    methods(Access=public)
        function obj=boundary_conditions(type,value)
            obj.type=type;
            obj.value=value;
        end
        
        function Edges=fixed_edge_matrix_boolean(obj)
            Edges=zeros(4,1);
            if(strcmp(obj.type{1},'S'))
                Edges(1)=1;
                Edges(3)=1;
            elseif(strcmp(obj.type{1},'Q'))
                Edges(3)=1;
            end
            if(strcmp(obj.type{2},'S'))
                Edges(2)=1;
                Edges(4)=1;
            elseif(strcmp(obj.type{2},'Q'))
                Edges(4)=1;
            end
        end
        
        function Edges=fixed_edge_matrix_values(obj)
            Edges=zeros(4,1);
            if(strcmp(obj.type{1},'S'))
                Edges(1)=obj.value(1);
            elseif(strcmp(obj.type{1},'Q'))
                Edges(3)=obj.value(1);
            end
            if(strcmp(obj.type{2},'S'))
                Edges(2)=obj.value(2);
            elseif(strcmp(obj.type{2},'Q'))
                Edges(4)=obj.value(2);
            end
        end
    end
end