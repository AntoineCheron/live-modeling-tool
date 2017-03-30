classdef parameter_grid_list
    % List of parameters grids
    
    properties(Access=private)
        p       % List of parameters grid
        mult    % Multiplicative option in parameter generation
    end
    
    methods(Access=public)
        
        %% Constructor
        function obj=parameter_grid_list(p,mult)
            if(nargin<2); obj.mult=0; else obj.mult=mult; end
            if(nargin>0)
                obj.p=p;
            else
                obj.p={}; 
            end
        end
        
        %% Addition of parameter to list
        function obj=add(obj,p)
            obj.p={obj.p,p}; 
        end
        
        %% Get parameters grid list
        function p=get_parameter_grid_list(obj)
            p=obj.p;
        end
        
        %% Discretization
        function dis=discretize(obj)
            for i=1:length(obj.p)
                dis{i}=obj.p{i}.discretize_interval; 
            end
        end
        
        %% Cross product 
        function val=cross_product(obj)
            dis=obj.discretize;
            if(length(dis)==1)
                val=dis{1}';
            elseif(length(dis)==2)
                val=obj.cross_product_2param(dis,obj.mult);
            elseif(length(dis)==3)
                val=obj.cross_product_3param(dis);
            elseif(length(dis)==4)
                val=obj.cross_product_4param(dis);
            end
        end
        
        %% Decross product 
        function [param1,param2,param3,val_matrix]=decross_product(obj,ps,val_vector)
            if(length(ps(1,:))==2)
                [param1,param2,param3,val_matrix]=obj.decross_product_2param(ps,val_vector);
            elseif(length(ps(1,:))==3)
                [param1,param2,param3,val_matrix]=obj.decross_product_3param(ps,val_vector);
            end
        end
        
        
    end
    methods(Static)
%         function val=cross_product_1param(dis)
%             val=Nan(1,length(dis{1}));
%             for i=1:length(dis{1})
%                 val(i)=dis{1}(i);
%             end
%         end
        
        function val=cross_product_2param(dis,mult)
            comp=1;
            for i=1:length(dis{1})
                for j=1:length(dis{2})
                    if(mult==0)
                        val(comp,:)=[dis{1}(i),dis{2}(j)];
                    else
                        val(comp,:)=[dis{1}(i),dis{1}(i)*dis{2}(j)];
                    end
                    comp=comp+1;
                end
            end
        end
        
        function val=cross_product_2param_mult(dis)
            comp=1;
            for i=1:length(dis{1})
                for j=1:length(dis{2})
                    val(comp,:)=[dis{1}(i),dis{1}(i)*dis{2}(j)];
                    comp=comp+1;
                end
            end
        end
        
        function val=cross_product_3param(dis)
            comp=1;
            for i=1:length(dis{1})
                for j=1:length(dis{2})
                    for k=1:length(dis{3})
                        val(comp,:)=[dis{1}(i),dis{2}(j),dis{3}(k)];
                        comp=comp+1;
                    end
                end
            end
        end
        
        function val=cross_product_4param(dis)
            comp=1;
            for i=1:length(dis{1})
                for j=1:length(dis{2})
                    for k=1:length(dis{3})
                        for l=1:length(dis{4})
                            val(comp,:)=[dis{1}(i),dis{2}(j),dis{3}(k),dis{4}(l)];
                            comp=comp+1;
                        end
                    end
                end
            end
        end
        
        %% De-cross product to param1 and param2
        function [param1,param2,param3,val_matrix]=decross_product_2param(ps,val_vector)
            n2=length(ps);
            n=sqrt(n2);
            i=1; j=1;
            for comp=1:length(ps)
                param1(i)=ps(comp,1);
                param2(j)=ps(comp,2);
                val_matrix(i,j)=val_vector(comp);
                % fprintf('i=%d\tj=%d\tmu=%f\tlambda=%f\tval=%f\n', i, j, mu(i), lambda(j),val_vector(comp));
                if(mod(comp,n)==0)
                    i=i+1;
                    j=1;
                else
                    j=j+1;
                end
            end
            param3=nan;
        end
        
        %% De-cross product to param1, param2 & param3
        %#JM to debug
        function [param1,param2,param3,val_matrix]=decross_product_3param(ps,val_vector)
            n3=length(ps);
            n=n3^(1/3);
            n2=n^2;
            i=1; j=1; k=1;
            for comp=1:length(ps)
                param1(i)=ps(comp,1);
                param2(j)=ps(comp,2);
                param3(k)=ps(comp,3);
                val_matrix(i,j,k)=val_vector(comp);
                % fprintf('i=%d\tj=%d\tmu=%f\tlambda=%f\tval=%f\n', i, j, mu(i), lambda(j),val_vector(comp));
                if(mod(comp,n2)==0)
                    i=i+1;
                    j=1;
                    k=1;
                elseif(mod(comp,n)==0)
                    j=j+1;
                    k=1;
                else
                    k=k+1;
                end
            end
        end
        
    end
    
    
end






