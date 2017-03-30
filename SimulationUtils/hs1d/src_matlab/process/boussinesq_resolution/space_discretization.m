classdef space_discretization
% class defining the discretization in space. Discretization has not to be regular (linearly spaced)

    properties(Access=private)
        xmin                    % minimum x coordinates [m]
        xmax                    % maximum x coordinates [m]
        discretization_type     % linear, log, square, custom
        x                       % (Nx X 1) edges of the space discretization: on these x coordinates Q is computed [m]
%         dx                      % (Nx X 1) space steps array [m]
        x_S
        w_resampled
        soil_depth_resampled
        angle_resampled
    end
    
    properties(Access=public)
        Nx                      % number of discretization time steps
        A
        B
        C 
        Omega
        Omega2
    end
    
    methods(Access=public)
        function obj=space_discretization(xmin,xmax,Nx,discretization_type,xcustom)
            if(nargin<4) discretization_type='lin'; end
            if(nargin<5) xcustom=nan; end 
            obj.xmin=xmin;
            obj.xmax=xmax;
            obj.Nx=Nx;
            obj.discretization_type=discretization_type;
            switch obj.discretization_type
                case 'lin'
                    obj.x=obj.xmin:(obj.xmax-obj.xmin)/obj.Nx:obj.xmax;
                case 'log'
                    if(obj.xmin==0) 
                        log_x_min=-2;
                        obj.xmin=10^log_x_min;
                    else
                        log_x_min=log10(obj.xmin);
                    end
                    log_x_max=log10(obj.xmax);
                    obj.x=logspace(log_x_min,log_x_max,obj.Nx+1);
                case 'custom'
                    obj.x=xcustom;
                    obj.Nx=length(obj.x)-1;
                    obj.xmin=min(obj.x);
                    obj.xmax=max(obj.x);
                case 'square'
                    obj.x=sqrt(obj.xmin):(sqrt(obj.xmax)-sqrt(obj.xmin))/obj.Nx:sqrt(obj.xmax);
                    obj.x=(obj.x).^2;
                otherwise
                    fprintf(strcat('no discretization type corresponding to ',obj.discretization_type,'  \n'));
            end
            obj.x=obj.x';
%             obj.dx=obj.compute_dx;
            obj=obj.compute_x_centered;
        end
        
        function x=get_edges_coordinates(obj)
            x=obj.x;
        end
        
        function x_S=get_center_coordinates(obj)
            x_S=obj.x_S;
        end
        
        function xmin=get_xmin(obj)
            xmin=obj.xmin;
        end
        
        function xmax=get_xmax(obj)
            xmax=obj.xmax;
        end
        
        function Nx=get_Nx(obj)
            Nx=obj.Nx;
        end
        
        function discretization_type=get_discretization_type(obj)
            discretization_type=obj.discretization_type;
        end
                
        function obj=resample_hs1D_spatial_variables(obj,x,w,soil_depth,angle)
%             obj.w_resampled=interpn(x,w,obj.x_S);
%             obj.soil_depth_resampled=interpn(x,soil_depth,obj.x_S);
%             obj.angle_resampled=interpn(x,angle,obj.x);
%             smooth_width_function = fit(x, w,  'smoothingspline', 'SmoothingParam', 0.0001);
% #JM to change
            smooth_width_function = fit(x, w,  'smoothingspline', 'SmoothingParam', 0.9);
            obj.w_resampled=smooth_width_function(obj.x_S);
%             smooth_slope_function = fit(x, angle,  'smoothingspline', 'SmoothingParam', 0.0001);
% #JM to change
            smooth_slope_function = fit(x, angle,  'smoothingspline', 'SmoothingParam', 0.9);
            obj.angle_resampled=smooth_slope_function(obj.x);
            obj.soil_depth_resampled=interpn(x,soil_depth,obj.x_S);
        end
        
        function [x_S,w_resampled,soil_depth_resampled,angle_resampled,x]=get_resampled_variables(obj)
            x_S=obj.x_S;
            w_resampled=obj.w_resampled;
            soil_depth_resampled=obj.soil_depth_resampled;
            angle_resampled=obj.angle_resampled;
            x=obj.x;
        end
                
        function obj=set_matrix_properties(obj)
            obj.A=obj.first_derivative_downstream;
            obj.B=obj.first_derivative_upstream;
            obj.Omega=obj.weight_matrix;
            obj.Omega2=obj.weight_matrix_bis;
%             obj.C=obj.first_derivative_centered;
        end
        
    end
    
    % #JM ideally change public to private ...
    methods(Access=public)
        
        % length(x_S)=Nx
        function obj=compute_x_centered(obj)
            obj.x_S=(obj.x(2:end)+obj.x(1:end-1))/2;
        end
        
        % length(dx)=Nx
        function dx=compute_dx(obj)
            dx=obj.x(2:end)-obj.x(1:end-1);
            %             dx=obj.x(2:end)-obj.x(1:end-1);
            %             dx(end+1)=dx(end);
%             xtemp=[obj.x(1)-(obj.x(2)-obj.x(1));obj.x;obj.x(end)+(obj.x(end)-obj.x(end-1))];
%             dx=(xtemp(3:end)-xtemp(1:end-2))/2;
        end
        
        % length(dx_S)=Nx-1
        function dx_S=compute_dx_centered(obj)
            dx_S=obj.x_S(2:end)-obj.x_S(1:end-1);
        end
        
%% stencil matrix (for basic derivation centered, downstream, upstream) + mean matrix (Omega)        
        %        [ 0   ...  0   0 ...  0 ]
        %        [ -1   1   0   0 ...  0 ]
        %   B  = [ 0   -1   1   0 ...  0 ]
        %        [ ...   ...    ...   ...]
        %        [ 0  ... ...  0  -1   1 ]
        %        [ 0   ...  0   0 ...  0 ]
        % size B (Nx+1) x (Nx)
        
        function B=first_derivative_upstream(obj)
            A=obj.first_derivative_downstream;
            A=A(1:end-1,1:end-1);
            A(A>0)=1;
            A(A<0)=-1;
            
            B=[zeros(1,length(obj.x)-1);A;zeros(1,length(obj.x)-1)];
            dx_S=obj.compute_dx_centered;
            dx_S=[0;dx_S;0];
            B=diag(1./dx_S)*B;
            B=sparse(B);
        end
        
        %        [ -1  1   0  0   0  ...  0 ]
        %   A  = [ 0  -1   1  0  0  ...   0 ]
        %        [ ...      ...      ...  0 ]
        %        [ 0    ...   ...  0  -1  1 ]
        % size A= (Nx) x (Nx+1)
        function A=first_derivative_downstream(obj)
            A=-diag(ones(length(obj.x)-1,1));
            A=[A,zeros(length(obj.x)-1,1)];
            Abis=diag(ones(length(obj.x)-1,1));
            Abis=[zeros(length(obj.x)-1,1),Abis];
            A=A+Abis;
            dx=obj.compute_dx;
            A=diag(1./dx)*A;
            A=sparse(A);
        end
        
        %        [ 0   1  0   0  ...  0 ]
        %   C  = [-1   0  1  0  ...   0 ]
        %        [ ...    ...    ...  1 ]
        %        [ 0  ... ...  0  -1  0 ]
        function C=first_derivative_centered(obj)
            C=-diag(ones(length(obj.x)-1,1),-1);
            C=C+diag(ones(length(obj.x)-1,1),1);
            C=diag(1./obj.dx)*C;
            C=sparse(C);
        end
                
        %        [0.5 0.5  0   0 ...  0 ]
        % Omega= [ 0  0.5 0.5  0 ...  0 ]
        %        [ ...    ...    ... 0.5]
        %        [ 0  ... ...  0  0  0.5]
        function Omega=weight_matrix(obj)
            Omega=obj.first_derivative_upstream;
            Omega(Omega>0)=0.5;
            Omega(Omega<0)=0.5;
            Omega=sparse(Omega);
        end
        
        function Omega2=weight_matrix_bis(obj)
            Omega2=obj.first_derivative_upstream;
            Omega2(Omega2<0)=0;
            Omega2=sparse(Omega2);
        end
    end

end