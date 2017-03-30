    % Etude de convergence spatiale (comme avant)
% Temporelle 

classdef runs
% class storing the parametrization to run a boussinesq simulation over a given hillslope of on a parametrized hillslope
    properties(Access=public)
        hs1D
        boussinesq_simulation
        simulation_results
    end
    
    methods(Static)
        function obj=run_direct(hs1D)
            obj=runs;
            if(nargin<1) 
                [obj,discretization]=obj.choose_structure_parameters;
                % change hydraulic parameters
                obj=obj.choose_hydraulic_parameters;
            else
                [obj,discretization]=obj.choose_structure_parameters(hs1D);
                % set hydraulic parameters got in hs1D
                [f,k]=hs1D.get_hydraulic_properties;
                obj.hs1D=obj.change_hillslope_hydraulic_properties(f,k*3600);
            end
            % initial conditions
            percentage_loaded=0.5;            % how the hillslope is filled at time t=0. percentage_loaded is the ratio: 0.5 -> filled at 50 % of its maximum capacity.
            % set the source (rain, infiltration, recharge) options
            [t,source_term]=obj.set_source_options;
            % set the boundary conditions
            boundary_cond=obj.set_boundary_conditions(discretization);
            % set the solver options
            solver_options=obj.set_solver_options;
            % initialize a boussinesq simulation
            obj=obj.initialize_boussinesq_simulation(discretization,source_term,boundary_cond,percentage_loaded,t.get_tmin);
            % Run the DAE
            [t,x_center,x_edge,S,Q,QS,obj.boussinesq_simulation]=obj.boussinesq_simulation.implicit_scheme_solver(t,solver_options);
            obj.simulation_results=simulation_results(t,x_center,x_edge,S,Q,QS);
            % Postprocessing options
%             [time_travel,time_ref]=obj.postprocess_results(t);  
%             obj.postprocess_results(t); 
%             obj.save_results;
        end
        
        function obj=run_simulation(hs1D,source_terms,percentage_loaded,solver_options,Sinitial)
            if(nargin<5)
                Sinitial=nan;
            end
            obj=runs;
            [obj,discretization]=obj.choose_structure_parameters(hs1D);
            obj=obj.choose_hydraulic_parameters(hs1D.f,hs1D.k);
            t=source_terms.time;
            % set the boundary conditions
            boundary_cond=obj.set_boundary_conditions(discretization);
            % initialize a boussinesq simulation
            obj=obj.initialize_boussinesq_simulation(discretization,source_terms,boundary_cond,percentage_loaded,t.get_tmin,Sinitial);
            % Run the DAE
            [t,x_center,x_edge,S,Q,QS,obj.boussinesq_simulation]=obj.boussinesq_simulation.implicit_scheme_solver(t,solver_options);
            obj.simulation_results=simulation_results(t,x_center,x_edge,S,Q,QS);
        end
    end
    
    methods(Access=public)
        function obj=runs
            obj.hs1D=[];
        end
        
        function [obj,discretization]=choose_structure_parameters(obj,hs1D)
            % Parametrization
            % geomorphologic properties
            if(nargin<2)
                x=[0,100];
                angle=0.05;
                w=[500,500];
%                 x=linspace(-100,100,100); %x=[0,200];          % in m
%                 angle=0.05.*(x>0)-0.05.*(x<0); %angle=0.05;    % in % 
%                 w=100*ones(size(x)); %w=[100,900];            % in m
                soil_depth=1;       % in m
                obj.hs1D=obj.set_hillslope_geomorphologic_properties(x,angle,w,soil_depth);
            else
                obj.hs1D=hs1D;
                [x,w,soil_depth,angle]=obj.hs1D.get_spatial_properties;
                obj.hs1D=obj.set_hillslope_geomorphologic_properties(x,angle,w,soil_depth);
            end
            % spatial discretization
            Nx=100;                         % Number of points in the discretization scheme so Nx points to compute S & qS and (Nx+1) points to compute Q
            discretization_type='lin';      % Type of discretization can be 'custom' or 'log'. if 'custom' type is chosen then you need to provide your own discretization
            % Instantiation
            % space discretization object
            discretization=obj.set_space_discretization(Nx,discretization_type);
        end
        
        function obj=choose_hydraulic_parameters(obj,f,k,d,i)
            if(nargin<5)
                i=nan; if(nargin<4) 
                d=2; if(nargin<3)
                k=1; if(nargin<2)
                f=0.3;
                     end
                     end
                      end
            end
            if(~isnan(i))
                [x,w,~,~]=obj.hs1D.get_spatial_properties;
                obj.hs1D=obj.set_hillslope_geomorphologic_properties(x,i,w,d);
            end
            obj.hs1D=obj.change_hillslope_hydraulic_properties(f,k);
        end
        
        function boundary_cond=set_boundary_conditions(obj,discretization)
            [f,~]=obj.hs1D.get_hydraulic_properties;
            [~,w_resampled,soil_depth_resampled,~]=discretization.get_resampled_variables;
            % boundary conditions
            boundary_type={'S','Q'};                                     % boundary type of the problem for xmin & xmax if 'S' Dirichlet condition, S is fixed, if 'Q' Neuman condition, Q is fixed
            % if saturated at the river bank boundary_value=[f*w(1)*soil_depth(1),0]; else : boundary_value=[0,0];
            boundary_value=[f*w_resampled(1)*soil_depth_resampled(1),0];           % boundary value at the edges for the value prescripted in boundary_type
%             boundary_value=[0,0]; 
            % boundary conditions
            boundary_cond=boundary_conditions(boundary_type,boundary_value);
        end
        
        function [t,source_term]=set_source_options(obj)
            recharge_type='square';         % type of recharge can also be 'periodical', 'data_based', 'square' or 'random'
            % source object for the source terms
            if(strcmp(recharge_type,'data_based'))
                Pluvio_file=input('Enter the name of the Pluvio file: \n');
                Pluvio_directory=which(Pluvio_file);
                [t,source_term]=obj.set_data_based_source_terms(recharge_type,Pluvio_directory);
            else
                % choose source terms
                time_unity_type='day';          % unity of time properties
                tmin=0;                         % minimum time for the simulation
                tmax=35;                        % maximum time for the simulation
                Nt=35*24*10;                         % Number of time steps where the solution will be exported (so here 10/day)
                recharge_rate=30;              % in mm/d mean recharge rate
                period=5;                       % period in days
                t=time_properties(tmin,tmax,Nt,time_unity_type);
                source_term=obj.set_source_terms(t,recharge_rate,recharge_type,period);
            end
        end
        
        function odeset_struct=set_solver_options(obj,odeset_struct)
            if(nargin<2)
                % Solver options ode15s used by default
                maxstep=nan;                    % maximum time steps allowed to do the integration in hours
                RelTol=3e-14;                   % Relative Tolerance for the solution
                AbsTol=1e-10;
                odeset_struct=odeset('maxstep',maxstep,'RelTol',RelTol,'AbsTol',AbsTol);
            end
            % solver options
%             solver_options=ode_options;
%             solver_options=solver_options.instantiate(maxstep,RelTol,AbsTol);
        end
        
        function obj=initialize_boussinesq_simulation(obj,discretization,source_term,boundary_cond,percentage_loaded,t_initial,Sinitial)
            if(nargin<7)
                Sinitial=nan;
            end
            obj.boussinesq_simulation=boussinesq_simulation;
            obj.boussinesq_simulation=obj.boussinesq_simulation.set_hydraulic_properties(obj.hs1D);
            obj.boussinesq_simulation=obj.boussinesq_simulation.simulation_parametrization(discretization,source_term,boundary_cond,percentage_loaded,t_initial,Sinitial);
        end
        
%         function [time_travel,time_ref]=postprocess_results(obj,t)
%         function [t_decrease,QS_decrease,dQSdt_decrease,Q_decrease,dQdt_decrease]=postprocess_results(obj,t)
        function postprocess_results(obj,t)
            S_max=obj.compute_S_max;
            [Watershed_Area,w,Watershed_area_spatialized]=obj.compute_watershed_area;
            Flux_in_total=obj.compute_flux_in_total(t);
            Flux_in=obj.compute_flux_in(t);
            obj.simulation_results.plot_results(S_max,Watershed_Area,w,Flux_in_total);
            yy=[obj.simulation_results.S;obj.simulation_results.Q;obj.simulation_results.QS];
%             [mass_balance,Flux_in_spatialized,Darcy_Flux_difference_spatialized,Seepage_Flux_spatialized,dS]=obj.boussinesq_simulation.compute_raw_mass_balance;
            dS=obj.boussinesq_simulation.compute_raw_storage_variation(t,yy);
            [Mass_balance_tot_mmd,Flux_in_tot_mmd,Darcy_Flux_tot_mmd,Seepage_tot_mmd,Storage_Var_mmd]=obj.simulation_results.mass_changes_study(Watershed_area_spatialized,Flux_in,dS);
            obj.simulation_results.plot_mass_balance(Flux_in_tot_mmd,Darcy_Flux_tot_mmd,Seepage_tot_mmd,Storage_Var_mmd,Mass_balance_tot_mmd);
%             obj.simulation_results.generate_storage_video(S_max,w,Watershed_Area,Watershed_area_spatialized,Flux_in,dS);
%             [time_travel,time_ref,N_in]=obj.simulation_results.sample_travel_time(Flux_in);
%             weight=round(N_in/min(N_in));
%             time_travel=repelem(time_travel,weight);
%             time_ref=repelem(time_ref,weight);
            % specify time after rain to take the curve
%             time_after_rain=0.1; % in days
%             [t_decrease,QS_decrease,dQSdt_decrease,Q_decrease,dQdt_decrease]=obj.simulation_results.get_decrease_parts(Flux_in,time_after_rain);
        end
        
        function save_parametrization(obj,foldername)
            [x_S,w_resampled,soil_depth_resampled,angle_resampled,x_Q]=obj.boussinesq_simulation.discretization.get_resampled_variables;
            angle_resampled=interpn(x_Q,angle_resampled,x_S);
            [f,k]=obj.hs1D.get_hydraulic_properties;
            filename=strcat(foldername,'\parameters.para');
            f=f*ones(size(x_S)); k=k*ones(size(x_S));
            M=nan(length(x_S),6); M(:,1)=x_S; M(:,2)=w_resampled; M(:,3)=soil_depth_resampled; M(:,4)=angle_resampled; M(:,5)=f; M(:,6)=k;
            fid = fopen(filename, 'w');
            string_char=sprintf('x\tw\td\ti\tf\tk\n');
            fprintf(fid, string_char);
            fclose(fid);
            dlmwrite(filename,M, '-append', 'precision', '%E','delimiter','\t');
        end
        
        function hs1D=set_hillslope_geomorphologic_properties(obj,x,angle,w,soil_depth)
            if(nargin<5)
                soil_depth=ones(size(x));
            elseif(length(soil_depth)==1)
                soil_depth=soil_depth*ones(size(x));
            elseif(length(soil_depth)~=length(x))
                fprintf('no consistency betweeen x and soil depth information \n');
                fprintf('uniform depth of 1 m will be assumed \n');
                soil_depth=ones(size(x));
            end
            if(nargin<4)
                w=ones(size(x));
            elseif(length(w)~=length(x))
                fprintf('no consistency betweeen x and width function w \n');
                fprintf('uniform hillslopes of width 1 will be assumed \n');
                w=ones(size(x));
            end
            if(nargin<3)
                angle=ones(size(x));
            elseif(length(angle)==1)
                angle=angle*ones(size(x));
            elseif(length(angle)~=length(x))
                fprintf('no consistency betweeen x and angle information i \n');
                fprintf('uniform slope of 0.05 \n');
                angle=0.05*ones(size(x));
            end
            switch length(x)
                case 1
                    fprintf('not enough information to characterize the length of your hillslope \n you need to inform at least xmin and xmax \n');
                case 2
                    % resample your parameters linearly
                    xint=(x(1):(x(2)-x(1))/100:x(2))';
                    wint=interpn(x,w,xint,'linear');
                    soil_depthint=interpn(x,soil_depth,xint,'linear');
                    angleint=interpn(x,angle,xint,'linear');
                case 0
                    fprintf('x cannot be empty! \n');
                otherwise
                    xint=x;
                    wint=w;
                    soil_depthint=soil_depth;
                    angleint=angle;
            end
            hs1D=hillslope1D(1);
            hs1D=hs1D.set_spatial_parameters(xint,wint,angleint,soil_depthint);
        end
        
        function hs1D=change_hillslope_hydraulic_properties(obj,f,k)
            hs1D=obj.hs1D.set_hydraulic_parameters(k,f);
            fprintf('WARNING: hydraulic properties of the hillslope have been changed: \n');
            fprintf(['porosity becomes equal to ' num2str(f) ' \n']);
            fprintf(['hydraulic conductivity becomes equal to ' num2str(k) ' m/h \n']);
        end
        
        function discretization=set_space_discretization(obj,Nx,type,xcustom)
            if(nargin<4) xcustom=-1; end
            if(nargin<3) type='lin'; end
            if(nargin<2) Nx=100; end
            [x,w,soil_depth,angle]=obj.hs1D.get_spatial_properties;
            xmin=x(1); xmax=x(end);
%             type='custom';
%             xcustom=[x(1):1:x(2)-1,x(2):10:x(end)];
            discretization=space_discretization(xmin,xmax,Nx,type,xcustom);
            discretization=discretization.resample_hs1D_spatial_variables(x,w,soil_depth,angle);
            discretization=discretization.set_matrix_properties;
        end
        
        function source_terms=set_source_terms(obj,t,recharge_rate,type,period)
            if(nargin<4) type='steady'; end
            if(nargin<5) period=10; end
            
            source_terms=source(type,period);
            source_terms=source_terms.set_recharge_chronicle(t,recharge_rate);
        end
        
        function [t,source_terms]=set_data_based_source_terms(obj,type,Pluvio_directory)
            if(~strcmp(type,'data_based')) 
                fprintf('Select data_based type for source terms if you want to upload your own recharge/pluvio chronicle \n'); 
                return 
            end
            source_terms=source(type);
            [t,source_terms]=source_terms.upload_recharge_chronicle(Pluvio_directory);
        end
        
        
        function S_max=compute_S_max(obj)
            [f,~]=obj.hs1D.get_hydraulic_properties;
            [~,w,soil_depth,~]=obj.boussinesq_simulation.discretization.get_resampled_variables;
            S_max=f*w.*soil_depth;
        end
        
        function [Watershed_area,w,Watershed_area_spatialized]=compute_watershed_area(obj)
            [~,w,~]=obj.boussinesq_simulation.discretization.get_resampled_variables;
            dx=obj.boussinesq_simulation.discretization.compute_dx;
            Watershed_area_spatialized=w.*dx;
            Watershed_area=sum(Watershed_area_spatialized);
        end
        
        function N_in_total=compute_flux_in_total(obj,t)
            N_in_total=obj.boussinesq_simulation.compute_recharge_total(t);
        end
        
        function N_in=compute_flux_in(obj,t)
            N_in=obj.boussinesq_simulation.compute_recharge(t);
        end
        
        function ETR_OUT=compute_ETR(obj,t,S)
            ETR_OUT=obj.boussinesq_simulation.compute_ETR_OUT(t,S)
        end
        
        function Sfin=get_final_storage(obj)
            Sfin=obj.simulation_results.S;
            Sfin=Sfin(:,end);
        end
        
        function state_values_fin=get_final_state_values(obj)
            Sfin=obj.simulation_results.S;
            Sfin=Sfin(:,end);
            Qfin=obj.simulation_results.Q;
            Qfin=Qfin(:,end);
            QSfin=obj.simulation_results.QS;
            QSfin=QSfin(:,end);
            state_values_fin=[Sfin;Qfin;QSfin];
        end
        
        function [S_max,Watershed_area,w,Flux_in_total,Watershed_area_spatialized]=get_key_parameters(obj,t)
            S_max=obj.compute_S_max;
            [Watershed_area,w,Watershed_area_spatialized]=obj.compute_watershed_area;
            Flux_in_total=obj.compute_flux_in_total(t);
        end
        
        function [S_max,Watershed_area,w,Flux_in_total,Watershed_area_spatialized]=save_key_parameters(obj,t,folder_output)
            [S_max,Watershed_area,w,Flux_in_total,Watershed_area_spatialized]=obj.get_key_parameters(t);
            key_param=[];
            key_param.S_max=S_max;
            key_param.Watershed_area=Watershed_area;
            key_param.w=w;
            key_param.Flux_in_total=Flux_in_total;
            key_param.Watershed_area_spatialized=Watershed_area_spatialized;
            x=obj.boussinesq_simulation.discretization.get_center_coordinates;
            key_param.x=x;
            filename_spatialized_key_param=strcat(folder_output,'\key_spatialized_param.param');
            fid = fopen(filename_spatialized_key_param, 'w');
            fprintf(fid, 'First column: x, second column: Smax [m2], third column: w [m], fourth column: watershed_area_spatialized [m2] \n');
            fclose(fid);            
            S=nan(length(x),4);
            S(:,1)=x;
            S(:,2)=S_max;
            S(:,3)=w;
            S(:,4)=Watershed_area_spatialized;
            dlmwrite(filename_spatialized_key_param,S, '-append', 'precision', '%E','delimiter','\t');
            
            save(strcat(folder_output,'\key_parameters.mat'),'key_param');
        end
        
        function save_sim_run(obj,folder_output,option_space_limited)
            hs1D=obj.hs1D;
            bouss_sim=obj.boussinesq_simulation;
            save(strcat(folder_output,'\hs1D.mat'),'hs1D');
            if(nargin>=3 && option_space_limited==2)
                bouss_sim.sol_simulated=[];
                save(strcat(folder_output,'\boussinesq_simulation.mat'),'bouss_sim');
            elseif(nargin>=3 && option_space_limited==1)
                filename_sol_simulated=strcat(folder_output,'\sol_simulated.code');
                fid = fopen(filename_sol_simulated, 'w');
                fprintf(fid, 'First row: t (sol_simulated.x), Next rows: state vector outputs from ode15s (sol_simulated.y) \n');
                fclose(fid);
                sol_simulated=nan(size(bouss_sim.sol_simulated.y,1)+1,size(bouss_sim.sol_simulated.y,2));
                sol_simulated(1,:)=bouss_sim.sol_simulated.x;
                sol_simulated(2:end,:)=bouss_sim.sol_simulated.y;
                dlmwrite(filename_sol_simulated,sol_simulated, '-append', 'precision', '%E','delimiter','\t');
                bouss_sim.sol_simulated=[];
                save(strcat(folder_output,'\boussinesq_simulation.mat'),'bouss_sim');
            else
                save(strcat(folder_output,'\boussinesq_simulation.mat'),'bouss_sim');
            end
        end
    end
    
    methods(Static)
        function [x,S_max,width_function,Watershed_area_spatialized]=load_key_parameters(key_param_file_path)
            fid=fopen(key_param_file_path);
            if(fid>0)
                data_key_param=dlmread(key_param_file_path,'\t',1,0);
                x=data_key_param(:,1);
                S_max=data_key_param(:,2);
                width_function=data_key_param(:,3);
                Watershed_area_spatialized=data_key_param(:,4);
            else
                fprintf('Cannot open key_param file where lies key parameters \n');
                x=nan;
                S_max=nan;
                width_function=nan;
                Watershed_area_spatialized=nan;
            end
        end
    end

end