classdef reservoir_simulation
    properties(Access=public)
        s                   % aperture surface at the outlet of the reservoir
        S                   % Surface that "see" the recharge
        hmax                % maximum height of the reservoir after it is spilling out the reservoir (SEOF)
        source_terms        % source object containing the time properties and the recharge related to it
        initial_conditions
        sol_simulated
    end
    
    methods(Access=public)
        % constructor
        function obj=reservoir_simulation(s,S)
            obj.s=s;
            obj.S=S;
        end
        
         % set the parameters for the resolution of the Boussinesq DAE
        function [t,x_center,x_edge,S,Q,QS,obj]=implicit_scheme_solver(obj,t_properties)
            % initial vector for the DAE
            y0=100;

            % Function to integrate
            f=@(t,y,yp)(obj.odefun(t,y,yp));
            
             % time range to solve
            [time_range,unit]=t_properties.get_properties;
            
            % Initial slope of the DAE
            yp0=0;
                        
            % get consistent initial conditions thanks to decic
            [y0new,yp0new] = decic(f, time_range(1), y0, [], yp0, []);
            options=odeset('InitialSlope',yp0new);

            obj.sol_simulated=ode15i(f,time_range,y0new,yp0new,options);
            
            t=time_range;
            if(time_range(end)>obj.sol_simulated.x(end))
                tmax_pos=find(time_range-obj.sol_simulated.x(end)<=0,1,'last');
                t=time_range(1:tmax_pos);
            end
            S=deval(obj.sol_simulated,t);
            Q=S(block_size+1:2*block_size+1,:);
            QS=S(2*block_size+2:3*block_size+1,:);
            S=S(1:block_size,:);
            x_edge=obj.discretization.get_edges_coordinates;
            x_center=obj.discretization.get_center_coordinates;
        end
        
        function f=odefun(obj,t,y,yp)
            R=obj.source_terms.compute_recharge_rate(t);
            g=9.81;
            f=yp+obj.s/obj.S*sqrt(yp.^2+2*g*y)-R;
%             f=yp^2+2*g*y-(obj.S/obj.s*(yp-R))^2;
        end
    end
end