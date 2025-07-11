function dx = quad2D_dynamics(t, x, p)
    % Extract states
    y = x(1); z = x(2); phi = x(3);
    y_dot = x(4); z_dot = x(5); phi_dot = x(6);

    % Errors
    e_y = p.y_des - y; e_y_dot = p.y_dot_des - y_dot; %p.ydot_des = 0 since reference is set-point
    e_z = p.z_des - z; e_z_dot = p.z_dot_des - z_dot; %p.zdot_des = 0 since reference is set-point

    % Gains
    Kp_y = 6; Kd_y = 10;
    Kp_z = 20; Kd_z = 8;
    Kp_phi = 50; Kd_phi = 15;

    % Altitude hold controller
    T = p.m * (Kd_z*e_z_dot + Kp_z*e_z + p.g);% Control zddot_des = 0 as zdes is constant

    
    %y-Position and phi controller

    % commanded yddot_c
    y_ddot_c = (p.y_ddot_des + Kp_y*e_y + Kd_y*e_y_dot); %p.yddot_des = 0 since reference is constant set-point
    %Desired phi angle generated by outer-loop
    phi_des = -(1/p.g)*y_ddot_c; %phi_c
    phi_dot_des = -(1/p.g)*(-Kp_y*e_y_dot + Kd_y*phi*p.g); %Consider 0 for simplicity

    % Roll controller
    tau_phi = Kp_phi*(phi_des - phi) + Kd_phi*(phi_dot_des - phi_dot); %   phi_dot_des = 0 (Assumption)         %(phidot_des-phi_dot),phidotdot,Ixx 
    


    % Disturbances
    if t > 5 && t < 20
        dy = double(1 * sin(2*pi*2/10 * t) > 0);  % 1 if positive, 0 otherwise
    else
        dy = 0;
    end

    % Accelerations
    y_ddot = -sin(phi)*T/p.m ;
    z_ddot = cos(phi)*T/p.m - p.g;
    phi_ddot = tau_phi / p.Ixx;


    % State derivative
    dx = [y_dot+ dy/p.m; z_dot; phi_dot; y_ddot; z_ddot; phi_ddot];
end
