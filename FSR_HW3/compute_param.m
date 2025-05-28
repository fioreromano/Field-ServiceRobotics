function [C, Q, R_b, M] = compute_param(eta_b, eta_dot_b, Ib)
    phi = eta_b(1);
    theta = eta_b(2);
    psi = eta_b(3);
    Q_dot=[0 0 -cos(theta)*eta_dot_b(2);
           0 -sin(phi)*eta_dot_b(1) -sin(theta)*sin(phi)*eta_dot_b(2)+cos(theta)*cos(phi)*eta_dot_b(1);
           0 -cos(phi)*eta_dot_b(1) -sin(theta)*cos(phi)*eta_dot_b(2)-cos(theta)*sin(phi)*eta_dot_b(1)];
    
    Q = [1 0 -sin(theta);
        0 cos(phi) cos(theta)*sin(phi);
        0 -sin(phi) cos(theta)*cos(phi)];

    S = skew(Q*eta_dot_b');

    C = Q'*S*Ib*Q + Q'*Ib*Q_dot;

    R_b = [cos(theta)*cos(psi) sin(phi)*sin(theta)*cos(psi)-cos(phi)*sin(psi) cos(phi)*sin(theta)*cos(psi)+sin(phi)*sin(psi);
          cos(theta)*sin(psi)  sin(phi)*sin(theta)*sin(psi)-cos(phi)*cos(psi) cos(phi)*sin(theta)*sin(psi)-sin(phi)*cos(psi);
          -sin(theta)               sin(phi)*cos(theta)                                           cos(phi)*cos(theta)];
    
    M = Q'*Ib*Q;
    

end

function S=skew(w)

    if(numel(w)~= 1)
        S = [0    -w(3)  w(2);
             w(3)  0    -w(1);
            -w(2)  w(1)  0];
    else
        S= zeros(3);
    end

end