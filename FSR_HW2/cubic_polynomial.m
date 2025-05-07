function [s, s_d, s_dd] = cubic_polynomial(T, t)
  
    % calcolo i coef considerando che s_i = 0; s_f = 1; e s_d_i = 0, s_d_f
    % = 0;
    a_0 = 0; 
    a_1 = 0;
    a_2 = 3 / T^2;
    a_3 = -2 / T^3;

    s = a_3 * t.^3 + a_2 * t.^2 + a_1 * t + a_0;
    s_d = 3 * a_3 * t.^2 + 2 * a_2 * t + a_1; 
    s_dd = 6 * a_3 * t + 2 * a_2; 
  
end