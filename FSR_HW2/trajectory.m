function [v_t, w_t, x_s, y_s, theta] = trajectory(s, s_d, s_dd, qf, qi, alpha_x, alpha_y, beta_x, beta_y)
    
    x_s = s.^3.*qf(1) - (s-1).^3.*qi(1) + alpha_x.*s.^2.*(s-1) + beta_x.*s.*(s-1).^2;
    y_s = s.^3.*qf(2) - (s-1).^3.*qi(2) + alpha_y.*s.^2.*(s-1) + beta_y.*s.*(s-1).^2;

    dx_ds = 3 .* s.^2 .* qf(1) - 3 .* (s - 1).^2 .* qi(1) + alpha_x .* s .* (3 .* s - 2) + beta_x .* (3 .* s - 1) .* (s - 1);
    dy_ds = 3 .* s.^2 .* qf(2) - 3 .* (s - 1).^2 .* qi(2) + alpha_y .* s .* (3 .* s - 2) + beta_y .* (3 .* s - 1) .* (s - 1);

    ddx_ds = 6 .* s .* qf(1) - 6 .* (s - 1) .* qi(1) + alpha_x .* (6 .* s - 2) + beta_x .* (6 .* s - 4);
    ddy_ds = 6 .* s .* qf(2) - 6 .* (s - 1) .* qi(2) + alpha_y .* (6 .* s - 2) + beta_y .* (6 .* s - 4);

    % Calcolo dell'angolo theta della traiettoria
    theta = atan2(dy_ds, dx_ds);

    % Calcolo della velocità tangenziale
    v_tilde = sqrt(dx_ds.^2 + dy_ds.^2);

    % Calcolo della velocità angolare
    w_tilde = (ddy_ds .* dx_ds - ddx_ds .* dy_ds) ./ ((dx_ds).^2 + (dy_ds).^2);

    % Calcolo di v_t e w_t
    v_t = v_tilde .* s_d;
    w_t = w_tilde .* s_d;
end