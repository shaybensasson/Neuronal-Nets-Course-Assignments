function nmh = dnmh_over_dt_fun(x, alpha_x, beta_x)
    %dn/dt or dm/dt or dh/dt - Equations 7,15,16
    nmh = alpha_x *(1-x) - beta_x * x;
end

