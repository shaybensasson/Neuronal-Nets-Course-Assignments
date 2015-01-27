function dydt = hodgkinHuxleyDE( ...
    p,tspan,t, ...
    prevStepData)
    %% Get External Current
    percentile = t/tspan(end) *100;
    percentile = percentile-mod(percentile,1); %floor
    I = feedExternalCurrent(p, percentile);
    global LastPercentileToTrace;
    curMod = mod(percentile,10);
    if (curMod==0 && LastPercentileToTrace < percentile)
        LastPercentileToTrace = percentile;
        fprintf('[perc:%d] feeding %d mA current ...\n', percentile, I);
    end
    
    %% Calc Derivatives
    V = prevStepData(1);
    n = prevStepData(2);
    m = prevStepData(3);
    h = prevStepData(4);
    
    %n gates - Potassium (K+)
    alpha_n = alpha_n_fun(V); %Equation 12
    beta_n = beta_n_fun(V); %Equation 13
        
    %m gates - Sodium (Na+)
    alpha_m = alpha_m_fun(V); %Equation 20
    beta_m = beta_m_fun(V); %Equation 21
    
    %h gates - Sodium (Na+)
    alpha_h = alpha_h_fun(V); %Equation 23
    beta_h = beta_h_fun(V); %Equation 24
            
    %calc the currents
    g_K = p.gbar_K*(n^4); %Equation 6
    I_K = g_K * (V-p.E_K);  %Equation 4
    
    g_Na = (m^3) * h * p.gbar_Na; %Equation 14
    I_Na = g_Na * (V-p.E_Na); %Equation 3
    
    I_L = p.gbar_L *(V-p.E_l); %Equation 5
    I_ion = I + -1*(I_K + I_Na + I_L); %Equation 2
    %NOTE: 
    %in order for depolarization will be a positive event -> positive ions 
    % fall into a neuron that is gonna be more positive.
    % So whenever we have voltage or equilibrium voltage, we're gonna flip the
    % sign on that.
    
    
    %---calculate the derivatives on a given point ---%
    dvdt = I_ion/p.C;
    dndt = dnmh_over_dt_fun(n, alpha_n, beta_n); %Equation 7
    dmdt = dnmh_over_dt_fun(m, alpha_m, beta_m); %Equation 15
    dhdt = dnmh_over_dt_fun(h, alpha_h, beta_h); %Equation 16
    
    dydt = [dvdt ; dndt ; dmdt ; dhdt];
    
    %{
    function dn(varname)
        try
            disp(sprintf('%s=%g', varname, evalin('caller',varname)));
        catch
            disp(sprintf('%s=N/A', varname));
        end
        
    end

    list = who('*');
    for n = 1:1:size(list,1)
        var = list{n};
        dn(var);
    end
    
    p
    pause    
    %}
end

