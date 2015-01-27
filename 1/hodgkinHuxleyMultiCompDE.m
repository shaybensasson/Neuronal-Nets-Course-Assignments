function dydt = hodgkinHuxleyMultiCompDE( ...
    p,tspan,t, ...
    conMat, prevStepData, ...
    feedExternalCurrent)

    %% Get External Current
    totalElements = numel(tspan);
    index = t/p.tf*totalElements;
    index= min(index+1, totalElements); %1 based indexes
    index = index-mod(index,1); %remove remainder
    tPrev = max(t-p.dt,0);
    iPrev = max(index,1);
    Iprev = feedExternalCurrent(tPrev, iPrev); %previous step current
    if (p.TRACE==1)
        fprintf('[t:%.2f,#%d] feeding %.2f mA current ...\n', tPrev, iPrev, Iprev);
    end
    
    %% Calc Derivatives
    prevStepData = reshape(prevStepData,[p.comps 4]);
    V = prevStepData(:,1);
    n = prevStepData(:,2);
    m = prevStepData(:,3);
    h = prevStepData(:,4);
    
    %% Functions
    % K channel
    %alpha_n = @(v) iif (v ~= 10, 0.01.*(-v + 10)./(exp((-v + 10)./10) - 1), 0.1);
    alpha_n = @(v) 0.01.*(-v + 10)./(exp((-v + 10)./10) - 1);
    beta_n  = @(v) 0.125*exp(-v./80);

    % Na channel (activating)
    %alpha_m = @(v) iif (v ~= 25, 0.1.*(-v + 25)./(exp((-v + 25)./10) - 1), 1);
    alpha_m = @(v) 0.1.*(-v + 25)./(exp((-v + 25)./10) - 1);
    beta_m  = @(v) 4*exp(-v./18);

    % Na channel (inactivating)
    alpha_h = @(v) 0.07*exp(-v./20);
    beta_h  = @(v) 1./(exp((-v + 30)./10) + 1);
    
    % HH channel currents
    hh = @(Vm, g_Na, g_K, g_l) ...
        ( ...
            g_Na .* (Vm - p.E_Na) ...
            + g_K  .* (Vm - p.E_K ) ...
            + g_l  .* (Vm - p.E_l ) ...
        );
    
    g_Na = p.gbar_Na*(m.^3).*h;
    g_K  = p.gbar_K.*(n.^4);
    g_l  = p.gbar_l;

    dndt = alpha_n(V).*(1 - n) - beta_n(V).*n;
    dmdt = alpha_m(V).*(1 - m) - beta_m(V).*m;
    dhdt = alpha_h(V).*(1 - h) - beta_h(V).*h;

    dVdt = -hh(V, g_Na, g_K, g_l);
    dVdt = dVdt - conMat * V ./ p.Ra;
    
    %{
    conMat * V is the equivalent of evaluating (Vm[i] - Vm[j])/Ra 
    for every neighboring combination of i and j without the need 
    for a loop or special boundary cases!
    %}
            
        
    %provide external stimulus
    dVdt = dVdt ./ p.Cm;
    dVdt(p.elec) = dVdt(p.elec) + Iprev;

    dydt = [dVdt ; dndt ; dmdt ; dhdt];
end

