function V = alpha_m_fun( V )
    %NOTE: 
    %in order for depolarization will be a positive event -> positive ions 
    % fall into a neuron that is gonna be more positive.
    % So whenever we have voltage or equilibrium voltage, we're gonna flip the
    % sign on that.
    V = -1*(V);
    
    %Equation 20
    V =.1*( (V+25) / (exp((V+25)/10)-1) );   
end

