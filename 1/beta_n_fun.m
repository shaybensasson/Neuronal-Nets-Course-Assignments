function V = beta_n_fun( V )
    
    %NOTE: 
    %in order for depolarization will be a positive event -> positive ions 
    % fall into a neuron that is gonna be more positive.
    % So whenever we have voltage or equilibrium voltage, we're gonna flip the
    % sign on that.
    V = -1*(V);
    
    %Equation 13
    V =.125*exp(V/80);
end

