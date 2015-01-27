function V = beta_h_fun( V )
    %NOTE: 
    %in order for depolarization will be a positive event -> positive ions 
    % fall into a neuron that is gonna be more positive.
    % So whenever we have voltage or equilibrium voltage, we're gonna flip the
    % sign on that.
    V = -1*(V);

    %Equation 24
    V = 1/(exp((V+30)/10)+1);
end

