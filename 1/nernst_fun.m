function mV = nernst_fun(temperature, valence, Cout, Cin)
%compute Nernst for a given ion
%NOTICE: temprature in in Celcious
    R = 8.31 ; % joules/mol - gas constant - Rydberg's Constant
    T = temperature + 273; %convert to Kelvin
    F = 9.649e+4; %Faraday's constant
    
    %output in mV
    mV = 1000 * ((R*T) / (valence*F) * log(Cout/Cin)) ;
end