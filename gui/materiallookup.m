function mat = materiallookup(material)
% Dictionary of materials. Add to as required (+ add to dropdowns in gui)
    % [a b c] where f(T) = a + bT + cT^2
    % k: thermal conductivity [W m-1 K-1]
    % c_p: specific heat capacity [J kg-1 K-1]
    % rho: density [kg m-3]
    switch material
        case "Copper"
            k = [385 0 0];
            c_p = [399.5814286 -.0585714 0];
            rho = [8940 0 0];
        case "Inconel 718"
            k = [9.5164989 .0216787 -.0000039];
            c_p = [363.8195515 .1233661 .0000527];
            rho = [8190 0 0];
        case "Haynes 25"
            k = [9.9905357 .0205437 -.000003];
            c_p = [396.5228931 .2075422 .0000134];
            rho = [9070 0 0];
    end
    
    mat = [k; c_p; rho];
end