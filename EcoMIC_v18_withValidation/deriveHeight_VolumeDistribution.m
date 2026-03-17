% This function returns Height in function of PNA and GA
function [Ht,Vd] = deriveHeight_VolumeDistribution(PNA,GA,BW,varargin)
days = PNA;
if ~isempty(varargin)
    if string(varargin{1}) == "v" % female
        if GA <= 26
            Ht =  7.419743e-04*(days.^2)+ 5.619694e-02*(days) + 3.343850e01;    
            Vd =  0.9*BW;
        elseif GA <= 29
            Ht =   3.61874e-08*(days.^3)- 1.08198e-04*(days.^2)+ 1.25024e-01*(days) + 3.53636e01;  
            Vd =   0.85*BW;
        elseif GA <= 32
            Ht =   3.727585e-08*(days.^3) - 1.089800e-04*(days.^2)+ 1.215347e-01*(days)+ 3.891467e01; 
            Vd =   0.85*BW;
        elseif GA <= 36
            Ht =  - 3.401428e-09*(days.^4) + 3.318122e-06*(days.^3)- 1.173090e-03*(days.^2) + ...
                     2.320859e-01*(days)+ 4.316323e01;    
            Vd =   0.8*BW;
        else 
            Ht =   1.777345e-08*(days.^3)- 6.086589e-05*(days.^2)+ 8.644990e-02*(days)+ 5.148881e01;    
            Vd =   0.75*BW;
        end

    else % male
        if GA <= 26
            Ht =  8.905395e-04*(days.^2)+ 4.139884e-02*(days)+ 3.375000e01;
            Vd =  0.9*BW;
        elseif GA <= 29
            Ht =  3.14905e-08*(days.^3)- 1.05086e-04*(days.^2)+ 1.29879e-01*(days)+ 3.52325e01;
            Vd =  0.85*BW;
        elseif GA <= 32
            Ht =   3.662092e-08*(days.^3)- 1.099299e-04*(days.^2)+1.254900e-01*(days)+ 3.901344e01;
            Vd =   0.85*BW;
        elseif GA <= 36
            Ht =   -3.501005e-09*(days.^4)+ 3.405542e-06*(days.^3)- 1.218821e-03*(days.^2)+ 2.462250e-01*(days)+ 4.342661e01;
            Vd =   0.8*BW;
        else
            Ht =   2.873236e-08*(days.^3)- 8.491953e-05*(days.^2) + 9.965299e-02*(days)+ 5.173054e01;
            Vd =   0.75*BW;
        end
    end
else
    if GA <= 26
        Ht =  8.905395e-04*(days.^2)+ 4.139884e-02*(days)+ 3.375000e01;
        Vd =  0.9*BW;
    elseif GA <= 29
        Ht =  3.14905e-08*(days.^3)- 1.05086e-04*(days.^2)+ 1.29879e-01*(days)+ 3.52325e01;
        Vd =  0.85*BW;
    elseif GA <= 32
        Ht =   3.662092e-08*(days.^3)- 1.099299e-04*(days.^2)+1.254900e-01*(days)+ 3.901344e01;
        Vd =   0.85*BW;
    elseif GA <= 36
        Ht =   -3.501005e-09*(days.^4)+ 3.405542e-06*(days.^3)- 1.218821e-03*(days.^2)+ 2.462250e-01*(days)+ 4.342661e01;
        Vd =   0.8*BW;
    else
        Ht =   2.873236e-08*(days.^3)- 8.491953e-05*(days.^2) + 9.965299e-02*(days)+ 5.173054e01;
        Vd =   0.75*BW;
    end
end
end

