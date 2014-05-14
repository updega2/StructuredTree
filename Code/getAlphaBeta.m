function [Alpha, Beta] = getAlphaBeta(a_Gamma,a_Zeta,z_Lrr)

Alpha = (1.0+a_Gamma^(a_Zeta/2.0))^(-1.0/a_Zeta);

Beta = Alpha*sqrt(a_Gamma);

end



