%The rate of the Cox protein
function f = functionforCox(delta, zeta, gamma, C, epsilon,DCox,Cox,eta,theta)
f = delta - zeta*gamma*C + epsilon-DCox*Cox-eta*Cox-theta*Cox;
end

%d[Cox]/dt=P_Pe+P_pBAD-D_Cox*[Cox]-I_(C_2)*M_(C_2)*[C]-T_(Cox_4)*[Cox]-I_Cox*T_(Cox_4)*[Cox]
%P_Pe = delta
%P_pBAD = epsilon
%I_(C_2) = zeta
%M_(C_2) = gamma
%T_(Cox_4) = eta
%I_Cox = theta
%D_Cox = DCox
