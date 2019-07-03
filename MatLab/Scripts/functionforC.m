%The rate of the C protein
function f = functionforC(alpha,beta,Tet,DC, C,gamma)
f = alpha - beta*Tet-DC*C-gamma*C;
end

%d[C]/dt=P_pTet-I_TetR*[TetR]-D_C*[C]-M_(C_2)*[C]  
%P_pTet = alpha
%I_TetR = beta
%M_(C_2) = gamma
%D_C = DC

