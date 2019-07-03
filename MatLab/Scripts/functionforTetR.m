%The rate of the TetR protein
function f = functionforTetR(epsilon,iota,Ara,DTet,Tet)
f = epsilon + iota*Ara-DTet*Tet;
end

%d[TetR]/dt=P_pBAD+I_Ara*[Arabinose]-D_TetR*[TetR]
%P_pBAD = epsilon
%I_Ara = iota
%D_TetR = DTet

