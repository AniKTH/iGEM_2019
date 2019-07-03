% Modeling the time series for all molecules
% f = functionforC(alpha,beta,Tet,DC, C,gamma)
% f = functionforCox(delta, zeta, gamma, C, epsilon,DCox,Cox,eta,theta)
% f = functionforTetR(epsilon,iota,Ara,DTet,Tet)

% Paramaters
alpha = 10;  %P_pTet = alpha, strength of the pTet promoter
beta = 1;   %I_TetR = beta, Inhibition by TetR 
gamma = 2;  %M_(C_2) = gamma, Dimerization of C
delta = 4;  %P_Pe = delta, strength of the Pe promoter
epsilon = 3;%P_pBAD = epsilon, strength of the pBAD promoter
zeta = 1.5;   %I_(C_2) = zeta, Inhibition by the C2 dimer
eta =1;       %T_(Cox_4) = eta, Tetramization of Cox
theta =2;     %I_Cox = theta, Inhibition by the Cox4 tetramer
iota =6;      %I_Ara = iota, Inducment from Arabinose
DC = 0.2;     %D_C = DC, Degredation rate of C
DCox = 0.2;   %D_Cox = DCox, Degredation rate of Cox
DTet = 0.2;   %D_TetR = DTet, Degredation of TetR
Ara = 0;      %Arabinose concentration

% Time steps
t1 = 0;
dt = 0.01;
t2 = 10;
t = [t1:dt:t2];
size = length(t);

% Initial concentration of molecules
C = zeros(size,4);
C(1,1) = 0.5;
C(1,2) = 1.5;
C(1,3) = 5;
C(1,4) = 10;
Cox = zeros(size,4);
Cox(1,1) = 0.5;
Cox(1,2) = 1.5;
Cox(1,3) = 5;
Cox(1,4) = 10;
Tet = zeros(size,4);
%{%}
Tet(1,1) = 0.5;
Tet(1,2) = 1.5;
Tet(1,3) = 5;
Tet(1,4) = 10;

% Making a container for the f(x) values
fC = zeros(size,4); % Makes the container for the f(x) values the same size of the C matrix
fC(1,1) = functionforC(alpha,beta,Tet(1,1),DC, C(1,1),gamma);
fC(1,2) = functionforC(alpha,beta,Tet(1,2),DC, C(1,2),gamma);
fC(1,3) = functionforC(alpha,beta,Tet(1,3),DC, C(1,3),gamma);
fC(1,4) = functionforC(alpha,beta,Tet(1,4),DC, C(1,4),gamma);
fCox = zeros(size,4); % Makes the container for the f(x) values the same size of the Cox matrix
fCox(1,1) = functionforCox(delta, zeta, gamma, C(1,1), epsilon,DCox,Cox(1,1),eta,theta);
fCox(1,2) = functionforCox(delta, zeta, gamma, C(1,2), epsilon,DCox,Cox(1,2),eta,theta);
fCox(1,3) = functionforCox(delta, zeta, gamma, C(1,3), epsilon,DCox,Cox(1,3),eta,theta);
fCox(1,4) = functionforCox(delta, zeta, gamma, C(1,4), epsilon,DCox,Cox(1,4),eta,theta);
fTet = zeros(size,4); % Makes the container for the f(x) values the same size of the TetR matrix
fTet(1,1) = functionforTetR(epsilon,iota,Ara,DTet,Tet(1,1));
fTet(1,2) = functionforTetR(epsilon,iota,Ara,DTet,Tet(1,2));
fTet(1,3) = functionforTetR(epsilon,iota,Ara,DTet,Tet(1,3));
fTet(1,4) = functionforTetR(epsilon,iota,Ara,DTet,Tet(1,4));

% f = functionforC(alpha,beta,Tet,DC, C,gamma)
% f = functionforCox(delta, zeta, gamma, C, epsilon,DCox,Cox,eta,theta)
% f = functionforTetR(epsilon,iota,Ara,DTet,Tet)
%{
For calculating single initial concentrations
g = 1; %Resets the row determination
for v = dt:dt:t2 %...for every time step...
    C(g+1) = C(g) + dt*functionforC(alpha,beta,Tet(g),DC, C(g),gamma)); %...calculate and save the new concentration of C in the next row, using Euler.
    Cox(g+1) = Cox(g) + dt*functionforCox(delta, zeta, gamma, C(g), epsilon,DCox,Cox(g),eta,theta); %...and calculate and save the new concentration of Cox in the next row, using Euler
    Tet(g+1) = Tet(g) + dt*functionforTetR(epsilon,eta,Ara,DTet,Tet(g)); %...and calculate and save the new concentration of TetR in the next row, using Euler
    g = g+1; % The next row becomes the current.
end
%}

%For calculating several initial concentrations
for y = 1:4 % For every column in the concentration matrices (for every initial concentration)...
    g = 1; %Resets the row determination
    for v = dt:dt:t2 %...for every time step...
        C(g+1,y) = C(g,y) + dt*functionforC(alpha,beta,Tet(g,y),DC, C(g,y),gamma); %...calculate and save the new concentration of C in the next row, using Euler.
        Cox(g+1,y) = Cox(g,y) + dt*functionforCox(delta, zeta, gamma, C(g,y), epsilon,DCox,Cox(g,y),eta,theta); %...and calculate and save the new concentration of Cox in the next row, using Euler
        Tet(g+1,y) = Tet(g,y) + dt*functionforTetR(epsilon,iota,Ara,DTet,Tet(g,y)); %...and calculate and save the new concentration of TetR in the next row, using Euler
        fC(g+1,y) = functionforC(alpha,beta,Tet(g,y),DC, C(g,y),gamma);
        fCox(g+1,y) = functionforCox(delta, zeta, gamma, C(g,y), epsilon,DCox,Cox(g,y),eta,theta);
        fTet(g+1,y) = functionforTetR(epsilon,iota,Ara,DTet,Tet(g,y));
        g = g+1; % The next row becomes the current.
    end
end
%{
% Producing the time series
subplot(2,2,1)
plot(t, C(:,1),'k')
hold on
plot(t, Cox(:,1),'b')
plot(t, Tet(:,1),'g')
xlabel('Time [time unit]')
ylabel('Concentration [concentration unit]')
title(legend('C','Cox',"TetR"), 'Molecules')
title('Time series')

subplot(2,2,2)
plot(t, C(:,2),'k')
hold on
plot(t, Cox(:,2),'b')
plot(t, Tet(:,2),'g')
xlabel('Time [time unit]')
ylabel('Concentration [concentration unit]')
title(legend('C','Cox',"TetR"), 'Molecules')
title('Time series')

subplot(2,2,3)
plot(t, C(:,3),'k')
hold on
plot(t, Cox(:,3),'b')
plot(t, Tet(:,3),'g')
xlabel('Time [time unit]')
ylabel('Concentration [concentration unit]')
title(legend('C','Cox',"TetR"), 'Molecules')
title('Time series')

subplot(2,2,4)
plot(t, C(:,4),'k')
hold on
plot(t, Cox(:,4),'b')
plot(t, Tet(:,4),'g')
title(legend('C','Cox',"TetR"), 'Molecules')
xlabel('Time [time unit]')
ylabel('Concentration [concentration unit]')
title('Time series')

%Producing a phase portrait
plot(Cox(:,1),C(:,1))
hold on
plot(Cox(:,2),C(:,2))
plot(Cox(:,3),C(:,3))
plot(Cox(:,4),C(:,4))
xlabel('[Cox]')
ylabel('[C]')
title(legend('0.5','1.5',"5",'10'), 'Initial concentration')
title('Phase portrait')
%}
%Producing a stability plot
plot(C(:,1),fC(:,1),'k')
hold on
plot(Cox(:,1),fCox(:,1),'b')
plot(Tet(:,1),fTet(:,1),'g')
xlabel('[Molecule]')
ylabel('F(x)')
title(legend('C','Cox',"TetR"), 'Molecules')
title('Stability')
%{%}