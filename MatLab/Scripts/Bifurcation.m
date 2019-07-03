% Modeling the bifurcation plots
tic % Start the timer

% Paramaters
alpha = 10;  %P_pTet = alpha, strength of the pTet promoter
%beta = 1;   %I_TetR = beta, Inhibition by TetR 
gamma = 2;  %M_(C_2) = gamma, Dimerization of C
delta = 4;  %P_Pe = delta, strength of the Pe promoter
epsilon = 3;%P_pBAD = epsilon, strength of the pBAD promoter
zeta = 1.5;   %I_(C_2) = zeta, Inhibition by the C2 dimer
eta =1;       %T_(Cox_4) = eta, Tetramization of Cox
theta =2;     %I_Cox = theta, Inhibition by the Cox4 tetramer
iota =6;      %I_Ara = iota, Inducment from Arabinose
DC = 0.2;     %D_C = DC, Degredation rate of C%
DCox = 0.2;   %D_Cox = DCox, Degredation rate of Cox
DTet = 0.2;   %D_TetR = DTet, Degredation of TetR
Ara = 0;      %Arabinose concentration

% Tested parameter
p1 = 0; %Lower limit to the altered parameter
p2 = 10; %Upper limit to the atered parameter
dp = 0.001; %Stepwise change in the parameter
beta = p1:dp:p2; %Matrix with the different values of gamma

% Initial concentration of molecules
%Varying the concentration of C
C1 = -15; %Lower limit to the concentration of C 
C2 = 50; %Upper limit to the concentration of C 
dC = 0.001; %Stepwise change in the concentration of C 
C = C1:dC:C2; %Matrix with the different concentrations of C

%Varying the concentration of Cox
Cox1 = -15; %Lower limit to the concentration of Cox
Cox2 = 50; %Upper limit to the concentration of Cox
dCox = 0.001; %Stepwise change in the concentration of Cox
Cox = Cox1:dCox:Cox2; %Matrix with the different concentrations of Cox

%Varying the concentration of TetR
Tet1 = -15; %Lower limit to the concentration of TetR
Tet2 = 50; %Upper limit to the concentration of TetR
dTet = 0.001; %Stepwise change in the concentration of TetR
Tet = Tet1:dTet:Tet2; %Matrix with the different concentrations of TetR

% The sizes of the matrixes, important for the containers and the loops
column = length(beta);
size = length(C);

% Creating containers for the steady state values and its corresponding
% gamma value
if size < column %The container should be as large as possible
    Cx = zeros(column,1);
    Coxx = zeros(column,1);
    Tetx = zeros(column,1);
    Cy = zeros(column,1);
    Coxy = zeros(column,1);
    Tety = zeros(column,1);
else
    Cx = zeros(size,1);
    Coxx = zeros(size,1);
    Tetx = zeros(size,1);
    Cy = zeros(size,1);
    Coxy = zeros(size,1);
    Tety = zeros(size,1);
end

% Row determining variable for the containers
gC = 1;
gCox = 1;
gTet = 1;

% Calculating f(x)
for j = 1:column % For every value of gamma...
    for i = 1:size % ...for every concentration value of C...
        fC = functionforC(alpha,beta(j),Tet(i),DC, C(i),gamma);%Calculate the f(x) using the corresponding concentration value of C
        fCox = functionforCox(delta, zeta, gamma, C(i), epsilon,DCox,Cox(i),eta,theta);
        fTet = functionforTetR(epsilon,iota,Ara,DTet,Tet(i));
        if abs(fC) < 0.001  %...if f(x) is close enough to zero...
            %if abs(C(i-1)-C(i)) > 0.0001 %... and the concentration values arent too close...
            Cx(gC) = beta(j); %...save the Arabinose value in the x container...
            Cy(gC) = C(i); % ...and the steady state value in the y container.
            gC = gC+1; % Moves the next collected steady state and gamma value to the next row.
            %end
        end
        if  abs(fCox) < 0.001  %...if f(x) is close enough to zero...
            Coxx(gCox) = beta(j); %...save the Arabinose value in the x container...
            Coxy(gCox) = Cox(i); % ...and the steady state value in the y container.
            gCox = gCox+1; % Moves the next collected steady state and gamma value to the next row.
        end
        if  abs(fTet) < 0.001  %...if f(x) is close enough to zero...
            Tetx(gTet) = beta(j); %...save the Arabinose value in the x container...
            Tety(gTet) = Tet(i); % ...and the steady state value in the y container.
            gTet = gTet+1; % Moves the next collected steady state and gamma value to the next row.
        end
    end
    j %For Analis peace of mind
end

plot(Cx,Cy,'k.') %Plot the values in the y container (the steady state values) over the values in the x container (the corresponding gamma values)
hold on
plot(Coxx,Coxy,'b.') %Plot the values in the y container (the steady state values) over the values in the x container (the corresponding gamma values)
plot(Tetx,Tety,'g.') %Plot the values in the y container (the steady state values) over the values in the x container (the corresponding gamma values)
xlabel('Beta') %Labels the x-axis
ylabel('[Concentration]') %Labels the y-axis
grid on
xlim([-1 11])
ylim([-1 17])
title('Bifurcation: I_T_e_t_R') %Gives the title of the graph
title(legend('C','Cox',"TetR"), 'Molecules')
toc %Stop the timer