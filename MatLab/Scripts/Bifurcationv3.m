clear all
clc
tic % Start the timer
% Paramaters
alpha = 10;  %P_pTet = alpha, strength of the pTet promoter
beta = 14;   %I_TetR = beta, Inhibition by TetR 
gamma = 2;  %M_(C_2) = gamma, Dimerization of C
delta = 4;  %P_Pe = delta, strength of the Pe promoter
epsilon = 3;%P_pBAD = epsilon, strength of the pBAD promoter
zeta = 1.5;   %I_(C_2) = zeta, Inhibition by the C2 dimer
eta =1;       %T_(Cox_4) = eta, Tetramization of Cox
theta =2;     %I_Cox = theta, Inhibition by the Cox4 tetramer
iota =6;      %I_Ara = iota, Inducment from Arabinose
%DC = 0.2;     %D_C = DC, Degredation rate of C
DCox = 0.2;   %D_Cox = DCox, Degredation rate of Cox
DTet = 0.2;   %D_TetR = DTet, Degredation of TetR
Ara = 1;      %Arabinose concentration

% Tested parameter
p1 = 0; %Lower limit to the altered parameter
p2 = 14; %Upper limit to the atered parameter
dp = 0.01; %Stepwise change in the parameter
parameter = p1:dp:p2; 
column = length(parameter);

% Time steps
t1 = 0;
dt = 0.001;
t2 = 60 ;
t = [t1:dt:t2];
size = length(t);

% Initial concentration of molecules
C = zeros(size,column);
C(1,:) = -5;

Cox = zeros(size,column);
Cox(1,:) = -5;

Tet = zeros(size,column);
Tet(1,:) = -5;

% Making a container for the f(x) values
fC = zeros(size,column); % Makes the container for the f(x) values the same size of the C matrix

fCox = zeros(size,column); % Makes the container for the f(x) values the same size of the Cox matrix
fTet = zeros(size,column); % Makes the container for the f(x) values the same size of the TetR matrix
for j = 1:column
    fC(1,j) = functionforC(alpha,beta,-5,parameter(j), -5,gamma);
    fCox(1,j) = functionforCox(delta, zeta, gamma, -5, epsilon,DCox,-5,eta,theta);
    fTet(1,j) = functionforTetR(epsilon,iota,Ara,DTet,-5);
end

Cx = zeros(size,1);
Coxx = zeros(size,1);
Tetx = zeros(size,1);
Cy = zeros(size,1);
Coxy = zeros(size,1);
Tety = zeros(size,1);
% Row determining variable for the containers
gC = 1;
gCox = 1;
gTet = 1;


%For calculating several initial concentrations
for j = 1:column
    g = 1; %Resets the row determination
    for v = dt:dt:t2 %...for every time step...
        C(g+1,j) = C(g,j) + dt*functionforC(alpha,beta,Tet(g,j),parameter(j), C(g,j),gamma); %...calculate and save the new concentration of C in the next row, using Euler.
        Cox(g+1,j) = Cox(g,j) + dt*functionforCox(delta, zeta, gamma, C(g,j), epsilon,DCox,Cox(g,j),eta,theta); %...and calculate and save the new concentration of Cox in the next row, using Euler
        Tet(g+1,j) = Tet(g,j) + dt*functionforTetR(epsilon,iota,Ara,DTet,Tet(g,j)); %...and calculate and save the new concentration of TetR in the next row, using Euler
        
        fC(g+1,j) = functionforC(alpha,beta,Tet(g,j),parameter(j), C(g,j),gamma);
        fCox(g+1,j) = functionforCox(delta, zeta, gamma, C(g,j), epsilon,DCox,Cox(g,j),eta,theta);
        fTet(g+1,j) = functionforTetR(epsilon,iota,Ara,DTet,Tet(g,j));
        
        g = g+1; % The next row becomes the current.
    end
end

for j = 1:column
    for i = 1:size
        if abs(fC(i,j)) < 0.07
            Cx(gC) = parameter(j); %...save the Arabinose value in the x container...
            Cy(gC) = C(i,j); % ...and the steady state value in the y container.    
            gC = gC+1; % Moves the next collected steady state and gamma value to the next row.
        end
        if abs(fCox(i,j)) < 0.01
            Coxx(gCox) = parameter(j); %...save the Arabinose value in the x container...
            Coxy(gCox) = Cox(i,j); % ...and the steady state value in the y container.    
            gCox = gCox+1; % Moves the next collected steady state and gamma value to the next row.
        end
        if abs(fTet(i,j)) < 0.01
            Tetx(gTet) = parameter(j); %...save the Arabinose value in the x container...
            Tety(gTet) = Tet(i,j); % ...and the steady state value in the y container.    
            gTet = gTet+1; % Moves the next collected steady state and gamma value to the next row.
        end
    end
end

subplot(1,3,3)
plot(Cx,Cy,'k.') %Plot the values in the y container (the steady state values) over the values in the x container (the corresponding gamma values)
hold on
plot(Coxx,Coxy,'b.') %Plot the values in the y container (the steady state values) over the values in the x container (the corresponding gamma values)
plot(Tetx,Tety,'g.') %Plot the values in the y container (the steady state values) over the values in the x container (the corresponding gamma values)
xlabel('D_C') %Labels the x-axis
ylabel('Protein concentration') %Labels the y-axis
grid on
xlim([-1 15])
title('Bifurcation: D_C') %Gives the title of the graph
title(legend('C','Cox',"TetR"), 'Protein')
ylim([-2025, 2000])

subplot(1,3,1)
plot(t, C(:,1001),'k')
hold on
plot(t, Cox(:,1001),'b')
plot(t, Tet(:,1001),'g')
xlabel('Time')
ylabel('Protein concentration')
title(legend('C','Cox',"TetR"), 'Protein')
title('Time series')
grid on

subplot(1,3,2)
plot(C(:,1001),fC(:,1001),'k')
hold on
plot(Cox(:,1001),fCox(:,1001),'b')
plot(Tet(:,1001),fTet(:,1001),'g')
xlabel('Protein concentration')
ylabel('f(x)')
title(legend('C','Cox',"TetR"), 'Protein')
title('Stability')
grid on
toc %Stop the timer