% Script principale da lanciare per risolvere il problema bi-dimensionale
% di aletta rettangolare. Tale script chiama le funzione pos, gauss,
% coverg, jacobi e seidel, che quindi devono essere incluse nella cartella di
% lavoro.

clear all
close all
clc

%Per poter ricostruire la matrice delle temperature dell'intero sistema
%load('T_shape1')

%Aletta
L = 0.01; %m (metà altezza aletta)
Hal = 0.2; %m %Ad 1.5 dovrebbe esserci l'indipendenza di Q da L, ma eta non può che decrescere
W = 0.3; %m

k = 236; % W/ m K
h1 = 45; %W / m^2 K
h2 = 50; %W / m^2 K

%Chip
Hc = 0.03; %m
u = 2.5e6; %W/m^3

%Ambiente
Too = 23; % °C

nx = 35;
ny = 5;

x = linspace(0,Hal,nx);
dx = x(2) - x(1);

y = linspace(0,L,ny);
dy = y(2) - y(1);

%Biot Locale
Bi1_dy = h1 * dy / k;

Bi2_dx = h2 * dx / k;

nodi = nx * ny;

A = zeros(nodi);
b = zeros(nodi,1);

% Matrice delle Posizioni (M)
M = zeros(ny,nx);

ik = 1;

for ic=1:nx
    for ir=1:ny
        M(ir,ic) = ik;
        ik = ik+1;
    end
end

for in = 1:nodi

        [irig,icol] = pos(M,in);
 
       % Elementi interni (pura conduzione)
        if irig>1  && irig<ny && icol>1  && icol<nx %A
            
                A(in,in) = -2*(dx/dy + dy/dx);
            
                A(in,M(irig-1,icol)) = dx/dy;
                A(in,M(irig+1,icol)) = dx/dy;
                
                A(in,M(irig,icol-1)) = dy/dx;
                A(in,M(irig,icol+1)) = dy/dx;
        end

      
        % Elementi a sinistra B
        if icol==1 && irig>1 && irig<ny

                A(in,in) = -2*(dx/dy + dy/dx);
            
                A(in,M(irig-1,icol)) = dx/dy;
                A(in,M(irig+1,icol)) = dx/dy;
                
                A(in,M(irig,icol+1)) = 2 * dy/dx;
                
                b(in) = - (2*u*dy*Hc)/(k);

        end

        % Elementi a destra I
        if icol==nx && irig ~=1 && irig ~=ny

                A(in,in) = -2*(dy/dx + dx/dy + Bi1_dy);

                A(in,M(irig,icol-1)) = 2*dy/dx;
 
                A(in,M(irig-1,icol)) = dx/dy;
                A(in,M(irig+1,icol)) = dx/dy; 
                
                b(in) = -2*Bi1_dy*Too;
        end

     
        % Elementi sopra G
        if irig==1 && icol ~=1 && icol ~=nx

               A(in,in) = -2*(dy/dx + dx/dy + Bi2_dx);
   
               A(in,M(irig,icol-1)) = dy/dx;
               A(in,M(irig,icol+1)) = dy/dx;

               A(in,M(irig+1,icol))= 2*dx/dy;
  
               b(in) = -2*Bi2_dx*Too;
        end
        
        % Elementi sotto H
        if irig==ny && icol ~=1 && icol ~=nx

               A(in,in) = -2*(dy/dx + dx/dy);

               A(in,M(irig,icol-1)) = dy/dx;
               A(in,M(irig,icol+1)) = dy/dx;

               A(in,M(irig-1,icol)) = 2*dx/dy;
        end
end

        % Spigolo destro alto D (1,nx)=(1,end)
         A(M(1,end),M(1,end)) = -(dy/dx + dx/dy + Bi2_dx + Bi1_dy);
         A(M(1,end),M(1,end-1)) = dy/dx;
         A(M(1,end),M(2,end)) = dx/dy;
         b(M(1,end)) = - (Bi2_dx + Bi1_dy) * Too;

        % Spigolo destro basso E (ny,nx)=(end,end) %%%%Controlla
         A(M(end,end),M(end,end)) = -(dy/dx + dx/dy + Bi1_dy);
         A(M(end,end),M(end-1,end)) = dx/dy;
         A(M(end,end),M(end,end-1)) = dy/dx;
         b(M(end,end)) = - (Bi1_dy) * Too;
         
        % Spigolo sinistro alto C (1,1)=(1,1)
         A(M(1,1),M(1,1)) = -(dy/dx + dx/dy + Bi2_dx);
         A(M(1,1),M(1,2)) = dy/dx;
         A(M(1,1),M(2,1)) = dx/dy;
         b(M(1,1)) = - (Bi2_dx) * Too - (u*dy*Hc)/(k);
         
        % Spigolo sinistro basso F (ny,1)=(end,1)
         A(M(end,1),M(end,1)) = -(dy/dx + dx/dy);
         A(M(end,1),M(end,2)) = dy/dx;
         A(M(end,1),M(end-1,1)) = dx/dy;
         b(M(end,1)) = - (u*dy*Hc)/(k);

 
 %Risoluzione del sistema A T = b:


 % Opzione 1): Metodo diretto scelto da MATLAB in base alle caratteristiche della
 % matrice A
T = A \ b;


%% Grafici

%Mappa di temperatura
%Passiamo da un vettore temperatura ad una matrice temperatura (con nodi
%assegnati)

T_shape = reshape(T,ny,nx);

%Funzione per salvare la matrice come file .txt
dlmwrite('TemperatureAletta.txt',T_shape);

T_shape2 = zeros(size(M));
for ir = 1:size(M,1)
    for ic = 1:size(M,2)
        T_shape2 (ir,ic) = T(M(ir,ic));
    end
end

Tm = zeros(size(M));
for ir = 1:size(Tm,1)
    Tm(ir,:) = T_shape(size(Tm,1)-ir+1,:);
end

figure(3)
surf(x,y,Tm);
xlabel('x [m]','Fontsize',20);
ylabel('y [m]','Fontsize',20);
zlabel('T [°C]','Fontsize',20);
%colormap hot
set(gca,'Fontsize',20);
shading interp

figure(4)
surface(x,y,Tm);
xlabel('x [m]','Fontsize',20);
ylabel('y [m]','Fontsize',20);
%colormap hot
hbi = colorbar;
title(hbi,'T [°C]');
set(gca,'Fontsize',20);
shading interp


figure(6)
[co,ho] = contour(x,y,Tm);
legend('T [°C]');
v = 30:0.5:150;
clabel(co,ho,v, 'LabelSpacing', 200, 'FontSize', 20);
%colormap hot
xlabel('x [m]','Fontsize',20);
ylabel('y [m]','Fontsize',20);
set(gca,'Fontsize',20);


figure(8)
plot(x,Tm(ny,:), x, Tm(1,:),'r--','LineWidth',2);
xlabel('x [m]','Fontsize',22)
ylabel('T [°C]','Fontsize',22)
hold on
grid on
set(gca,'FontSize',20)
legend('y = L/2', 'y = 0')

set(gca,'Fontsize',20);

Tm5050 = Tm;
save('Tm5050');


%% Soluzione analitica

%In realtà non ci vuole h medio
m = sqrt(((h2+h1)/(2)*2*(2*L+W))/(k*2*L*W));

%m = sqrt((h2*2*(Hal+W))/(k*2*L*(Hal+W)));

A = (u*Hc)/(k*m);

c1 = (A) / (exp(2*m*Hal)*((1+(h1)/(k*m))/(1-(h1)/(k*m)))-1);
c2 = A * (1 + (1) / (exp(2*m*Hal)*((1+(h1)/(k*m))/(1-(h1)/(k*m)))-1));

T_analitica = c1*exp(m*x) + c2*exp(-m*x) + Too;



THC = T_analitica(1);
save('THC');

save('T_analitica');

figure(9)
plot(x,Tm(ny,:), x, Tm(1,:), x,T_analitica,'r--','LineWidth',2);
xlabel('x [m]','Fontsize',22)
ylabel('T [°C]','Fontsize',22)
hold on
grid on
set(gca,'FontSize',20)
legend('y = W/2', 'y = 0','Soluzione Analitica')
%plot(x,T_analitica);
hold off

%Soluzione analitica con punta adiabatica:



%% Potenza termica

%Calcolo potenza dissipata dall'aletta: la potenza scambiata è la somma
%della potenza termica scambiata per convezione dalle due pareti suepriori
%e dalla parete laterale

% Q scambiata dalle due pereti laterali (h2)
Q_singola_h2 = ones(nx,1);

for iq = 2:nx-1
    Q_singola_h2(iq) = h2*dx*(T_shape(1,iq) - Too);
end

% Q ai bordi le calcolo separatamente perchè il volumetto è tagliato
Q_singola_h2(1) = h2*dx*(T_shape(1,1) - Too)/2;
Q_singola_h2(end) = h2*dx*(T_shape(1,end) - Too)/2;

% Due volte per la simmetria del problema
Q_totale_h2 = 2 * (sum(Q_singola_h2)); %W/m

% Q scambiata dalla parete superiore (h1)
Q_singola_h1 = ones(ny,1);

for iq1 = 2:ny-1
    Q_singola_h1(iq1) = h1*dy*(T_shape(iq1,end)-Too);
end

Q_singola_h1(1) = h1*dy*(T_shape(1,end) - Too)/2;
Q_singola_h1(end) = h1*dy*(T_shape(end,end) - Too)/2;

Q_totale_h1 = 2 * (sum(Q_singola_h1)); %W/m

% Fondamentalmente nella somma il contributo di h1 è trascurabile in quanto
% la superficie di scambio è molto minore di quella relativa a h2
Q_aletta = Q_totale_h1 + Q_totale_h2;

%Supponiamo che l'aletta abbia una profondità pari a 10 volte la lunghezza
%L; in tale caso la Q_totale è pari a:

%Q_no_aletta = h * W * (Ts - Too); %W/m

%% Rendimento 

%La temperatura a contatto con i chip
Ts = max(max(T_shape)); %°C

% Ricaviamo eta dalla definizione (!nb si sono considerati solo i lati che 
% scambiano e non le aree = è potenza su unità di lunghezza)
eta_numerico = Q_aletta / ((Ts-Too)*(((2*Hal*h2)+2*L*h1)));

% Ricaviamo eta con una relazione presente sul manuale
%m = sqrt((h2*2*(2*L+W))/(k*2*L*W));
% Con ottima approssimazione possiamo dire che la punta dell'aletta è
% adiabatica
eta_analitico = tanh(m*Hal)/(m*Hal);

%% Parete alettata
%Vi è un dissipatore costituito da 4 alette uguali, con passo tra le alette
%pari al loro spessore 

%Dati del problema
N_alette = 4; 
Spessore = 2*L; %Per noi L è la metà della larghezza dell'aletta
Passo = Spessore;

%Possiamo considerare con buona approssimazione che la temperatura del chip
%è uniforme su tutto il chip: se prendiamo il minimo abbiamo il caso 
%peggiore, cioè quello con uno scambio termico inferiore, ma nello
%specifico essendo già stato utilizzato Ts si userà il massimo del massimo

Gpna = h1 * ((N_alette-1)*Passo);
Gpa = N_alette * eta_numerico * (2*Hal*h2 + 2*L*h1);

Q_parete = (Ts-Too) * (Gpna + Gpa);

eta_parete = Q_parete / (4*((Ts-Too)*((2*Hal*h2+2*L*h1))) + (N_alette-1)*h1*Passo*(Ts-Too));

%!nb: in verità quando non c'è l'aletta la potenza termica scambiata
%dovrebbe essere ricavata diversamente, la temperatura non è più Ts, stiamo
%approssimando

%% Parete alettata (equazioni semplificate)

%Attenzione, si fa notare che con area dell'aletta si intende la superficie
%di scambio, che nello specifico è in realtà solo il perimetro esterno

Ae = ((2*L+2*Hal)*N_alette + (N_alette-1)*2*Passo);
% A rigore andrebbe usata h1, ma siccome la superficie con cui lo scambio
% convettivo dovrebbe avvenire con h1 è molto minore della superficie per
% cui lo scambio dovrebbe avvenire con h2, è una buona ipotesi
% semplificativa
Geq = h2 * (Ae + (eta_numerico-1)*N_alette*(2*L+2*Hal));
Q_parete_1 = (Ts-Too) * Geq
eta_parete_1 = (Ae + (eta_numerico-1)*N_alette*(2*L+2*Hal))/(Ae)

%% Parete alettata ipotizzando che il passo possa variare
% Mi serve una funzione per iterare un ciclo in cui il parametro da variare
% è il passo tra le alette

% C'è un problema quando il passo è 0, in quel caso, l'aletta è singola, le
% equazioni non tengono conto di questa cosa

multiplier = linspace(1,10,10); 

for i = 1:10
    Passo_var(i) = Spessore * multiplier(i);    
    
    [X(i),Y(i)] = parete(N_alette, Spessore, Passo_var(i), L,h1,h2,k,Ts,Too,Hal,W,eta_numerico);
end

%Plot della potenza termica smaltita al variare del passo tra le alette
figure(10)
plot(Passo_var,X,'LineWidth',2);
xlabel('Passo [m]','Fontsize',12)
ylabel('Potenza smaltita [W]','Fontsize',12)
hold on
grid on
set(gca,'FontSize',12)
%legend('y = W/2', 'y = 0','Soluzione Analitica')
%plot(x,T_analitica);
hold off



%Non viene richiesto di ricavare eta della parete al variare del passo, ma,
%una volta fatto perchè all'aumentare del passo eta della parete aumenta
%invece di diminuire

%{

% Efficienza singola aletta
m_sing = sqrt((h2*2*(Spessore+W))/(k*Spessore*W));
eta_singola_aletta = tanh(m_sing*Hal)/(m_sing*Hal);

% Potenza ideale singola aletta
Qid = ((Ts-Too)*(((2*Hal*h2)+Spessore*h1)));

Gpa = eta_singola_aletta * Qid/(Ts-Too) * N_alette;

Apna = Passo * N_alette;
Gpna = h1*Apna;

% Conduttanza equivalente
Geq = Gpa + Gpna;

% Potenza termica totale scambiata
Qtot = Geq * (Ts-Too);
% Efficienza parete alettata (che considera tutta la parete, sia quella 
% alettata che quella non alettata)
eta_parete = Qtot / (N_alette * Qid + 4 * h1 * Passo * (Ts-Too));
%}


