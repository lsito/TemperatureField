% Script principale da lanciare per risolvere il problema bi-dimensionale
% nello specifico si vuole ottenere il campo di temperatura in un chip

% Attenzione: si è usata metà altezza del chip!

clear all
close all
clc

%Importo il file dell'andamento delle temperature nell'aletta
%per poter usare una condizione al contorno del primo tipo
Temperatures = importdata('TemperatureAletta.txt');

%Mi serve solo la prima colonna (quella a contatto con il chip)
Ts = Temperatures(1:end,1);

%Chip
L = 0.01; %m (metà altezza chip)
Hc = 0.03; %m
W = 0.3; %m
k = 149; % W/ m K
u = 2.5e6; %W/m^3

nx = 5;
ny = 5;

x = linspace(0,Hc,nx);
dx = x(2) - x(1);

y = linspace(0,L,ny);
dy = y(2) - y(1);

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
                
                b(in) = -u*dx*dy/k;
        end
      
        % Spigolo (1,1) B
        if icol==1 && irig==1 %B
                
            
                A(in,in) = -(dx/dy + dy/dx);
            
                A(in,M(irig+1,icol)) = dx/dy;
                
                A(in,M(irig,icol+1)) = dy/dx;
                
                b(in) = -(u*dx*dy)/(2*k);
        end

        % Spigolo (1,ny) C
        if icol==1 && irig==ny
            
                A(in,in) = -(dx/dy + dy/dx);
            
                A(in,M(irig-1,icol)) = dx/dy;
                
                A(in,M(irig,icol+1)) = dy/dx;
                
                b(in) = -(u*dx*dy)/(2*k); 
        end

        % Tutto il lato a contatto con l'aletta
        if icol==nx

               A(in,in) = 1;
               b(in) = Ts(irig);

        end
        
        % Lato adiabatico F
        if icol==1 && irig ~=1 && irig ~=ny

                A(in,in) = -2*(dx/dy + dy/dx);
            
                A(in,M(irig-1,icol)) = dx/dy;
                A(in,M(irig+1,icol)) = dx/dy;
                
                A(in,M(irig,icol+1)) = 2*dy/dx;
                
                b(in) = -(u*dx*dy)/(k); 
        end
        
        % Lato adiabatico G
        if irig==1 && icol ~=1 && icol ~=nx

                A(in,in) = -2*(dx/dy + dy/dx);
            
                A(in,M(irig,icol-1)) = dy/dx;
                A(in,M(irig,icol+1)) = dy/dx;
                
                A(in,M(irig+1,icol)) = 2*dx/dy;
                
                b(in) = -(u*dx*dy)/(2*k); 
        end
        
        % Lato adiabatico H
        if irig==ny && icol ~=1 && icol ~=nx

                A(in,in) = -2*(dx/dy + dy/dx);
            
                A(in,M(irig,icol-1)) = dy/dx;
                A(in,M(irig,icol+1)) = dy/dx;
                
                A(in,M(irig-1,icol)) = 2*dx/dy;
                
                b(in) = -(u*dx*dy)/(2*k); 
        end
        
end
 
 %Risoluzione del sistema A T = b:


 % Opzione 1): Metodo diretto scelto da MATLAB in base alle caratteristiche della
 % matrice A
T = A \ b;

%  % Opzione 2): Metodo diretto di Gauss
% [Tg,Am,Bm] = gauss(A, b);   % 
% 
% figure(1)
% spy(A)
% title('Matrice dei coefficienti (A)','Fontsize',20)
% 
% figure(2)
% spy(Am)
% title('Matrice dei coefficienti modificata in seguito all''applicazione del metodo di eliminazione di Gauss (A)','Fontsize',20)
% 
% if converg(A)==1
%     To = 0*ones(nodi,1);
%     tol = 1e-4;
%      % Opzione 3): Metodo iteratico di Jacobi
%     [Tj,iterj,residj]=jacobi(A,b,To,tol);
%     % Opzione 4): Metodo iteratico di Gauss Seidel
%     [Tse,iterse,residse]=seidel(A,b,To,tol);
%     
%     TT = zeros(nodi,3);   % matrice che confronta le soluzione
%     TT(:,1) = T;
%     TT(:,2) = Tj;
%     TT(:,3) = Tse;
% end

%% Grafici

%Mappa di temperatura
%Passiamo da un vettore temperatura ad una matrice temperatura (con nodi
%assegnati)

T_shape = reshape(T,ny,nx);
T_shape1 = T_shape;

%Per poter ricostruire la matrice delle temperature dell'intero sistema
save('T_shape1');

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

TmGen = Tm;
save('TmGen');

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
legend('y = W/2', 'y = 0')

set(gca,'Fontsize',20);

%% Soluzione analitica
load('THC','THC') %Devo caricare nel workspace solo il valore THc

c2 = (u * Hc^2) / (2*k) + THC;
 
T_analitica_chip = - (u*x.^2) / (2*k) + c2;
 
save('T_analitica_chip')

figure(9)
plot(x,Tm(ny,:), x, Tm(1,:), x, T_analitica_chip,'r--','LineWidth',2);
xlabel('x [m]','Fontsize',22)
ylabel('T [°C]','Fontsize',22)
hold on
grid on
set(gca,'FontSize',20)
legend('y = W/2', 'y = 0', 'soluzione analitica')

 
 
return
%Calcolo potenza dissipata dall'aletta

Q_singola = ones(nx,1);


for iq = 2:nx-1
    Q_singola(iq) = h*dx*(T_shape(1,iq) - Too);
end

Q_singola(1) = h*dx*(T_shape(1,1) - Too)/2;
Q_singola(end) = h*dx*(T_shape(1,end) - Too)/2;

Q_totale = 2 * (sum(Q_singola)) %W/m

%Supponiamo che l'aletta abbia una profondità pari a 10 volte la lunghezza
%L; in tale caso la Q_totale è pari a:

Q_no_aletta = h * W * (Ts - Too) %W/m
