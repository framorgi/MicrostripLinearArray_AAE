
clear all;
close all;
clc;
%% COSTANTI
c=3*10^8;
v0=3e8;
mu0=4*3.14e-7;
e0=8.85e-12;
%% SPECIFICHE DI PROGETTO 
%Specifiche dell'antenna:
%-Frequenza centrale
f0= 2.45e9; %GHz
% da cui lunghezza onda freespace lambda0
lambda0= c/f0
% da cui k0 freespace wavenumber
k0=2*pi/lambda0;
%-Banda
Bw= 100e6;  %MHz
%-Guadagno
Gain=1; %dB 

%%Specifiche del substrato :
Er=4.4;
h=0.0016; %m

% Note: -metallizzazione PEC di spessore infinitesimo

%% Design della geometria 
%Calcolo della larghezza W
W=(v0/(2*f0))*sqrt(2/(Er+1))

%Calcolo Ereff dalla Balanis (14-1)3
Ereff= (Er+1)/2+((Er-1)/2)*((1+12*h/W)^(-1/2))

%Calcolo lunghezza  L
%Se considero gli effetti ai bordi, Leff= L+2deltaL (fringing)
%Dove deltaL è approssimata dalla Balanis (14-2)ed è funzione del rapporto fra
%W e l'altezza h del substrato 

deltaL=h*0.412*((Ereff+0.3)*(W/h+0.264))/((Ereff-0.258)*(W/h+0.8));

%Ricavo L dalla Balanis (14-7)  

L=(1/(2*f0*sqrt(Ereff)*sqrt(mu0*e0)))-2*deltaL



%% Z-Match
%% Lambda/4 Transformer
% References Balanis, DocB

%  Calcolo la resistenza in ingresso (considero il bordo) della patch
%  Considero una linea feed di microstriscia  dimensionata a 50 ohm
%  Dimensiono un trasformatore lamba/4 per adattare i due sistemi

%CALCOLO DI Rin0 (Resistenza di ingresso al bordo)
%della patch usando la Balanis 14-20

%slotrate è il rapporto fra W e lamba0--> serve solo per scegliere quale
%funzione a tratti scegliere
slotrate=W/lambda0



% dalla Balanis (14-12)Svolgendo l'integrale I1 ricavo la  prima conduttanza dell modello in figura (14.8b)
%La seconda conduttanza è uguale, quindi il parallelo è dato da 2*G1
%

x=k0*(W);
i1=-2+cos(x)+(x*sinint(x))+(sin(x)/x);
G1=i1/(120*pi*pi);          %Conductance


% dalla Balanis (14-17) ricavo Rin calcolata al bordo della patch e la
% chiamo Rin0
%G12 induttanza mutua fra i due slot, ed è data dalla Balanis (14-18a)


a=@(th)(((sin((x./2).*cos(th))./cos(th)).^2).*(besselj(0,(k0.*L.*sin(th)))).*(sin(th)).^3);

a1=integral(a,0,pi);

G12=a1/(120*pi*pi);     %in siemens


Rin0= 1/(2*(G1+G12))


%CALCOLO DELLA MICROSTRISCIA
%Definisco una  linea di microstriscia di dimensione Wm tale che sia a Zc=
%50 ohm

Wm=0.00303;%( per avere una microstriscia a 50 ohm)


% CALCOLO DEL TRASFORMATORE LAMBDA/4
%
%Ricavo l'impedenza Z0 del trasformatore 
Z0=sqrt(50*111)

%per ricavare la lunghezza d'onda lungo la microstriscia calcolo la
%relativa Eeff, tenendo conto che non sono più nella patch ma nella niea di
%microstriscia, quindi la larghezza è Wm
Ereff_mstrip= (Er+1)/2+((Er-1)/2)*((1+12*h/Wm)^(-1/2))
lambdag= lambda0/sqrt(Ereff_mstrip)
%il trasformatore ha L pari a lambdag/4
L_transf=lambdag/4

%NOTABENE   Wm è la larghezza della microstrip a 50 ohm, mentre la W del
%trasformatore la calcolo con un tool imponendo Z0

%% Inset feed
%  Considero una linea feed di microstriscia  dimensionata a 50 ohm
%  Calcolo la resistenza in ingresso della patch scegliendo un punto y0 in

%  cui l'impedenza matcha 50 ohm

% invertendo la Balanis (14-20a), imponendo l'impedenza R(y0)=50 ohm, ricavo 
%
 y0= L/pi*acos(sqrt(50/111))
 
%% PLOT TABLE
% 

%  Patch = table(W,L,h,f0,Er,Ereff);
 fig=uifigure('Name','Design della Patch');
 fig.Position = [500 500 740 210];
 uit1=uitable(fig,'Data',[W L h f0 Er Ereff],'ColumnName',{'W'; 'L'; 'h'; 'f0'; 'Er'; 'Ereff'});
 uit1.Position = [20 130 680 60];
 uit2=uitable(fig,'Data',[Rin0 y0 ],'ColumnName',{'Rin0'; 'y0 a  Rin=50 ohm'});
  uit2.Position = [20 70 680 60];
   uit2=uitable(fig,'Data',[Z0 L_transf Ereff_mstrip],'ColumnName',{'Z0 del trasformatore'; 'L del trasformatore( lambdag/4)';'Ereff microstriscia'});
  uit2.Position = [20 10 680 60];


%% ARRAY
% 
%
% Coded by Chan Sokthai (sokthai@msn.com) -  Department of Engineering, University of Fukui

% element numbers
N = 25;
% element spacing
d = 5;
% theta zero direction
% 90 degree for braodside, 0 degree for endfire.
theta_zero = 0;
An = 1;
j = sqrt(-1);
AF = zeros(1,360);
for theta=1:360
    
    % change degree to radian
    deg2rad(theta) = (theta*pi)/180;
    
    %array factor calculation
    for n=0:N-1
        AF(theta) = AF(theta) + An*exp(j*n*2*pi*d*(cos(deg2rad(theta))-cos(theta_zero*pi/180))) ;
    end
    AF(theta) = abs(AF(theta));
    
end
% plot the array factor
% figure;
% polar(deg2rad,AF);

%%