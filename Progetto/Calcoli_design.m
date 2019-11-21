
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

%Valore W corretto
W=0.04
%Calcolo Ereff dalla Balanis (14-1)3
Ereff= (Er+1)/2+((Er-1)/2)*((1+12*h/W)^(-1/2))
lambdag_patch= lambda0/sqrt(Ereff)
%Calcolo lunghezza  L
%Se considero gli effetti ai bordi, Leff= L+2deltaL (fringing)
%Dove deltaL è approssimata dalla Balanis (14-2)ed è funzione del rapporto fra
%W e l'altezza h del substrato 

deltaL=h*0.412*((Ereff+0.3)*(W/h+0.264))/((Ereff-0.258)*(W/h+0.8));

%Ricavo L dalla Balanis (14-7)  

L=(1/(2*f0*sqrt(Ereff)*sqrt(mu0*e0)))-2*deltaL

%Valore L corretto
L=0.0288

%% Z-Match
%% Lambda/4 Transformer
% References Balanis, DocB

%  Calcolo la resistenza in ingresso (considero il bordo) della patch
%  Considero una linea feed di microstriscia  dimensionata a 50 ohm
%  Dimensiono un trasformatore lamba/4 per adattare i due sistemi

%CALCOLO DI Rin0 (Resistenza di ingresso al bordo)
%della patch usando la Balanis 14-20




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
%https://www.emtalk.com/mscalc.php?er=4.4&h=1.6&h_units_list=hmm&f=2.45&Zo=74.4983&EL=0&Operation=Synthesize&Wa=1.4527986584&W_units_list=Wmm&La=0.172456338325&L_units_list=Lmm

Wm=0.00305;%( per avere una microstriscia a 50 ohm)

Ereff_mstrip= (Er+1)/2+((Er-1)/2)*((1+12*h/Wm)^(-1/2))

lambdag_mstrip= lambda0/sqrt(Ereff_mstrip)

%CALCOLO DEL TRASFORMATORE LAMBDA/4
%Za valore di impedenza da adattare
Za=111;
%Ricavo l'impedenza Zt del trasformatore 
Zt=sqrt(50*Za)
W_transf= 0.0013 %( per avere un trasformatore  a Zt ohm)
Ereff_transf= (Er+1)/2+((Er-1)/2)*((1+12*h/W_transf)^(-1/2))
lambdag_transf= lambda0/sqrt(Ereff_transf)
%il trasformatore ha L pari a lambdag/4
L_transf=lambdag_transf/4


%%TRASFORMATORE 2 LAMBDA/4
%CALCOLO DEL TRASFORMATORE LAMBDA/4
%Za valore di impedenza da adattare
Za=111;
%Ricavo l'impedenza Zt del trasformatore 
Zt2=sqrt(50*100)
W_transf2= 0.00165 %( per avere un trasformatore  a Zt ohm)
Ereff_transf2= (Er+1)/2+((Er-1)/2)*((1+12*h/W_transf)^(-1/2))
lambdag_transf2= lambda0/sqrt(Ereff_transf)
%il trasformatore ha L pari a lambdag/4
L_transf2=lambdag_transf2/4


%% Inset feed
%  Considero una linea feed di microstriscia  dimensionata a 50 ohm
%  Calcolo la resistenza in ingresso della patch scegliendo un punto y0 in

%  cui l'impedenza matcha 50 ohm

% invertendo la Balanis (14-20a), imponendo l'impedenza R(y0)=50 ohm, ricavo 
%
y0= L/pi*acos(sqrt(50/Rin0))

%% PLOT TABLE
% 

%  Patch = table(W,L,h,f0,Er,Ereff);
 fig=uifigure('Name','Design della Patch');
 fig.Position = [500 500 740 240];
 uit1=uitable(fig,'Data',[W L h f0 lambdag_patch Er Ereff],'ColumnName',{'W'; 'L'; 'h'; 'f0';'lambdag'; 'Er'; 'Ereff'});
 uit1.Position = [20 160 680 40];
 uit2=uitable(fig,'Data',[Rin0 y0 ],'ColumnName',{'Rin0'; 'Inset y0 a  Rin=50 ohm'});
 uit2.Position = [20 120 680 40];
 uit3=uitable(fig,'Data',[ Ereff_mstrip Wm lambdag_mstrip],'ColumnName',{'Ereff microstriscia';'W microstriscia'; 'lambdag microstriscia'});
 uit3.Position = [20 80 680 40];
 uit4=uitable(fig,'Data',[Zt lambdag_transf W_transf L_transf Ereff_transf],'ColumnName',{'Zt del trasformatore';'lambdag trasformatore';'W trasformatore ';'L del trasformatore( lambdag/4)';'Ereff trasf'});
 uit4.Position = [20 40 680 40];

fig=uifigure('Name','Design del FEED');
 fig.Position = [500 500 740 240];
  uit5=uitable(fig,'Data',[Zt2 lambdag_transf2 W_transf2 L_transf2 Ereff_transf2],'ColumnName',{'Zt del trasformatore';'lambdag trasformatore';'W trasformatore ';'L del trasformatore( lambdag/4)';'Ereff trasf'});
 uit5.Position = [20 40 680 40];
%% ARRAY
clear all;
close all;
clc;
%
% Coded by Chan Sokthai (sokthai@msn.com) -  Department of Engineering, University of Fukui
f0= 2.45e9
c=3*10^8;
lambda0= c/f0

% element numbers
N =4;
% element spacing
d = 0.5;
% theta zero direction
% 90 degree for braodside, 0 degree for endfire.
theta_zero = 90;
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
 figure;
 polar(deg2rad,AF);

%%