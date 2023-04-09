%CASO DE ESTUDIO 2 - SISTEMA DE TRES VARIABLES DE ESTADO - MOTOR CC

clear all; close all; clc;

%PUNTO 1:
%Implementar un algoritmo de simulaci�n para inferir el comportamiento 
%de las variables inter�s mediante integraci�n Euler con deltat=10e-7 s.

%Se procede a aplicar la metodolog�a utilizada en el Caso 1, puntos 1 y 2:

%Par�metros necesarios para simular:
T=5; At=1e-7; Kmax=T/At; t=linspace(0,T,Kmax);Va=linspace(0,0,Kmax);Tl=linspace(0,0,Kmax);
Ia=zeros(1,Kmax); wr=zeros(1,Kmax); Zeta=zeros(1,Kmax);

%Par�metros f�sicos:
Laa=366e-6;
J=5e-9;
Ra=55.6; 
Bm=0; 
Ki=6.49e-3;
Km=6.53e-3;

%Condiciones iniciales:
Ia(1)=0; Zeta(1)=0; wr(1)=0; Va(1)=0; Tl(1)=0; %Vectores de inter�s

%Modelo Lineal del sistema:

%Se define las variables de estado:
...x1:ia
...x2:wr
...x3:zeta
...modelo objetivo: xp=A.x+B.u
...Se identifica como entradas u al torque de carga Tl y a la tensi�n
...u1=Va
...u2=Tl
...de alimentaci�n Va
%A partir de las edos del modelo otorgado y aplicando linealizaci�n por
%Taylor

A=[-Ra/Laa -Km/Laa 0 ; Ki/J   -Bm/J  0 ; 0       1    0 ];

B1=[(1/Laa) 0   0]';         %Matriz B para entrada Va
B2=[ 0    -1/J  0]';         %Matriz B para entrada TL
 

Xop=[0     0      0]';       %punto de operacion 


%Simulaci�n utilizando Euler:
Ia1(1)=0;                    %Matrices de simulacion comienzan siendo 0
wr1(1)=0;
zeta1(1)=0;
     
Vin=0;                       %Par�metros de entrada iniciales
TL=0;

ii=0;
 for i=1:Kmax-1                             %itera desde 1 hasta kmax-1 (5s)
     ii=ii+At;                              %ii:tiempo
     
     if(ii>=0.1)                            %Alimento al motor con 12v
         Vin=12;                            
     end
     
     if(ii>=0.15)                           %Para t>=0.15s
         TL = i/10000000000;                %Tl es aumenta linealmente con el tiempo
     end                                    
     
     Va(i)=Vin;                             %Se llenan los vectores de entrada
     Tl(i)=TL;                              
     
     if (Tl(i)~=0 && wr(i)<0.001 && wr(i)>-0.001) %Se busca y guarda el TLMAX
         TLMAX=Tl(i);
     end
     
%      if (Tl(i)==TLMAX) %CORREGIR, SE DEBE RESIMULAR CON TLMAX DE ENTRADA Y OBSERVAR EL PICO DE CORRIENTE
%          iMAX=Ia(i);
%      end
     
     %sistema modelado por las ecuaciones diferenciales:
     Iap = -Ra/Laa*Ia(i) - Km/Laa*wr(i) + 1/Laa*Va(i) ;%EDO de ipunto
     wrp =   Ki/J *Ia(i) -   Bm/J*wr(i) -   1/J*Tl(i) ;%EDO de wrpunto
     zetap=     wr(i);                                 %EDO de zetapunto
     
     Ia(i+1)=Ia(i)+Iap*At;                             %Se aplica m�todo num�rico Euler
     wr(i+1)=wr(i)+wrp*At; 
     Zeta(i+1)=Zeta(i)+zetap*At;
     
     %Sistema modelado por su ecuacion matricial-vectorial de estado
     x=[Ia(i) wr(i) Zeta(i)]';
     xp=A *(x-Xop) + B1*Va(i) + B2*Tl(i);      %se calcula el vector xp
     x=x+xp*At;                                %vector de estados actualiza su valor
     Ia1(i+1)=x(1);                            %Matrices necesarias para graficar en el tiempo
     wr1(i+1)=x(2);
     zeta1(i+1)=x(3);
     
     
     
 end

 
 figure(1)
 subplot(5,1,1);
 plot(t,Va,'b');title(' Va');
 subplot(5,1,2);
 plot(t,Ia1,'r');title(' ia');
 subplot(5,1,3);
 plot(t,Tl,'b');title(' Tl');
 subplot(5,1,4);
 plot(t,wr1,'b');title(' wr');
 subplot(5,1,5);
 plot(t,zeta1,'b');title(' zeta');
 

 disp('EL TORQUE M�XIMO QUE PUEDE SOPORTAR EL MOTOR CON Va:12V ES: ')
 disp(TLMAX)

 
 %PUNTO 2
 %C�lculo de la corriente m�xima = pico de corriente cuando TL=TLMAX
 %Se realiza la misma simulacion que en el caso anterior pero con TL fijo
 
 %Par�metros necesarios para simular:
T=5; At=1e-7; Kmax=T/At; t=linspace(0,T,Kmax);Va=linspace(0,0,Kmax);Tl=linspace(0,0,Kmax);
Ia=zeros(1,Kmax); wr=zeros(1,Kmax); Zeta=zeros(1,Kmax);

%Par�metros f�sicos:
Laa=366e-6;
J=5e-9;
Ra=55.6; 
Bm=0; 
Ki=6.49e-3;
Km=6.53e-3;

%Condiciones iniciales:
Ia(1)=0; Zeta(1)=0; wr(1)=0; Va(1)=0; Tl(1)=0; %Vectores de inter�s

%Modelo Lineal del sistema:

%Se define las variables de estado:
...x1:ia
...x2:wr
...x3:zeta
...modelo objetivo: xp=A.x+B.u
...Se identifica como entradas u al torque de carga Tl y a la tensi�n
...u1=Va
...u2=Tl
...de alimentaci�n Va
%A partir de las edos del modelo otorgado y aplicando linealizaci�n por
%Taylor

A=[-Ra/Laa -Km/Laa 0 ; Ki/J   -Bm/J  0 ; 0       1    0 ];

B1=[(1/Laa) 0   0]';         %Matriz B para entrada Va
B2=[ 0    -1/J  0]';         %Matriz B para entrada TL
 

Xop=[0     0      0]';       %punto de operacion 


%Simulaci�n utilizando Euler:
Ia1(1)=0;                    %Matrices de simulacion comienzan siendo 0
wr1(1)=0;
zeta1(1)=0;
     
Vin=0;                       %Par�metros de entrada iniciales
TL=0;

ii=0;
 for i=1:Kmax-1                             %itera desde 1 hasta kmax-1 (5s)
     ii=ii+At;                              %ii:tiempo
     
     if(ii>=0.1)                            %Alimento al motor con 12v
         Vin=12;                            
     end
     
     if(ii>=0.1)                           %Para t>=0.15s
         TL = TLMAX; 
     end                                    
     
     Va(i)=Vin;                             %Se llenan los vectores de entrada
     Tl(i)=TL;                              
     
     
     
     if (i>2 && Ia(i)>Ia(i-1))
          iMAX=Ia(i);
     end
     
     %sistema modelado por las ecuaciones diferenciales:
     Iap = -Ra/Laa*Ia(i) - Km/Laa*wr(i) + 1/Laa*Va(i) ;%EDO de ipunto
     wrp =   Ki/J *Ia(i) -   Bm/J*wr(i) -   1/J*Tl(i) ;%EDO de wrpunto
     zetap=     wr(i);                                 %EDO de zetapunto
     
     Ia(i+1)=Ia(i)+Iap*At;                             %Se aplica m�todo num�rico Euler
     wr(i+1)=wr(i)+wrp*At; 
     Zeta(i+1)=Zeta(i)+zetap*At;
     
     %Sistema modelado por su ecuacion matricial-vectorial de estado
     x=[Ia(i) wr(i) Zeta(i)]';
     xp=A *(x-Xop) + B1*Va(i) + B2*Tl(i);      %se calcula el vector xp
     x=x+xp*At;                                %vector de estados actualiza su valor
     Ia1(i+1)=x(1);                            %Matrices necesarias para graficar en el tiempo
     wr1(i+1)=x(2);
     zeta1(i+1)=x(3);
     
     
     
 end

 
 figure(2)
 subplot(5,1,1);
 plot(t,Va,'b');title(' Va');
 subplot(5,1,2);
 plot(t,Ia1,'r');title(' ia');
 subplot(5,1,3);
 plot(t,Tl,'b');title(' Tl');
 subplot(5,1,4);
 plot(t,wr1,'b');title(' wr');
 subplot(5,1,5);
 plot(t,zeta1,'b');title(' zeta');
 
 disp('LA CORRIENTE M�XIMA CUANDO TL=TLMAX ')
 disp(iMAX)
 
 



















