%SISTEMA DE DOS VARIABLES DE ESTADO - CIRCUITO RLC
clear all; close all; clc;

%PUNTO 1 y 2:
%Asignar valores a R=4,7K?, L=10?Hy, y C=100nF. 
...Obtener simulaciones que permitan estudiar la dinámica del sistema
...con una entrada de tensión escalón de 12V, que cada 
...1ms cambia de signo.
%Asignar valores a R=5,6K?, L=10?Hy, y C=100nF;
...repetir lo anterior para comparar el 
...resultado y verificar la correcta simulación.

%Parámetros necesarios para simular:
T=2e-3; At=1e-9; Kmax=T/At; t=linspace(0,T,Kmax); 
%T: tiempo de simulación
%At: Período de muestreo de la simulación
%Kmax: cantidad de muestras
I=zeros(1,Kmax);Vc=zeros(1,Kmax);u=linspace(0,0,Kmax); 
%Matrices/vetctores fila de Kmax lugares vacíos

%Parámetros físicos:
Vin=12;
%R=4.7e3; L=10e-6; C=100e-9;      %Para PUNTO 1
R=5.6e3; L=10e-6; C=100e-9;      %Para PUNTO 2

%Condiciones iniciales:
I(1)=0; Vc(1)=0; 
u(1)=Vin; %Se define escalón de entrada u 


%Modelo Lineal del sistema:

A=[-R/L -1/L ; 1/C  0]; 
B=[1/L; 0];
E=[R 0];

x=[I(1) Vc(1)]';             %vector de estado
Xop=[0 0]';                  %punto de operacion 

Il(1)=0;                     % vectores de variables de interés
Vcl(1)=0;                    ...inicializados en 0
Vr(1)=0;



ii=0;
 for i=1:Kmax-1                             %itera desde 1 hasta kmax-1 (2ms)
     ii=ii+At;                              %ii:tiempo
     
     if(ii>=1e-3)                           %Se cumplió el ms, entonces se invierte Vin
         ii=0;
         Vin = Vin*-1;                      
     end                                    
     
     u(i)=Vin;                              %asignamos el valor de entrada u
     
     %sistema modelado por las ecuaciones diferenciales:
     Ip =-(R/L)*I(i)-(1/L)*Vc(i)+(1/L)*u(i);%EDO de ipunto
     Vcp = 1/C * I(i);                      %EDO de vpunto
     I(i+1)=I(i)+Ip*At;                     %Se aplica método numérico Euler
     Vc(i+1)=Vc(i)+Vcp*At;                  
     
     %sistema modelado por su ecuación matricial-vectorial:           
     xp=A*(x-Xop)+B*u(i);                   %se calcula el vector xp
     x=x+xp*At;                             %vector de estados actualiza su valor
     Y=E*x;                                 %vector de salida actualiza su valor
     Vr(i+1)=Y(1);                          %Tension en la R
     Il(i+1)=x(1);                          %Corriente 
     Vcl(i+1)=x(2);                         %Tension del capacitor
    
 end
 
 figure(1)
 subplot(2,2,1);
 plot(t,u,'b');title(' tension de entrada Vin');
 subplot(2,2,2);
 plot(t,Il,'r');title('corriente i');
 subplot(2,2,3);
 plot(t,Vcl,'b');title('tension Vc');
 subplot(2,2,4);
 plot(t,Vr,'b');title(' tension de salida VR');
 
 
 