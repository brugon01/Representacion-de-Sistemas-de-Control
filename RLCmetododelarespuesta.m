%SISTEMA DE DOS VARIABLES DE ESTADO - CIRCUITO RLC
clear all; close all; clc;

%PUNTO 3: M�TODO DE LA RESPUESTA AL ESCAL�N

%En el archivo Curvas_Medidas_RLC.xls (datos en la hoja 1 y etiquetas en la
%hoja 2)encontrar�n las series de datos que deber�an emplear para deducir 
...los valores de R, L y C del circuito. Emplear el m�todo de la respuesta 
...al escal�n, tomando como salida la tensi�n en el capacitor.

%Importaci�n y Graficaci�n a partir de los datos:

Data=readmatrix('Curvas_Medidas_RLC');      %Data guarda en columnas los datos otorgados
Corriente=Data(:,2);                        %Se define a la corriente como la columna 2 de los datos
Vcapacitor=Data(:,3);                       %Se define a VC como la columna 3 de los datos
Time=Data(:,1);                             %Se define a Time como la comna 1 de los datos

figure('Name','Gr�ficas obtenidas desde los datos')   %Gr�fica de corriente y Vc vs tiempo
subplot(2,1,1)
plot(Time,Corriente,'b');title('Corriente');
subplot(2,1,2)
plot(Time,Vcapacitor,'r');title('Voltaje en el capacitor')

%Aplicaci�n del m�todo de la respuesta al escal�n - Chen
%Salida: Voltaje sobre el capacitor.

K=12;                                      %Amplitud escal�n
z=145;                                     %Ubicaci�n de la muestra inicial
d=100;                                     %distancia entre muestras

%Se toman tres muestras de la salida en funci�n del tiempo:
y1=Vcapacitor(z);
t1=Time(z);

y2=Vcapacitor(z+d);
t2=Time(z+d);

y3=Vcapacitor(z+d+d);
t3=Time(z+d+d);

%Se define el valor de los Kn:
k1=y1/K - 1;
k2=y2/K - 1;
k3=y3/K - 1;

%Se define el valor del par�metro b
b= 4*k1^3*k3 - 3*k1^2*k2^2 - 4*k2^3+k3^2 + 6*k1*k2*k3;

%Se calculan los alfa
a1=(k1*k2 + k3 - sqrt(b)) / (2 * (k1^2 + k2));
a2=(k1*k2 + k3 + sqrt(b)) / (2 * (k1^2 + k2));

%Se calcula el par�metro Beta
Beta=(k1+a2) / (a1-a2);

%C�lculo de constantes de tiempo

T1= -t1/log(a1);
T2= -t1/log(a2);

T3= Beta*(T1 - T2) + T1;     


%Definici�n de la funci�n de transferencia:
s=tf('s');

%G=((12)*(T3*s+1))/((T1*s +1)*(T2*s +1)) 
G=12/((T1*s +1)*(T2*s +1)) %Elimino el cero para luego despejar par�metros 


figure('Name','Modelo obtenido con m�todo de Chen');
hold on
step(G*exp(-s*0.01),'y');   %Se aplica retardo para asemejarla. 12 es la amplitud del escal�n
plot(Time,Vcapacitor,'r');
title('Tensi�n en el capacitor');
legend({'Por Chen','Observada'},'Location','northwest','Orientation','horizontal')
hold off

%Funci�n de Transferencia Vc/Vin= 1 / [L*C*s^2 + R*C*s + 1]
%A partir de la FdT:

LC=7.533e-06

RC=0.00744

%Iterando el valor de R comercial hasta que las ecuaciones se asemejen:

R=220
C=RC/R
L=LC/C

Gp=1 / [L*C*s^2 + R*C*s + 1]

figure('Name','Gr�fica con los valores de R, L y C encontrados');
hold on
step(G*exp(-s*0.01),'y');   
plot(Time,Vcapacitor,'r');
step(12*Gp*exp(-s*0.01),'b');            %Azul tapa a la amarilla
title('Tensi�n en el capacitor');
legend({'Encontrada','Medida','Con valores encontrados de R L y C'},'Location','northwest','Orientation','horizontal')
hold off

%PUNTO 4 
%Una vez determinados los par�metros R, L y C, 
...emplear la serie de corriente desde 0.05seg en adelante 
...para validar el resultado

%Se repite la misma metodolog�a de simulacion que en los puntos 1 y 2:

%Par�metros necesarios para simular:
T=0.1; At=1e-4; Kmax=T/At; t=linspace(0,T,Kmax); 
%T: tiempo de simulaci�n
%At: Per�odo de muestreo de la simulaci�n
%Kmax: cantidad de muestras
I=zeros(1,Kmax);Vc=zeros(1,Kmax);u=linspace(0,0,Kmax); 
%Matrices/vetctores fila de Kmax lugares vac�os
%Par�metros f�sicos:
Vin=0;
%R L y C ya est�n cargados
%Condiciones iniciales:
I(1)=0; Vc(1)=0; 
u(1)=Vin; %Se define escal�n de entrada u 


%Modelo Lineal del sistema:

A=[-R/L -1/L ; 1/C  0]; 
B=[1/L; 0];
E=[R 0];

x=[I(1) Vc(1)]';             %vector de estado
Xop=[0 0]';                  %punto de operacion 

Il(1)=0;                     % vectores de variables de inter�s
Vcl(1)=0;                    ...inicializados en 0


ii=0;
 for i=1:Kmax-1                             %itera desde 1 hasta kmax-1 (2ms)
     ii=ii+At;                              %ii:tiempo
     
     if(ii>=0)                              %Vin comienza siendo 0v
         Vin = 0;                      
     end 
     if(ii>=0.01)                           %Vin pasa a ser 12v
         Vin = 12;                      
     end 
     if(ii>=0.05)                           %Se invierte Vin a -12v
         Vin = Vin*-1;                      
     end                                    
     
     u(i)=Vin;                              %asignamos el valor de entrada u
     
     %sistema modelado por las ecuaciones diferenciales:
     Ip =-(R/L)*I(i)-(1/L)*Vc(i)+(1/L)*u(i);%EDO de ipunto
     Vcp = 1/C * I(i);                      %EDO de vpunto
     I(i+1)=I(i)+Ip*At;                     %Se aplica m�todo num�rico Euler
     Vc(i+1)=Vc(i)+Vcp*At;                  
     
     %sistema modelado por su ecuaci�n matricial-vectorial:           
     xp=A*(x-Xop)+B*u(i);                   %se calcula el vector xp
     x=x+xp*At;                             %vector de estados actualiza su valor
     Y=E*x;                                 %vector de salida actualiza su valor
     Il(i+1)=x(1);                          %Corriente 
     Vcl(i+1)=x(2);                         %Tension del capacitor
    
 end
 
 figure('Name','Validaci�n')
 
 subplot(3,1,1);
 plot(t,u,'b');title(' tension de entrada Vin');
 subplot(3,1,2);
 hold on
 plot(t,Il,'b');plot(Time,Corriente,'r');title('corriente i');
 legend({'M�todo de la Rta al escal�n','Real'},'Location','northeast','Orientation','horizontal')
 hold off
 subplot(3,1,3);
 hold on
 plot(t,Vcl,'b');plot(Time,Vcapacitor,'r');title('tension Vc');
 legend({'M�todo de la Rta al escal�n','Real'},'Location','northeast','Orientation','horizontal')
 hold off

 disp('END')













