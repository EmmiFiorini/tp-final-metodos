%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%% Ejercicio Nro. 3 TP#FINAL Métodos Numéricos
%%%
%%% Aproximación Lineal por Mínimos Cuadrados:
%%%
%%%         y = Ax + B
%%%
%%% function [A,B,CC] = Ajuste_Lineal_MC(X,Y)
%%%
%%% Parámetros de Entrada:
%%%       X = vector Nx1 con las abscisas de los pares de puntos
%%%       Y = vector Nx1 con las ordenadas de los pares de puntos
%%%
%%% Parámetros de Salida:
%%%
%%%      A = pendiente del Ajuste Lineal
%%%      B = Ordenada al origen
%%%      CC = Coeficiente de Correlación
%%%
%%% Dr. Ing. Franco Pessana
%%% FICEN
%%% Universidad Favaloro
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
clc

% puntos a aproximar
X=[0;1.5;3;4.5;6;7.5;9;10.5;12;13.5;15;16.5;18;19.5];
Y=[19.7;21.6;29.1;27.4;32.1;36;40.2;47.9;54.2;65.2;70.1;77.8;93.9;105.4];
x_eval= linspace(0, 20, 1000);

%% Aproximacion Lineal a)

[A_1,B_1,CC] = Ajuste_Lineal_MC(X,Y); 
Y_AL_a = Funcion_1(x_eval,A_1,B_1);
F1=figure(1);
set(F1,'position',[140 20 1100 650],'Menubar','none',...
        'NumberTitle','off','name', 'Ejercicio 3a TP#FINAL. M�todos Num�ricos');
plot(X,Y,'*r',x_eval,Y_AL_a, 'b') % Nube de puntos originales
xlabel('D [mm]')
ylabel('P(D) [mmHg]')
legend ('Pk', 'Aproximación lineal Pk')
grid

%% Aproximacion Exponencial; b)

Y_log = log(Y);
[A_2,B_2,CC] = Ajuste_Lineal_MC(X,Y_log);
p0 = exp(B_2);
Y_AL_b =  p0 * exp(A_2 * x_eval); % Recta de Ajuste evaluada en los puntos xk
F2=figure(2);
set(F2,'position',[140 20 1100 650],'Menubar','none',...
        'NumberTitle','off','name', 'Ejercicio 3b TP#FINAL. M�todos Num�ricos');
plot(X,Y,'*r',x_eval,Y_AL_b, 'b') % Nube de puntos originales
xlabel('D [mm]')
ylabel('P(D) [mmHg]')
legend ('Pk', 'Aproximación exponencial Pk')
grid

%% Lagrange; c)

C_3 = Interp_Lagrange(X,Y');
Y_AL_c = @(X) Eval_Polinomio(X,C_3);% Armado del polinomio
F3=figure(3);
set(F3,'position',[140 20 1100 650],'Menubar','none',...
        'NumberTitle','off','name', 'Ejercicio 3c TP#FINAL. M�todos Num�ricos');
plot(X,Y,'*r',x_eval, Y_AL_c(x_eval),'b')
xlabel('D [mm]')
ylabel('P(D) [mmHg]')
legend ('Pk', 'Polinomio Interpolador de Lagrange Pk')
grid 

%% Interpolador Newton; d)

C_4 = Interp_Newton(X,Y);
Y_AL_d = @(X)Eval_Polinomio(X,C_4);% Armado del polinomio interpolador de Newton
F4=figure(4);
set(F4,'position',[140 20 1100 650],'Menubar','none',...
        'NumberTitle','off','name', 'Ejercicio 3d TP#FINAL. M�todos Num�ricos');
plot(X,Y,'*r',x_eval, Y_AL_d(x_eval),'b')
xlabel('D [mm]')
ylabel('P(D) [mmHg]')
legend ('Pk', 'Polinomio Interpolador de Newton Pk')
grid 

%% Spline Cubicas; e)


MC = Spline_Cubica(X,Y);
N = size(MC,1);
L = 14; % Cantidad de puntos entre intervalos de X discretos
x = zeros(N,L+1); % Matriz de abscisas (101 puntos por cada intervalo)
for k=1:N
    d = (X(k+1)-X(k))/L; % Armado de abscisas temporales para los 
    x(k,:) = X(k):d:X(k+1); % Polinomios spline cúbicos
end    
F5=figure(5);
set(F5,'position',[140 20 1100 650],'Menubar','none',...
        'NumberTitle','off','name', 'Ejercicio 3e TP#FINAL. Métodos Numéricos');
hold on
plot(X,Y,'*r')
xlabel('D [mm]')
ylabel('P(D) [mmHg]')
grid 
Y_AL_e=zeros(N,length(x));
for k=1:N
    Y_AL_e(k,:) = Eval_Polinomio2(x(k,:),MC(k,:),X(k));
    plot(x(k,:),Y_AL_e(k,:))
end
legend ('Pk', 'Evolucion Spline Cubicas')
hold off

%% Interpolacion de Richardson; f)

%i.
D_m = 12.5;
Delta = 1e-13;
Tol = 1e-13;
[D,err,relerr,n] = Extrapolacion_Richardson_O2n(@(X) A_1*X+B_1,D_m,Delta,Tol);% referencia a esa funci�n
F6=figure(6);
set(F6,'position',[140 20 1100 650],'Menubar','none',...
        'NumberTitle','off','name', 'Ejercicio 3f (i) TP#FINAL. M�todos Num�ricos');
Elasticidad1 = D(end,end);
tangente = Elasticidad1*(X - D_m) + A_1*D_m+B_1;
plot(X,Y,'*r',x_eval,Y_AL_a, 'k',X, tangente, 'y--') % Nube de puntos originales
xlabel('D [mm]')
ylabel('P(D) [mmHg]');
legend('Pk','Aproximacion Lineal' ,'Tangente en D_M');
grid 

 %ii.
[D,err,relerr,n] = Extrapolacion_Richardson_O2n(@(X)  p0 * exp(A_2 * X),D_m,Delta,Tol);% referencia a esa función
F7=figure(7);
set(F7,'position',[140 20 1100 650],'Menubar','none',...
        'NumberTitle','off','name', 'Ejercicio 3f (ii) TP#FINAL. M�todos Num�ricos');
Elasticidad_2 = D(end,end);
tangente = Elasticidad_2*(X - D_m) +  p0 * exp(A_2 * D_m);
plot(X,Y,'*r',x_eval,Y_AL_b, 'k', X, tangente, 'm--') % Nube de puntos originales
xlabel('D [mm]')
ylabel('P(D) [mmHg]');
legend('Pk','Aproximacion Exponencial' ,'Tangente en D_M');
grid 

%iii.
[D,err,relerr,n] = Extrapolacion_Richardson_O2n(Y_AL_c,D_m,Delta,Tol);% referencia a esa función
F8=figure(8);
set(F8,'position',[140 20 1100 650],'Menubar','none',...
        'NumberTitle','off','name', 'Ejercicio 3f (iii) TP#FINAL. M�todos Num�ricos');
Elasticidad_3 = D(end,end);
tangente = Elasticidad_3*(X - D_m) + Y_AL_c(D_m);
plot(X,Y,'*r',x_eval,Y_AL_c(x_eval), 'k',X, tangente, 'm--') % Nube de puntos originales
xlabel('D [mm]')
ylabel('P(D) [mmHg]');
legend('Pk','Aproximacion Polinomial' ,'Tangente en D_M');
grid 

%iv
[D,err,relerr,n] = Extrapolacion_Richardson_O2n(Y_AL_d,D_m,Delta,Tol);% referencia a esa función
F9=figure(9);
set(F9,'position',[140 20 1100 650],'Menubar','none',...
        'NumberTitle','off','name', 'Ejercicio 3f (iv) TP#FINAL. M�todos Num�ricos');
Elasticidad_4 = D(end,end);
tangente = Elasticidad_4*(X - D_m) + Y_AL_d(D_m);
plot(X,Y,'*r',x_eval,Y_AL_d(x_eval), 'k',X, tangente, 'm--') % Nube de puntos originales
xlabel('D [mm]')
ylabel('P(D) [mmHg]');
legend('Pk','Aproximacion Polinomial' ,'Tangente en D_M');
grid

%v
[D,err,relerr,n] = Extrapolacion_Richardson_O2n(Y_AL_e,D_m,Delta,Tol);% referencia a esa función
F10=figure(10);
set(F10,'position',[140 20 1100 650],'Menubar','none',...
        'NumberTitle','off','name', 'Ejercicio 3f (v) TP#FINAL. M�todos Num�ricos');
Elasticidad_5 = D(end,end);
tangente = Elasticidad_5*(X - D_m) + Y_AL_e(D_m);
plot(X,Y,'*r',x_eval,Y_AL_e(x_eval), 'k',X, tangente, 'm--') % Nube de puntos originales
xlabel('D [mm]')
ylabel('P(D) [mmHg]');
legend('Pk','Aproximacion Polinomial' ,'Tangente en D_M');
grid


