%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%% Ejercicio Nro. 3 TP#FINAL Metodos Numericos
%%%
%%% Juana Kallis, Emma fiorini y Agustina Vidaurreta
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
clc

% puntos a aproximar
X=[0;1.5;3;4.5;6;7.5;9;10.5;12;13.5;15;16.5;18;19.5];
Y=[19.7;21.6;29.1;27.4;32.1;36;40.2;47.9;54.2;65.2;70.1;77.8;93.9;105.4];
x_eval= linspace(0, 20, 1000);
N = length(X);
F1 = figure(1);
set(F1, 'Position', [100, 100, 1200, 700],'Menubar','none',...
   'NumberTitle','off','name', 'Ejercicio 3a-e TP#FINAL- M閠odos Num閞icos');
%% Aproximacion Lineal a)

[A_1,B_1,CC] = Ajuste_Lineal_MC(X,Y); 
Y_AL_a = Funcion_1(x_eval,A_1,B_1);
subplot(2,3,1)
plot(X,Y,'*r',x_eval,Y_AL_a, 'b') % Nube de puntos originales
xlabel('D [mm]')
ylabel('P(D) [mmHg]')
legend ('Pk', 'Aproximaci贸n lineal Pk')
title('Ejercicio 3a');
grid

%% Aproximacion Exponencial; b)

Y_log = log(Y);
[A_2,B_2,CC] = Ajuste_Lineal_MC(X,Y_log);
p0 = exp(B_2);
Y_AL_b =  p0 * exp(A_2 * x_eval); % Recta de Ajuste evaluada en los puntos xk
subplot(2,3,2)
plot(X,Y,'*r',x_eval,Y_AL_b, 'b') % Nube de puntos originales
xlabel('D [mm]')
ylabel('P(D) [mmHg]')
legend ('Pk', 'Aproximacion exponencial Pk')
title('Ejercicio 3b');
grid

%% Lagrange; c)

C_3 = Interp_Lagrange(X,Y');
Y_AL_c = @(X) Eval_Polinomio(X,C_3);% Armado del polinomio
subplot(2,3,3)
plot(X,Y,'*r',x_eval, Y_AL_c(x_eval),'b')
xlabel('D [mm]')
ylabel('P(D) [mmHg]')
legend ('Pk', 'Polinomio Interpolador de Lagrange Pk')
title('Ejercicio 3c');
grid 

%% Interpolador Newton; d)

C_4 = Interp_Newton(X,Y);
Y_AL_d = @(X)Eval_Polinomio(X,C_4);% Armado del polinomio interpolador de Newton
subplot(2,3,4)
plot(X,Y,'*r',x_eval, Y_AL_d(x_eval),'b')
xlabel('D [mm]')
ylabel('P(D) [mmHg]')
legend ('Pk', 'Polinomio Interpolador de Newton Pk')
title('Ejercicio 3d');
grid 

%% Spline Cubicas; e)
MC = Spline_Cubica(X,Y);
N = size(MC,1);
L = 14; % Cantidad de puntos entre intervalos de X discretos
x = zeros(N,L+1); % Matriz de abscisas (101 puntos por cada intervalo)
for k=1:N
    d = (X(k+1)-X(k))/L; % Armado de abscisas temporales para los 
    x(k,:) = X(k):d:X(k+1); % Polinomios spline c煤bicos
end    
subplot(2,3,5)
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
title('Ejercicio 3e');
hold off

%% Interpolacion de Richardson; f)

%i.
D_m = 12.5;
Delta = 1e-13;
Tol = 1e-13;
F2 = figure(2);
set(F2,'position',[140 20 1300 750],'Menubar','none',...
        'NumberTitle','off','name', 'Ejercicio 3f TP#FINAL. M閠odos Num閞icos');
subplot(2,3,1)
[D,err,relerr,n] = Extrapolacion_Richardson_O2n(@(X) A_1*X+B_1,D_m,Delta,Tol);% referencia a esa funci髇
Elasticidad1 = D(end,end);
tangente = Elasticidad1*(X - D_m) + A_1*D_m+B_1;
plot(X,Y,'*r',x_eval,Y_AL_a, 'k',X, tangente, 'y--') % Nube de puntos originales
xlabel('D [mm]')
ylabel('P(D) [mmHg]');
legend('Pk','Aproximacion Lineal' ,'Tangente en D_M');
title('Ejercicio 3f (i)');
grid 

 %ii.
subplot(2,3,2)
[D,err,relerr,n] = Extrapolacion_Richardson_O2n(@(X)  p0 * exp(A_2 * X),D_m,Delta,Tol);% referencia a esa funci贸n
Elasticidad_2 = D(end,end);
tangente = Elasticidad_2*(X - D_m) +  p0 * exp(A_2 * D_m);
plot(X,Y,'*r',x_eval,Y_AL_b, 'k', X, tangente, 'm--') % Nube de puntos originales
xlabel('D [mm]')
ylabel('P(D) [mmHg]');
legend('Pk','Aproximacion Exponencial' ,'Tangente en D_M');
title('Ejercicio 3f (ii)');
grid 

%iii.
subplot(2,3,3)
[D,err,relerr,n] = Extrapolacion_Richardson_O2n(@(x) Y_AL_c(x),D_m,Delta,Tol);% referencia a esa funci贸n
Elasticidad_3 = D(end,end);
tangente = Elasticidad_3*(X - D_m) + Y_AL_c(D_m);
plot(X,Y,'*r',x_eval,Y_AL_c(x_eval), 'k',X, tangente, 'm--') % Nube de puntos originales
xlabel('D [mm]')
ylabel('P(D) [mmHg]');
legend('Pk','Aproximacion Polinomial' ,'Tangente en D_M');
title('Ejercicio 3f (iii)');
grid 

%iv
subplot(2,3,4)
[D,err,relerr,n] = Extrapolacion_Richardson_O2n(@(x) Y_AL_d(x),D_m,Delta,Tol);% referencia a esa funci贸n
Elasticidad_4 = D(end,end);
tangente = Elasticidad_4*(X - D_m) + Y_AL_d(D_m);
plot(X,Y,'*r',x_eval,Y_AL_d(x_eval), 'k',X, tangente, 'm--') % Nube de puntos originales
xlabel('D [mm]')
ylabel('P(D) [mmHg]');
legend('Pk','Aproximacion Polinomial' ,'Tangente en D_M');
title('Ejercicio 3f (iv)');
grid

%v
subplot(2,3,5)
Spline_9 = @(x) Eval_Polinomio2(x, MC(9,:), X(9)); %funcion
[D, err, relerr, n] = Extrapolacion_Richardson_O2n(Spline_9, D_m, Delta, Tol);
Elasticidad_5 = D(end,end);
tangente = Elasticidad_5 * (X - D_m) + Spline_9(D_m);

y_spline = zeros(size(x_eval));
for i = 1:length(x_eval)%evaluo la spline a mano punto por punto
    idx = find(X <= x_eval(i), 1, 'last');%busco tramo para cada x_eval(i)
    if idx == length(X), idx = idx - 1; end
    y_spline(i) = Eval_Polinomio2(x_eval(i), MC(idx,:), X(idx));
end
plot(X, Y, '*r', x_eval, y_spline, 'k', X, tangente, 'm--')
xlabel('D [mm]')
ylabel('P(D) [mmHg]')
legend('Pk', 'Spline c鷅ica', 'Tangente en D_M')
title('Ejercicio 3f (v)')
grid

%% g:
P_v = Funcion_3g(x_eval);
F3 = figure(3);
set(F3,'position',[100 100 1300 750],'Menubar','none',...
    'NumberTitle','off','name','Ejercicio 3g TP#FINAL. M閠odos Num閞icos');
%i:
subplot(2,3,1)
plot(X,Y,'*r',x_eval,Y_AL_a, 'b',x_eval,P_v,'m') % Nube de puntos originales
xlabel('D [mm]')
ylabel('P(D) [mmHg]')
legend ('Pk', 'Aproximacion lineal Pk','Funcion Verdadera')
title('Ejercicio 3g (i)');
grid

Dif=0;
for i=1:N
    Dif = Dif + (feval('Funcion_3g', X(k)) - A_1*X(i)+B_1)^2;
end
ECM_a=1/N *sqrt(Dif);
fprintf('El ECM de la aproximaci髇 lineal es: %.10f\n', ECM_a);

%ii:
subplot(2,3,2)
plot(X,Y,'*r',x_eval,Y_AL_b,'b',x_eval,P_v,'m') % Nube de puntos originales
xlabel('D [mm]')
ylabel('P(D) [mmHg]')
legend ('Pk', 'Aproximacion exponencial Pk','Funcion Verdadera')
title('Ejercicio 3g (ii)');
grid

Difb=0;
for i=1:N
    Difb = Difb + (feval('Funcion_3g', X(k)) - exp(B_2)*exp(X(i)*A_2))^2;
end
ECM_b=1/N *sqrt(Difb);
fprintf('El ECM de la aproximaci髇 exponencial es: %.10f\n', ECM_b);


%iii:
subplot(2,3,3)
plot(X,Y,'*r',x_eval,Y_AL_c(x_eval), 'b',x_eval,P_v,'m') % Nube de puntos originales
xlabel('D [mm]')
ylabel('P(D) [mmHg]')
legend ('Pk', 'Polinomio Interpolador de Lagrange Pk','Funcion Verdadera')
title('Ejercicio 3g (iii)');
grid

Difc=0;
for i=1:N
    Difc = Difc + (feval('Funcion_3g', X(k)) - Y_AL_c(X(i)))^2;
end
ECM_c=(1/N) *sqrt(Difc);
fprintf('El ECM del Polinomio Interpolador de Lagrange es: %.10f\n', ECM_c);

%iv:
subplot(2,3,4)
plot(X,Y,'*r',x_eval,Y_AL_d(x_eval), 'b',x_eval,P_v,'m') % Nube de puntos originales
xlabel('D [mm]')
ylabel('P(D) [mmHg]')
legend ('Pk', 'Polinomio Interpolador de Newton Pk','Funcion Verdadera')
title('Ejercicio 3g (iv)');
grid

Difc=0;
for i=1:N
    Difc = Difc + (feval('Funcion_3g', X(k)) - Y_AL_d(X(i)))^2;
end
ECM_d=(1/N) *sqrt(Difc);
fprintf('El ECM del Polinomio Interpolador de Newton es: %.10f\n', ECM_d);

%v
subplot(2,3,5)
plot(X,Y,'*r') % Nube de puntos originales
hold on;

Difc=0;
for k=1:N
    Y_AL_e(k,:) = Eval_Polinomio2(x(k,:),MC(k,:),X(k));
    plot(x(k,:),Y_AL_e(k,:),'b');
    y_interp = Eval_Polinomio2(X(k), MC(k,:), X(k)); 
    Difc = Difc + (feval('Funcion_3g', X(k)) - y_interp).^2;
end

hold on;
plot(x_eval,P_v,'m');
hold off;
xlabel('D [mm]')
ylabel('P(D) [mmHg]')
legend ('Pk', 'Evolucion Spline Cubicas','Funcion Verdadera')
title('Ejercicio 3g (v)');
grid

ECM_e=(1/N) *sqrt(Difc);
fprintf('El ECM de la interpolacion por Spline Cubicas es: %.10f\n', ECM_e);

ECMs = [ECM_a, ECM_b, ECM_c, ECM_d, ECM_e];
[minECM, idxMin] = min(ECMs);%digo cual es el minimo
nombresECM = {'Aproximaci髇 Lineal', 'Aproximaci髇 Exponencial','Polinomio Interpolador de Lagrange', ...
              'Polinomio Interpolador de Newton','Interpolaci髇 por Spline C鷅icas'};
fprintf('El ECM m韓imo es el del %s\n', nombresECM{idxMin});