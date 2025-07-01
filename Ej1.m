%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%% Ejercicio Nro. 1 TP#FINAL Metodos Numericos
%%%
%%% Juana Kallis, Emma fiorini y Agustina Vidaurreta
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;
clc;

M = 100;
m = 5;
Ang = pi/6;
g = 9.81;
cte_1 = ((m+M)* g * sin(Ang))/(M+(m*sin(Ang)*sin(Ang)));
cte_2 = -(m*g*sin(Ang)*cos(Ang))/(M+(m*sin(Ang)*sin(Ang)));

%% c. 
A = [cos(Ang) 1; (m+M) m*cos(Ang)];
B = [g*sin(Ang) 0]';

X = Triang_Gauss(A,B);

fprintf('La segunda derivada de x, calculada por Triangulacion gaussiana es :%.5f\n',X(1));
fprintf('La segunda derivada de q, calculada por Triangulacion gaussiana es :%.5f\n',X(2)); 
fprintf('La segunda derivada de x, calculada analiticamente es :%.5f\n',cte_2);
fprintf('La segunda derivada de q, calculada analiticamente es :%.5f\n',cte_1); 

%% d.
t0 = 0;
tf = 4;
dt = 0.001;
M = (tf-t0)/ dt; % cantidad de puntos
t_eval = linspace(t0, tf, M+1);% incluye t=0

a_x = @(t) cte_2;
a_q = @(t) cte_1; 

v_q = zeros(size(t_eval));% inicializo las velocidades
v_x = zeros(size(t_eval));

for i = 2:2:length(t_eval) % integral para cada puntos, creo una nube de puntos
    ti = t_eval(i);
    Mi = i / 2;  % cantidad de subintervalos
    v_q(i) = Regla_Simpson_Compuesta(a_q, t0, ti, Mi);
    v_x(i) = Regla_Simpson_Compuesta(a_x, t0, ti, Mi);
end

[Aq, Bq, Cq] = Ajuste_Lineal_MC(t_eval', v_q');% hago ajuste lineal
[Ax, Bx, Cx] = Ajuste_Lineal_MC(t_eval', v_x');

F1 = figure(1);%grafico;
set(F1, 'Position', [100, 100, 1200, 700],'Menubar','none',...
   'NumberTitle','off','name', 'Ejercicio 1 TP#FINAL- Métodos Numéricos');

y_vq= @(t) Aq*t + Bq;% funci�n de velocidad ajustada
y_vx= @(t) Ax*t + Bx;

subplot(2,2,1);
hold on;
plot(t_eval, y_vq(t_eval), 'b');
xlabel('t [s]');
ylabel('v_q [m/s]');
title('Velocidad en q vs tiempo');
grid on;
subplot(2,2,2);
hold on;
plot(t_eval, y_vx(t_eval), 'r');
xlabel('t [s]');
ylabel('v_x [m/s]');
title('Velocidad en x vs tiempo');
grid on;
%% e.

x_M = zeros(size(t_eval));
x_m = zeros(size(t_eval));

for i = 2:2:length(t_eval) %integral para cada puntos, creo una nube de puntos
    ti = t_eval(i);
    Mi = i / 2;  % cantidad de subintervalos
    x_M(i) = Regla_Simpson_Compuesta(y_vq, t0, ti, Mi);
    x_m(i) = Regla_Simpson_Compuesta(y_vx, t0, ti, Mi);
end

Cm = polyfit(t_eval, x_m, 2);% llevo a un ajuste cuadratico
y_xm = polyval(Cm, t_eval); % evaluo el polinomio

CM = polyfit(t_eval, x_M, 2);% llevo a un ajuste cuadratico
y_xM = polyval(CM, t_eval); % evaluo el polinomio

subplot(2,2,3);%grafico;
hold on;
plot(t_eval, y_xM, 'b');
xlabel('t [s]');
ylabel('x_M [m]');
title('Posicion q vs tiempo');
grid on;
subplot(2,2,4);
hold on;
plot(t_eval, y_xm, 'r');
xlabel('t [s]');
ylabel('x_m [m]');
title('Posicion x vs tiempo');
grid on;
