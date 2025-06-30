%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%% Ejercicio Nro. 1 TP#FINAL MÃ©todos NumÃ©ricos
%%%
%%% Juana Kallis, Emma fiorini y Agustina Vidaurreta
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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


for i = 2:2:length(t_eval) %integral para cada puntos, creo una nube de puntos
    ti = t_eval(i);
    Mi = i / 2;  % cantidad de subintervalos
    v_q(i) = Regla_Simpson_Compuesta(a_q, t0, ti, Mi);
    v_x(i) = Regla_Simpson_Compuesta(a_x, t0, ti, Mi);
end

%hago ajuste lineal
[Aq, Bq, Cq] = Ajuste_Lineal_MC(t_eval', v_q');
[Ax, Bx, Cx] = Ajuste_Lineal_MC(t_eval', v_x');

F1 = figure(1);
set(F1, 'Position', [100, 100, 1200, 700],'Menubar','none',...
   'NumberTitle','off','name', 'Ejercicio 1 TP#FINAL- Métodos Numéricos');

%grafico;
subplot(1,2,1);
hold on;
y_vq = Aq * t_eval + Bq;
plot(t_eval, y_vq, 'b');
xlabel('t [s]');
ylabel('v_q [m/s]');
title('Velocidad en q vs tiempo');
grid on;

subplot(1,2,2);
hold on;
y_vx=Ax*t_eval + Bx;
plot(t_eval, y_vx, 'r');
xlabel('t [s]');
ylabel('v_x [m/s]');
title('Velocidad en x vs tiempo');
grid on;