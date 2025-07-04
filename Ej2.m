%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%% Ejercicio Nro. 2 TP#FINAL Metodos Numericos
%%%
%%% Juana Kallis, Emma fiorini y Agustina Vidaurreta
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% c
M = 100;
m = 5;
Ang = pi/6;
g = 9.81;
b_M = 5;
b_m= 10;
t0 = 0;
tf = 4;
dt = 0.001;
N = (tf-t0)/dt;
den = M + m * sin(Ang)^2;

A = [0 1 0 0; 0 -b_M/den 0 b_m*cos(Ang)/den; 0 0 0 1; 0 b_m*cos(Ang)/den 0 -(M + m)*b_M / (m * den)];
B = [0 -g*m*cos(Ang)*sin(Ang)/(M + m*sin(Ang)*sin(Ang)) 0 (m+M)*g*cos(Ang)*sin(Ang)/(M + m*sin(Ang)*sin(Ang))]';

Mod_Est = @(t, X) (A * X' + B)';  
X0 = [0; 0; 0; 0];

[T, X] = Ec_Dif_Runge_Kutta_O4_Sistemas(Mod_Est, t0, tf,X0, N);

x = X(:,1);
v_x = X(:,2);
q = X(:,3);
v_q = X(:,4);

F1 = figure(1);
set(F1, 'Position', [100, 100, 1200, 700],'Menubar','none',...
   'NumberTitle','off','name', 'Ejercicio 2c TP#FINAL- Métodos Numéricos');
subplot(2,3,1)
plot(T, v_q, 'r')
xlabel(' t [s]');
ylabel('v_q[m/s]');
title('Velocidad de q vs tiempo');
grid;
subplot (2,3,2)
plot(T, v_x, 'r')
xlabel(' t [s]');
ylabel('v_x[m/s]');
title('Velocidad en x vs tiempo');
grid;
subplot (2,3,5)
plot(T, x, 'r')
xlabel(' t [s]');
ylabel('x[m]');
title('Posicion en x vs tiempo');
grid;
subplot (2,3,4)
plot(T ,q,'r')
xlabel(' t [s]');
ylabel('q[m]');
title('Posicion en q vs tiempo');
grid;

%% d
a_x = zeros(size(v_x));
a_q = zeros(size(v_q));

for i = 2:N 
    a_x(i) = (v_x(i+1) - v_x(i-1)) / (2*dt);
    a_q(i) = (v_q(i+1) - v_q(i-1)) / (2*dt);
end

F2 = figure(2);
set(F2, 'Position', [100, 100, 1200, 500],'Menubar','none',...
   'NumberTitle','off','name', 'Ejercicio 2d TP#FINAL - Aceleraciones');

subplot(1,2,1)
plot(T, a_x, 'b')
xlabel('t [s]');
ylabel('a_M [m/s^2]');
title('Aceleración de la cuña M');
grid;
subplot(1,2,2)
plot(T, a_q, 'm')
xlabel('t [s]');
ylabel('a_m [m/s^2]');
title('Aceleración del bloque m');
grid;