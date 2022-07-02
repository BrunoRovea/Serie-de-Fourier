clear; clc; close all;
tic;
%=========================Vari�veis de entrada============================%
T0 = 0.1; %Per�odo fundamental do sinal
N = 1000; %N�mero de termos da s�rie
Ne = 10; %Numero de harm�nicos no espectro em frequ�ncia



%=========================Caracter�stica do sinal=========================%
f0 = 1/T0; %Frequ�ncia fudamental em Hz
W0 = 2*pi*f0; %Frequ�ncia angular fundamental em Rad/s



%===================C�lculo dos coeficientes da S�rie=====================%
syms t;syms n;
x1 = (58/3 + 80*t);
I1 = [-1/60 2/60];
x2 = (74/3 - 80*t);
I2 = [2/60 5/60];

%Resolvo a integral simb�lica definida
aa0  = (1/T0)*int(x1,t, I1)*t;
aan = (2/T0)*int(x1*cos(n*W0*t),t, I1);
bbn = (2/T0)*int(x1*sin(n*W0*t),t, I1);
ddn = 1/T0*int(x1*exp(-1i*n*W0*t),t, I1);

aa0  = aa0 + (1/T0)*int(x2,t, I2);
aan = aan + (2/T0)*int(x2*cos(n*W0*t),t, I2);
bbn = bbn + (2/T0)*int(x2*sin(n*W0*t),t, I2);
ddn = 1/T0*int(x1*exp(-1i*n*W0*t),t, I2);


%Converte-se as fun��es de vari�veis simb�licas para fun��es tipo handle
%Fun��es que dependem de n
a0 = matlabFunction(aa0);
a0 = a0(1);
an = matlabFunction(aan);
bn = matlabFunction(bbn);
dn = matlabFunction(ddn);

%============C�lculo dos coeficientes e determina��o da s�rie=============%
x = @(t) a0;% Ou d0
a = zeros(1, N);
b = zeros(1, N);
c = zeros(1, N);
d = zeros(1, N);
phi = zeros(1, N);
for n = 1: N
   a(n) = an(n);
   b(n) = bn(n);
   c(n) = sqrt(a(n)^2 + b(n)^2);
   d(n) = dn(n);
   
   % Corrigindo um erro de arredondamento do matlab
   if((abs(a(n)) < 1e-14))
       a(n) = 0;
   end
   if((abs(b(n)) < 1e-14))
       b(n) = 0;
   end
   if((abs(c(n)) < 1e-14))
       c(n) = 0;
   end
   if((abs(d(n)) < 1e-14))
       d(n) = 0 + 0j;
   end
   
   phi(n) = atan(b(n)/a(n));
   % Corrigindo o poss�vel erro da divis�o 0/0 no c�lculo de phi
   if((a(n) == 0) && (b(n) == 0))
      phi(n) = 0; 
   end
   
   %Determina��o da fun��o
   x = @(t) x(t) + a(n)*cos(n*W0*t) + b(n)*sin(n*W0*t);%Forma trigonom�trica
   %x = @(t) x(t) + c(n)*cos(n*W0*t - phi(n));% Forma compacta
   %x = @(t) x(t) + d(n)*exp(1i*n*W0*t);% Forma complexa
   %x = @(t) x(t) + conj(d(n))*exp(-1i*n*W0*t);% n assume valor de -n nessa linha
                                               % por isso do sinal negativo
end



%================================Gr�ficos=================================%
n = 1:1:Ne; %Abscissa do espectro
t = linspace(-2*T0, 2*T0, N); %Vetor tempo
% Consertar Eixos

%Sinal no Tempo
figure(1);
plot(t, x(t));
title('Sinal x(t) no dom�nio do tempo');
xlabel('Tempo[s]'); ylabel('x(t)[v]');
% legend('1000 pontos');
grid on;


%Espectro em frequ�ncia

%an
figure(2);
a(N+1) = a0;
plot(0, a0, 'Xr');
hold on;
plot(n, an(n),'Xr');
title('Espectro em frequ�ncia do sinal x(t)');
xlabel('nW0'); ylabel('Coeficientes an');
% legend('1000 pontos');
if(min(a) == max(a))
   axis([0 Ne+0.8 -1 1]);
end
if(min(a) ~= max(a))
   axis([0 Ne+0.8 min(a)*1.5 max(a)*1.3]);
end
grid on;

%bn
figure(3);
plot(n, bn(n),'Xr');
title('Espectro em frequ�ncia do sinal x(t)');
xlabel('nW0'); ylabel('Coeficientes bn');
% legend('1000 pontos');
if(min(b) == max(b))
   axis([0 Ne+0.8 -1 1]);
end
if(min(b) ~= max(b))
    axis([0 Ne+0.8 min(b)*1.5 max(b)*1.3]);
end
grid on;

%cn
figure(4);
c(N+1) = a0;
plot(0, a0, 'Xr');
hold on;
plot(n, c(n),'Xr');
title('Espectro em frequ�ncia do sinal x(t)');
xlabel('nW0'); ylabel('Coeficientes cn');
% legend('1000 pontos');
if(min(c) == max(c))
   axis([0 Ne+0.8 -1 1]);
end
if(min(c) ~= max(c))
    axis([0 Ne+0.8 min(c)*1.5 max(c)*1.3]);
end
grid on;

%|dn|
figure(5);
plot(0, a0, 'Xr');
hold on;
plot(n, abs(d(n)),'Xr');
title('Espectro em frequ�ncia do sinal x(t)');
xlabel('nW0'); ylabel('M�dulo dos coefcientes dn');
% legend('1000 pontos');
if(min(abs(d(n))) == max(abs(d(n))))
   axis([0 Ne+0.8 -1 1]);
end
%if(min(abs(d(n))) ~= max(abs(d(n))))
%    axis([0 Ne+0.8 min(abs(d(n)))*1.5 a0*1.3]);
%end
grid on;

%phi
figure(6);
d(N+1) = a0;
plot(n, phi(n),'Xr');
title('Fase do sinal x(t)');
xlabel('nW0'); ylabel('phi');
% legend('1000 pontos');
if(min(phi) == max(phi))
   axis([0 Ne+0.8 -1 1]);
end
if(min(phi) ~= max(phi))
    axis([0 Ne+0.8 min(phi)*1.5 max(phi)*1.3]);
end
grid on;
%clear aa0 aan an ans bbn bn ccn cn ddn dn I1 I2 x1 x2 Ne N;
t = toc;