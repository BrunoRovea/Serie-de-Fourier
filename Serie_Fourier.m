clear; clc; close all;
tic;
%=========================Variï¿½veis de entrada============================%
T0 = 0.1; %Perï¿½odo fundamental do sinal
N = 10000; %Nï¿½mero de termos da sï¿½rie
Ne = 10; %Numero de harmï¿½nicos no espectro em frequï¿½ncia



%=========================Caracterï¿½stica do sinal=========================%
f0 = 1/T0; %Frequï¿½ncia fudamental em Hz
W0 = 2*pi*f0; %Frequï¿½ncia angular fundamental em Rad/s



%===================Cï¿½lculo dos coeficientes da Sï¿½rie=====================%
%Definindo o sinal no tempo
syms t;syms n;
x1 = (1);
I1 = [0, T0/2];
x2 = (0);
I2 = [T0/2 T0];

%Resolvo a integral simbï¿½lica definida
aa0  = (1/T0)*int(x1,t, I1)*t;
aan = (2/T0)*int(x1*cos(n*W0*t),t, I1);
bbn = (2/T0)*int(x1*sin(n*W0*t),t, I1);
ddn = 1/T0*int(x1*exp(-1i*n*W0*t),t, I1);

aa0  = aa0 + (1/T0)*int(x2,t, I2);
aan = aan + (2/T0)*int(x2*cos(n*W0*t),t, I2);
bbn = bbn + (2/T0)*int(x2*sin(n*W0*t),t, I2);
ddn = 1/T0*int(x1*exp(-1i*n*W0*t),t, I2);


%Converte-se as funï¿½ï¿½es de variï¿½veis simbï¿½licas para funï¿½ï¿½es tipo handle
%Funï¿½ï¿½es que dependem de n
a0 = matlabFunction(aa0);
a0 = a0(1);
an = matlabFunction(aan);
bn = matlabFunction(bbn);
dn = matlabFunction(ddn);

%============Cï¿½lculo dos coeficientes e determinaï¿½ï¿½o da sï¿½rie=============%
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
   % Corrigindo o possï¿½vel erro da divisï¿½o 0/0 no cï¿½lculo de phi
   if((a(n) == 0) && (b(n) == 0))
      phi(n) = 0; 
   end
   
   %Determinaï¿½ï¿½o da funï¿½ï¿½o
   x = @(t) x(t) + a(n)*cos(n*W0*t) + b(n)*sin(n*W0*t);%Forma trigonomï¿½trica
   %x = @(t) x(t) + c(n)*cos(n*W0*t - phi(n));% Forma compacta
   %x = @(t) x(t) + d(n)*exp(1i*n*W0*t);% Forma complexa
   %x = @(t) x(t) + conj(d(n))*exp(-1i*n*W0*t);% n assume valor de -n nessa linha
                                               % por isso do sinal negativo
end



%================================Grï¿½ficos=================================%
n = 1:1:Ne; %Abscissa do espectro
t = linspace(-2*T0, 2*T0, N); %Vetor tempo
% Consertar Eixos

%Sinal no Tempo
figure(1);
plot(t, x(t));
title('Sinal x(t) no domï¿½nio do tempo');
xlabel('Tempo[s]'); ylabel('x(t)[v]');
% legend('1000 pontos');
grid on;


%Espectro em frequï¿½ncia

%an
figure(2);
a(N+1) = a0;
plot(0, a0, 'Xr');
hold on;
plot(n, an(n),'Xr');
title('Espectro em frequï¿½ncia do sinal x(t)');
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
title('Espectro em frequï¿½ncia do sinal x(t)');
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
title('Espectro em frequï¿½ncia do sinal x(t)');
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
title('Espectro em frequï¿½ncia do sinal x(t)');
xlabel('nW0'); ylabel('Mï¿½dulo dos coefcientes dn');
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
