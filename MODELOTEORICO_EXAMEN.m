%% Modelo teorico de 3GDL --- Original
close all, clc, clear all;
mnom=100; %Masa nominal de entreplanta [ton]
mi=mnom*1000; %Masa [kg]
Ixx=5.6969E-05; %Momento de inercia [m4]
E=2.1e+11; %Modulo de elasticidad del acero [N/m2]
h=3.00; %Altura [m]
Kxx=24*E*Ixx/h^3; %Rigidez del muelle [N/m]
M=[mi 0 0;0 mi 0;0 0 20*mi]; %Matriz de masas
K=[Kxx -Kxx 0;-Kxx 2*Kxx -Kxx;0 -Kxx 2*Kxx]; %Matriz de rigidez
z=[0.02 0.02 0.03]; %Vector de amortiguamiento
nn=10; %Multiplo para resolucion

[fn phi V mk kk ck dw H]=Matransferencia(K,M,z,nn);

%Ploteo
npt=2^nn;
fmax=1.4*max(fn); df=fmax/npt;
f=0:df:fmax-df; %Vector de frecuencias [Hz]
H11=zeros(1,npt); H12=zeros(1,npt); H13=zeros(1,npt);
H21=zeros(1,npt); H22=zeros(1,npt); H23=zeros(1,npt);
H31=zeros(1,npt); H32=zeros(1,npt); H33=zeros(1,npt);
for mm=1:npt
    H11(1,mm)=H(1,1,mm);
    H12(1,mm)=H(1,2,mm);
    H13(1,mm)=H(1,3,mm);
    H21(1,mm)=H(2,1,mm);
    H22(1,mm)=H(2,2,mm);
    H23(1,mm)=H(2,3,mm);
    H31(1,mm)=H(3,1,mm);
    H32(1,mm)=H(3,2,mm);
    H33(1,mm)=H(3,3,mm);
end

% figure()
% loglog(f,abs(H11),'displayname','H11'); hold on;
% loglog(f,abs(H21),'displayname','H21');
% loglog(f,abs(H31),'displayname','H31'); legend('show');
% xlabel('Frecuencia [Hz]');ylabel('|H(w)|'); %xlim([0 100])
% hold off;
% 
% figure()
% loglog(f,abs(H12),'displayname','H12'); hold on;
% loglog(f,abs(H22),'displayname','H22');
% loglog(f,abs(H32),'displayname','H32'); legend('show');
% xlabel('Frecuencia [Hz]');ylabel('|H(w)|'); %xlim([0 100])
% hold off;
% 
% figure()
% loglog(f,abs(H13),'displayname','H13'); hold on;
% loglog(f,abs(H23),'displayname','H23');
% loglog(f,abs(H33),'displayname','H33'); legend('show');
% xlabel('Frecuencia [Hz]');ylabel('|H(w)|'); %xlim([0 100])
% hold off;

tablafn=table({'Modo 1';'Modo 2';'Modo 3'},fn,'variablename',{'Modos ', ' Frecuencias [Hz]'});
display(tablafn);

tablaphi=table({'GDL 1';'GDL 2';'GDL 3'},phi,'variablename',{'GDL', ['MODO1 ' ' MODO2 ' ' MODO 3']});
display(tablaphi);

tabla_dinak=table({'Modo 1';'Modo 2';'Modo 3'},[mk' kk' ck'],'variablename',{'Prop Dinamica', ['Mk ' ' Kk ' ' Ck']});
display(tabla_dinak);

%% Parte 1 --- Caracterizacion del modelo --- 2GDL (3er GDL empotrado)
%Modelo teorico de 2GDL
close all, clc, clear all;
mnom=100; %Masa nominal de entreplanta [ton]
mi=mnom*1000; %Masa [kg]
Ixx=5.6969E-05; %Momento de inercia [m4]
E=2.1e+11; %Modulo de elasticidad del acero [N/m2]
h=3.00; %Altura [m]
Kxx=24*E*Ixx/h^3; %Rigidez del muelle [N/m]
M=[mi 0;0 mi]; %Matriz de masas
K=[Kxx -Kxx;-Kxx 2*Kxx]; %Matriz de rigidez
z=[0.02 0.02]; %Vector de amortiguamiento
nn=10; %Multiplo para resolucion

[fn phi V mk kk ck]=Matransferencia(K,M,z,nn);

tablafn=table({'Modo 1';'Modo 2'},fn,'variablename',{'Modos ', ' Frecuencias [Hz]'});
display(tablafn);

tablaphi=table({'GDL 1';'GDL 2'},phi,'variablename',{'GDL', ['MODO1 ' ' MODO2 ']});
display(tablaphi);

% tabla_dinak=table({'Modo 1';'Modo 2'},[mk' kk' ck'],'variablename',{'Prop Dinamica', ['Mk ' ' Kk ' ' Ck']});
% display(tabla_dinak);

%% Modelo teorico de 3GDL----(Empotramiento---> Aumentando la masa) ---- valido
close all, clc, clear all;
mnom=100; %Masa nominal de entreplanta [ton]
mi=mnom*1000; %Masa [kg]
Ixx=5.6969E-05; %Momento de inercia [m4]
E=2.1e+11; %Modulo de elasticidad del acero [N/m2]
h=3.00; %Altura [m]
Kxx=24*E*Ixx/h^3; %Rigidez del muelle [N/m]
M=[mi 0 0;0 mi 0;0 0 100*mi]; %Matriz de masas
K=[Kxx -Kxx 0;-Kxx 2*Kxx -Kxx;0 -Kxx 2*Kxx]; %Matriz de rigidez
z=[0.02 0.02 0.03]; %Vector de amortiguamiento
nn=10; %Multiplo para resolucion

[fn phi V mk kk ck dw H]=Matransferencia(K,M,z,nn);

%Ploteo
npt=2^nn;
fmax=1.4*max(fn); df=fmax/npt;
f=0:df:fmax-df; %Vector de frecuencias [Hz]
H11=zeros(1,npt); H12=zeros(1,npt); H13=zeros(1,npt);
H21=zeros(1,npt); H22=zeros(1,npt); H23=zeros(1,npt);
H31=zeros(1,npt); H32=zeros(1,npt); H33=zeros(1,npt);
for mm=1:npt
    H11(1,mm)=H(1,1,mm);
    H12(1,mm)=H(1,2,mm);
    H13(1,mm)=H(1,3,mm);
    H21(1,mm)=H(2,1,mm);
    H22(1,mm)=H(2,2,mm);
    H23(1,mm)=H(2,3,mm);
    H31(1,mm)=H(3,1,mm);
    H32(1,mm)=H(3,2,mm);
    H33(1,mm)=H(3,3,mm);
end

figure()
loglog(f,abs(H11),'displayname','H11'); hold on;
loglog(f,abs(H21),'displayname','H21');
loglog(f,abs(H31),'displayname','H31'); legend('show');
xlabel('Frecuencia [Hz]');ylabel('|H(w)|'); %xlim([0 100])
hold off;

figure()
loglog(f,abs(H12),'displayname','H12'); hold on;
loglog(f,abs(H22),'displayname','H22');
loglog(f,abs(H32),'displayname','H32'); legend('show');
xlabel('Frecuencia [Hz]');ylabel('|H(w)|'); %xlim([0 100])
hold off;

figure()
loglog(f,abs(H13),'displayname','H13'); hold on;
loglog(f,abs(H23),'displayname','H23');
loglog(f,abs(H33),'displayname','H33'); legend('show');
xlabel('Frecuencia [Hz]');ylabel('|H(w)|'); %xlim([0 100])
hold off;

tablafn=table({'Modo 1';'Modo 2';'Modo 3'},fn,'variablename',{'Modos ', ' Frecuencias [Hz]'});
display(tablafn);

tablaphi=table({'GDL 1';'GDL 2';'GDL 3'},phi,'variablename',{'GDL', ['MODO1 ' ' MODO2 ' ' MODO 3']});
display(tablaphi);

tabla_dinak=table({'Modo 1';'Modo 2';'Modo 3'},[mk' kk' ck'],'variablename',{'Prop Dinamica', ['Mk ' ' Kk ' ' Ck']});
display(tabla_dinak);

%% Modelo teorico de 3GDL----(Empotramiento---> Disminuyendo la masa)
close all, clc, clear all;
mnom=100; %Masa nominal de entreplanta [ton]
mi=mnom*1000; %Masa [kg]
Ixx=5.6969E-05; %Momento de inercia [m4]
E=2.1e+11; %Modulo de elasticidad del acero [N/m2]
h=3.00; %Altura [m]
Kxx=24*E*Ixx/h^3; %Rigidez del muelle [N/m]
M=[mi 0 0;0 mi 0;0 0 0.001*mi]; %Matriz de masas
K=[Kxx -Kxx 0;-Kxx 2*Kxx -Kxx;0 -Kxx 2*Kxx]; %Matriz de rigidez
z=[0.02 0.02 0.03]; %Vector de amortiguamiento
nn=12; %Multiplo para resolucion

[fn phi V mk kk ck dw H]=Matransferencia(K,M,z,nn);

%Ploteo
npt=2^nn;
fmax=1.4*max(fn); df=fmax/npt;
f=0:df:fmax-df; %Vector de frecuencias [Hz]
H11=zeros(1,npt); H12=zeros(1,npt); H13=zeros(1,npt);
H21=zeros(1,npt); H22=zeros(1,npt); H23=zeros(1,npt);
H31=zeros(1,npt); H32=zeros(1,npt); H33=zeros(1,npt);
for mm=1:npt
    H11(1,mm)=H(1,1,mm);
    H12(1,mm)=H(1,2,mm);
    H13(1,mm)=H(1,3,mm);
    H21(1,mm)=H(2,1,mm);
    H22(1,mm)=H(2,2,mm);
    H23(1,mm)=H(2,3,mm);
    H31(1,mm)=H(3,1,mm);
    H32(1,mm)=H(3,2,mm);
    H33(1,mm)=H(3,3,mm);
end
figure()
loglog(f,abs(H11),'displayname','H11'); hold on;
loglog(f,abs(H21),'displayname','H21');
loglog(f,abs(H31),'displayname','H31'); legend('show');
xlabel('Frecuencia [Hz]');ylabel('|H(w)|'); %xlim([0 100])
hold off;

figure()
loglog(f,abs(H12),'displayname','H12'); hold on;
loglog(f,abs(H22),'displayname','H22');
loglog(f,abs(H32),'displayname','H32'); legend('show');
xlabel('Frecuencia [Hz]');ylabel('|H(w)|'); %xlim([0 100])
hold off;

figure()
loglog(f,abs(H13),'displayname','H13'); hold on;
loglog(f,abs(H23),'displayname','H23');
loglog(f,abs(H33),'displayname','H33'); legend('show');
xlabel('Frecuencia [Hz]');ylabel('|H(w)|'); %xlim([0 100])
hold off;

% tablafn=table({'Modo 1';'Modo 2';'Modo 3'},fn,'variablename',{'Modos ', ' Frecuencias [Hz]'});
% display(tablafn);
% 
% tablaphi=table({'GDL 1';'GDL 2';'GDL 3'},phi,'variablename',{'GDL', ['MODO1 ' ' MODO2 ' ' MODO 3']});
% display(tablaphi);
% 
% tabla_dinak=table({'Modo 1';'Modo 2';'Modo 3'},[mk' kk' ck'],'variablename',{'Prop Dinamica', ['Mk ' ' Kk ' ' Ck']});
% display(tabla_dinak);


%% Parte 2 - pto a
%Usando la contribucion de los dos primeros modos
close all; clear all; clc;
f1=0.65; f2=1.82;
zheta1=0.02;
w1=2*pi*f1; w2=2*pi*f2;

%Ploteo de la funcion
x=@(t) (sin(w1*t)+(0.2*sin(w2*t)))*exp(-zheta1*w1*t);
figure()
fplot(x,[0 120]);
xlabel('Tiempo [s]'); ylabel('y(t)'); grid on;

%% Parte 2 --- Parte b
%funcion-original
close all; clear all; clc;
%f1=0.65; f2=1.82; zheta1=0;
f1=0.362; f2=2.349; zheta1=0;
w1=2*pi*f1; w2=2*pi*f2;
n=12; npt=2^n;
xx=@(t) sin(w1*t)+(0.2.*sin(w2*t))*exp(zheta1*w1*t);
%i.- Condiciones originales
Tmax=1000/f1;
dt=Tmax/npt; fs=1/dt; df=1/Tmax;
xx=@(t) sin(w1*t)+(0.2.*sin(w2*t));
tt=[0:dt:Tmax-dt]; y=feval(xx,tt);
f=[0:df:fs-df];
Yw=fft(y,npt)/npt; %Espectro de fourier

%ii.- eliminando leakage
fc=50; Tmaxii=fc/f1;
dtii=Tmaxii/npt; fsii=1/dtii; dfii=1/Tmaxii;
fii=[0:dfii:fsii-dfii];
tii=[0:dtii:Tmaxii-dtii]; yii=feval(xx,tii);
Ywii=fft(yii,npt)/npt;

%iii.- Eliminando Aliasing
fciii=fc; Tmaxiii=fciii/f1;
dtiii=Tmaxiii/npt; fsiii=1/dtiii; dfiii=1/Tmaxiii;
fiii=[0:dfiii:fsiii-dfiii];
tiii=[0:dtiii:Tmaxiii-dtiii]; yiii=feval(xx,tiii);
Ywiii=fft(yiii,npt)/npt;

%Ploteo parte b
figure()
plot(f(1:length(f)/2),abs(Yw(1:length(f)/2)*2),'displayname','Condicion Original'); hold on;
plot(fii,abs(Ywii)*2,'displayname','Minimizando Leakage');
xlim([0 2*f2]); xlabel('Frecuencias [Hz]'); ylabel('|H(w)| [m]'); grid on;
legend('show'); hold off;

% figure()
% plot(tt,y,'displayname','Señal original'); hold on;
% plot(tii,yii,'displayname','Minimizando Leakage');
% xlim([0 20]); xlabel('Tiempo [s]'); ylabel('y(t) [m]'); grid on;
% legend('show'); hold off;

figure()
plot(tt,y,'displayname','Señal original'); hold on;
plot(tiii,yiii,'displayname','Sin aliasing'); fplot(xx,[0 20]);
xlim([0 20]); xlabel('Tiempo [s]'); ylabel('y(t) [m]'); grid on;
legend('show'); hold off;

%% Parte 2 --- Parte c
%funcion-original
close all; clear all; clc;
%f1=1.0143; f2=2.6556; zheta1=0.02;
f1=0.362; f2=2.349; zheta1=0.02;
w1=2*pi*f1; w2=2*pi*f2;
n=12; npt=2^n;
xx=@(t) (sin(w1*t)+(0.2*sin(w2*t))).*exp(-zheta1*w1*t);
%i.- Condiciones originales
Tmax=1000/f1;
dt=Tmax/npt; fs=1/dt; df=1/Tmax;
tt=[0:dt:Tmax-dt]; y=feval(xx,tt);
f=[0:df:fs-df];
Yw=fft(y,npt)/npt; %Espectro de fourier

%Ploteo de funcion original
figure()
fplot(xx,[0 50]); xlabel('Tiempo [s]'); ylabel('y(t) [m]'); grid on;
legend('Original');

%ii.- eliminando leakage
fc=50; Tmaxii=fc/f1;
dtii=Tmaxii/npt; fsii=1/dtii; dfii=1/Tmaxii;
fii=[0:dfii:fsii-dfii];
tii=[0:dtii:Tmaxii-dtii]; yii=feval(xx,tii);
Ywii=fft(yii,npt)/npt;

%iii.- Eliminando Aliasing
fciii=fc; Tmaxiii=fciii/f1;
dtiii=Tmaxiii/npt; fsiii=1/dtiii; dfiii=1/Tmaxiii;
fiii=[0:dfiii:fsiii-dfiii];
tiii=[0:dtiii:Tmaxiii-dtiii]; yiii=feval(xx,tiii);
Ywiii=fft(yiii,npt)/npt;

%Ploteo parte c
figure()
plot(f,abs(Yw),'displayname','Condicion Original'); hold on;
plot(fii,abs(Ywii),'displayname','Minimizando Leakage');
xlim([0 2*f2]); xlabel('Frecuencias [Hz]'); ylabel('|H(w)| [m]'); grid on;
legend('show'); hold off;

figure()
plot(tt,y,'displayname','Señal original'); hold on;
plot(tiii,yiii,'displayname','Sin aliasing'); fplot(xx,[0 20]);
xlim([0 20]); xlabel('Tiempo [s]'); ylabel('y(t) [m]'); grid on;
legend('show'); hold off;


%% Parte 3----> Modelo teorico de 3GDL
close all, clc, clear all;
mnom=100; %Masa nominal de entreplanta [ton]
mi=mnom*1000; %Masa [kg]
Ixx=5.6969E-05; %Momento de inercia [m4]
E=2.1e+11; %Modulo de elasticidad del acero [N/m2]
h=3.00; %Altura [m]
Kxx=24*E*Ixx/h^3; %Rigidez del muelle [N/m]
M=[mi 0 0;0 mi 0;0 0 20*mi]; %Matriz de masas
K=[Kxx -Kxx 0;-Kxx 2*Kxx -Kxx;0 -Kxx 2*Kxx]; %Matriz de rigidez
z=[0.02 0.02 0.03]; %Vector de amortiguamiento
nn=11; %Multiplo para resolucion

[fn phi V mk kk ck dw Hdesp Hvel Hacc]=Matransferencia(K,M,z,nn);

npt=2^nn;
fmax=1.4*max(fn); df=fmax/npt;
f=0:df:fmax-df; %Vector de frecuencias [Hz]

%Matriz funcion de transferencia de desplazamiento
H11d=zeros(1,npt); H12d=zeros(1,npt); H13d=zeros(1,npt);
H21d=zeros(1,npt); H22d=zeros(1,npt); H23d=zeros(1,npt);
H31d=zeros(1,npt); H32d=zeros(1,npt); H33d=zeros(1,npt);
for mm=1:npt
    H11d(1,mm)=Hdesp(1,1,mm);
    H12d(1,mm)=Hdesp(1,2,mm);
    H13d(1,mm)=Hdesp(1,3,mm);
    H21d(1,mm)=Hdesp(2,1,mm);
    H22d(1,mm)=Hdesp(2,2,mm);
    H23d(1,mm)=Hdesp(2,3,mm);
    H31d(1,mm)=Hdesp(3,1,mm);
    H32d(1,mm)=Hdesp(3,2,mm);
    H33d(1,mm)=Hdesp(3,3,mm);
end

%Matriz funcion de transferencia de velocidad
H11v=zeros(1,npt); H12v=zeros(1,npt); H13v=zeros(1,npt);
H21v=zeros(1,npt); H22v=zeros(1,npt); H23v=zeros(1,npt);
H31v=zeros(1,npt); H32v=zeros(1,npt); H33v=zeros(1,npt);
for mm=1:npt
    H11v(1,mm)=Hvel(1,1,mm);
    H12v(1,mm)=Hvel(1,2,mm);
    H13v(1,mm)=Hvel(1,3,mm);
    H21v(1,mm)=Hvel(2,1,mm);
    H22v(1,mm)=Hvel(2,2,mm);
    H23v(1,mm)=Hvel(2,3,mm);
    H31v(1,mm)=Hvel(3,1,mm);
    H32v(1,mm)=Hvel(3,2,mm);
    H33v(1,mm)=Hvel(3,3,mm);
end


%Matriz funcion de transferencia de aceleracion
H11=zeros(1,npt); H12=zeros(1,npt); H13=zeros(1,npt);
H21=zeros(1,npt); H22=zeros(1,npt); H23=zeros(1,npt);
H31=zeros(1,npt); H32=zeros(1,npt); H33=zeros(1,npt);
for mm=1:npt
    H11(1,mm)=Hacc(1,1,mm);
    H12(1,mm)=Hacc(1,2,mm);
    H13(1,mm)=Hacc(1,3,mm);
    H21(1,mm)=Hacc(2,1,mm);
    H22(1,mm)=Hacc(2,2,mm);
    H23(1,mm)=Hacc(2,3,mm);
    H31(1,mm)=Hacc(3,1,mm);
    H32(1,mm)=Hacc(3,2,mm);
    H33(1,mm)=Hacc(3,3,mm);
end

%Grafica de matriz de funcion de transferencia Hij(w) -- desplazamiento

% figure()
% loglog(f,abs(H11d),'displayname','H11'); hold on;
% loglog(f,abs(H21d),'displayname','H21');
% loglog(f,abs(H31d),'displayname','H31'); legend('show');
% xlabel('Frecuencia [Hz]');ylabel('|H(w)| despl'); %xlim([1e-1 1e+1])
% hold off;
% 
% figure()
% semilogy(f,abs(H12d),'displayname','H12'); hold on;
% semilogy(f,abs(H22d),'displayname','H22');
% semilogy(f,abs(H32d),'displayname','H32'); legend('show');
% xlabel('Frecuencia [Hz]');ylabel('|H(w) despl|'); %xlim([0 100])
% hold off;
% 
% figure()
% semilogy(f,abs(H13d),'displayname','H13'); hold on;
% semilogy(f,abs(H23d),'displayname','H23');
% semilogy(f,abs(H33d),'displayname','H33'); legend('show');
% xlabel('Frecuencia [Hz]');ylabel('|H(w) despl|'); %xlim([0 100])
% hold off;

%Ploteo diagramas de nyquist - desplazamiento
% figure()
% subplot(3,2,1)
% plot(real(H11d),imag(H11d),'displayname','H11 despl'); legend('show');
% xlabel('Real(Hw)');ylabel('Im(Hw)');
% xlim([min(real(H11d)) max(real(H11d))]);
% subplot(3,2,2)
% plot(real(H12d),imag(H12d),'displayname','H12 despl'); legend('show');
% xlim([min(real(H12d)) max(real(H12d))]);
% xlabel('Real(Hw)');ylabel('Im(Hw)');
% subplot(3,2,3)
% plot(real(H13d),imag(H13d),'displayname','H13 despl'); legend('show'); 
% xlim([min(real(H13d)) max(real(H13d))]);
% xlabel('Real(Hw)');ylabel('Im(Hw)');
% subplot(3,2,4)
% plot(real(H22d),imag(H22d),'displayname','H22 despl'); legend('show'); 
% xlim([min(real(H22d)) max(real(H22d))]);
% xlabel('Real(Hw)');ylabel('Im(Hw)');
% subplot(3,2,5)
% plot(real(H23d),imag(H23d),'displayname','H23 despl'); legend('show'); 
% xlim([min(real(H23d)) max(real(H23d))]);
% xlabel('Real(Hw)');ylabel('Im(Hw)');
% subplot(3,2,6)
% plot(real(H33d),imag(H33d),'displayname','H33 despl'); legend('show'); 
% xlim([min(real(H33d)) max(real(H33d))]);
% xlabel('Real(Hw)');ylabel('Im(Hw)');

%Grafica de matriz de funcion de transferencia Hij(w) -- velocidad

% figure()
% loglog(f,abs(H11v),'displayname','H11'); hold on;
% loglog(f,abs(H21v),'displayname','H21');
% loglog(f,abs(H31v),'displayname','H31'); legend('show');
% xlabel('Frecuencia [Hz]');ylabel('|H(w)| velocidad'); %xlim([1e-1 1e+1])
% hold off;
% 
% figure()
% semilogy(f,abs(H12v),'displayname','H12'); hold on;
% semilogy(f,abs(H22v),'displayname','H22');
% semilogy(f,abs(H32v),'displayname','H32'); legend('show');
% xlabel('Frecuencia [Hz]');ylabel('|H(w)| velocidad'); %xlim([0 100])
% hold off;
% 
% figure()
% semilogy(f,abs(H13v),'displayname','H13'); hold on;
% semilogy(f,abs(H23v),'displayname','H23');
% semilogy(f,abs(H33v),'displayname','H33'); legend('show');
% xlabel('Frecuencia [Hz]');ylabel('|H(w)| velocidad'); %xlim([0 100])
% hold off;

%Ploteo diagramas de nyquist - velocidad
% figure()
% subplot(3,2,1)
% plot(real(H11v),imag(H11v),'displayname','H11 veloc'); legend('show');
% xlabel('Real(Hw)');ylabel('Im(Hw)');
% xlim([min(real(H11v)) max(real(H11v))]);
% subplot(3,2,2)
% plot(real(H12v),imag(H12v),'displayname','H12 veloc'); legend('show');
% xlim([min(real(H12v)) max(real(H12v))]);
% xlabel('Real(Hw)');ylabel('Im(Hw)');
% subplot(3,2,3)
% plot(real(H13v),imag(H13v),'displayname','H13 veloc'); legend('show'); 
% xlim([min(real(H13v)) max(real(H13v))]);
% xlabel('Real(Hw)');ylabel('Im(Hw)');
% subplot(3,2,4)
% plot(real(H22v),imag(H22v),'displayname','H22 veloc'); legend('show'); 
% xlim([min(real(H22v)) max(real(H22v))]);
% xlabel('Real(Hw)');ylabel('Im(Hw)');
% subplot(3,2,5)
% plot(real(H23v),imag(H23v),'displayname','H23 veloc'); legend('show'); 
% xlim([min(real(H23v)) max(real(H23v))]);
% xlabel('Real(Hw)');ylabel('Im(Hw)');
% subplot(3,2,6)
% plot(real(H33v),imag(H33v),'displayname','H33 veloc'); legend('show'); 
% xlim([min(real(H33v)) max(real(H33v))]);
% xlabel('Real(Hw)');ylabel('Im(Hw)');

%Grafica de matriz de funcion de transferencia Hij(w) -- aceleracion

figure()
semilogy(f,abs(H11),'displayname','H11'); hold on;
semilogy(f,abs(H21),'displayname','H21');
semilogy(f,abs(H31),'displayname','H31'); legend('show');
xlabel('Frecuencia [Hz]');ylabel('|H(w)| aceleración'); %xlim([1e-1 1e+1])
hold off;
% 
figure()
%semilogy(f,abs(H12),'displayname','H12'); hold on;
semilogy(f,abs(H22),'displayname','H22');
% semilogy(f,abs(H32),'displayname','H32'); legend('show');
xlabel('Frecuencia [Hz]');ylabel('|H(w)| aceleración'); %xlim([0 100])
hold off;
% 
figure()
semilogy(f,abs(H13),'displayname','H13'); hold on;
semilogy(f,abs(H23),'displayname','H23');
semilogy(f,abs(H33),'displayname','H33'); legend('show');
xlabel('Frecuencia [Hz]');ylabel('|H(w)| aceleración'); %xlim([0 100])
hold off;

%Ploteo diagramas de nyquist - aceleracion
% figure()
% subplot(3,2,1)
% plot(real(H11),imag(H11),'displayname','H11 acc'); legend('show');
% xlabel('Real(Hw)');ylabel('Im(Hw)');
% xlim([min(real(H11)) max(real(H11))]);
% subplot(3,2,2)
% plot(real(H12),imag(H12),'displayname','H12 acc'); legend('show');
% xlim([min(real(H12)) max(real(H12))]);
% xlabel('Real(Hw)');ylabel('Im(Hw)');
% subplot(3,2,3)
% plot(real(H13),imag(H13v),'displayname','H13 acc'); legend('show'); 
% xlim([min(real(H13)) max(real(H13))]);
% xlabel('Real(Hw)');ylabel('Im(Hw)');
% subplot(3,2,4)
% plot(real(H22),imag(H22),'displayname','H22 acc'); legend('show'); 
% xlim([min(real(H22)) max(real(H22))]);
% xlabel('Real(Hw)');ylabel('Im(Hw)');
% subplot(3,2,5)
% plot(real(H23),imag(H23),'displayname','H23 acc'); legend('show'); 
% xlim([min(real(H23)) max(real(H23))]);
% xlabel('Real(Hw)');ylabel('Im(Hw)');
% subplot(3,2,6)
% plot(real(H33),imag(H33),'displayname','H33 acc'); legend('show'); 
% xlim([min(real(H33)) max(real(H33))]);
% xlabel('Real(Hw)');ylabel('Im(Hw)');

%Tablas resumen
tablafn=table({'Modo 1';'Modo 2';'Modo 3'},fn,'variablename',{'Modos ', ' Frecuencias [Hz]'});
display(tablafn);

tablaphi=table({'GDL 1';'GDL 2';'GDL 3'},phi,'variablename',{'GDL', ['MODO1 ' ' MODO2 ' ' MODO 3']});
display(tablaphi);

tabla_dinak=table({'Modo 1';'Modo 2';'Modo 3'},[mk' kk' ck'],'variablename',{'Prop Dinamica', ['Mk ' ' Kk ' ' Ck']});
display(tabla_dinak);

aa=[max(abs(H11)) max(abs(H12)) max(abs(H13)) max(abs(H21)) max(abs(H22)) max(abs(H23)) max(abs(H31)) max(abs(H32)) max(abs(H33))];
max_aa=max(aa); display(max_aa);
ruido=(0.1/100)*max(aa); display(ruido);

%Graficas a mostrar
%4.1 H11 --- Desplazamiento LogLog
figure()
loglog(f,abs(H11d),'displayname','H11'); hold on;
legend('show'); xlabel('Frecuencia [Hz]');ylabel('|H(w)| despl'); xlim([1e-1 1e+1])
hold off;

%4.2 H21 --- Velocidad Nyquist
figure()
plot(real(H21v),imag(H21v),'displayname','H21 veloc','linestyle','none','marker','.'); legend('show');
xlim([min(real(H21v)) max(real(H21v))]);
xlabel('Real(Hw)');ylabel('Im(Hw)'); title('Diagrama Nyquist');

%4.3 H33 ---- Velocidad | Parte Real - Imag
figure()
plot(f,real(H33v),'displayname','H33 Real'); hold on;
plot(f,imag(H33v),'displayname','H33 Imag');
legend('show'); xlabel('Frecuencia [Hz]');ylabel('|H(w)| velocidad'); %xlim([0 100])
hold off;

%4.4 H22 ---- Aceleracion Loglog
figure()
loglog(f,abs(H22),'displayname','H22');
legend('show'); xlabel('Frecuencia [Hz]');ylabel('|H(w)| aceleración'); xlim([1e-2 1e1])
hold off;
