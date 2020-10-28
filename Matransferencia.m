function [fn Vn V mk kk ck dw Hdesp Hvel Hacc]=Matransferencia(K,M,zheta,nn)
n=size(K); n=n(1); %Dimension de la matriz de rigidez y de masas
[V D]=eig(K,M); %Autovector y autovalores
wn=sqrt(diag(D)); %Vector de frecuencias modales [rad/s]
fn=wn/(2*pi); %Vector de frecuencias modales [Hz]
npt=2^nn; %Numero de puntos
%Definicion del vector de frecuencias
a=1.4*max(fn); wmax=2*pi*a; dw=wmax/npt;
w=0:dw:wmax-dw; %Vector de frecuencias [rad/s]
f=w./(2*pi); %Vector de frecuencias [Hz]
betta=zeros(npt,n); %Matriz de relacion de frecuencias (w/wn)
%Definicion y los valores caracteristicos y normalizacion del vector modal
mk=zeros(1,n); %Vector para normalizar el vector modal PHI*M*PHI
kk=zeros(1,n); %Vector de rigidez caracteristica
Vn=zeros(n); %Matriz modal normalizada  PHI*M*PHI=1
for h=1:n
V(:,h)=V(:,h)/max(abs(V(:,h)));
mk(1,h)=V(:,h)'*M*V(:,h); %masa caracteristica
Vn(:,h)=V(:,h)./sqrt(mk(1,h)); %Normalizacion final del vector modal
kk(1,h)=Vn(:,h)'*K*Vn(:,h); %masa caracteristica
mk(1,h)=Vn(:,h)'*M*Vn(:,h); %recalculo masa caracteristica para comprobar
ck(1,h)=zheta(1,h)*2*mk(1,h)*wn(h); %Coeficiente de amortiguamiento
betta(:,h)=w./wn(h);
end
%Generacion funciones de transferencia
h2=zeros(n); h2d=zeros(n); h2v=zeros(n);
Hacc=zeros(n,n,npt); %Matriz de funciones de transferencia de aceleraciones
Hvel=zeros(n,n,npt); %Matriz de funciones de transferencia de velocidad
Hdesp=zeros(n,n,npt); %Matriz de funciones de transferencia de desplazamiento
 for h=1:npt
    opa_acc=zeros(n); %Actualizar matriz en cada paso de tiempo
    opa_vel=zeros(n);
    opa_desp=zeros(n);
    for hh=1:n
     h1d=(1/(wn(hh)^2))/(1+(2*1i*zheta(1,hh)*betta(h,hh))-(betta(h,hh)^2)); %Hw de desplazamiento
     h2d=(Vn(:,hh)*Vn(:,hh)').*h1d;
     h1v=(1i*betta(h,hh)/wn(hh))/(1+(2*1i*zheta(1,hh)*betta(h,hh))-(betta(h,hh)^2)); %Hw de velocidad
     h2v=(Vn(:,hh)*Vn(:,hh)').*h1v;
     h1=(-(betta(h,hh)^2))/(1+(2*1i*zheta(1,hh)*betta(h,hh))-(betta(h,hh)^2)); %Hw de aceleracion
     h2=(Vn(:,hh)*Vn(:,hh)').*h1;
     opa_desp=opa_desp+h2d;
     opa_vel=opa_vel+h2v;
     opa_acc=opa_acc+h2;
    end
     Hacc(:,:,h)=opa_acc;
     Hvel(:,:,h)=opa_vel;
     Hdesp(:,:,h)=opa_desp;
 end
end