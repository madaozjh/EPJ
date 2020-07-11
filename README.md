# EPJ
matlab details for EPJ Quantum technology
%% calculate the time-evolution of the probability from different intial states
clear;
J=1;
mft=0;
norm=0;
z=10000;
tau=1;
h=0.01; %the step of regkuta
for cc=0:0
F(1:floor(z/tau))=0;
phi(1:floor(z/tau))=0;
a(1:floor(z/tau))=0;
b(1:floor(z/tau))=0;
A(1:floor(z/tau))=0;
B(1:floor(z/tau))=0;
e(1:floor(z/tau))=0;
f(1:floor(z/tau))=0;
mft=0;
norm=0;
tau0=0;
c=11+2*cc;
d=3+2*cc;%the hopping strength between different structure
N=c; % the number of structures
H=zeros(N,N); % the hamiltonian


H(1:c,1:c)=0;
for i=1:d-1
    H(i,i+1)=0.9+0.2*rand(1);
    H(i+1,i)=0.9+0.2*rand(1);
end
for i=d+1:c-1
H(i,i+1)=0.9+0.2*rand(1);
H(i+1,i)=0.9+0.2*rand(1);
end
H(d+1,(d+1)/2)=0.9+0.2*rand(1);
H((d+1)/2,d+1)=0.9+0.2*rand(1);
H=-H;
rho(1:N,1:z)=0;
rho(1,1)=1;
Rho(1:N,1:z)=0;
Rho(d,1)=1;
%% 
for i=2:1:z;
rho(1:c,i)=expm(-1i*H*h*i)*rho(1:c,1);
Rho(1:c,i)=expm(-1i*H*h*i)*Rho(1:c,1);
end
for i=1:floor(z/tau)
    a(i)=real(rho(d,i));
    b(i)=imag(rho(d,i));
    A(i)=real(Rho(d,i));
    B(i)=imag(Rho(d,i));
end
e(1)=a(1)/tau/h;
f(1)=b(1)/tau/h;
for i=2:1:floor(z/tau)%计算F（tau）
    e(i)=a(i)/tau/h;
    f(i)=b(i)/tau/h;
    for j=1:i-1
    e(i)=e(i)-e(j)*A(i-j)+f(j)*B(i-j);
    f(i)=f(i)-e(j)*B(i-j)-f(j)*A(i-j);
    end
    F(i)=e(i)*e(i)+f(i)*conj(f(i));
    if F(i)<=0.0001 && F(i-1)>0.0001;
        tau0=i;
        %break
    end
    plot(i*h*tau,F(i),'r.')
    hold on
    %plot(i*h*tau,rho(d,i*tau)*conj(rho(d,i*tau)),'ro')
    %hold on
end
for i=1:floor(z/tau)
mft=mft+F(i)*i*tau*h*tau*h;
norm=norm+F(i)*h*tau
end
mft=mft/norm;
Mft(cc+1)=mft;
end
return
fid1=fopen('D:\2\G2.txt','wt');
for i=1:1:21
 p1=Mft(i);  
fprintf(fid1,'%11.8f %13.8f\r\n',i*h,p1);
end
