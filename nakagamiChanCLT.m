

clc;
clear all;
%nakagamiFading Channel
 m =3.15;mu=m/2; % value of fading parameter
 N=10^6; % Number of samples generated for verification of PDF
 Nr=16;
 Nt=N;
 sn=[-1 1]; %variable to define sign of the generated coefficients
 %% Nak−m function starts
 n=sn((rand(Nr,Nt)>0.5)+1).*sqrt(gamrnd(mu,1/mu,Nr,Nt)/2) ...
+j*sn((rand(Nr,Nt)>0.5)+1).*sqrt(gamrnd(mu,1/mu,Nr,Nt)/2);
%%Instructions for obtaining uniformly distributed phase
phi=2*pi*rand(Nr,Nt);
n=abs(n);
H=(n).*cos(phi)+j*(n).*sin(phi);

H=sum(abs(H));
