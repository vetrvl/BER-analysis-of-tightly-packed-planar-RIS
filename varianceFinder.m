clc
clear
close all
A=1;

allpos_cos=cos(linspace(-pi/2,pi/2,1000));
allpos_phi=(linspace(-pi/2,pi/2,1000));
Npath = 100;
for dist_var = [1/8 ]
    
    SNRrang1=-45:2:-10;
    
    iterations =1e4;
    
    lambda = 1;
    d_H = dist_var*lambda; %Horizontal antenna spacing
    d_V = dist_var*lambda; %Vertical antenna spacing
    N = floor((A*lambda^2)/(d_H*d_V));
%     M_H = floor(sqrt(N));
%     M_V = floor(sqrt(N));
M_H = N;
M_V = 1;
U = zeros(3,N); %Matrix containing the position of the antennas
    % U = ones(3,M);
    i = @(m) mod(m-1,M_H); %Horizontal index
    j = @(m) floor((m-1)/M_H); %Vertical index
    
    for m = 1:N
        U(:,m) = [0; i(m)*d_H; j(m)*d_V]; %Position of the mth element
    end
    for nn = 1:N
        for mm = 1:N
            RR(nn,mm)=(sinc(2*norm(U(:,nn)-U(:,mm)))); % RIS R;
        end
    end
    %             hcorr=(mvnrnd(zeros(1,N),RR,bits_total)+1j*mvnrnd(zeros(1,N),RR,bits_total))./sqrt(2);
    for iteration=1:iterations

%         a =zeros(N,Npath);
        for ii = 1:Npath
            varphi = allpos_phi(randi(1000));
            theta=asin(rand);
            k = (-2*pi/lambda) * [cos(varphi)*cos(theta); cos(theta)*sin(varphi); sin(theta)];
            a(:,ii) = transpose(exp(1i* k'*U));
            cl(:,ii) =((A)*((randn+1j*randn)./sqrt(2)))*a(:,ii);
        end
        gg(iteration,:,:) = squeeze(sum(cl'))./sqrt(Npath);
    end
end
hcorr_npath = squeeze(gg);
hcorr = hcorr_npath;
%         n=(sigma).*(randn(1,bits_total) + 1i*randn(1,bits_total));

h_napth = var(sum(abs(hcorr')))
% RR=eye(N);
cha = (mvnrnd(zeros(1,length(RR)),A*RR,10000)+1j*mvnrnd(zeros(1,length(RR)),A*RR,10000))./sqrt(2);
%mean(sum(abs(cha')))
h_mvn=var(sum(abs(cha')))

% cc=sqrt(A)*(randn(N,10000)+1j*randn(N,10000))./sqrt(2);
% mean(sum(abs(cc)));
% %var(sum(abs(cha')))
theo = A*(4-pi)*sum(sum(RR))/4
histogram(sum(abs(hcorr')))
hold on
histogram(sum(abs(cha')))
