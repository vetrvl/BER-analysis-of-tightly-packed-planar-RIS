clc;
clear ;
 close all;
SNR_range1=-25:2:50;
% SNR_range1=40:2:50;
SNR_range=SNR_range1;
for iiii =1
    for dist_var = [0.0625 0.0625*2 0.0625*4 0.0625*8]
        dist_var
        lambda = 1;
        d_H = dist_var*lambda; %Horizontal antenna spacing
        d_V = dist_var*lambda; %Vertical antenna spacing
        %Define the antenna geometry
        M_H = 4;
        M_V = 4;
        N = M_H*M_V; %Total number of antennas
        U = zeros(3,N); %Matrix containing the position of the antennas
        % U = ones(3,M);
        i = @(m) mod(m-1,M_H); %Horizontal index
        j = @(m) floor((m-1)/M_H); %Vertical index
        for m = 1:N
            U(:,m) = [0; i(m)*d_H; j(m)*d_V]; %Position of the mth element
        end
        for nn = 1:N
            for mm = 1:N
                RR(nn,mm)=(sinc(2*norm(U(:,nn)-U(:,mm)))).^2; % RIS R;
            end
        end
        if iiii==4
            var1 = d_H*d_V*((4-pi)/4)*sum(sum(RR));
            mu1 =N*sqrt(d_H*d_V*pi/4);
        else
            b = iiii;
            [mean_quant, var_quant] = quantmeanvarfinder(b);
            var1 = d_H*d_V*var_quant*sum(sum(RR));
            mu1 =N*sqrt(d_H*d_V)*mean_quant;
        end
        %     var1 = d_H*d_V*((4-pi)/4)*sum(sum(RR));    %correlated case
        SNR_count = 1;
        for SNR=SNR_range1
            SNR% test SNR
            % number of reflecting elements in LIS
            M=2;                              % Signaling order
            bpcu=log2(M);                     % bits per channel use (SE)
            bits_total=bpcu*2e5;               % total number of bits % total number of bits (increase this number of higher SNR)
            decoded=zeros(1,bits_total);
            bit_seq=randi(2,1,bits_total)-1;   % bit sequence
            sigma=sqrt((1)/(2*(10^(SNR/10)))); % SNR=Es/NP (Es=1)  noise samples ~ CN(0,NP)  NP: noise power = 2*sigma^2
            ss=pskmod(0:M-1,M,0,'Gray');       % Gray-encoded M-PSK signal constellation
            ss=ss./sqrt(mean(abs(ss).^2));
            kk=1;
            bit_errors = 0;
            for iteration=1:bits_total/bpcu
                bits=bit_seq(kk:kk+bpcu-1);
                % Bin2Dec & Symbols
                dec=sum(bits(1:log2(M)).*[2.^(log2(M)-1:-1:0)])+1;
                x=ss(dec);   % virtual PSK symbol exp(j^wm)
                % h=(randn(N,1)+1i*randn(N,1))/sqrt(2);     % Channels between LIS and Destination
                % Determination of Phases
                %  Phi=exp(-1i*angle(h));
                %  If no optimization at LIS (BLIND scheme)
                % Phi=ones(N,1);
                % Received Signal
                n=(sigma).*(randn + 1i*randn);
                % A=sum(h.*Phi);
                %              A = sqrt(N*d_H*d_V*0.6318)*randn+N*sqrt(d_H*d_V*pi/4)*0.6366;
                A = sqrt(var1)*randn+mu1;
                % or simply A=sum(abs(h)) for intelligent transmission
                %           A=sum(h)      for blind transmission
                r = A * x + n;
                % Maximum Likelihood Detector
                metrics=zeros(1,M);
                ss1=metrics;
                count=1;
                for count1=1:M
                    metrics(count)=norm(r - A*ss(count1))^2;
                    ss1(count)=count1;
                    count=count+1;
                end
                [a1,a2]=min(metrics);
                % a1 % should be close to zero at high SNR
                % Dec2Bin
                dec1=zeros(1,log2(M));
                for pp=1:log2(M)
                    dec1(pp)=mod(floor((ss1(a2)-1)/(2^(log2(M)-pp))),2);
                end
                decoded(kk:kk+log2(M)-1)=dec1;
                kk=kk+log2(M);
            end
            bit_errors=sum(xor(bit_seq,decoded));
            BER(SNR_count)=bit_errors/bits_total;
            SNR_count = SNR_count+1;
        end
        
        % numerical integration
        count=1;
        % syms t
        %     for SNR=SNR_range
        %         sigma=sqrt((1)./(2*(10.^(SNR/10))));
        %         NP=2*sigma.^2;
        %         func =  @(t) (1/pi).*((1+2.*(var1./(NP.*((sin(t)).^2)))).^(-0.5)).*exp((-(mu1.^2)./(NP.*((sin(t)).^2)))./(1+2.*(var1./(NP.*((sin(t)).^2)))));
        %         Pb2(count)= integral(func,0,pi/2);
        %         count=count+1;
        %     end
        %     if iiii ==4
        %         strgg = strcat("Ideal phase shift Simulation");
        %         strgg1 = strcat("Ideal phase shift Theory");
        %     else
        %         strgg = strcat (num2str(iiii)," bit phase shift Simulation");
        %         strgg1 = strcat (num2str(iiii)," bit phase shift Theory ");
        %     end
        %     semilogy(SNR_range,Pb2,'-','LineWidth',1.5,'Displayname',strgg1);grid on;hold on;
        %     hold on
        semilogy(SNR_range1,BER,'-','LineWidth',1.5);
        hold on
    end
end
ylim([1e-5 1e0]);
xlim([-25 50]);
grid on
 xlabel('SNR [dB]');
  ylabel('BER');
%  legend show
