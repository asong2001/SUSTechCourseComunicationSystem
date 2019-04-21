clear all
close all
clc
datasize = 100000;
EbNo = 0:2:20;
M = 16;
x = randsrc(2,datasize/2,[0:15]);
x1 = pskmod(x,M,pi/M);
h = randn(4,datasize/2)+1j*randn(4,datasize/2);
h = h./sqrt(2);
for indx = 1:length(EbNo)
    sigma1 = sqrt(1/(4*10.^(EbNo(indx)/10)));
    n = sigma1*(randn(2,datasize/2)+1j*randn(2,datasize/2));
    y = x1+n;%SISO AWGN
    y1 = x1+n./h(1:2,:);%SISO Rayleigh
    x2 = pskdemod(y,M,pi/M);
    x3 = pskdemod(y1,M,pi/M);
    sigma2 = sqrt(1/(2*10.^(EbNo(indx)/10)));
    n = sigma2*(randn(4,datasize/2)+1j*randn(4,datasize/2));
    n1(1,:) = (conj(h(1,:)).*n(1,:)+h(2,:).*conj(n(2,:)))./(sum(abs(h(1:2,:)).^2));
    n1(2,:) = (conj(h(2,:)).*n(1,:)-h(1,:).*conj(n(2,:)))./(sum(abs(h(1:2,:)).^2));
    y = x1+n1;
    x4 = pskdemod(y,M,pi/M);
    n2(1,:) = (conj(h(1,:)).*n(1,:)+h(2,:).*conj(n(2,:))+...
               conj(h(3,:)).*n(3,:)+h(4,:).*conj(n(4,:)))./(sum(abs(h).^2));
    n2(2,:) = (conj(h(2,:)).*n(1,:)-h(1,:).*conj(n(2,:))+...
               conj(h(4,:)).*n(3,:)-h(3,:).*conj(n(4,:)))./(sum(abs(h).^2));
    y1 = x1+n2;
    x5 = pskdemod(y1,M,pi/M);
    [temp,ber1(indx)] = biterr(x,x2,log2(M));
    [temp,ber2(indx)] = biterr(x,x3,log2(M));
    [temp,ber3(indx)] = biterr(x,x4,log2(M));
    [temp,ber4(indx)] = biterr(x,x5,log2(M));
end
semilogy(EbNo,ber1,'-r*',EbNo,ber2,'-go',EbNo,ber3,'-bd',EbNo,ber4,'-k.')
grid on
axis([0 20 10^-5 1])
set(gca,'XTick',0:2:20);
ylabel('BER')
xlabel('EbNo(dB)')
legend('AWGN Channel','SISO Rayleigh channel','Alamouti 2x1','Alamouti 2x2')