clear;
close all;

datasize = 100000;
EbNo = 0:2:20;
M = 4;
x = randsrc(2,datasize/2,[0:3]);
x1 = pskmod(x,M,pi/4);% 调制后的信号
h = randn(4,datasize/2) + 1j*randn(4,datasize/2);
h = h./sqrt(2);

for index = 1:length(EbNo)
    sigmal = sqrt(1/(4*10.^(EbNo(index)/10)));
    n = sigmal*(randn(2,datasize/2) + 1j*randn(2,datasize/2));
    y = x1 + n; % SISO AWGN
    y1 = x1 + n./h(1:2,:);
    x2 = pskdemod(y,M,pi/4);
    x3 = pskdemod(y1,M,pi/4);
    sigmal2 = sqrt(1/(2*10.^(EbNo(index)/10)));

    n = sigmal2*(randn(4,datasize/2) + 1j*randn(4,datasize/2));
    n1(1,:) = (conj(h(1,:)).* n(1,:)+ h(2,:).* conj(n(2,:)))./ (sum(abs(h(1:2,:)).^2));
    n2(2,:) = (conj(h(2,:)).* n(1,:)- h(1,:).* conj(n(2,:)))./ (sum(abs(h(1:2,:)).^2));

    y = x1 + n1;
    x4 = pskdemod(y,M,pi/4);
    n2(1,:) = (conj(h(1,:)).* n(1,:) + h(2,:) .* conj(n(2,:))+ ...
                conj(h(3,:)).*n(3,:) + h(4,:) .* conj(n(4,:))) ./ (sum(abs(h).^2));
    n2(2,:) = (conj(h(2,:)).* n(1,:) - h(1,:) .* conj(n(2,:))+ ...
                conj(h(4,:)).*n(3,:) - h(3,:) .* conj(n(4,:))) ./ (sum(abs(h).^2));
    y1 = x1 + n2;
    x5 = pskdemod(y1,M,pi/4);

    [~,ber1(index)] = biterr(x,x2,log2(M));
    [~,ber2(index)] = biterr(x,x3,log2(M));
    [~,ber3(index)] = biterr(x,x4,log2(M));
    [~,ber4(index)] = biterr(x,x5,log2(M));
end
semilogy(EbNo, ber1, '-r*',EbNo, ber2,'-go',EbNo, ber3,'-bd',EbNo,ber4,'-k');
grid on
legend('AWGN','Rayleigh channel','2x1 Alamouti Code','2x2 Alamouti Code');
