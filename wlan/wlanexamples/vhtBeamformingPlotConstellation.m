function vhtBeamformingPlotConstellation(sym,titleStr)
% VHTBeamformingExample Featured example helper function
% Plot constallation

%   Copyright 2015 The MathWorks, Inc.

[Nsd,NSym,Nss] = size(sym);

figure;
h = plot(squeeze(reshape(sym(:,:,end:-1:1),Nsd*NSym,1,Nss)),'.');
str = cell(Nss,1);
for i = 1:Nss
    str{i} = ['Spatial stream ' num2str(i)];
end
legend(h(end:-1:1),str,'Location','Best')
title(titleStr)
xlabel('Real')
ylabel('Imag')
xlim([-2 2])
ylim([-2 2])

end