function [emdSTEM,EWamp] = PRISM02_multislice(emdSTEM, probe)

% Colin Ophus - 2020 Sept
% 02 - conventional multislice STEM for comparison, from same potential input

% Output detector settings
flagOutput3D = true;
drBins3D = 1 / 1000;  % spacing of bins in 3D output [rads]

% 4D output
flagOutput4D = true;
thicknessOutput = 10:50:500;
% thicknessOutput = 100;

emdSTEM.E0 = 300e3;       % Microscope voltage [Volts]
emdSTEM.probeSemiangleArray = probe / 1000;  % [rads]
emdSTEM.probeDefocusDF = 0;  % relative to the entrance surface of the potentials [Angstroms]

% % Probe positions
% dxy = emdSTEM.cellDim(1:2) / 500;
% xR = [0 1]*emdSTEM.cellDim(1);
% yR = [0 1]*emdSTEM.cellDim(2);
% emdSTEM.xp = (xR(1)+dxy/2):dxy:(xR(2)-dxy/2);
% emdSTEM.yp = (yR(1)+dxy/2):dxy:(yR(2)-dxy/2);
% Single probe position in center of UC
emdSTEM.xp = emdSTEM.cellDim(1) * 0.5;
emdSTEM.yp = emdSTEM.cellDim(2) * 0.5;

% Calculate wavelength and electron interaction parameter
m = 9.109383*10^-31;
e = 1.602177*10^-19;
c =  299792458;
h = 6.62607*10^-34;
emdSTEM.lambda = h/sqrt(2*m*e*emdSTEM.E0) ...
    /sqrt(1 + e*emdSTEM.E0/2/m/c^2) * 10^10; % wavelength in A
emdSTEM.sigma = (2*pi/emdSTEM.lambda/emdSTEM.E0) ...
    *(m*c^2+e*emdSTEM.E0)/(2*m*c^2+e*emdSTEM.E0);

% Fourier coordinates
qx = makeFourierCoords(emdSTEM.imageSize(1),emdSTEM.pixelSize(1));
qy = makeFourierCoords(emdSTEM.imageSize(2),emdSTEM.pixelSize(2));
[emdSTEM.qya,emdSTEM.qxa] = meshgrid(qy,qx);
emdSTEM.q2 = emdSTEM.qxa.^2 + emdSTEM.qya.^2;
emdSTEM.q1 = sqrt(emdSTEM.q2);

% Initial probe amplitude
qProbe = emdSTEM.probeSemiangleArray / emdSTEM.lambda;
dqx = qx(2) - qx(1);
dqy = qy(2) - qy(1);
emdSTEM.PsiInit = min(max( ...
    (qProbe*emdSTEM.q1 - emdSTEM.q2) ./ ...
    sqrt(dqx^2*emdSTEM.qxa.^2 + dqy^2*emdSTEM.qya.^2) ...
    + 0.5, 0), 1);
emdSTEM.PsiInit(1,1) = 1;
emdSTEM.PsiInit(:) = emdSTEM.PsiInit / sqrt(sum(abs(emdSTEM.PsiInit(:)).^2));

% Probe defocus
emdSTEM.probeDefocusC1 = -1 * emdSTEM.probeDefocusDF;
chi = (pi*emdSTEM.lambda*emdSTEM.probeDefocusC1)*emdSTEM.q2;
PsiInit = emdSTEM.PsiInit .* exp(-1i*chi);

% Probe shift basis
qxShift = -2i*pi*emdSTEM.qxa;
qyShift = -2i*pi*emdSTEM.qya;
sub = abs(PsiInit(:)) > 0;
qxShiftSub = qxShift(sub);
qyShiftSub = qyShift(sub);

% AA aperture reduced coordinate indices and mask
emdSTEM.xAA = [(1:(emdSTEM.imageSize(1)/4)) ...
    ((1-(emdSTEM.imageSize(1)/4)):0)+emdSTEM.imageSize(1)];
emdSTEM.yAA = [(1:(emdSTEM.imageSize(2)/4)) ...
    ((1-(emdSTEM.imageSize(2)/4)):0)+emdSTEM.imageSize(2)];
emdSTEM.maskAA = false(emdSTEM.imageSize);
emdSTEM.maskAA(emdSTEM.xAA,emdSTEM.yAA) = true;
emdSTEM.maskAAinv = ~emdSTEM.maskAA;

% propagator
emdSTEM.prop = exp( emdSTEM.q2(emdSTEM.xAA,emdSTEM.yAA) * ...
    (-1i*pi*emdSTEM.lambda*emdSTEM.sliceThickness));

% init output
PsiOutput = zeros(emdSTEM.imageSize / 2);

% Detector coordinates of 3D output
if flagOutput3D == true
    qxa = emdSTEM.qxa(emdSTEM.xAA,emdSTEM.yAA);
    qya = emdSTEM.qya(emdSTEM.xAA,emdSTEM.yAA);
    q1 = sqrt(qxa.^2 + qya.^2);
    
    emdSTEM.qMax = min(max(qxa(:,1)),max(qya(1,:)));
    alphaMax = emdSTEM.qMax * emdSTEM.lambda;
    emdSTEM.detectorAngles = (drBins3D/2):drBins3D:(alphaMax-drBins3D/2);
    numDetBins = length(emdSTEM.detectorAngles);
    
    % detector bin indices
    dqDet = drBins3D / emdSTEM.lambda;
    qDetInds = max(q1 / dqDet + 0.5,1);
    qF = floor(qDetInds);
    dq = qDetInds - qF;
    
    % lower index detector bins
    qDet1sub = qF <= numDetBins;
    qDetBins1 = qF(qDet1sub);
    qDetWeights1 = 1 - dq(qDet1sub);
    
    qDet2sub = qF <= numDetBins-1;
    qDetBins2 = qF(qDet2sub)+1;
    qDetWeights2 = dq(qDet2sub);
    
    % init 3D output array
    emdSTEM.output3D = zeros(length(emdSTEM.xp),length(emdSTEM.yp),numDetBins);
end


if flagOutput4D == true
    % planes to output thickness
    z = (1:emdSTEM.numPlanes)*emdSTEM.sliceThickness;
    indOutput4D = zeros(length(thicknessOutput),1);
    for a0 = 1:length(thicknessOutput)
        [~,indOutput4D(a0)] = min(abs(z - thicknessOutput(a0)));
    end
    
    emdSTEM.output4D = zeros( ...
        length(emdSTEM.xAA),...
        length(emdSTEM.yAA),...
        length(emdSTEM.xp),...
        length(emdSTEM.yp),...
        length(thicknessOutput),...
        'single');
    emdSTEM.output4Dinds = indOutput4D;
    emdSTEM.output4Dthicknesses = indOutput4D * emdSTEM.sliceThickness;
end



% Main loop
tic
for a0 = 1:emdSTEM.numFP
    trans = exp((1i*emdSTEM.sigma) * emdSTEM.pot(:,:,:,a0));
    
    for ax = 1:length(emdSTEM.xp)
        
        
        for ay = 1:length(emdSTEM.yp)
            
            % Shift probe
            Psi = PsiInit;
            Psi(sub) = Psi(sub) .* exp( ....
                qxShiftSub*emdSTEM.xp(ax) + qyShiftSub*emdSTEM.yp(ay));
            
            % Propagate through foil
            for a2 = 1:emdSTEM.numPlanes
                Psi(:) = fft2(ifft2(Psi) .* trans(:,:,a2));
                if a2 < emdSTEM.numPlanes
                    Psi(emdSTEM.xAA,emdSTEM.yAA) = ...
                        Psi(emdSTEM.xAA,emdSTEM.yAA) .* emdSTEM.prop;
                    Psi(emdSTEM.maskAAinv) = 0;
                end
                
                if flagOutput4D == true
                    [val,ind] = min(abs(indOutput4D - a2));
                    if val == 0
                        emdSTEM.output4D(:,:,ax,ay,ind) = ...
                            emdSTEM.output4D(:,:,ax,ay,ind) + ...
                            abs(Psi(emdSTEM.xAA,emdSTEM.yAA)).^2;
                    end
                end
                
            end
            
        end
        PsiOutput(:) = abs(Psi(emdSTEM.xAA,emdSTEM.yAA)).^2;
        
        if flagOutput3D == true
            emdSTEM.output3D(ax,ay,:) = squeeze(emdSTEM.output3D(ax,ay,:)) + ...
                accumarray(qDetBins1, ...
                PsiOutput(qDet1sub) .* qDetWeights1,[numDetBins 1]) + ...
                accumarray(qDetBins2, ...
                PsiOutput(qDet2sub) .* qDetWeights2,[numDetBins 1]);
        end
    end
    
end
if flagOutput3D == true
    emdSTEM.output3D(:) = emdSTEM.output3D / emdSTEM.numFP;
end
if flagOutput4D == true
    emdSTEM.output4D(:) = emdSTEM.output4D / emdSTEM.numFP;
end
toc

EWamp = sqrt(PsiOutput);

% fig1 = figure(1002)
%% clf
% plot(emdSTEM.detectorAngles * 1e3,...
%     squeeze(emdSTEM.output3D(ax,ay,:)),...
%     'linewidth',2,'color','r')
 % plot(emdSTEM.detectorAngles * 1e3,...
 %     squeeze(emdSTEM.output3D(ax,ay,:)) ./ emdSTEM.detectorAngles(:),...
 %     'linewidth',2,'color','r')

%saveas(fig1,'trainset_input1','png')

%fig2 =  figure(1001)
% clf
% EWamp = sqrt(PsiOutput);
% imagesc(fftshift(EWamp))
% % imagesc(fftshift(abs(Psi)))
% % imagesc(fftshift(abs(emdSTEM.PsiInit)))
% % imagesc((abs(ifft2(Psi))))
% axis equal off
% %colormap(jetBlack)
% colorbar
% set(gca,'position',[0 0 0.85 1])

%saveas(fig2,'trainset_input2','png')

end
