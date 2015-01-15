function simSFEwMOC(varargin)
if nargin >= 1
    a= varargin(1);  ruff = a{1};
    if nargin >= 2
        a= varargin(2);  Trial = a{1};
    end
end

doPlots = 1;
Data = [];               % data prealocate
FigN =100;
%~~~ Level and other parameters ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Lsfe = 0;                % target level
SNRs = [25];             % in dB re 1
Attens = [0.3];
Nsegs = 3500;            % I think this # should be = to # of OHCs in length
Nsubj = 1;
Nreps = 1;

%~~~ Frequency Specifications ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
fmax = 10;
fmin = 0.5;
df = 0.083;              % 83 Hz resolution
r0 = 0.05;				 % roughness magnitude (irrelevant)

%~~~ Parmaeters for Fibers ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Nfibers = 100;
fiber_freqs = fmin:fmax/Nfibers:fmax+1;
fiber_freqs = fiber_freqs(1:Nfibers);

x = detail(linscale(0,1,Nsegs));
CFs = cochlear_map(x,'guinea pig')/1000;

% specify propertiees of BM filters
Ntau = (1.78)*CFs.^0.44;                % SGO 2002 gpig Nsfe/2
Qerb = 4*CFs.^0.35;
tau = Ntau./CFs;                        % delay in ms
db = (1./Qerb)/sqrt(2*pi);              % for Gaussian filter

Nfreqs = round((fmax-fmin)/df)+1;
freqs = linscale(fmin,fmax,Nfreqs);

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Generating irregularities function if not prvided 
if nargin == 0
        ruff = r0*randn(size(x));               % roughness
        Trial = 190;
end

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Generating irregularities function, MOC fibers and irragularity
% function with MOC fibers affecting it

for cnt1 = 1:Nfibers
    f = fiber_freqs(cnt1)
    mFiber(:,cnt1) = generate_moc_fiberx(f);
end
allOHCs = sum(mFiber,2);
indxF = find(allOHCs>0);
disp(sprintf('%d fibers affecting %d OHCs',Nfibers, length(indxF)));
innervation_pattern = allOHCs;
invtdOHCs = allOHCs(indxF).^(-1).*Attens;
allOHCs(indxF) = invtdOHCs;
indxNF = find(allOHCs==0);
allOHCs(indxNF) = 1;
ruff_moc = ruff.*allOHCs';

figure(FigN); cla;
subplot(4,1,1)
plot(CFs, ruff);
xlim([0.5 12])
subplot(4,1,3)
plot(CFs, ruff_moc);
xlim([0.5 12])
subplot(4,1,2)
scatter(CFs(find(innervation_pattern)), innervation_pattern(find(innervation_pattern)), 'rx');
xlim([0.5 12])
subplot(4,1,4)
plot(CFs, ruff - ruff_moc);
xlim([0.5 12])
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





for nsnr=1:numel(Attens)
    SNR = SNRs(1);
    Atten = Attens(nsnr);
    disp(sprintf('Doing Atten %d of %d',nsnr,numel(Attens)));
    for n=1:Nfreqs
        f = freqs(n);
        b = f./CFs;
        T = (exp(-(log(b)./(2*db)).^2).*exp(-2*pi*i*log(b).*Ntau)).^2;
        Tmoc = (ruff_moc.*T)./(ruff);
        % no pobint in integrating more than necessary...
        ok = find(b>0.5&b<2);
        indx = round(mean(ok));

        R(:,n) = numint(ruff(ok).*T(ok),x(ok));
        Rmoc(:,n) = numint(ruff_moc(ok).*T(ok),x(ok));

        if n == floor(Nfreqs/2)
            figure(FigN+2)
            subplot(3,1,1); cla
            plot(CFs,dB(T)); hold on
            plot(CFs,dB(Tmoc),'-r');
            subplot(3,1,2); cla;
            plot(CFs,dB(T) - dB(Tmoc),'-r');
%             plot(CFs,unwrap(angle(T')));hold on
%             scatter(CFs(ok),unwrap(angle(T(ok)')),'rs');
            subplot(3,1,3); cla
            plot(CFs,ruff - ruff_moc,'b'); hold on;
            plot(CFs(ok),ruff(ok) - ruff_moc(ok),'r');
        end
    end


    Lavg = mean(mean(dB(R)));
    Lavg_moc = mean(mean(dB(Rmoc)));
    R = R * adB(Lsfe-Lavg);
    Rmoc1 = Rmoc* adB(Lsfe-Lavg);
    %     Rmoc2 = Rmoc* adB(Lsfe-Lavg_moc);

    % add noise...
    Noise = adB(Lsfe-SNR)*(randn(size(R))+i*randn(size(R)));
    Rnn = R;                                % no noise
    R = Rnn + Noise;
    Rmoc = Rmoc + Noise;
    
%~~~~~   
    
    [RHO,PVAL] = corr(dB(abs(R')), dB(abs(Rmoc')));
    
    % load the data...
    clear D;
    for m=1:Nsubj
        D(m).f = freqs;
        D(m).Psfe = R(m,:);
        D(m).Psfe_moc = Rmoc(m,:);
        D(m).Pnoise = Noise(m,:);
    end
    if(doPlots)
        figure(FigN+1); 
        subplot(2,1,1); cla;
        plot(freqs, dB(abs(R(m,:))), '-b');hold all;
        plot(freqs, dB(abs(Rmoc1(m,:))), '-r');hold all;
        text(4, 15,[' Rho = ' num2str(RHO)],'FontSize',14)
        text(4, 11,[' p = ' num2str(PVAL)],'FontSize',14)
        scatter(freqs, dB(abs(Noise(m,:))), '+k');hold all;
        ylim([-20 20])
        xlim([0.5 10])
        %             subplot(3,1,2)
        %             plot(freqs, unwrap(angle(R(m,:))), '-b');hold all;
        %             plot(freqs, unwrap(angle(Rmoc1(m,:))), '-r');hold all;
        %             xlim([0.5 10])
        subplot(2,1,2);cla; hold all;
        scatter(freqs, dB(abs(Rmoc1(m,:))) - dB(abs(R(m,:))), 'o');
        xlim([0.5 10])
        ylim([-10 10])
    end
    Data{1} = D;
    saveas(FigN+1, num2str(Trial), 'png');
    saveas(FigN+1, num2str(Trial), 'fig')

end








