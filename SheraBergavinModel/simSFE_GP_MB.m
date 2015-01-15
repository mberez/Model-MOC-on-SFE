function simSFE_GP_MB(ruff,Trial)

plotit = true;
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
fiber_freqs = fmin:fmax/Nfibers:40+1;
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
% Generating irregularities function, MOC fibers and irragularity
% function with MOC fibers affecting it
%         ruff = r0*randn(size(x));               % roughness
for cnt1 = 1:Nfibers
    f = fiber_freqs(cnt1)
    mFiber(:,cnt1) = generate_moc_fiber1(f);
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

    for nr=1:Nreps
        disp(sprintf('Doing repetition %d of %d',nr,Nreps));


        for m=1:Nsubj
            for n=1:Nfreqs
                f = freqs(n);
                b = f./CFs;
                T = (exp(-(log(b)./(2*db)).^2).*exp(-2*pi*i*log(b).*Ntau)).^2;
                % no pobint in integrating more than necessary...
                ok = find(b>0.5&b<2);
                indx = round(mean(ok));

                R(m,n) = numint(ruff(ok).*T(ok),x(ok));
                Rmoc(m,n) = numint(ruff_moc(ok).*T(ok),x(ok));

                %                 figure(FigN+1);
                %                 subplot(3,1,1);
                %                 scatter(x(indx), dB(abs(R(:,n))), 'bo');hold all;
                %                 scatter(x(indx), dB(abs(Rmoc(:,n))), 'ro');
                %                 xlim([0 1]);

                if n == floor(Nfreqs/2)
                    figure(FigN+2)
                    subplot(3,1,1); cla
                    plot(CFs,dB(T)); hold on
                    scatter(CFs(ok),dB(T(ok)),'rs');
                    subplot(3,1,2); cla;
                    %                 plot(x,ruff,'b');
                    %                 hold all
                    %                 plot(x(ok),ruff(ok),'g');
                    plot(CFs,unwrap(angle(T')));hold on
                    scatter(CFs(ok),unwrap(angle(T(ok)')),'rs');
                    subplot(3,1,3); cla
                    plot(CFs,ruff - ruff_moc,'b'); hold on;
                    plot(CFs(ok),ruff(ok) - ruff_moc(ok),'r');
                end
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
        % load the data...
        clear D;
        for m=1:Nsubj
            D(m).f = freqs;
            D(m).Psfe = R(m,:);
            D(m).Psfe_moc = Rmoc(m,:);
            D(m).Pnoise = Noise(m,:);
        end
        if(plotit && Nreps < 2)
            figure(FigN+1);
            subplot(2,1,1)
            plot(freqs, dB(abs(R(m,:))), '-b');hold all;
            plot(freqs, dB(abs(Rmoc1(m,:))), '-r');hold all;
            %     plot(freqs, dB(abs(Rmoc2(m,:))), '-m');hold all;
            scatter(freqs, dB(abs(Noise(m,:))), '+k');hold all;
            xlim([0.5 10])
            %             subplot(3,1,2)
            %             plot(freqs, unwrap(angle(R(m,:))), '-b');hold all;
            %             plot(freqs, unwrap(angle(Rmoc1(m,:))), '-r');hold all;
            %             xlim([0.5 10])
            %             subplot(3,1,3)
            %             plot(CFs, ruff - ruff_moc);
            %             xlim([0.5 12])
            subplot(2,1,2);hold all;
            scatter(freqs, dB(abs(Rmoc1(m,:))) - dB(abs(R(m,:))), 'o');
            xlim([0.5 10])
        end
        Data{nr} = D;
                saveas(FigN+1, num2str(Trial), 'png');
                saveas(FigN+1, num2str(Trial), 'fig')
    end
end








