clear all; 
saveit = ~true;
optionalfilelabel = 'Rstapes0';         % leave empty if none
include_CFab = false;                   % for Maria GPig version 2
scaling_case = ~true;

%~~~~ Multiple Reflections ~~NO~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
use_Rstapes = true;
Rstapes = 0;                            % magnitude of Rstapes

%~~~~ Multiple Sources ~~NO~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% use "two point-source" model
use_twosources = ~true;
% CFs at source locations...
srcCF2 = 3.05;
srcCF1 = 2.95;

use_linear_freqs = true;
plotit = true;

%~~~~ Seed for SNRs ~~~~~NO~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% rng('default');
% rng('shuffle');
fix_seed = ~true;
if (fix_seed)                           % make it reproducibly random
  seed = 2^19+23^5;
  rng(seed);
end

use_same_seed_for_all_SNRs = true;
% if (use_same_seed_for_all_SNRs)
%   seed = rng;                           % stash it away
% end

% determine overall level...
Lsfe = 0;                               % target level
SNRs = [25];                            % in dB re 1

Nsegs = 3500;
Nsubj = 10;
Nreps = 1;

if (plotit & Nreps>1)
  disp('Too many reps: Not plotting intermediates');
  plotit = false;
end

x = detail(linscale(0,1,Nsegs));
CFs = cochlear_map(x,'guinea pig')/1000;

% specify propertiees of BM filters
Ntau = (1.78)*CFs.^0.44;                % SGO 2002 gpig Nsfe/2
if (include_CFab)
  CFab = 3;
  ok = find(CFs<=CFab);                          
  % note that CFs go from high to low 
  Ntau(ok) = Ntau(ok(1))*(CFs(ok)/CFab).^0.6; % steeper below CFa|b

  % round out the seam a little bit...
  ok = find(~(CFs>CFab*0.9&CFs<CFab*1.1));
  Ntau = exp(spline(log(CFs(ok)),log(Ntau(ok)),log(CFs)));
end

Qerb = 4*CFs.^0.35;

if (scaling_case)
  % scaling case (everywhere the same as 1 kHz CF)...
  Ntau = (1.78)*ones(size(CFs));
  Qerb = 4*ones(size(CFs));
end

tau = Ntau./CFs;			% delay in ms
db = (1./Qerb)/sqrt(2*pi);              % for Gaussian filter

% measurement frequencies...
fmax = 12; 
fmin = 0.5;
ppo = 65;                               % points/octave (Schairer et al. used 65)

df = 0.083;                             % 83 Hz resolution
Nfreqs = round((fmax-fmin)/df)+1;
freqs = linscale(fmin,fmax,Nfreqs);

Nsfe_act = interp1(CFs,2*Ntau,freqs);

r0 = 0.05;				% roughness magnitude (irrelevant)

% preallocate...
exitR = zeros(Nsubj,Nfreqs);
if (use_Rstapes)
  Rorig = zeros(Nsubj,Nfreqs);
end

Data = [];

legendStr = {};

for nsnr=1:numel(SNRs)
  SNR = SNRs(nsnr);
  disp(sprintf('Doing SNR %d of %d',nsnr,numel(SNRs)));
%   if (use_same_seed_for_all_SNRs)
%     rng(seed);
%   end
  for nr=1:Nreps
    disp(sprintf('Doing repetition %d of %d',nr,Nreps));

    for m=1:Nsubj
      ruff = r0*randn(size(x));		% roughness
      if (use_twosources)
        ruff = zeros(size(x));
        [tmpCF,idx2] = closest(CFs,srcCF2);
        [tmpCF,idx1] = closest(CFs,srcCF1);
        ruff(idx2) = r0;
        % ruff(idx1) = r0*(0.95-0.1i);
        ruff(idx1) = r0*exp(-i*pi/5);
      end
      
      for n=1:Nfreqs
        f = freqs(n);
        b = f./CFs;
        
        T = (exp(-(log(b)./(2*db)).^2).*exp(-2*pi*i*log(b).*Ntau)).^2;
        
        % no point in integrating more than necessary...
        ok = find(b>0.5&b<2);
        R(m,n) = numint(ruff(ok).*T(ok),x(ok));

        if (~true && m==1 && mod(n,15)==1)
          ok = find(b>0.5&b<2);
          figure(1)
          subplot(211)
          plot(x,dB(T));
          hold on
          plot(x(ok),dB(T(ok)),'r:');
          subplot(212)
          plot(x,cycs(T));
          hold on
        end
      end
      if (use_Rstapes)                  % add multiple internal reflections
        Rorig(m,:) = R(m,:) ./ max(abs(R(m,:)));
        % give it random starting phase
        Rst = Rstapes.*exp(1i*2*pi*rand(1));
        % make Rstapes interesting...
        % Rstapes = Rstapes * sin(2*pi*((freqs-freqs(1))/freqs(end)));
        R(m,:) = Rorig(m,:) ./ (1 - Rst.*Rorig(m,:));
      end
    end
    
    Lavg = mean(mean(dB(R)));
    if (use_Rstapes)
      LavgOrig = mean(mean(dB(Rorig)));
      LevelDiff = Lavg-LavgOrig;
      % use LavgOrig
      Lavg = LavgOrig;
    end

    R = R * adB(Lsfe-Lavg);
    % add noise...
    Noise = adB(Lsfe-SNR)*(randn(size(R))+i*randn(size(R)));
    Rnn = R;                                % no noise
    R = Rnn + Noise;
    if (use_Rstapes)
      Rorig = Rorig * adB(Lsfe-Lavg) + Noise;
    end

    % load the data...
    clear D;
    for m=1:Nsubj
      D(m).f = freqs;
      D(m).Psfe = R(m,:);
      D(m).Pnoise = Noise(m,:);
    end
    Data{nr} = D;
  end
  
  if (saveit)
    cwd = pwd;
    dir = ['DATA_SNR',num2str(SNR),'_',datestr(clock,'dd-mmm-yyyy_HH-MM')];
    if (~isempty(optionalfilelabel))
      dir = strcat(dir,'_',optionalfilelabel);
    end
    
    mkdir(dir);
    cd(dir);
    save DATA_simNsfe.mat                 % save entire workspace
    cd(cwd);
  end
  
end








