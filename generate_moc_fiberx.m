function mFiber = generate_moc_fiber1(f)
doPlots = 1;
nOHCs = 3500;
fmax = 12; 
fmin = 0.5;
x = detail(linscale(0,1,nOHCs));
CFs = cochlear_map(x,'guinea pig')/1000;
b = f./CFs;  

% index that corresponds to CF
indx = find(b> 0.999 & b < 1.001); indx = indx(end);

% number of arborizations:
% Narb = ceil(rand*3);
Narb = 1;
% Span can be from 2.5 to 25 % 
SpanType = ceil(rand*25);
while SpanType <= 2.5
    SpanType = ceil(rand*25);
    Span = nOHCs*SpanType*0.01;
end
if SpanType <= 10 && SpanType > 2.5
    Span = nOHCs*SpanType*0.01;
end
if SpanType > 10 
    Span = nOHCs*SpanType*0.01; 
end   

% number of OHCs per arbor
cnt = 1;
for narb = 1:Narb
    while cnt < 5
    cnt = ceil(rand*20);
    end
    Nohc(narb,:) = cnt;
    cnt = 1;
end

%~~~ Constructing the Innervation Function ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
mFiber = zeros(numel(CFs), 1);
ApicalOffset = round(nOHCs*0.05);
ApicalBoundary = indx+ApicalOffset;
BasalBoundary = ApicalBoundary - Span; 
if BasalBoundary < 0
while BasalBoundary < 0
Span = nOHCs*3*0.01;
BasalBoundary = ApicalBoundary - Span; 
end
end

mFiber(ApicalBoundary) =1;
mFiber(BasalBoundary) =1;

for narb = 1:Narb
    for nohc = 1:Nohc-2
    cnt = ceil(rand*Span);
    indx2 = ApicalBoundary - cnt;
    mFiber(indx2) =1;
    end
end

if doPlots == 1
figure(2); cla;
scatter(CFs(find(mFiber == 1)), mFiber(find(mFiber == 1)), 'rx'); hold on;
plot([CFs(BasalBoundary) CFs(BasalBoundary)], [-2 2], '-b');
plot([CFs(ApicalBoundary) CFs(ApicalBoundary)], [-2 2], '-b');
plot([CFs(indx) CFs(indx)], [-2 2], '-r');
xlim([fmin fmax]) 
text(6,-1,[' nOHC ' num2str(Nohc)],'FontSize',18)
text(6,-1.25,[' Span ' num2str(SpanType) '%'],'FontSize',18)
end
end