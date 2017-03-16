function [zlayt,zlayb,vlay,varargout] = layerise(Z,V,dzmin,varargin)
% zz = 2e4:1e3:3e4;
% vv = 5e3 + 5e2*(zz-1e4*cos((zz-2e4)*2*pi/1e4))/1e4;
% Z = [0 1000 2000 2001 2002 19999,zz]';
% V = [3500 3500 3501 4000 4500 5000,vv]';

% profile on
% dzmin = 20;

% profile on

Zdo = Z(:);Vdo = V(:);

dz = diff(Zdo);
dv = diff(Vdo);
discgrad = 0.15; % gradients of more than 0.15 m/s per m are discontinuities
constgrad = 0.002; % gradients of less than 1 m/s per 100m are constant V
grad = dv./dz;

%% discontinuities
dind = find(abs(grad) >= discgrad); % find discontinuities
a = diff(dind);
b = find([a;inf]>1);
c = diff([0;b]);% length of sequences with discs
i1 = cumsum(c); % end points of sequences with discs
i0 = i1-c+1; % start points of sequences with discs
Ndisc = length(i0); % number of discontinuities

Zdisc = zeros(Ndisc,1);
Wdisc = zeros(Ndisc,1);
Vdto = zeros(Ndisc,1);
Vdbo = zeros(Ndisc,1);
kill = [];
for id = 1:Ndisc
    Zdisc(id) = mean(Zdo(dind(i0(id)):1+dind(i1(id))));
    Wdisc(id) = Zdo(dind(i1(id))+1) - Zdo(dind(i0(id)));
    Vdto(id) = Vdo(dind(i0(id)));
    Vdbo(id) = Vdo(dind(i1(id))+1);
    kill = [kill;find(Zdo>=Zdisc(id)-Wdisc(id)/2 & Zdo<=Zdisc(id)+Wdisc(id)/2)];
end
[Zdisc,Vdto,Vdbo,Wdisc]

%% make discontinuous V and Z
% kill within zone of discontinuity
Zdo(kill) = []; 
Vdo(kill) = [];
% add in discontinuities explicitly and sort into place
Zdo = [Zdo;Zdisc;Zdisc];
Vdo = [Vdo;Vdto;Vdbo];
[Zdo,isort] = sort(Zdo);
Vdo = Vdo(isort);
% add on surface vel if needed
if min(Zdo)~=0, Zdo = [0;Zdo]; Vdo = [Vdto(1);Vdo];end 
Nz = length(Zdo);

figure(1); clf, hold on
plot(V,Z,'-ko')
plot(Vdo,Zdo,'-ob')
set(gca,'ydir','reverse','ylim',[0, max(Zdo)],'xlim',[0.9*min(Vdo) 1.1*max(Vdo)])


%% Now can just work on making constant layers

% extablish easy constant layers
dz = diff(Zdo);
dv = diff(Vdo);
grad = dv./dz;
cind = find(abs(grad) <= constgrad); % find constant layers
a = diff(cind);
b = find([a;inf]>1);
c = diff([0;b]);% length of constant sequences
i1 = cumsum(c); % end points of constant sequences
i0 = i1-c+1; % start points of constant sequences
Nconst = length(i0); % number of constant layers

zlayt = zeros(Nconst,1);
zlayb = zeros(Nconst,1);
vlay = zeros(Nconst,1);
for id = 1:Nconst
    zlayt(id) = Zdo(cind(i0(id)));
    zlayb(id) = Zdo(cind(i1(id))+1);
    vlay(id) = mean(Vdo(cind(i0(id)):1+cind(i1(id))));
end

% any seds/surface?
if zlayt(1)~=0
    zsb = zlayt(1);
    zlayt = [0;zlayt];
    zlayb = [zsb;zlayb];
    vlay = [mean(Vdo(Zdo<=zsb));vlay];
end

% insert discontinuities as zero-thickness layers
[zlayt,isort] = sort([zlayt;Zdisc]);
zlayb = sort([zlayb;Zdisc]);
vtemp = [vlay;mean([Vdto,Vdbo],2)];
vlay = vtemp(isort);

% find gaps not yet filled by layers
[zlayt zlayb vlay]
gp = find([zlayt;max(Z)] - [0;zlayb] > 0);
zlayt_temp = [zlayt;max(Z)];
zlayb_temp = [0;zlayb];
for igap = 1:length(gp) % loop through regions between discontinuities
gpt = zlayb_temp(gp(igap));
gpb = zlayt_temp(gp(igap));
Zgp = Zdo(Zdo<=gpb & Zdo>=gpt);
Vgp = Vdo(Zdo<=gpb & Zdo>=gpt);
% kill wrong side of discontinuities if it's there
if sum(Zgp==gpt)>1, Zgp(1) = [];Vgp(1) = []; end
if sum(Zgp==gpb)>1, Zgp(end) = [];Vgp(end) = []; end

nsplit = ceil((gpb-gpt)/dzmin); 
zsplit = linspace(gpt,gpb,nsplit+1)';
vsplit = zeros(nsplit,1);
for isplit = 1:nsplit
    vsplit(isplit) = mean(Vgp(Zgp>=zsplit(isplit) & Zgp<=zsplit(isplit+1)));
end

zlayt = [zlayt;zsplit(1:end-1)];
zlayb = [zlayb;zsplit(2:end)];
vlay = [vlay;vsplit];

end

% re-sort layer top & bottoms
[zlayt,isort] = sort(zlayt);
zlayb = sort(zlayb);
vlay = vlay(isort);
% remove zero thickness layers (i.e. discontinuities)
wlay = zlayb-zlayt;
zlayt(wlay==0)=[];
zlayb(wlay==0)=[];
vlay(wlay==0)=[];
wlay(wlay==0)=[];

[zlayt zlayb vlay];



nlay = length(vlay);

zlayp = reshape([zlayt';zlayb'],2*nlay,1);
vlayp = reshape([vlay';vlay'],2*nlay,1);
plot(vlayp,zlayp,'-go')

% % 
% % 
% % % while any(diff(dind)==1) % collate discontinuities
% % %     ii = dind(find(diff(dind)==1,1,'first')) ;
% % %     dz(ii) = mean([dz(ii),dz(ii+1)]); dz(ii+1) = [];
% % %     dv(ii) = dv(ii) + dv(ii+1); dv(ii+1) = [];
% % %     grad = dv./dz;
% % %     dind = find(abs(grad) >= discgrad);
% % % end
% % 
% % cind = find(abs(grad) <= constgrad); % find constant layers
% % % while any(diff(cind)==1) % collate constant layers
% % %     ii = cind(find(diff(cind)==1,1,'first'));
% % %     dz(ii) = dz(ii) + dz(ii+1); dz(ii+1) = [];
% % %     dv(ii) = mean([dv(ii),dv(ii+1)]); dv(ii+1) = [];
% % %     grad = dv./dz;
% % %     cind = find(abs(grad) <= constgrad);
% % % end
% % gind = setdiff([1:length(grad)]',[dind;cind]);
% % Zdo = [0;cumsum(dz)];
% % Vdo = Vdo(1)+[0;cumsum(dv)];
% %     
% % figure(1); clf, hold on
% % plot(V,Z,'-ko')
% % plot(Vdo,Zdo,'-b')
% % % zlayp = reshape([zlayt;zlayb],2*nlay,1);
% % % vlayp = reshape([vlay;vlay],2*nlay,1);
% % % plot(vlayp,zlayp,'-g')
% % % plot(vlay,zlayp,'-g')
% % set(gca,'ydir','reverse','ylim',[0, max(Z)],'xlim',[3e3 max(V)+100])
% % %     zlayt
% % %     zlayb
% % %     vlay
% % % profile viewer
% % return
% %     
% %     
% % %% split Z into regions according to existing constant regions or dzmin 
% % Zdo = Z(1);
% % Vdo = V(1);
% % cumz = 0;
% % cumv = 0;
% % for iz = 1:Nz-1
% %     if dz(iz)>=dzmin || dz(iz)==0
% %         Zdo = nan;
% %         if cumz==0
% %             Zdo = [Zdo Z(iz+1)];
% %             Vdo = [Vdo V(iz+1)];
% %         else
% %             Zdo = [Zdo Zdo(end)+cumz Z(iz+1)];
% %             Vdo = [Vdo Vdo(end)+cumv V(iz+1)];
% %         end
% %         cumz = 0; cumv=0;
% %     elseif (iz>1) && (sign(dv(iz))~=sign(dv(iz-1)))
% %         Zdo = [Zdo Zdo(end)+cumz];
% %         Vdo = [Vdo Vdo(end)+cumv];
% %         newcumz = dz(iz);
% %         newcumv = dv(iz);
% %     else
% %         newcumv = cumv+dv(iz);
% %         newcumz = cumz+dz(iz);
% %         if newcumz>=dzmin
% %             fz = (newcumz-dzmin)/dz(iz);
% %             Zdo = [Zdo Zdo(end)+dzmin];
% %             Vdo = [Vdo Vdo(end)+cumv+fz*dv(iz)];
% %             newcumz = newcumz-dzmin;
% %             newcumv = (1-fz)*dv(iz);
% %         end
% %     end
% %     cumz = newcumz;
% %     cumv = newcumv;
% % end
% % % pin at the end
% % if cumz~=0
% %     Zdo = [Zdo Z(end)];
% %     Vdo = [Vdo V(end)];
% % end
% % % Z = Zdo;
% % % V = Vdo;
% %             
% % zto = Zdo(1:end-1);
% % zbo = Zdo(2:end);
% % dz = zbo-zto;
% % 
% % Nz = length(Zdo);
% % 
% % %% discretise regions into layers depending on gradient
% % zlayt = [];
% % zlayb = [];
% % vlay = [];
% % for iz = 1:Nz-1
% %     dv = diff(Vdo(iz:iz+1))./(zbo(iz)-zto(iz));
% %     if dv==0 % constant V layer
% % %         disp('const') 
% %         zlayt = [zlayt zto(iz)];
% %         zlayb = [zlayb zbo(iz)];
% %         vlay  = [vlay Vdo(iz)];
% %     elseif isinf(dv) % discontinuity
% %         continue
% %     else % gradient - must split it up
% %         nsplit = max([ceil((zbo(iz)-zto(iz))/dzmin),2]); 
% %         zsplit = linspace(zto(iz),zbo(iz),nsplit+1);
% %         vsplit = linspace(Vdo(iz),Vdo(iz+1),nsplit+1);
% %         
% %         zlayt = [zlayt zsplit(1:end-1)];
% %         zlayb = [zlayb zsplit(2:end)];
% %         vlay = [vlay midpts(vsplit)];
% %     end
% % end       
% % 
% % %% add on the end as a constant layer zero thickness;
% % zlayt = [zlayt Z(end)];
% % zlayb = [zlayb Z(end)];
% % vlay = [vlay V(end)];
% % 
% % nlay = length(vlay);
% % 
%% Do to other variables, using linterp
for iv = 1:length(varargin)
    ival = varargin{iv}; ival = ival(:);
    oval = nan(size(vlay));
    zzz = [min(Z):dzmin/100:max(Z)]; zzz=zzz(:);
    try
        val = linterp(Z,ival,zzz);
    catch
        t
    end
    for ilay = 1:length(vlay)
        oval(ilay) = mean(val(zzz>zlayt(ilay) & zzz<zlayb(ilay)));
    end
    % pin on the end
%     oval(end) = val(end);
    
    varargout{iv} = oval;
end
% % % profile viewer
% % 
% % % close all
% % % figure(1); clf, hold on
% % % plot(V,Z,'-b')
% % % plot(Vdo,Zdo,'-ko')
% % % zlayp = reshape([zlayt;zlayb],2*nlay,1);
% % % vlayp = reshape([vlay;vlay],2*nlay,1);
% % % plot(vlayp,zlayp,'-g')
% % % plot(vlay,zlayp,'-g')
% % % set(gca,'ydir','reverse','ylim',[0, max(Z)],'xlim',[3e3 max(V)+100])
% % %     zlayt
% % %     zlayb
% % %     vlay
% %     
% %     
end
    