function [zlayt,zlayb,vlay,varargout] = layerise(Z,V,dzmin,varargin)
% zz = 2e4:1e3:3e4;
% vv = 5e3 + 5e2*(zz-1e4*cos((zz-2e4)*2*pi/1e4))/1e4;
% Z = [0 2000 2000 19999,zz];
% V = [3500 3500 4500 5000,vv];
% 
% dzmin = 1000;

% profile on

Z = Z(:);
dz = diff(Z);
dv = diff(V);
Nz = length(Z);

%% split Z into regions according to existing constant regions or dzmin 
Zdo = Z(1);
Vdo = V(1);
cumz = 0;
cumv = 0;
for iz = 1:Nz-1
    if dz(iz)>=dzmin || dz(iz)==0
        Zdo = nan;
        if cumz==0
            Zdo = [Zdo Z(iz+1)];
            Vdo = [Vdo V(iz+1)];
        else
            Zdo = [Zdo Zdo(end)+cumz Z(iz+1)];
            Vdo = [Vdo Vdo(end)+cumv V(iz+1)];
        end
        cumz = 0; cumv=0;
    elseif (iz>1) && (sign(dv(iz))~=sign(dv(iz-1)))
        Zdo = [Zdo Zdo(end)+cumz];
        Vdo = [Vdo Vdo(end)+cumv];
        newcumz = dz(iz);
        newcumv = dv(iz);
    else
        newcumv = cumv+dv(iz);
        newcumz = cumz+dz(iz);
        if newcumz>=dzmin
            fz = (newcumz-dzmin)/dz(iz);
            Zdo = [Zdo Zdo(end)+dzmin];
            Vdo = [Vdo Vdo(end)+cumv+fz*dv(iz)];
            newcumz = newcumz-dzmin;
            newcumv = (1-fz)*dv(iz);
        end
    end
    cumz = newcumz;
    cumv = newcumv;
end
% pin at the end
if cumz~=0
    Zdo = [Zdo Z(end)];
    Vdo = [Vdo V(end)];
end
% Z = Zdo;
% V = Vdo;
            
zto = Zdo(1:end-1);
zbo = Zdo(2:end);
dz = zbo-zto;

Nz = length(Zdo);

%% discretise regions into layers depending on gradient
zlayt = [];
zlayb = [];
vlay = [];
for iz = 1:Nz-1
    dv = diff(Vdo(iz:iz+1))./(zbo(iz)-zto(iz));
    if dv==0 % constant V layer
%         disp('const') 
        zlayt = [zlayt zto(iz)];
        zlayb = [zlayb zbo(iz)];
        vlay  = [vlay Vdo(iz)];
    elseif isinf(dv) % discontinuity
        continue
    else % gradient - must split it up
        nsplit = max([ceil((zbo(iz)-zto(iz))/dzmin),2]); 
        zsplit = linspace(zto(iz),zbo(iz),nsplit+1);
        vsplit = linspace(Vdo(iz),Vdo(iz+1),nsplit+1);
        
        zlayt = [zlayt zsplit(1:end-1)];
        zlayb = [zlayb zsplit(2:end)];
        vlay = [vlay midpts(vsplit)];
    end
end       

%% add on the end as a constant layer zero thickness;
zlayt = [zlayt Z(end)];
zlayb = [zlayb Z(end)];
vlay = [vlay V(end)];

nlay = length(vlay);

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
    oval(end) = val(end);
    
    varargout{iv} = oval;
end
% profile viewer

% close all
% figure(1); clf, hold on
% plot(V,Z,'-b')
% plot(Vdo,Zdo,'-ko')
% zlayp = reshape([zlayt;zlayb],2*nlay,1);
% vlayp = reshape([vlay;vlay],2*nlay,1);
% plot(vlayp,zlayp,'-g')
% plot(vlay,zlayp,'-g')
% set(gca,'ydir','reverse','ylim',[0, max(Z)],'xlim',[3e3 max(V)+100])
%     zlayt
%     zlayb
%     vlay
    
    
end
    