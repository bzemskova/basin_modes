function params = gendata( om, k, rdir, flags)

params = gendata_params();
Lx_basin = params.Lxbasin;
basin_alpha = params.alpha;
filename = params.filename;

rname = sprintf('run_om%.8f',om);
prec='real*4';
ieee='b';
fs = 12; fn = 'times';

g = params.g;
f = params.f; 
np = [8 4]; % processors [nx, ny]

% sponge cells on north, east and west boundaries
nsponge = 15;

% vertical grid parameters
Lz = params.Lz; %4000; % [m]

% this is a crude estimate for the KW neglecting the shelf
cp = sqrt(g*Lz);
lam = 2*pi/k;

Ld = cp/f; % offshore length scale

% vertical grid, exponential high res near surface
Lz_tot = params.Lz_tot;
nzc = 30;
mindz = 25; maxdz = 500;
dz = smooth([ones(1,floor(700/mindz))*mindz ...
    logspace(log10(mindz),log10(maxdz),nzc)],5)';
dz = smooth(dz,5)';
tmp = [cumsum(dz)];
ind = find(tmp<Lz_tot);
dze = Lz_tot-max(abs(tmp(ind)));  
% make sure that depth is Lz 
dz = [dz(ind) dze];
zf = -[0 cumsum(dz)]; % this is RF
z = 0.5*(zf(1:end-1)+zf(2:end));
nzc = length(dz);


%horizontal grid parameters
dx_outer = lam/10; 
% set the outer grid resolution to resolve the forcing wave
dy_outer = dx_outer;
dx_inner = 5e3; 
dy_inner = 5e3; 

Lx = 2.*lam + 2*nsponge*dx_outer;
Ly = 2.5*Ld + basin_alpha*Lx_basin; % in m

LF = Lx_basin; %width of the basin

high_res_pad = 800e3;
x0inner = 0.7*Lx - LF/2 - high_res_pad;
x1inner = 0.7*Lx + LF/2 + high_res_pad;

% now make horizontal grids
% x
dx = [dx_outer*ones(1,ceil(x0inner/dx_outer)) ...
      dx_inner*ones(1,ceil((x1inner-x0inner)/dx_inner))...
      dx_outer*ones(1,ceil((Lx-x1inner)/dx_outer))];
% add cells on right to make domain divisible by # of processors
n_add = np(1) - mod(length(dx),np(1));
dx = [dx, dx_outer*ones(1,n_add)];

% smooth dx
dx = smooth(smooth(dx,5),5)';
xg = [0 cumsum(dx)];
xc = 0.5*(xg(2:end)+xg(1:end-1));

% y
dy = [dy_inner*ones(1,floor((0.25*Ly)/dy_inner)) ...
      dy_outer*ones(1,floor((0.75*Ly)/dy_outer))];
% add cells offshore to make domain divisible by # of processors
n_add = np(2) - mod(length(dy),np(2));
dy = [dy, dy_outer*ones(1,n_add)];

% smooth dy
dy = smooth(smooth(dy,5),5)';
yg = [0 cumsum(dy)];
yc = 0.5*(yg(2:end)+yg(1:end-1));

nxc = length(xc); nyc = length(yc);
Lx = max(xc); Ly = max(yc);

%% Topography


[~,mindx] = min(dx);

[~,y_inlet]=min(abs(yc-basin_alpha*Lx_basin));
x_basin0 = mindx+ceil((200e3)/dx_inner);
x_basin1 = mindx+ceil((200e3+Lx_basin)/dx_inner);

Lz_basin = Lz*params.z_ratio;

basin_width = xc(x_basin1)-xc(x_basin0);

hsh = 250; 
ysh = 50e3; 
dysl = 25e3;
prof_y = -hsh -0.5*(Lz_basin-hsh)*(1+tanh((yc(1:y_inlet)-ysh)/dysl));

hsh_x = 250; ysh_x = 50e3; dysl_x = 25e3;
prof = 0*xc(x_basin0:x_basin1);
x = xc(x_basin0:x_basin1);


for i=1:length(x)
   if i<length(x)/2
       prof(i) = -hsh_x ...
           -0.5*(Lz_basin-hsh_x)*(1+tanh((x(i)-xc(x_basin0)-ysh_x)/dysl_x));
   else
       prof(i) = -hsh_x ...
        -0.5*(Lz_basin-hsh_x)*(1-tanh(((x(i)-basin_width-xc(x_basin0))+ysh_x)/dysl_x));
    end
end

prof(1) = 0; % vertical wall
prof(end) = 0;
PROF2 = repmat(prof',1,y_inlet);
PROF2(:,1) = 0;

PROF = PROF2;
for i=1:(x_basin1-x_basin0)+1
PROF(i,:) = -(PROF2(i,:).*prof_y)/Lz_basin;
end


topo = -Lz*ones(nxc,nyc);
topo(:,1:y_inlet)=0;
topo(x_basin0:x_basin1,1:y_inlet)=PROF; 
topo(xc<(xc(x_basin0)+basin_width*params.inlet_frac),y_inlet) = 0;
topo(xc>xc(x_basin1)-basin_width*params.inlet_frac,y_inlet) = 0;

PROF = topo;

%% now stratification 
rho0 = 999.8; g = 9.81; alpha = 2e-4;

r1 = 992; r2 = 995;
r0 = (r1+r2)/2; dr = r2-r1;
N2back = 0.0035^2; %(2*pi/(0.5*3600))^2;
mupyc = 400;
Zpyc = -400;
r = r2 - 0.5*dr*(1+tanh((z-Zpyc)/mupyc)) - z*N2back*r0/g;
n2 = 0.5*(dr/r0)*(g/mupyc)*sech((z-Zpyc)/mupyc).^2 + N2back;

t =  (1-r/rho0)/alpha+5;


%% forcing region

x0mask = dx_outer*nsponge;
x1mask = x0mask + lam;

maskxy = zeros(nxc,nyc);
maskxy(xc>x0mask&xc<x1mask,:) = 1;
maskxy = conv2(maskxy,ones(4,1)/4,'same');
maskxy = conv2(maskxy,ones(4,1)/4,'same');
maskxy = conv2(maskxy,ones(4,1)/4,'same');


mask = repmat(maskxy,[1 1 nzc]);

%% rbcs forcing: modal structures to be read into rbcs_fields_load

vg = exp(-yc/Ld);
UMODE = repmat(vg,[nxc 1 nzc]);


%% output some parameters
params.nxc = nxc;
params.nyc = nyc;
params.nzc = nzc;
params.np = np;


%% strat done, grids done, topo done

% initial fields
T = permute(repmat(t',[1 nxc nyc]),[2 3 1]);

U = 0*T;
V = 0*T;

% boundary sponge region fields

Uzonal = zeros(nyc,nzc,2); % [ny nz nt]
Vzonal = zeros(nyc,nzc,2);
Tzonal = repmat(squeeze(T(1,:,:)),[1 1 2]);

Umerid = zeros(nxc,nzc,2); % [nz nx nt]
Vmerid = zeros(nxc,nzc,2);
Tmerid = repmat(squeeze(T(:,end,:)),[1 1 2]);


%% write some files

fid = fopen(fullfile(rdir,'Uinit.bin'),'w',ieee);
fwrite(fid,U,prec);
fclose(fid);

fid = fopen(fullfile(rdir,'Vinit.bin'),'w',ieee);
fwrite(fid,V,prec);
fclose(fid);

fid = fopen(fullfile(rdir,'Tinit.bin'),'w',ieee);
fwrite(fid,T,prec);
fclose(fid);

fid=fopen(fullfile(rdir,'delZ.bin'),'w',ieee);
fwrite(fid,dz,prec);
fclose(fid);

fid=fopen(fullfile(rdir,'delX.bin'),'w',ieee);
fwrite(fid,dx,prec);
fclose(fid);

fid=fopen(fullfile(rdir,'delY.bin'),'w',ieee);
fwrite(fid,dy,prec);
fclose(fid);

fid=fopen(fullfile(rdir,'Umerid.bin'),'w',ieee);
fwrite(fid,Umerid,prec);
fclose(fid);

fid=fopen(fullfile(rdir,'Vmerid.bin'),'w',ieee);
fwrite(fid,Vmerid,prec);
fclose(fid);

fid=fopen(fullfile(rdir,'Tmerid.bin'),'w',ieee);
fwrite(fid,Tmerid,prec);
fclose(fid);

fid=fopen(fullfile(rdir,'Uzonal.bin'),'w',ieee);
fwrite(fid,Uzonal,prec);
fclose(fid);

fid=fopen(fullfile(rdir,'Vzonal.bin'),'w',ieee);
fwrite(fid,Vzonal,prec);
fclose(fid);

fid=fopen(fullfile(rdir,'Tzonal.bin'),'w',ieee);
fwrite(fid,Tzonal,prec);
fclose(fid);

fid=fopen(fullfile(rdir,'mask.bin'),'w',ieee);
fwrite(fid,mask,prec);
fclose(fid);

fid=fopen(fullfile(rdir,'umode.bin'),'w',ieee);
fwrite(fid,UMODE,prec);
fclose(fid);


% Write topography every time
fid = fopen(fullfile(rdir,'topog.bin'),'w',ieee);
fwrite(fid,PROF,prec);
fclose(fid);

% Save (x,y) coordinates of the basin
save(filename,'y_inlet','x_basin0','x_basin1')