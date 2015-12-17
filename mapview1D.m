function [vint, XI,h]=mapview1D(v,center,v1, settings)
%function [vint, XI]=mapview1D(v,center,v1, settings)
%----------------------
% Input:
% v: 3D matrix with density
% center:  coordinates of the center position of the plot
% v1: vector to define the direction to plot
% settings: structure standard fields seen below
% default.boundaries1D=[-0.5 0.5];
% default.showaxis=1;
% default.normalize=1;
% default.color='r';
% default.units='fractional';
% default.unitcell=[1 1 1 90 90 90];
% default.Xlabel='X'; default.Ylabel='Y'; default.Zlabel='Z';



%------------------------
default.boundaries1D=[-0.5 0.5];
default.showaxis=1;
default.normalize=1;
default.color='r';
default.linestyle='-';
default.units='fractional';
default.unitcell=[1 1 1 90 90 90];
default.Xlabel='X'; default.Ylabel='Y'; default.Zlabel='Z';

if nargin>3
    if (isfield(settings,'boundaries'));                                 default.boundaries=settings.boundaries;  end
    if (isfield(settings,'boundaries1D'));                               default.boundaries1D=settings.boundaries1D;  end
    if (isfield(settings,'units') && isfield(settings,'unitcell'));      default.units=settings.units;end
    if (isfield(settings,'color'));                                      default.color=settings.color;end
    if (isfield(settings,'style'));                                      default.linestyle=settings.linestyle;end
    if (isfield(settings,'showaxis'));                                   default.showaxis=settings.showaxis;end
    if (isfield(settings,'normalize'));                                  default.normalize=settings.normalize;end
elseif nargin < 3
    error('not enough input arguments')
else
    
end

settings=default;
gridsize=size(v);

if iscolumn(v1); v1=v1'; end

%vectors spanning the plane are normalized
v1=v1/norm(v1);

XI=settings.boundaries1D(1):0.001:settings.boundaries1D(2);


xd=XI.*v1(1)+center(1);
yd=XI.*v1(2)+center(2);
zd=XI.*v1(3)+center(3);

%find the size of the box needed
b=[min(xd) max(xd); min(yd) max(yd); min(zd) max(zd)];
r={floor(b(1,1)*gridsize(1)):ceil(b(1,2)*gridsize(1))+1;...
    floor(b(2,1)*gridsize(2)):ceil(b(2,2)*gridsize(2))+1;...
    floor(b(3,1)*gridsize(3)):ceil(b(3,2)*gridsize(3))+1};
for k=1:3
    if length(r{k})==1
        r{k}=[r{k}-1 r{k} r{k}+1];
    end
end
%take account of cyclic boundaries
ri={mod(r{1}-1,gridsize(1))+1; mod(r{2}-1,gridsize(2))+1; mod(r{3}-1,gridsize(3))+1};

vi=v(ri{1},ri{2},ri{3});
%coordinates for all point in the new box

x=(r{1}-1)./gridsize(1);
y=(r{2}-1)./gridsize(2);
z=(r{3}-1)./gridsize(3);
% dx=1./gridsize(1); dy=1./gridsize(2); dz=1./gridsize(3);

[X,Y,Z]=ndgrid(x,y,z);
P = [2 1 3];
X = permute(X, P);
Y = permute(Y, P);
Z = permute(Z, P);
vi = permute(vi, P);


%get data for 1D plotting

vint = interp3(X,Y,Z,vi,xd,yd,zd);
if strcmp(settings.units,'absolute')
    %get aspect ratio for plot
    a=settings.unitcell(1); b=settings.unitcell(2); c=settings.unitcell(3); alpha=settings.unitcell(4); beta=settings.unitcell(5); gamma=settings.unitcell(6);
    cosalphastar=(cosd(beta)*cosd(gamma)-cosd(alpha))/(sind(beta)*sind(gamma));
    cstar=a*b*sind(gamma)/(a*b*c*(1-cosd(alpha)^2-cosd(beta)^2-cosd(gamma)^2+2*cosd(alpha)*cosd(beta)*cosd(gamma)).^.5);
    A=[a 0 0; b*cosd(gamma) b*sind(gamma) 0 ; c*cosd(beta) -c*sind(beta)*cosalphastar 1/cstar];
    G=A^2;
    vnorm(1)=(v1*G*v1').^.5;
end
if settings.normalize
    maxdens=max(max(vint));
    vint=vint/maxdens;
end


switch settings.units
    case 'absolute'
        XI=XI*vnorm(1);
    case 'fractional'
        XI=XI;
end

hold on
h=plot(XI,vint,'Color',settings.color,'linestyle',settings.linestyle);
switch settings.units
    case 'absolute'
        xlim([settings.boundaries1D(1:2)*vnorm(1)])
    case 'fractional'
        xlim(settings.boundaries1D)
end

end