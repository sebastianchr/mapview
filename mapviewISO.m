function [fv11,info]=mapviewISO(v,settings)
% function fv1=mapviewISO(v,settings)
% Program for plotting isosurfaces of a 3D density file as 
% Author: Sebastian Chrisensen (Aarhus University)
%
%Description:
%Compulsary input: density array, v. The viewing direction defined by the 
%vector, r0. 
%additional input: Additional settings can be specified in settings
%structure. The various settings are described below.
%Settings:
% boundaries: ([0 1 0 1 0 1]) Define the boudaries of the box to be plotted. [x_min x_max y_min y_max z_min z_max] 
% isosurface: (0.5) Specify level of isosurface. >0: absolute units. <0: relative
% to maximal density.
% showaxis: (1) 0 or 1 to hide or show axis
% units: 'fractional' (default) or 'absolute'. Absolute units require size of
% unitcell to be defined.
% unitcell: ([1 1 1]) Currently only works correctly for orthogonal cells
% color: Color of isosurface. letter input e.g. 'r' 'g' 'b' etc
% grid: hide or show grid: 0 or 1 
% interpolation_faactor: (1) 
% point_of_view: ([1 1 1]) Define the initial viewing direction.



%------------------------
%Default settings
default.boundaries=[0 1 0 1 0 1];
default.isosurface=0.5;
default.normalize=0;
default.showaxis=1;
default.units='fractional';
default.unitcell=[1 1 1];
default.color='r';
default.grid=0;
default.interpolation_factor=1;
default.point_of_view=[1 1 1];
if nargin==2
    if isa(settings,'struct')
        
        if 0==(isfield(settings,'boundaries'));                             settings.boundaries=default.boundaries;  end
        if 0==(isfield(settings,'units') && isfield(settings,'unitcell'));  settings.units=default.units;end
        if 0==(isfield(settings,'unitcell'));                               settings.unitcell=default.unitcell;end
        if 0==(isfield(settings,'isosurface'));                             settings.isosurface=default.isosurface;end
        if 0==(isfield(settings,'showaxis'));                               settings.showaxis=default.showaxis;end
        if 0==(isfield(settings,'normalize'));                              settings.normalize=default.normalize;end
        if 0==(isfield(settings,'color'));                                  settings.color=default.color;end
        if 0==(isfield(settings,'grid'));                                   settings.grid=default.grid;end
        if 0==(isfield(settings,'interpolation_factor'));                   settings.interpolation_factor=default.interpolation_factor;end
        if 0==(isfield(settings,'point_of_view'));                                    settings.point_of_view=default.point_of_view;end
    
    else
        error('settings not a structure')
    end
elseif nargin == 1
    
    settings=default;
else
    error('not enough input arguments')
end
info=[];

verbose=0;

%##### Define borders  ################
gridsize=size(v);
b=reshape(settings.boundaries,2,3)';
r={floor(b(1,1)*gridsize(1)):ceil(b(1,2)*gridsize(1))+1;...
   floor(b(2,1)*gridsize(2)):ceil(b(2,2)*gridsize(2))+1;...
   floor(b(3,1)*gridsize(3)):ceil(b(3,2)*gridsize(3))+1};

ri={mod(r{1}-1,gridsize(1))+1; mod(r{2}-1,gridsize(2))+1; mod(r{3}-1,gridsize(3))+1};

vi=v(ri{1},ri{2},ri{3});
  

mindens=min(vi(:));
maxdens=max(vi(:));

info.min=mindens;
info.max=maxdens;

if settings.normalize
    vi=vi/maxdens;
end

switch settings.units
    case 'absolute'
        x=(r{1}-1)./gridsize(1).*settings.unitcell(1);
        y=(r{2}-1)./gridsize(2).*settings.unitcell(2);
        z=(r{3}-1)./gridsize(3).*settings.unitcell(3);
        
        dx=1./gridsize(1).*settings.unitcell(1); dy=1./gridsize(2).*settings.unitcell(2); dz=1./gridsize(3).*settings.unitcell(3);
    case 'fractional'
        x=(r{1}-1)./gridsize(1);
        y=(r{2}-1)./gridsize(2);
        z=(r{3}-1)./gridsize(3);
        dx=1./gridsize(1); dy=1./gridsize(2); dz=1./gridsize(3);
end

r0=settings.point_of_view;
if norm(r0(1:2))~=0
    AZ=acosd(-r0(2)./norm(r0(1:2)));
    EL=acosd((r0(1).^2.+r0(2).^2)./(norm(r0)*norm(r0(1:2))));
else
    AZ=0;
    EL=90;
    if verbose; disp([AZ EL] ); end
end


[X,Y,Z]=ndgrid(x,y,z);
   P = [2 1 3];
   X = permute(X, P);
   Y = permute(Y, P);
   Z = permute(Z, P);
   vi = permute(vi, P);
   
   if settings.interpolation_factor~=1
        dx=x(2)-x(1); dy=y(2)-y(1); dz=z(2)-z(1);
        xi=min(x):dx/settings.interpolation_factor:max(x);
        yi=min(y):dy/settings.interpolation_factor:max(y);
        zi=min(z):dz/settings.interpolation_factor:max(z);       
        [Xi,Yi,Zi]=meshgrid(xi,yi,zi); 
        vii=interp3(X,Y,Z,vi,Xi,Yi,Zi,'cubic');
        X=Xi; Y=Yi; Z=Zi; vi=vii;
   end
   
   
fv11=isosurface(X,Y,Z,vi,settings.isosurface(1));
fv12=patch(fv11,'FaceColor',settings.color(1),'EdgeColor','none');
isonormals(X,Y,Z,vi,fv12);

if numel(settings.isosurface)>1
    hold on
fv21=isosurface(X,Y,Z,vi,settings.isosurface(2));
fv22=patch(fv21,'FaceColor',settings.color(2),'EdgeColor','none');
isonormals(X,Y,Z,vi,fv22);
end




if settings.showaxis
    box on
else
    box off
    axis off
end


if settings.grid
    grid on
else
    grid off
end

daspect([settings.unitcell(1)^-1 settings.unitcell(2)^-1 settings.unitcell(3)^-1])
axis(settings.boundaries)
view(AZ,EL);

camlight;
lighting gouraud
%camproj perspective
end

