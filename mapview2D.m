function [output,h,info]=mapview2D(varargin)
%function [output,h,info]=mapview2D(v,center,v1, v2, settings)
% Program for taking pdf or other 3D density from Jana and plot them
% S. Chrisensen to make arbitrary slices
% Input:
% v: N1xN2xN3-array with density
% center: coordinates (1x3-arrays) for the center of the plot.
% v1, v2: two vectors (1x3-arrays) that span the plane to be plotted. Sofar vectors must
% be orthogonal
% settings: structure standard fields seen below
%
% ------  Default settings  ------------------
% default.boundaries2D=[-0.5 0.5 -0.5 0.5]; %1x4-array [xmin xmax ymin ymax]
% default.units='fractional';   %'fractional' or 'absolute'
% default.unitcell=[1 1 1 90 90 90]; %1x6-array. The unitcell in the 
%                                    %format: [a b c alpha beta gamma]
% default.showaxis=1;           % 0 or 1
% default.contour=0;            % 0: no contours, 1: Draw positive and 
%                               %negative contours with individual style, 
%                               %2: draw all contours with same style
% default.contour_color=[0.99 0.99 0.99]; %1x3 array with RGB-code (MATLAB style) 
% default.surf=1;               % 0: do not make surface plot or 
%                               % 1: make surface plot
% default.colormap='jet';       %any string with name of colormaps 
%                               %recognized by MATLAB
% default.coloraxis='auto';     %'auto' or 1x2-array with desired minimum 
%                               %and maximum values for the colorscale
% default.normalize=1;          %0: do not normalize density. All other 
%                               %values: The density is normalized to the 
%                               %given value
% default.interpolate=0;        %0: do not interpolate, any positive non-
%                               %zero scalar: enlarges the grid by this 
%                               %factor in each direction. NOT TESTED YET
% default.contourvalues='auto'; %'auto': Draws 10 equally spaced contours. 
%                               %1xN-array: Draws N contourlines at the 
%                               %specified values.
% default.contour_style='-';    % Use standard MATLAB notation to specify 
%                               %the linestyle for contour lines.
% default.contour_style_negative='--';  %The linestyle for negative contours
% default.linewidth=1;          %Positive number. Thickness of contour lines.

v=varargin{1};
center=varargin{2};
v1=varargin{3};
v2=varargin{4};



default.boundaries2D=[-0.5 0.5 -0.5 0.5]; %1x4-array [xmin xmax ymin ymax]
default.units='fractional'; %'fractional' or 'absolute'
default.unitcell=[1 1 1 90 90 90]; %1x6-array. The unitcell in the format [a b c alpha beta gamma]
default.showaxis=1; % 0 or 1
default.contour=0;  % 0: no contours, 1: Draw positive and negative contours with individual style, 2: draw all contours with same style
default.contour_color=[0.01 0.01 0.01]; %1x3 array with RGB-code (MATLAB style) 
default.surf=1;         % 0 or 1;
default.colormap='jet';  %any string with name of colormaps recognized by MATLAB
default.coloraxis='auto'; %'auto' or 1x2-array with desired minimum and maximum values for the colorscale
default.normalize=1;      %0: do not normalize density. All other values: The density is normalized to the given value
default.interpolate=0;    %0: do not interpolate, any positive non-zero scalar: enlarges the grid by this factor in each direction. NOT TESTED YET
default.contourvalues='auto'; %'auto': Draws 10 equally spaced contours. 1xN-array: Draws N contourlines at the specified values.
default.contour_style='-';    % Use standard notation to specify the linestyle for contour lines.
default.contour_style_negative='--';  %The linestyle for negative contours
default.linewidth=1; %Positive number. Thickness of contour lines.

if nargin>4
    
%     settings=struct(varargin{5:end})
    settings=varargin{5};
   % if (isfield(settings,'boundaries'));                                 default.boundaries=settings.boundaries;  end
    if (isfield(settings,'boundaries2D'));                               default.boundaries2D=settings.boundaries2D;  end
    if (isfield(settings,'units') && isfield(settings,'unitcell'));      default.units=settings.units;end
    if  isfield(settings,'unitcell');                                    default.unitcell=settings.unitcell;end
    if (isfield(settings,'coloraxis'));                                  default.coloraxis=settings.coloraxis;end
    if (isfield(settings,'colormap'));                                   default.colormap=settings.colormap;end
    if (isfield(settings,'showaxis'));                                   default.showaxis=settings.showaxis;end
    if (isfield(settings,'contour'));                                    default.contour=settings.contour;end
    if (isfield(settings,'contour_color'));                              default.contour_color=settings.contour_color;end
    if (isfield(settings,'contour_color_negative'));                     default.contour_color_negative=settings.contour_color_negative;
    elseif (isfield(settings,'contour_color'));                          default.contour_color_negative=settings.contour_color;
    else                                                                 default.contour_color_negative=default.contour_color;
    end
    if (isfield(settings,'contour_style'));                              default.contour_style=settings.contour_style;end
    if (isfield(settings,'contour_style_negative'));                     default.contour_style_negative=settings.contour_style_negative; end  
    if (isfield(settings,'contourvalues'));                              default.contourvalues=settings.contourvalues;end
        if (isfield(settings,'surf'));                                       default.surf=settings.surf; end
    if (isfield(settings,'normalize'));                                  default.normalize=settings.normalize;end
    if (isfield(settings,'interpolate'));                                default.interpolate=settings.interpolate;end
    if (isfield(settings,'linewidth'));                                  default.linewidth=settings.linewidth;end
  
    
elseif nargin < 4
    error('not enough input arguments')
else
    
end
previous_hold=ishold;
%clean current figure
% clf    %Dette er problematisk hvis figuren indeholder subplots da disse
% ogs� bliver fjernet. Kan ikke huske hvorfor clf oprindeligt er
% inkluderet
settings=default;
info=[];

%check if contour color has been set to white. Of some reason pure white is plotted as black when the figure is saved
if isequal(settings.contour_color, [1 1 1])
    settings.contour_color=[0.99 0.99 0.99];
    disp('Changing the contour color to: [0.99 0.99 0.99]')
end

unitcell=settings.unitcell;


% change graphics renderer
set(gcf, 'Renderer','zbuffer')

if iscolumn(v1); v1=v1'; end
if iscolumn(v2); v2=v2'; end

%check for orthonormal axis
if round2(v1'*v2,0.01)~=0
    error('Axis are not orthogonal. Program not adapted for this')
end

%vectors spanning the plane are normalized
v1=v1/norm(v1);
v2=v2/norm(v2);

gridsize=size(v);

b2d=settings.boundaries2D;
xi=[b2d(1):0.001:b2d(2)];
yi=[b2d(3):0.001:b2d(4)];
[XI, YI]=meshgrid(xi,yi);


xd=XI.*v1(1)+YI.*v2(1)+center(1);
yd=XI.*v1(2)+YI.*v2(2)+center(2);
zd=XI.*v1(3)+YI.*v2(3)+center(3);

%find the size of the box needed
b=[min(min(xd)) max(max(xd)); min(min(yd)) max(max(yd)); min(min(zd)) max(max(zd))];
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
dx=1./gridsize(1); dy=1./gridsize(2); dz=1./gridsize(3);




[X,Y,Z]=ndgrid(x,y,z);
P = [2 1 3];
X = permute(X, P);
Y = permute(Y, P);
Z = permute(Z, P);
vi = permute(vi, P);


%get data for 2D plotting

vint = interp3(X,Y,Z,vi,xd,yd,zd);

mindens=min(vint(:));
maxdens=max(vint(:));

info.min=mindens;
info.max=maxdens;


disp(['The minimum value in the map is: ' num2str(mindens)]);
disp(['The maximum value in the map is: ' num2str(maxdens)]);

if settings.normalize~=0
    vint=vint/maxdens*settings.normalize;
    
end

if isstr(settings.coloraxis)
    if settings.coloraxis=='auto'
        settings.coloraxis=[min(vint(:)), max(vint(:))];
    end
end
%set aspect ratio for plot
a=settings.unitcell(1); b=settings.unitcell(2); c=settings.unitcell(3); alpha=settings.unitcell(4); beta=settings.unitcell(5); gamma=settings.unitcell(6);
cosalphastar=(cosd(beta)*cosd(gamma)-cosd(alpha))/(sind(beta)*sind(gamma));
cstar=a*b*sind(gamma)/(a*b*c*(1-cosd(alpha)^2-cosd(beta)^2-cosd(gamma)^2+2*cosd(alpha)*cosd(beta)*cosd(gamma)).^.5);
%her er noget galt med A(3,3)
A=[a 0 0; b*cosd(gamma) b*sind(gamma) 0 ; c*cosd(beta) -c*sind(beta)*cosalphastar 1/cstar];
G=A^2;
vnorm(1)=(v1*G*v1').^.5; vnorm(2)=(v2*G*v2').^.5;


switch settings.units
    case 'absolute'
        xi=xi*vnorm(1);
        yi=yi*vnorm(2);
    case 'fractional'
        xi=xi;
        yi=yi;
end


hold on
if settings.surf==1
    try
        map=colormap(settings.colormap);
    catch err1
        try 
            map=colormap(settings.colormap(2:end));
            if strcmp(settings.colormap(1),'i')
            map=map(end:-1:1,:);
            end
        catch err2
            error('unknown colormap')
        end
        
        % trying if the given colormap name is part of the pmkmp-package
        % http://www.mathworks.com/matlabcentral/fileexchange/28982-perceptually-improved-colormaps
%         map=pmkmp(64, settings.colormap);
        
    end
    colormap(map)
%     h=surf(xi,yi,vint);

    h=imagesc(xi([1 end]),yi([1 end]),vint);     
set(gca,'ydir','normal')
    shading interp
caxis(settings.coloraxis)

end

if settings.contour==1 % positive and negative contourlevels are drawn
    if ischar(settings.contourvalues)
        if strcmp(settings.contourvalues,'auto')
            vmax=max(vint(:));
            vmin=min(vint(:));
            cstep=max(vint(:))/10;
            if vmin<0
                pos_contour=[0:cstep:vmax]*0.999;
                neg_contour=[0:cstep:-vmin]*0.999;
                settings.contourvalues=[-neg_contour(end:-1:2) pos_contour]-0.001;
                
            else
%                 settings.contourvalues=[vmin:cstep:vmax]-0.0001;
                    settings.contourvalues=[vmin:cstep:vmax]*0.999; %NOTE reducing the extrema-contours slightly, since nothing is drawn it is only a single point
            end
            
        else
            error('settings.contourvalues is a string but not "auto" as expected')
        end
    end
    if length(settings.contourvalues)==1
        settings.contourvalues=[settings.contourvalues settings.contourvalues];
    end
    settings.contourvalues=sort(settings.contourvalues);
    neg_ind=settings.contourvalues<0;
    pos_ind=settings.contourvalues>0;
    
    
    
    if settings.interpolate~=0
        disp('WARNING: You are using an untested feature')
        vint=interp2(vint,settings.interpolate);
        xi=interp2(xi,settings.interpolate);
        yi=interp2(yi,settings.interpolate);
    end
    
    if sum(neg_ind)~=0
    %plot negative contours
            negative_contours=settings.contourvalues(neg_ind);
        if sum(neg_ind)==1
            negative_contours=repmat(negative_contours,2,1);
        end
    
        [C,h]=contour(xi,yi,vint,negative_contours, '--','LineWidth',settings.linewidth,'LineColor',settings.contour_color);
        %          [C,h]=contourspline(xi,yi,vint, settings.contourvalues(neg_ind));
        % apparently the structure of h is different when using
        % contourspline thus the Zoffset part below will need to be adapted
        % to use contourspline
        set(h, 'linewidth', settings.linewidth)
        set(h, 'EdgeColor', settings.contour_color_negative)
        set(h, 'LineStyle', settings.contour_style_negative)
        
        hh = get(h,'Children');    %# get handles to patch objects
        Zoffset=max(vint(:));
        for i=1:numel(hh)
            zdata = Zoffset*ones(size( get(hh(i),'XData') ));
            set(hh(i), 'ZData',zdata)
        end
        hold on
        %make small white dots in the corner of the plot to avoid the size
        %of the plot is reduced when saved without axis
        c=[min(xi(:)) max(xi(:)) min(yi(:)) max(yi(:))];
        plot([c(1) c(1) c(2) c(2) c(1)], [c(3) c(4) c(4) c(3) c(3)],'.','Color',[1 1 1 ],'MarkerSize',0.1)
    end
    if sum(pos_ind)~=0
        %plot negative contours
        positive_contours=settings.contourvalues(pos_ind);
        if sum(pos_ind)==1
            positive_contours=repmat(positive_contours,2,1);
        end

        [C,h]=contour(xi,yi,vint,positive_contours, '-','LineWidth',settings.linewidth,'LineColor',settings.contour_color);
        
        set(h, 'linewidth', settings.linewidth)
        set(h, 'EdgeColor', settings.contour_color)
        set(h, 'LineStyle', settings.contour_style)
        hh = get(h,'Children');    %# get handles to patch objects
        Zoffset=max(vint(:));
        for i=1:numel(hh)
            zdata = Zoffset*ones(size( get(hh(i),'XData') ));
            set(hh(i), 'ZData',zdata)
        end
        %         set(hh, {'ZData'}, cellfun(@(x) Zoffset*ones(size(x)),get(hh,{'XData'}), 'UniformOutput',false))
        %should be equivalent to the above for-loop
        
        hold on
        c=[min(xi(:)) max(xi(:)) min(yi(:)) max(yi(:))];
        plot([c(1) c(1) c(2) c(2) c(1)], [c(3) c(4) c(4) c(3) c(3)],'.','Color',[1 1 1 ],'MarkerSize',0.1)
        
    end
    output=C;
elseif settings.contour==2  %all contour levels are drawn alike.
    
    if ischar(settings.contourvalues)
        if strcmp(settings.contourvalues,'auto')
            settings.contourvalues=linspace(min(vint(:)),max(vint(:))*0.999,10);
        else
            error('settings.contourvalues is a string but not "auto" as expected')
        end
    end
    
    if length(settings.contourvalues)==1
        settings.contourvalues=[settings.contourvalues settings.contourvalues];
    end
    
    settings.contourvalues=sort(settings.contourvalues);
    [C,h]=contour(xi,yi,vint, settings.contourvalues, '-','LineWidth',settings.linewidth,'LineColor',settings.contour_color);
    set(h, 'linewidth', settings.linewidth)
    hh = get(h,'Children');    %# get handles to patch objects
    Zoffset=max(vint(:));
    for i=1:numel(hh)
        zdata = Zoffset*ones(size( get(hh(i),'XData') ));
        set(hh(i), 'ZData',zdata)
    end
    %make small white dots in the corner of the plot to avoid the size
    %of the plot is reduced when saved without axis
    
    hold on
    c=[min(xi(:)) max(xi(:)) min(yi(:)) max(yi(:))];
    plot([c(1) c(1) c(2) c(2) c(1)], [c(3) c(4) c(4) c(3) c(3)],'.','Color',[1 1 1 ],'MarkerSize',0.1)
    output=C;
else
    output=1;
end

if settings.surf==0 && settings.contour==0
    error('You must give at least one plotting method: surf or contour')
end


switch settings.units
    case 'absolute'
        axis([settings.boundaries2D(1:2)*vnorm(1) settings.boundaries2D(3:4)*vnorm(2)])
    case 'fractional'
        axis(settings.boundaries2D)
end
view([0 90])


if settings.showaxis==1
    xlabel([ '[' num2str(v1(1)) ' , ' num2str(v1(2)) ' , ' num2str(v1(3)) ']' ],'Fontsize',14);
    ylabel([ '[' num2str(v2(1)) ' , ' num2str(v2(2)) ' , ' num2str(v2(3)) ']' ],'Fontsize',14);
    axis on
    set(gca,'Fontsize',14)
else
    axis off
end


%reset old hold status
if previous_hold==0
    hold off
elseif previous_hold==1
    hold on
end

end

function z = round2(x,y)
%ROUND2 rounds number to nearest multiple of arbitrary precision.
%   Z = ROUND2(X,Y) rounds X to nearest multiple of Y.
%
%Example 1: round PI to 2 decimal places
%   >> round2(pi,0.01)
%   ans =
%         3.14

%% defensive programming
error(nargchk(2,2,nargin))
error(nargoutchk(0,1,nargout))
if numel(y)>1
  error('Y must be scalar')
end

%%
z = round(x/y)*y;
end 