function drawXtal(coords,Mm,offset)
format compact

if nargin < 3
    offset = [0 0 0];
end

% define variable parameters:
bgColor = 1; % color of figure background
r0 = 0.6;    % base radius of spheres

h = figure('Color','black','ToolBar','none','MenuBar','none');
if bgColor == 1
    set(h,'Color','white')
end
set(gca,'Position',[0 0 1 1],'visible','off','DataAspectRatio',[1 1 1]);
cameratoolbar(h);


Natom = size(coords,1);
cla
light;
light('Position',[-1 -1 -2]);
[x,y,z] = sphere(20);

if Natom > 500
    [x,y,z] = sphere(12);
end
if Natom > 1000
    [x,y,z] = sphere(6);
end
colormap('default');
cmap = colormap();
Zmin = min(coords(:,1));
Zmax = max(coords(:,1));
if (Zmax == Zmin), Zmax = Zmin+1; end

for j=1:size(coords,1)
    if (1)
        color = cmap(1+round(63*(coords(j,1)-Zmin)/(Zmax-Zmin)),:); 
        r = r0+0.5*(coords(j,1)-Zmin)/(Zmax-Zmin);
    else
        switch(coords(j,1))
            case  13, color = [0.8 0 0]; r = 0.6; % Al
            case  22, color = [1.0 0 0]; r = 0.6; % Ti
            case  6, color = [0.7 0.7 0.7]; r = 1.0;  % C
            case  8, color = [0.3 0.3 1.0]; r = 0.7;  % O
            case  38, color = [1.0 0.3 1.0]; r = 0.8; % Sr
            case  14, color = [0.8 0.8 1.0]; r = 0.8; % Si
            case  7, color = [0.8 0.8 0.8]; r = 0.8; % N
            otherwise, color = cmap(1+round(63*(coords(j,1)-Zmin)/(Zmax-Zmin)),:); r = 0.6+0.5*(coords(j,1)-Zmin)/(Zmax-Zmin);
        end
    end
    vec = coords(j,2:4);
    surface('XData',vec(1) + r*x,'YData',vec(2) + r*y,...
        'ZData',vec(3) + r*z,'FaceColor',color,...
        'EdgeColor','none','FaceLighting','gouraud')
end


% draw the unit cell boundaries:
if bgColor == 1
    lineColor = [0 0 0];
else
    lineColor = [1 1 1];
end

if nargin > 1
    % l = [0 0 0;1 0 0;1 1 0;0 1 0;0 0 0;0 0 1;1 0 1;1 1 1;0 1 1;0 0 1]*(Mm);
    l = (Mm*[0 0 0;1 0 0;1 1 0;0 1 0;0 0 0;0 0 1;1 0 1;1 1 1;0 1 1;0 0 1].').';
    l = l+repmat(offset,size(l,1),1);
    line(l(:,1),l(:,2),l(:,3),'Color',lineColor);        

    l = (Mm*[0 1 0;0 1 1].').';
    l = l+repmat(offset,size(l,1),1);
    line(l(:,1),l(:,2),l(:,3),'Color',lineColor);        

    l = (Mm*[1 1 0;1 1 1].').';
    l = l+repmat(offset,size(l,1),1);
    line(l(:,1),l(:,2),l(:,3),'Color',lineColor);        

    l = (Mm*[1 0 0;1 0 1].').';
    l = l+repmat(offset,size(l,1),1);
    line(l(:,1),l(:,2),l(:,3),'Color',lineColor);        

end
