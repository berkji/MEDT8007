function hS = show_xdc_geir(Th,Mval)
% Mval = 0 -> delay as color
% Mval = 1 -> physical number mod Mval as color
% Mval = 2 -> random number to each physical element as color
% Mval = 3 -> Apodization as color

if nargin < 2, 
    Mval = 10000;
end

datar = xdc_get(Th,'rect'); % try to read rectangles
datat = xdc_get(Th,'tri'); % try to read triangles

if ~isempty(datar), % got rectangles
    disp('Read rectangular data for plotting....');
    [N,M] = size(datar);
        X = [datar(11,:); datar(14,:); datar(17,:); datar(20,:)]*1000;
        Y = [datar(12,:); datar(15,:); datar(18,:); datar(21,:)]*1000;
        Z = [datar(13,:); datar(16,:); datar(19,:); datar(22,:)]*1000;
        if Mval == 0,
            disp('Plots aperture with delay...');
            C = 1e9*[1; 1; 1; 1]*datar(23,:); % delay value in index 23 (converted to [ns])
%            Z = -C;    
        elseif Mval == 2,
            disp('Plots aperture with radom element numbers...');
            MySign = 1;
            if 0
                C = [1; 1; 1; 1]*datar(1,:);
                Ctmp = zeros(size(C));
                for no = 1:max(datar(1,:)),
                    Ctmp(C == no) = MySign*mod(C(C == no),-Mval);   
                    MySign = -MySign;
                end
                C = Ctmp;
            else
                C = [1; 1; 1; 1]*datar(1,:);
                Ctmp = zeros(size(C));
                nPhMax = max(C(:));
                noPhRand = rand(size(C));
                for no = 1:nPhMax,
                    Ctmp(C == no) = round(-Mval*noPhRand(no));
                end
                C = Ctmp;
            end
        elseif Mval == 1,
            disp('Plots aperture with physical element number...');
            C = [1; 1; 1; 1]*mod(datar(1,:),Mval); % physical element number
        elseif Mval == 3,
            disp('Plots aperture with apodization...');
            C = 1e9*[1; 1; 1; 1]*datar(5,:); % delay value in index 23 (converted to [ns])
        end    
elseif ~isempty(datat) % got triangles?
    [N,M] = size(datat);
    disp('Read triangular data for plotting....');    
    X = [datat(7,:); datat(10,:); datat(13,:)]*1000;
    Y = [datat(8,:); datat(11,:); datat(14,:)]*1000;
    Z = [datat(9,:); datat(12,:); datat(15,:)]*1000;
    C = [1; 1; 1]*mod(datat(1,:),Mval); % physical element number
end

MyAxis = [min(X(:)) max(X(:)) min(Y(:)) max(Y(:))];

if 1
    hS = fill3(X,Y,Z,C);
%    rotate(hS,[0 0 1],-45);
else % draw one at the time...
    MyAxis = [min(X(:)) max(X(:)) min(Y(:)) max(Y(:))];
    clf; hold on;
    for i = 0:max(C(:));
        iThis = (C(1,:) == i)
        hS = fill3(X(:,iThis),Y(:,iThis),Z(:,iThis),C(:,iThis)); drawnow; axis(MyAxis); caxis([0 max(C(:))]);
        pause
    end
    hold off;
end
    
if Mval == 0, colorbar; end


% hc = colorbar
colormap(jet);
view(2)
xlabel('x [mm]')
ylabel('y [mm]')
zlabel('z [mm]')
grid on;
%axis('image')
axis(MyAxis)
hold off;
drawnow;
figure(gcf);