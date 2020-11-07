% Saves all open figures with their title as filename. Figures stored in 
% a "figures" directory and a directory named with the day and time the
% image was stored.
%
% function [] = saveallfigs(figDir, format, noName)
%
% inputs:
%           figDir      The directory in which to save the figures
%                       use input [] to create a folder with name equal to
%                       the date and time the figures was saved
%           format      Defines which format the figures are stored in
%                       0 = .fig, 1 = .eps(res 300), 2 = .png (res140)
%                       3 = .swf, 4 = .jpg (res150), -1 = save as fig, png
%                       and flash
%           noName      0 -> file name = figure title
%                       1 -> file name = figure number
%                       2 -> file name = saved file name
%
% Example: 
%           saveallfigs([], 
%
% Created by Tore Bjåstad 09.02.2005


function [] = saveallfigs(figDir, format, noName)
if (nargin==1)
    format = 0;
    noName = 0;
end

if (nargin==0)
    clk = clock;
    save_dir = char;
    for ii = 1:6
        temp = clk(ii);
        if(temp<10)
            save_dir = [save_dir '0' num2str(round(temp)) '_'];
        else
            save_dir = [save_dir num2str(round(temp)) '_'];
        end
    end
    figDir = ['figures/' save_dir(1:end-1)];
    mkdir(figDir);
    format = 0;
    noName = 0;
end

if (nargin==2&isempty(figDir))    
    clk = clock;
    save_dir = char;
    for ii = 1:6
        temp = clk(ii);
        if(temp<10)
            save_dir = [save_dir '0' num2str(round(temp)) '_'];
        else
            save_dir = [save_dir num2str(round(temp)) '_'];
        end
    end
    figDir = ['figures/' save_dir(1:end-1)];
    mkdir(figDir);
    noName = 0;
end

if (nargin==2&~isempty(figDir))        
    noName = 0;
end

if (nargin==3&isempty(figDir))    
    clk = clock;
    save_dir = char;
    for ii = 1:6
        temp = clk(ii);
        if(temp<10)
            save_dir = [save_dir '0' num2str(round(temp)) '_'];
        else
            save_dir = [save_dir num2str(round(temp)) '_'];
        end
    end
    figDir = ['figures/' save_dir(1:end-1)];
    mkdir(figDir);
end


if format == -1
    saveallfigs(figDir, 0)
    saveallfigs(figDir, 2)
    saveallfigs(figDir, 3)
    return
end

if (exist(figDir)~=7)
    createdir = input(['Figure directory non existent.' char(10)...
        'Create directory: "' figDir '"? (y/n): '], 's');
    if (createdir=='y')
        mkdir(figDir);
    else
        return;
    end
end

if (sum(format==3)==1) | (sum(format==5)==1)
    filenameFlash = 'flash_slideshow';
end

h_figs = get(0, 'children');
h_figs = sort(h_figs);

groupsize= 1;
for jj=1:groupsize
for ii = jj:groupsize:length(h_figs)
    figure(h_figs(ii));
    cax = get(gcf, 'children');
    fname = get(get(gca, 'Title'), 'String');   
    if (length(fname)==0)
        fname = ['Figure_' num2str(h_figs(ii))];
    end
    if (size(fname, 1)~=1) %More than one line in title
        if iscell(fname)
            temp = fname;
            fname = '';
            for jj=1:size(temp, 1)
                fname = [fname temp{jj}];
            end
        else
            [n, m] = size(fname);
            temp_jj = '';
            for jj = 1:n
                kk = m;
                temp_kk = fname(jj, :);
                while (fname(jj, kk-1:kk)=='  ')
                    temp_kk = temp_kk(1:end-1);
                    kk = kk-1;
                end
                temp_jj = [temp_jj temp_kk];
            end
            fname = temp_jj;
        end
    end
    fname(regexp(fname, '\s')) = '_';
    fname(regexp(fname, '[.]')) = ',';
    fname(regexp(fname, '\n')) = '_';
    fname(regexp(fname, '\')) = '_';
    fname(regexp(fname, '/')) = '_';
    fname(regexp(fname, ':')) = '=';
    if (length(fname)>149)|(noName==1)
        fname = ['figure_' num2str(ii)];
    elseif noName == 2
        fig_filename = get(gcf, 'filename');
        fname_start = find(fig_filename=='\');
        fname_start = fname_start(end);
        fname = fig_filename(fname_start+1:end-4);
    end
    set(gcf, 'PaperPositionMode', 'auto')
    for kk = 1:length(format)
        switch format(kk)
            case 0
                saveas(gcf,[figDir,'/',fname(1:end), '.fig']);
            case 1
                print('-depsc','-r300',[figDir,'/',fname, '.eps']);
            case 2
                %             print('-dpng','-r140',[figDir,'/',fname, '.png']);
%                 print('-dpng','-r200',[figDir,'/',fname, '.png']);
                print('-dpng','-r450',[figDir,'/',fname, '.png']);
            case 3
                MakeFlashMovie(filenameFlash,'addframe',h_figs(ii),'-painters');
            case 4
                print('-djpeg90','-r150',[figDir,'/',fname, '.jpg']);
            case 5
                MakeFlashMovie(filenameFlash,'addframe',h_figs(ii),'-opengl');
            case 6
                print('-painters', '-dpdf','-r300', [figDir,'/',fname, '.pdf']);
        end
    end
end
end

if (sum(format==3)==1) | (sum(format==5)==1)
    copyfile([filenameFlash '.pdf'], [figDir '/' filenameFlash '.pdf']);
    MakeFlashMovie(filenameFlash,'finish', h_figs(ii), 8);
    movefile([filenameFlash '.swf'], [figDir '/' filenameFlash '.swf']);
end


disp('Done!!')
