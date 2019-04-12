function publish_fig(printing, filename, f, figtitle, xlab, ylab,w,fsz)
% this function takes in a figure and formats it for inclusion in a tex
% document. set printing arg to 1 to save figure for publication.

    figure(f) %bring focus to figure fighandle

    %set title and labels (note that latex can be used in the title and
    %label strings
    title(figtitle)
    xlabel(xlab,'Interpreter','Latex')
    ylabel(ylab,'Interpreter','Latex')
    
    h=w*f.Position(4)/f.Position(3);
    psz=1.1*[w,h];
    %sets font sizing for figure and axis
    set(findall(gcf,'type','text'),'FontSize',fsz)
    set(gca,'FontSize',fsz)
    set(gca,'FontName','Helvetica')

    set(gca,'LooseInset',get(gca,'TightInset'))    
    box off

    set(gcf,'InvertHardcopy','off');
    set(gcf, 'color', 'none');

    f.Units='inches';

    f.PaperPositionMode='manual';
    f.PaperUnits='inches';
    f.PaperSize=psz;
    f.PaperPosition=[(psz(1)-w)/2,(psz(2)-h)/2,w,h];
    
    %checks printing boolean and prints to either pdf or eps file in
    %current working directory
    if printing
        print(f,'-dpdf','-cmyk',filename)
    end
end


