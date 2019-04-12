clc
clear all
%close all

fid1 = load('RF.rpt');

rf = -fid1(:,2);
u = -fid1(:,3);

fps =24;


%% Make video
outputVideo = VideoWriter(fullfile(pwd,'DATAMOVIE.avi')); % generate an output file
outputVideo.FrameRate = fps; %set the frame rate
open(outputVideo); % open the output file
 
h=figure('Renderer','zbuffer','Color',[0 0 0],'Position', [100, 100, 1200, 900]);
    
for ii = 1:length(rf) % loop through the sorted list
    plot(u(1:ii),rf(1:ii),'rsq','MarkerSize',5,'MarkerFaceColor',[1,0,0])
    hold on
    plot(u(ii),rf(ii),'wo','MarkerSize',10)
    hold off
    set(findall(gcf,'type','text'),'FontSize',18)
    
    set(gca,'Color',[0 0 0]);
    set(gca,'XColor',[0.85 0.85 0.85]);
    set(gca,'YColor',[0.85 0.85 0.85]);
    set(gca,'FontSize',18)
    set(gca,'FontName','Helvetica')
    box off
    xlim([min(u) 1.02*max(u)])
    ylim([min(rf) 1.02*max(rf)])
    xlabel('$\Delta$(mm)','Color',[0.85,0.85,0.85],'Interpreter','latex')
    ylabel('$P$(N)','Color',[0.85,0.85,0.85],'Interpreter','latex')
%     set(get(gca,'YLabel'),'Rotation',0);
    pbaspect([1 .75 1])
    
    frame=getframe(h);
    writeVideo(outputVideo,frame); % write img to the video file
end
close(h)
 
%% Conclusion
close(outputVideo); % close the video file
disp('DONE')