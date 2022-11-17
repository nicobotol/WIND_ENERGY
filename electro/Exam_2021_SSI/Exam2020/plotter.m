function plotter(N,x,y1,y2,y1name,y2name,ylabelname,xlabelname,legendloc)
    figure(N)
    plots(N,1) = plot(x,y1,'b-','DisplayName',y1name);
        hold on
    plots(N,2) = plot(x,y2,'r-','DisplayName',y2name);
    yline(0);
        legend([plots(N,1), plots(N,2)],'Location',legendloc)
        ylabel(ylabelname)
        xlabel(xlabelname)
        ylim([min(cat(1,y1,y2))*1.1, max(cat(1,y1,y2))*1.1])
        hold off 
end