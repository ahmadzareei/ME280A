function meshImages
clear all;
close all;

NEth = [10,20];
NEr = [5,10,20];
ro = 0.25;
ri = 0.1;


for i = 1:length(NEth)
    for j=1:length(NEr)
        [p,e] =  mesh (ri, ro, NEr(j), NEth(i));
        clf; 
        patch('Faces',e,'Vertices',p,'FaceColor','w','Linewidth',2)
        title (sprintf('%d Elements in radial direction \n %d Elements in azimuthal direction',NEr(j), NEth(i)),'FontSize',14);
        set(gca, 'FontSize',16);
        axis square
        axis([-.25,.25,-.25,.25]);

        nam = sprintf('Nr%dNEth%d.eps',NEr(j), NEth(i));
        nam
        print(nam,'-depsc')
    end
end
