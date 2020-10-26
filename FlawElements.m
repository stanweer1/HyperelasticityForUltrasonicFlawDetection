clear
gs=1;           %element size
w=100;          %model width (number of element columns)
d=100;          %model depth (number of element rows )
dd=4;          %diameter of the defect
cn=w*d*10;      
hl=5;           %number of removable hole layer through thickness
xc=50*gs;			%x-coordinate of the hole center 
yc=50*gs;			%y-coordinate of the hole center 

%element map and coordinates 

for i=1:w
for j=1:d
    els(j,i)=(i-1)*d+j;
    gcx(j,i)=(i-1)*gs+gs/2;
    gcy(j,i)=(j-1)*gs+gs/2;
end
end

%exact hole coordinates

ns=round(pi*dd/1.8/gs);
ro=dd/2*ones(1,ns);
fi=(pi:pi/(ns-1):2*pi);

acx=xc+ro.*cos(fi);
acy=yc+ro.*sin(fi);

%approximate coordinates and numbers of removable hole boundary elements

for k=1:ns
    minimum=100;
    for i=1:w/gs
    for j=1:d/gs
    mindist(j,i)=sqrt((gcx(j,i)-acx(k))^2+(gcy(j,i)-acy(k))^2);
    if(mindist(j,i)<minimum)
        minimum=mindist(j,i);
        r=j;
        c=i;
    end
    end
    end
    aacx(k)=gcx(r,c);
    aacy(k)=gcy(r,c);
    elnr(k)=els(r,c);
    elmid(k)=els(yc,c);
end

for k=1:ns
    saacx(k)=aacx(k);
    saacy(k)=aacy(k)+2*(yc-aacy(k));
    selnr(k)=elnr(k)+2*(elmid(k)-elnr(k));
end

l=1;
xcoord=aacx(1);
minycoord=aacy(1);
maxycoord=saacy(1);
minelset=elnr(1);
maxelset=selnr(1);
for k=1:ns
    if(aacx(k)==xcoord)
         if(aacy(k)<minycoord)
         minelset=elnr(k);
         minycoord=aacy(k);
         end
         if(saacy(k)>maxycoord)
         maxelset=selnr(k);
         maxycoord=saacy(k);
         end 
    else
        elset(l,1)=minelset;
        elset(l,2)=maxelset;
        l=l+1;
        minelset=elnr(k);
        maxelset=selnr(k);
        minycoord=aacy(k);
        maxycoord=saacy(k);
        xcoord=aacx(k);
    end
end

% plot hole

figure
plot(aacx,aacy,'s','MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',3)
hold on
plot(saacx,saacy,'s','MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',3)
axis equal
grid on
axis([0 100 0 100])

% write hole data to file

fid = fopen(['hole_2_geometry_layer_',num2str(hl),'.txt'],'w');

text=char(['*ELSET, ELSET=hole_2_layer_',num2str(hl),', GENERATE']);
fprintf(fid,'%s\n',text); 
 
for i=1:length(elset(:,1))
    fprintf(fid,'%8.0f,%8.0f,%1.0f\n',elset(i,1)+cn*(hl-1),elset(i,2)+cn*(hl-1),1);    
end
fclose(fid)

