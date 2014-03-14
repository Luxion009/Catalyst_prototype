%soild用翼出力
%出力関数
function p_end=print_sld(zpos,ypos,chord,alpha,dihedral,sparpos,ed_foil,dwingnum,ribnum,pj_name,now_work)

mkdir("design_data")
cd("design_data")
mkdir(pj_name)
cd(pj_name)
mkdir("sld_out")
cd("sld_out")

zpos=zpos*1000;
ypos=ypos*1000;
dihedral=-dihedral;

fs=0;
%-----翼の前縁を探索
do	
	fs=fs+1;
until(ed_foil(fs,1)<=min(ed_foil(:,1)))

n_xdiv_l=80;
n_xdiv_t=60;

x_foil(1:n_xdiv_t+n_xdiv_l-1,1)=1./(n_xdiv_t+n_xdiv_l-1)^2.*([n_xdiv_t+n_xdiv_l-1:-1:1]').^2;	
x_foil(n_xdiv_t+n_xdiv_l,:)=0;
x_foil(n_xdiv_t+n_xdiv_l+1:2*n_xdiv_t+2*n_xdiv_l-1,1)=1./(n_xdiv_t+n_xdiv_l-1)^2.*([1:n_xdiv_t+n_xdiv_l-1]').^2;

y_foil(1:n_xdiv_l+n_xdiv_t,1)=interp1(ed_foil(1:fs,1),ed_foil(1:fs,2),x_foil(1:n_xdiv_l+n_xdiv_t,1),"extrap");
y_foil(n_xdiv_l+n_xdiv_t:2*n_xdiv_t+2*n_xdiv_l-1,1)=interp1(ed_foil(fs:rows(ed_foil),1),ed_foil(fs:rows(ed_foil),2),x_foil(n_xdiv_l+n_xdiv_t:2*n_xdiv_t+2*n_xdiv_l-1,1),"extrap");
y_foil(n_xdiv_l+n_xdiv_t,1)=0;


profx=x_foil;
profy=y_foil;
profz(1:rows(x_foil),1)=0;

prof=[profx profy profz];

%キャンバーライン
cambpoint=[0.01:0.01:0.99];
for i=1:length(cambpoint);
	cambline(i,2)=spline(ed_foil(fs-2:-1:2,1),ed_foil(fs-2:-1:2,2),cambpoint(1,i))-spline(ed_foil(fs+2:1:rows(ed_foil)-1,1),ed_foil(fs+2:1:rows(ed_foil)-1,2),cambpoint(1,i));
	cambline(i,1)=cambpoint(1,i);
	bufu=spline(ed_foil(fs-2:-1:2,1),ed_foil(fs-2:-1:2,2),cambpoint(1,i));
	bufl=spline(ed_foil(fs+2:1:rows(ed_foil)-1,1),ed_foil(fs+2:1:rows(ed_foil)-1,2),cambpoint(1,i));

	cambline(i,4)=(bufu+bufl)/2;
endfor

cambline(1:rows(cambline),3)=zpos;

spar_camb=spline(cambline(:,1),cambline(:,2),sparpos);
spar_camb_h=spar_camb/2;

sparmidp(1,2)=min(spline(ed_foil(fs-2:-1:2,1),ed_foil(fs-2:-1:2,2),sparpos),spline(ed_foil(fs+2:1:rows(ed_foil)-1,1),ed_foil(fs+2:1:rows(ed_foil)-1,2),sparpos))+spar_camb_h;
sparmidp(1,1)=sparpos;

%スパー中心に変換

profx_ref(:,1)=profx(:,1).-sparmidp(1,1);
profy_ref(:,1)=profy(:,1).-sparmidp(1,2);
profz_ref=profz;

%コード長考慮
profx_chord=profx_ref.*chord;
profy_chord=profy_ref.*chord;
profz_chord=profz_ref;

%回転変換
for i=1:rows(profx_chord)
	profx_rot(i,1)=profx_chord(i,1)*cos(-alpha/180*pi)-profy_chord(i,1)*sin(-alpha/180*pi);
	profy_rot(i,1)=profx_chord(i,1)*sin(-alpha/180*pi)+profy_chord(i,1)*cos(-alpha/180*pi);
endfor

profz_rot=profz_ref;

%上半角用回転変換
for i=1:rows(profx_chord)
	profz_rot(i,1)=profz_rot(i,1)*cos(-dihedral)-profy_rot(i,1)*sin(-dihedral);
	profy_rot(i,1)=profz_rot(i,1)*sin(-dihedral)+profy_rot(i,1)*cos(-dihedral);
endfor

profx_rot=profx_rot;

profy_rot=profy_rot+ypos;
profz_rot=profz_rot+zpos;


%disp("generating ps code")

openfile=strcat(dwingnum,'_',dec2base(ribnum,10),'.sldcrv');

fpw=fopen(openfile,'wt');
	
for i=1:rows(profx_chord)
	fprintf(fpw,"%0.2f,%0.2f,%0.2f\n",profz_rot(i,1),profy_rot(i,1),-profx_rot(i,1));
endfor

fclose(fpw);

% disp("output complete")
p_end=1;
cd(now_work)

endfunction