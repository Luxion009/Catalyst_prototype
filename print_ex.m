%出力関数
%コード長、取り付け角、桁位置、桁径、肉抜きデータ行列、翼型座標点データ、分割翼名、リブ番号、混合率、プロジェクトネーム、メインディレクトリ
function p_end=print_ex(chord,alpha,sparpos,spardia,cutdata,ed_foil,dwingnum,ribnum,mixture,pj_name,now_work)

mkdir("design_data")
cd("design_data")
mkdir(pj_name)
cd(pj_name)
mkdir("ps_out")
cd("ps_out")

fs=0;
%-----翼の前縁を探索
%disp("extending foildata")
do	
	fs=fs+1;
until(ed_foil(fs,1)<=min(ed_foil(:,1)))

n_xdiv_l=200;
n_xdiv_t=200;

x_foil(1:n_xdiv_t+n_xdiv_l-1,1)=1./(n_xdiv_t+n_xdiv_l-1)^2.*([n_xdiv_t+n_xdiv_l-1:-1:1]').^2;	
x_foil(n_xdiv_t+n_xdiv_l,:)=0;
x_foil(n_xdiv_t+n_xdiv_l+1:2*n_xdiv_t+2*n_xdiv_l-1,1)=1./(n_xdiv_t+n_xdiv_l-1)^2.*([1:n_xdiv_t+n_xdiv_l-1]').^2;

y_foil(1:n_xdiv_l+n_xdiv_t,1)=interp1(ed_foil(1:fs,1),ed_foil(1:fs,2),x_foil(1:n_xdiv_l+n_xdiv_t,1),"extrap");
y_foil(n_xdiv_l+n_xdiv_t:2*n_xdiv_t+2*n_xdiv_l-1,1)=interp1(ed_foil(fs:rows(ed_foil),1),ed_foil(fs:rows(ed_foil),2),x_foil(n_xdiv_l+n_xdiv_t:2*n_xdiv_t+2*n_xdiv_l-1,1),"extrap");
y_foil(n_xdiv_l+n_xdiv_t,1)=0;

%翼外形
profx=x_foil;
profy=y_foil;

%翼弦線
clx=[0 1];
cly=[0 0];

%プランク端
plx=[0.05 0.1 0.17 0.28 0.65];

%キャンバーライン
%disp("searching camber")
cambpoint=[0.01:0.01:0.99];
for i=1:length(cambpoint);
	cambline(i,2)=spline(ed_foil(fs-2:-1:2,1),ed_foil(fs-2:-1:2,2),cambpoint(1,i))-spline(ed_foil(fs+2:1:rows(ed_foil)-1,1),ed_foil(fs+2:1:rows(ed_foil)-1,2),cambpoint(1,i));
	cambline(i,1)=cambpoint(1,i);
	bufu=spline(ed_foil(fs-2:-1:2,1),ed_foil(fs-2:-1:2,2),cambpoint(1,i));
	bufl=spline(ed_foil(fs+2:1:rows(ed_foil)-1,1),ed_foil(fs+2:1:rows(ed_foil)-1,2),cambpoint(1,i));

	cambline(i,4)=(bufu+bufl)/2;
endfor

spar_camb=spline(cambline(:,1),cambline(:,2),sparpos);
spar_camb_h=spar_camb/2;

sparmidp(1,2)=min(spline(ed_foil(fs-2:-1:2,1),ed_foil(fs-2:-1:2,2),sparpos),spline(ed_foil(fs+2:1:rows(ed_foil)-1,1),ed_foil(fs+2:1:rows(ed_foil)-1,2),sparpos))+spar_camb_h;
sparmidp(1,1)=sparpos;

%スパー中心に変換
%disp("converting sparpos")
profx_ref(:,1)=profx(:,1).-sparmidp(1,1);
profy_ref(:,1)=profy(:,1).-sparmidp(1,2);
clx_ref=clx.-sparmidp(1,1);
cly_ref=cly.-sparmidp(1,2);
plx_ref=plx.-sparmidp(1,1);

%disp("converting chord")

%コード長考慮
profx_chord=profx_ref.*chord;
profy_chord=profy_ref.*chord;
clx_chord=clx_ref.*chord;
cly_chord=cly_ref.*chord;
plx_chord=plx_ref.*chord;

cntx=400;%ｘ方向基準
cnty=200;%ｙ方向基準

%disp("generating ps code")

openfile=strcat(dwingnum,'_',dec2base(ribnum,10),'.ps');

fpw=fopen(openfile,'wt');
	
fprintf(fpw,'%%!PS-Adobe-3.0\n');
fprintf(fpw,'/Times-Roman findfont 36 scalefont setfont\n')
fprintf(fpw,'/mm {2.834646 mul } def\n')
fprintf(fpw,'%d mm %d mm moveto\n',cntx-55,cnty+10);
if ribnum<10
	fprintf(fpw,'(%s  %s-0%d) show\n','CUTTY SARK',dwingnum,ribnum);
else
	fprintf(fpw,'(%s  %s-%d) show\n','CUTTY SARK',dwingnum,ribnum);
endif
fprintf(fpw,'/Times-Roman findfont 18 scalefont setfont\n')

fprintf(fpw,'%d mm %d mm moveto\n',cntx-70,cnty-15);
fprintf(fpw,'/Times-Roman findfont 18 scalefont setfont\n')
fprintf(fpw,'(chord:%.2f(mm)) show\n',chord);
fprintf(fpw,'%d mm %d mm moveto\n',cntx-70,cnty-22);
fprintf(fpw,'(spar dia:%.1f(mm)(min) %.1f(mm)(max) ) show\n',spardia(1,2),spardia(1,1));
fprintf(fpw,'%d mm %d mm moveto\n',cntx-70,cnty-29);
fprintf(fpw,'(alpha:%.2f(deg) mix:%.2f ) show\n',alpha,mixture);

fprintf(fpw,'newpath\n');
%用紙外枠
fprintf(fpw,'0 mm 0 mm moveto\n 0 mm 841 mm lineto\n 1189 mm 841 mm lineto\n 1189 mm 0 mm lineto\n 0 mm 0 mm lineto\n')
fprintf(fpw,'stroke\n');
%桁中心クロス
fprintf(fpw,'newpath\n');
fprintf(fpw,'%f mm 0 mm moveto\n 0 mm 841 mm rlineto\n',cntx)
fprintf(fpw,'stroke\n');
fprintf(fpw,'newpath\n');
fprintf(fpw,'0 mm %f mm moveto\n 1189 mm 0 mm rlineto\n',cnty)
fprintf(fpw,'stroke\n');
%桁円
fprintf(fpw,'newpath\n');
fprintf(fpw,' %f mm %f mm %f mm 0 360 arc\n',cntx,cnty,spardia(1,1)/2)
fprintf(fpw,' %f mm %f mm %f mm 0 360 arc\n',cntx,cnty,spardia(1,2)/2)
fprintf(fpw,'stroke\n');
%ストリンガ、プランク端指示
for i=1:5
	fprintf(fpw,'newpath\n');
	fprintf(fpw,'%f mm %f mm moveto\n 0 mm %f mm rlineto\n',plx_chord(1,i)+cntx,cly_chord(1,1)+cnty,140);
	fprintf(fpw,'%f mm %f mm moveto\n 0 mm %f mm rlineto\n',plx_chord(1,i)+cntx,cly_chord(1,1)+cnty,-100);
	fprintf(fpw,'stroke\n');
endfor
%よくげんせん
fprintf(fpw,'newpath\n');
fprintf(fpw,'%f mm %f mm moveto\n %f mm %f mm lineto\n',clx_chord(1,1)+cntx,cly_chord(1,1)+cnty,clx_chord(1,2)+cntx,cly_chord(1,2)+cnty);

fprintf(fpw,'stroke\n');
%肉抜き
if cutdata(1,1)==0

else
	for i=1:5
		fprintf(fpw,'newpath\n');
		fprintf(fpw,' %f mm %f mm %f mm 0 360 arc\n',cutdata(i,1)+cntx,cutdata(i,2)+cnty,cutdata(i,3)/2)
		fprintf(fpw,'stroke\n');
		fprintf(fpw,'newpath\n');
		fprintf(fpw,'%f mm %f mm moveto\n 0 mm 20 mm rlineto\n',cutdata(i,1)+cntx,cutdata(i,2)+cnty);
		fprintf(fpw,'%f mm %f mm moveto\n 0 mm -20 mm rlineto\n',cutdata(i,1)+cntx,cutdata(i,2)+cnty);
		fprintf(fpw,'stroke\n');
		fprintf(fpw,'newpath\n');
		fprintf(fpw,'%f mm %f mm moveto\n 20 mm 0 mm rlineto\n',cutdata(i,1)+cntx,cutdata(i,2)+cnty);
		fprintf(fpw,'%f mm %f mm moveto\n -20 mm 0 mm rlineto\n',cutdata(i,1)+cntx,cutdata(i,2)+cnty);
		fprintf(fpw,'stroke\n');
		fprintf(fpw,'/Times-Roman findfont 18 scalefont setfont\n')
		fprintf(fpw,'%d mm %d mm moveto\n',cutdata(i,1)+cntx-25,cutdata(i,2)+cnty-15);
		fprintf(fpw,'(dia:%.1f(mm)) show\n',cutdata(i,3));

	endfor
endif

%迎角線
fprintf(fpw,'newpath\n');
fprintf(fpw,'%f mm %f mm moveto\n %f mm %f mm lineto\n',0,-400*tan(alpha/180*pi)+cnty,1100,700*tan(alpha/180*pi)+cnty);
fprintf(fpw,'stroke\n');

fprintf(fpw,'newpath\n');
for i=1:rows(profx_chord)
	fprintf(fpw,'%f',profx_chord(i,1)+cntx)
	fprintf(fpw,' mm ')
	fprintf(fpw,'%f',profy_chord(i,1)+cnty)
	fprintf(fpw,' mm ');		

	if(i==1)
		fprintf(fpw,' moveto')
		fprintf(fpw,'\n');
	else
		fprintf(fpw,' lineto')
		fprintf(fpw,'\n');
	endif
endfor

fprintf(fpw,'stroke\n');	
fprintf(fpw,'showpage\n');

fclose(fpw);

disp("output complete")
p_end=1;
cd(now_work)

endfunction