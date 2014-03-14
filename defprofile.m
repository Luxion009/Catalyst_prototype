%分割翼端の翼形を決メル

%翼形を混合する

%-----翼型の読込
%-----使用翼型名を読込

clear foil foildata x_foil y_foil

fp=fopen('use_foila.txt','r');
use_foil{1}=fgetl(fp);
fclose(fp);

fp=fopen('use_foilb.txt','r');
use_foil{2}=fgetl(fp);
fclose(fp);

fp=fopen('use_foilc.txt','r');
use_foil{3}=fgetl(fp);
fclose(fp);

fp=fopen('use_foild.txt','r');
use_foil{4}=fgetl(fp);
fclose(fp);

fp=fopen('use_foile.txt','r');
use_foil{5}=fgetl(fp);
fclose(fp);

disp("input profiledata");

for i=1:5
	%-----翼型を読込
	fp=fopen(strcat(use_foil{i},'.dat'));
	buff=fgetl(fp);
	foil{i}=(fscanf(fp,'%f',[2,350]))';
	foil{i}=(foil{i});
	fclose(fp);
	fs(i)=0;

	%-----翼の前縁を探索
	do	
		fs(i)=fs(i)+1;
	until(foil{i}(fs(i),1)<=min(foil{i}(:,1)))

end


%中央翼の取り付け角を決定
%任意混合率における任意clを満たすαを探索する関数
% disp("generating central wing");

aoa_g=[1 2 3 4 5 6 7];
for n=1:7
	n
	a_cl(n,1)=xfoil(foil{1},aoa_g(n),chord_re(1,1),0).CL;
end
%clを基準に迎え角を取得する
disp("spline");
ang_cw=spline(a_cl(:,1),aoa_g(1,:),chord_cl(1,1));
%adefはmix,aoa,cl,cd,cm,xcpの順

mt_ang=ang_cw+ed_Hex_i_ang(1,1)
	
%-----翼型データを整理
%-----x方向分割数
n_xdiv_l=80;
n_xdiv_t=60;

disp("generating xy_foil");

x_foil(1:n_xdiv_t+n_xdiv_l-1,1)=1./(n_xdiv_t+n_xdiv_l-1)^2.*([n_xdiv_t+n_xdiv_l-1:-1:1]').^2;	
%x_foil(n_xdiv_t+1:n_xdiv_t+n_xdiv_l,1)=linspace(0.3,0,n_xdiv_l)';
x_foil(n_xdiv_t+n_xdiv_l,:)=0;
x_foil(n_xdiv_t+n_xdiv_l+1:2*n_xdiv_t+2*n_xdiv_l-1,1)=1./(n_xdiv_t+n_xdiv_l-1)^2.*([1:n_xdiv_t+n_xdiv_l-1]').^2;
%x_foil(n_xdiv_t+2*n_xdiv_l:2*n_xdiv_t+2*n_xdiv_l,1)=linspace(0.30001,1,n_xdiv_t)';

for i=1:5
	y_foil(1:n_xdiv_l+n_xdiv_t,i)=interp1(foil{i}(1:fs(i),1),foil{i}(1:fs(i),2),x_foil(1:n_xdiv_l+n_xdiv_t,1),"extrap");
	y_foil(n_xdiv_l+n_xdiv_t:2*n_xdiv_t+2*n_xdiv_l-1,i)=interp1(foil{i}(fs(i):rows(foil{i}),1),foil{i}(fs(i):rows(foil{i}),2),x_foil(n_xdiv_l+n_xdiv_t:2*n_xdiv_t+2*n_xdiv_l-1,1),"extrap");
	%y_foil(n_xdiv_l+n_xdiv_t,i)=0;
end

% foildata=[x_foil,y_foil(:,1)*a(1)+y_foil(:,2)*a(2)+y_foil(:,3)*a(3)+y_foil(:,4)*a(4)+y_foil(:,5)*a(5)];
% thn_f=spline(foildata(1:(rows(x_foil)+1)/2-1,1),foildata(1:(rows(x_foil)+1)/2-1,2),0.37)-spline(foildata((rows(x_foil)+1)/2:(rows(x_foil)+1)/2*2-1,1),foildata((rows(x_foil)+1)/2:(rows(x_foil)+1)/2*2-1,2),0.37)
% dlmwrite(strcat('foildata.dat'),foildata,' ');
% [pol xfoil_res err]=xfoil(foildata,4,500000,0);
% if thn_f>=0.13
% 	thc=0
% else
% 	thc=0.13-thn_f

% endif
% a_cd=pol.CD
% a_cl=pol.CL


function phi1=phi(x,x_foil,y_foil,d_cl,aoa,re,thn)
	printf("opt_cl=%f opt_thn=%f\n",d_cl,thn);
	printf("%f\t%f\t%f\t%f\t%f\n",x(1),x(2),x(3),x(4),x(5));
	thp=0.37;
	adata=[x_foil,y_foil(:,1)*x(1)+y_foil(:,2)*x(2)+y_foil(:,3)*x(3)+y_foil(:,4)*x(4)+y_foil(:,5)*x(5)];
	thn_f=spline(adata(1:(rows(x_foil)+1)/2-1,1),adata(1:(rows(x_foil)+1)/2-1,2),thp)-spline(adata((rows(x_foil)+1)/2:(rows(x_foil)+1)/2*2-1,1),adata((rows(x_foil)+1)/2:(rows(x_foil)+1)/2*2-1,2),thp)
	[pol foil err]=xfoil(adata,aoa,re,0);
	if thn_f>=thn
		thc=(thn-thn_f)/50
	else
		thc=(thn-thn_f)

	endif

	a_cd=pol.CD
	a_cl=pol.CL
	a_LD=a_cl/a_cd;

	if abs(d_cl-a_cl)<0.05
		aclc=(d_cl-a_cl)/10;
	else
		aclc=d_cl-a_cl;
	endif


	phi1=exp(25*abs(aclc))*100*abs(thc)*1000

	disp("------------------------------------")

endfunction

do
	
	for i=2:w_div_c
		do


			clear adata
			if(i==2)
				a=[1.0;0.0;0.0;0.0;0.0];
			endif


			printf("%d翼",i);
			chord_thn(i,1)=input("翼厚定義=");

			opt_cl=chord_cl(i*dw_divn,1);
			opt_aoa=mt_ang-ed_Hex_i_ang(i*dw_divn,1);
			opt_re=chord_re(i*dw_divn,1);
			opt_thn=chord_thn(i,1);
			disp("optimization")
			x=fminunc(@(x)phi(x,x_foil,y_foil,opt_cl,opt_aoa,opt_re,opt_thn),a);
			% [x,obj,info,iter,nf,lambda]=sqp(a,@(x)phi(x,x_foil,y_foil),[],[],[-1;-1;-1;-1],[1;1;1;1])
			a=[x(1);x(2);x(3);x(4);x(5)];
			foil_mixture(i,1:5)=x(1:5);

			adata=[x_foil,y_foil(:,1)*x(1)+y_foil(:,2)*x(2)+y_foil(:,3)*x(3)+y_foil(:,4)*x(4)+y_foil(:,5)*x(5)];
			% plot(adata(:,1),adata(:,2));
			% axis equal;

			thp=0.37;
			de_cl(i,1)=xfoil(adata,mt_ang-ed_Hex_i_ang(i*dw_divn,1),chord_re(i*dw_divn,1),0).CL;
			de_thn(i,1)=spline(adata(1:(rows(x_foil)+1)/2-1,1),adata(1:(rows(x_foil)+1)/2-1,2),thp)-spline(adata((rows(x_foil)+1)/2:(rows(x_foil)+1)/2*2-1,1),adata((rows(x_foil)+1)/2:(rows(x_foil)+1)/2*2-1,2),thp);
			disp("de_cl");
			disp(de_cl(i,1));
			disp("de_thn");
			disp(de_thn(i,1));

			real_chord(i,1)=zerofoil_def(c_cl(i*dw_divn),mt_ang-ed_Hex_i_ang(i*dw_divn),d_vel,adata);
			disp("real_chord")
			disp(real_chord(i,1));

			% thp=0.4;
			% thn_f=spline(adata(1:(rows(x_foil)+1)/2-1,1),adata(1:(rows(x_foil)+1)/2-1,2),thp)-spline(adata((rows(x_foil)+1)/2:(rows(x_foil)+1)/2*2-1,1),adata((rows(x_foil)+1)/2:(rows(x_foil)+1)/2*2-1,2),thp)

			opsign=input("1:ok 0:return");
		
		until(opsign==1)

	endfor

	printf("翼型データ\n");
	printf("中央翼\n");
	printf("必要cl:%f---定義cl:%f\n",chord_cl(1,1),xfoil(foil{1},ang_cw,chord_re(1,1),0).CL);
	printf("内翼\n");
	printf("必要cl:%f---定義cl:%f\n",chord_cl(50,1),de_cl(50,1));
	printf("外翼\n");
	printf("必要cl:%f---定義cl:%f\n",chord_cl(100,1),de_cl(100,1));
	printf("準最外翼\n");
	printf("必要cl:%f---定義cl:%f\n",chord_cl(150,1),de_cl(150,1));
	printf("最外翼\n");
	printf("必要cl:%f---定義cl:%f\n",chord_cl(200,1),de_cl(200,1));
	firure(20)
	plot(1:5,foil_mixture(:,1),"-;a;",1:5,foil_mixture(:,2),"-;b;",1:5,foil_mixture(:,3),"-;c;",1:5,foil_mixture(:,4),"-;d;",1:5,foil_mixture(:,5),"-;e;");

	prof_sig=input("1:ok 0:return");


until (prof_sig==1)