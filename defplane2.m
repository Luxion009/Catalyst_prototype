%循環を元に平面形を作る

%翼型解析データの読み込み
%読込フォルダへの移動

clear foil x_foil y_foil a_cl c_cl_o real_chord

figure(10);
plot(nang_line(1:we_div),chord(1:we_div),"-");

for i=1:w_div_c
	printf("%d翼\n",i)
	disp("Re数");
	disp(chord_re(dw_divn*i,1));
	disp("cl");
	disp(chord_cl(dw_divn*i,1));
	disp("誘導迎え角")
	disp(ed_Hex_i_ang(dw_divn*i,1));
	disp("---------------------")
endfor
buf=input("");

now_work=pwd;
%平面形作成開始

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
	y_foil(n_xdiv_l+n_xdiv_t,i)=0;
end

c_list=0.3:0.1:1.0;
ang_list=-3:1:7;
re_list=c_list.*d_vel/nu;

%各分割翼端のよく弦長を決定する
for i=1:w_div_c-1
	printf("define end_chord:%d\n",i);

	for rec=1:length(c_list)
		% a_cl(i,rec)=xfoil([x_foil y_foil(:,i)],mt_ang-ed_Hex_i_ang(i*dw_divn),re_list(rec),0).CL;
	
		[bufx foil err]=xfoil([x_foil y_foil(:,i)],mt_ang-ed_Hex_i_ang(i*dw_divn),re_list(rec),0);
		if(err==0)
			printf("%d,",rec);
			a_cl(i,rec)=bufx.CL;
		else
			printf("e,");
			if rec==length(c_list)
				a_cl(i,rec)=0;
			else
				a_cl(i,rec)=1000;
			endif

		endif
		%あるRE数におけるあるAOAのCL
	endfor

	for j=1:length(c_list)
		if a_cl(i,j)==1000
			a_cl(i,j)=(a_cl(i,j-1)+a_cl(i,j+1))/2;
		endif
	endfor

	%CLをもとにC＊CLを計算する
	for j=1:length(c_list)
		c_cl_o(i,j)=a_cl(i,j)*c_list(1,j);
	endfor

	real_chord(i,1)=interp1(c_cl_o(i,:),c_list,c_cl(i*dw_divn),"spline");
	disp("");


endfor

%分割翼端における風圧中心位置と翼厚を定義する
for i=1:w_div_c-1
	buffpol=xfoil([x_foil y_foil(:,i)],mt_ang-ed_Hex_i_ang(i*dw_divn),real_chord(i,1)*d_vel/nu,0);
	xcp(i,1)=0.25-buffpol.Cm/buffpol.CL;
	
	%ダイバーじぇんすを考慮して桁を前に
	xcp(i,1)=xcp(i,1)-0.01;

	thn_xcp(i,1)=(spline(x_foil(1:(rows(x_foil)+1)/2-1,1),y_foil(1:(rows(x_foil)+1)/2-1,i),xcp(i,1))-spline(x_foil((rows(x_foil)+1)/2:(rows(x_foil)+1)/2*2-1,1),y_foil((rows(x_foil)+1)/2:(rows(x_foil)+1)/2*2-1,i),xcp(i,1)));
	real_thn_xcp(i,1)=real_chord(i,1).*thn_xcp(i,1);
endfor

real_chord(5,1)=chord(we_div,1);

% real_chord
% xcp
% thn_xcp
% real_thn_xcp

buf=input("");

%分割翼における翼弦長変化率と桁位置変化率を決定する
for i=2:w_div_c
	real_chord_d(i,1)=(real_chord(i,1)-real_chord(i-1,1))/dw_divn;
	if i!=w_div_c
		spar_pos_d(i,1)=(xcp(i,1)-xcp(i-1,1))/dw_divn;
	endif
endfor
% chord((i-1)*dw_divn+1,1)=chord((i-1)*dw_divn,1)+chord_d(i,1);
% 	for j=2:dw_divn
% 		chord((i-1)*dw_divn+j,1)=chord((i-1)*dw_divn+j-1,1)+chord_d(i,1);
% 	endfor

disp("generating A wing");
%中央翼を定義
for m=1:50
	f_data(m,1:3)=[nang_line(m,1),real_chord(1,1),1];
	% f_data(m,4:13)=root_adata(1,:);
	f_mtang(m,1)=mt_ang;
	spar_pos(m,1)=xcp(1,1);
endfor

disp("generating taper wing");
for i=2:w_div_c-1
	for t1=dw_divn*(i-1)+1:dw_divn*i
		% disp(t1);
		f_data(t1,2)=f_data(t1-1,2)+real_chord_d(i,1);
		spar_pos(t1,1)=spar_pos(t1-1,1)+spar_pos_d(i,1);
		f_data(t1,1)=nang_line(t1,1);
		if t1==dw_divn*i
			f_data(t1,3)=1;
		else
			f_data(t1,3)=1-(1/dw_divn)*(t1-dw_divn*(i-1)+1);
		endif
		f_mtang(t1,1)=mt_ang;

	endfor	
endfor

errct=1;
disp("generating twist wing");
aoa_g=[-1 0 1 2 3 4 5];
for n=1:7
	n
	[bufx foil err]=xfoil([x_foil y_foil(:,4)],aoa_g(n),chord_re(we_div,1),0);
	if(err==0)
		a_cl(n,1)=bufx.CL;
	else
		errlistt(errct,1)=n;
		errct++;
	endif
end

if errct!=1
	for n=1:errct-1
		ncct=errlistt(n,1)
		a_cl(ncct,:)=interp1(aoa_g(1,ncct-1:2:ncct+1),a_cl(ncct-1:2:ncct+1,1),aoa_g(1,ncct));
	endfor
endif

a_cl(:,1)

%clを基準に迎え角を取得する
disp("spline");
ang_ew=spline(a_cl(:,1),aoa_g(1,:),chord_cl(we_div,1))

buf=input("")

f_mtang(we_div,1)=ang_ew+ed_Hex_i_ang(we_div,1);
mt_ang_d=(f_mtang(we_div,1)-mt_ang)/50;

buffpol=xfoil([x_foil y_foil(:,4)],ang_ew,real_chord(5,1)*d_vel/nu,0);
xcp(5,1)=0.25-buffpol.Cm/buffpol.CL;
%ダイバージェンスを考慮して桁を前に
xcp(5,1)=xcp(5,1)-0.02;
thn_xcp(5,1)=(spline(x_foil(1:(rows(x_foil)+1)/2-1,1),y_foil(1:(rows(x_foil)+1)/2-1,5),xcp(5,1))-spline(x_foil((rows(x_foil)+1)/2:(rows(x_foil)+1)/2*2-1,1),y_foil((rows(x_foil)+1)/2:(rows(x_foil)+1)/2*2-1,5),xcp(5,1)));
real_thn_xcp(5,1)=real_chord(5,1).*thn_xcp(5,1);

% disp(xcp(5,1))
% disp(thn_xcp(5,1))
% disp(real_thn_xcp(5,1))

real_chord
xcp
thn_xcp
real_thn_xcp


spar_pos_d(5,1)=(xcp(5,1)-xcp(4,1))/dw_divn;

for tw1=dw_divn*4:dw_divn*5
	% disp(tw1);
	f_data(tw1,2)=f_data(tw1-1,2)+real_chord_d(5,1);
	spar_pos(tw1,1)=spar_pos(tw1-1,1)+spar_pos_d(5,1);
	f_data(tw1,1)=nang_line(tw1,1);
	f_data(tw1,3)=0;
	f_mtang(tw1,1)=f_mtang(tw1-1,1)+mt_ang_d;
endfor

%ウイングレット部計算
disp("generating winglet");

disp("generating chord length");
for i=1:4
	i
	f_data(we_div+i*5,2)=zerofoil_def_wl(c_cl(we_div+i*5,1),f_mtang(we_div,1)-ed_Hex_i_ang(we_div+i*5,1),d_vel,[x_foil y_foil(:,4)]);
	% disp("fdata");
	% disp(f_data(we_div+i*5,2));
	wlf_d(i,1)=(f_data(we_div+i*5,2)-f_data(we_div+(i-1)*5,2))/5;

endfor



% wlc_list=0.1:0.1:0.8;
% wlre_list=wlc_list.*d_vel/nu;
% for i=1:4
% 	i

% 	for rec=1:length(wlc_list)
% 		printf("%d,",rec);
% 		wl_a_cl(i,rec)=xfoil([x_foil y_foil(:,4)],f_mtang(we_div,1)-ed_Hex_i_ang(we_div+i*5,1),wlre_list(rec),0).CL;
% 		%あるRE数におけるあるAOAのCL
% 	endfor

% 	%CLをもとにC＊CLを計算する
% 	for j=1:length(wlc_list)
% 		wl_c_cl_o(i,j)=wl_a_cl(i,j)*wlc_list(1,j);
% 	endfor

% 	f_data(we_div+i*5,2)=interp1(wl_c_cl_o(i,:),wlc_list,c_cl(we_div+i*5,1),"spline");

% 	wlf_d(i,1)=(f_data(we_div+i*5,2)-f_data(we_div+(i-1)*5,2))/5;


% endfor






for i=1:4
	for j=1:4
		we_div+(i-1)*5+j
		f_data(we_div+(i-1)*5+j,2)=f_data(we_div+(i-1)*5,2)+wlf_d(i,1)*j;
	endfor
endfor

for i=we_div+1:we_div+wld_c
	i
	f_data(i,1)=nang_line(i,1);
	f_data(i,3)=0;
	f_mtang(i,1)=f_mtang(we_div,1);
	spar_pos(i,1)=spar_pos(we_div,1);
endfor

%翼面積を求める
for i=1:we_div+wld_c
	s_d(i,1)=f_data(i,2)*line_e_d(i,1);
endfor

mw_s=sum(s_d)*2;

figure(10)
	plot(f_data(1:we_div,1),f_data(1:we_div,2),f_data(we_div+1:we_div+wld_c,1),f_data(we_div+1:we_div+wld_c,2));
	ylabel("chord");
	grid on;
	cd("design_data")
	cd(pj_name)
	print("re_chord.png",'-dpng','-r100')
	cd(now_work)

figure(11)
	plot(f_data(:,1),f_data(:,3));
	ylabel("mixture");
	grid on;
	cd("design_data")
	cd(pj_name)
	print("re_mixture.png",'-dpng','-r100')
	cd(now_work)

% figure(12)
% 	plot(f_data(:,1),mt_ang(:,1));
% 	ylabel("mt_ang");
% 	grid on;
% 	cd("design_data")
% 	cd(pj_name)
% 	print("re_mt_ang.png",'-dpng','-r100')
% 	cd(now_work)

figure(13)
	plot(f_data(1:we_div,1),f_data(1:we_div,2),f_data(we_div+1:we_div+wld_c,1),f_data(we_div+1:we_div+wld_c,2));
	ylabel("chord");
	axis("equal");
	grid on;
	cd("design_data")
	cd(pj_name)
	print("re_chord_equal.png",'-dpng','-r100')
	cd(now_work)

figure(14)
	plot(f_data(:,1),f_data(:,2),"-;chord;",f_data(:,1),f_data(:,3),"-;mixture;");
	xlabel('y[m]');
	grid on;
	cd("design_data")
	cd(pj_name)
	print("c_mix.png",'-dpng','-r100')
	cd(now_work)

% figure(15)
% 	plot(f_data(:,1),ed_gamma(:,1),"-;ed_gamma;",f_data(:,1),re_gamma(:,1),"-;re_gamma;");
% 	xlabel('y[m]');
% 	grid on;
% 	cd("design_data")
% 	cd(pj_name)
% 	print("ed_re_gamma.png",'-dpng','-r100')
% 	cd(now_work)

figure(15)
	plot(f_data(:,1),spar_pos(:,1),"-;sparpos;");
	xlabel("y[m]");
	grid on;
	cd("design_data")
	cd(pj_name)
	print("sparpos.png",'-dpng','-r100')
	cd(now_work)