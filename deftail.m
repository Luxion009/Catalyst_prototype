%運動解析を行い尾翼を決定


%アスペクト比を求める
aspect_ratio=(span+wl_span*2)^2/mw_s;

%主翼空力平均翼弦を求める
cbar=2/mw_s*sum(f_data(:,2).^2.*line_e_d(:,1));

%空力平均翼弦の位置（主翼前縁を原点とする座標系）を求める
cbar_pos=interp1(f_data(51:250,2),f_data(51:250,1),cbar,'linear');%補間関数は狭義単調増加でなくてはならない
%空力平均翼弦の桁位置(重心位置)を求める
cbar_sparpos=interp1(f_data(:,1),spar_pos(:,1),cbar_pos,'linear');
cog=cbar_sparpos;

%空力平均翼弦における桁位置から空力中心(0.25)までの距離(m)を求める
cbar_ac=cbar*(cbar_sparpos-0.25);

%主翼揚力傾斜を求める(0.5度で解析)
%L_cla=clalpha_analysis_m(f_data,Q_ij,mt_ang,line_e_d,lst_re,lst_i_ang,lst_comp_vel,d_vel,wld_c,we_div,a_density,adata)(1,1);

cla_analysis

alphaw=(cla_L_lst*9.8-L_lst*9.8)*2/a_density/d_vel^2/mw_s*180/pi/0.5;

%尾翼位置での吹き下ろし傾斜を求める(楕円翼と近似する)
% dep_dal=(mean(lst_i_ang(1:we_div/5)).*2.-mean(cla_lst_i_ang(1:we_div/5)).*2)*180/pi/0.5;
dep_dal=2*alphaw/pi/aspect_ratio;

% %水平尾翼空力データ読み込み
% cd("NACA0012")

% %-----翼型の読込
% fp=fopen('hsfoil.dat');
% hsfoil=fgetl(fp);
% hsfoil=(fscanf(fp,'%f',[2,250]))';
% fclose(fp);
% fshs=0;


% %-----尾翼空力解析データの
% %-----レイノルズリスト読み込み
% Relist={50000,'0.050';100000,'0.100'; 150000,'0.150'; 200000,'0.200'; 250000,'0.250'; 300000,'0.300';350000,'0.350'; 400000,'0.400';450000,'0.450'; 500000,'0.500';550000,'0.550';600000,'0.600';};;

% for m=1:12
% %mは順に50000,100000,150000,200000,250000,300000,350000,400000,450000,500000,550000,600000

% 	foilname='hsfoil';
% 	datafile_a='_T1_Re';
% 	Re=Relist{12+m};
% 	datafile_b='_M0.00_N9.0.txt';
% 	dataname=strcat(foilname,datafile_a,Re,datafile_b);

% 	fpr=fopen(dataname(1,:));
% 	do
% 		buffer=fscanf(fpr,'%c',[1,1]);
% 	until(buffer(1,1)=='-')
% 	buffer=fgetl(fpr);
% 	ahdata{1,m}=(fscanf(fpr,'%f',[10,100]))';
% 	fclose(fpr);
% endfor

% cd(now_work)

% %-----垂直尾翼翼型の読込
% %-----使用翼型を読込
% cd("NACA0010")

% %-----翼型の読込
% fp=fopen('vsfoil.dat');
% vsfoil=fgetl(fp);
% vsfoil=(fscanf(fp,'%f',[2,250]))';
% fclose(fp);
% fsvs=0;


% %-----尾翼空力解析データの
% %-----レイノルズリスト読み込み
% Relist={50000,'0.050';100000,'0.100'; 150000,'0.150'; 200000,'0.200'; 250000,'0.250'; 300000,'0.300';350000,'0.350'; 400000,'0.400';450000,'0.450'; 500000,'0.500';550000,'0.550';600000,'0.600';};;

% for m=1:12
% %mは順に50000,100000,150000,200000,250000,300000,350000,400000,450000,500000,550000,600000

% 	foilname='vsfoil';
% 	datafile_a='_T1_Re';
% 	Re=Relist{12+m};
% 	datafile_b='_M0.00_N9.0.txt';
% 	dataname=strcat(foilname,datafile_a,Re,datafile_b);

% 	fpr=fopen(dataname(1,:));
% 	do
% 		buffer=fscanf(fpr,'%c',[1,1]);
% 	until(buffer(1,1)=='-')
% 	buffer=fgetl(fpr);
% 	avdata{1,m}=(fscanf(fpr,'%f',[10,100]))';
% 	fclose(fpr);
% endfor

% cd(now_work)


	% %テールパイプ作成
	% disp("テールパイプ構成（全手動入力）");
	% t_manddata(1,1)=input("テール長(mm)=");
	% t_manddata(1,2)=input("テール径(mm)=");
	% t_manddata(1,3)=input("接合部長さ(mm)=");
	% t_manddata(1,4)=t_manddata(1,2);%マンドレル径、テール径と同じ
	% t_manddata(1,5)=input("テールマンドレル長さ(mm)=");

	% %積層構成入力

	% %初期構成決定
	% t_plydata(1:4,2)=24;
	% %24t
	% t_plydata(1:4,3)=90;
	% %積層角度
	% t_plydata(1:3:4,4)=90;
	% t_plydata(2:3,4)=45;
	% %積層方向
	% t_plydata(1:4,5)=0;
	% t_plydata(1:4,6)=t_manddata(1,1);

	% t_plydata(1,1)=t_manddata(cw,4)+C_T(lookup(C_T(:,1),t_plydata(1,2)),2)*2;
	% for j=2:4
	% 		t_plydata(j,1)=t_plydata(j-1,1)+C_T(lookup(C_T(:,1),t_plydata(j,2)),2)*2;
	% endfor
	% %外形計算

	% %手動入力
	% do
	% 	ts=input("１：部分積層増加　２：一周積層増加　３：積層減少 ---");

	% 	if ts==1
	% 		ra=input("積層角度=");
	% 		t_plydata=ply_ex(t_plydata,1,ra,C_T,t_manddata(1,1))

	% 	elseif ts==2
	% 		rs=input("層番号=");
	% 		t_plydata=ply_ins(t_plydata,rs,C_T,t_manddata(1,1));

	% 	elseif ts==3
	% 		rs=input("層番号=");
	% 		t_plydata(rs,:)=[];
	% 		for j=2:rows(t_plydata)
	% 			t_plydata(j,1)=t_plydata(j-1,1)+C_T(lookup(C_T(:,1),t_plydata(j,2)),2)*2;
	% 		end
	% 	endif

	% 	%積層構成を表示
	% 	ply_disp(t_plydata,000);
		
	% 	tsign=input("ok:1 return:0 ---")
			
	% until(tsign==0)

buf=input("deftail end")