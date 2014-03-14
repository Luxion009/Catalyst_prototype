for n=1:div_c %翼分割数
	%-----座標、積層データ、マンドレルデータ、翼番号、物性を送って各座標点での断面二次モーメントを計算する関数(kgf.m^2)
	%disp("generating S_I");
	S_I(dp_num*(n-1)+1:dp_num*n,1:2)=calc_i(M_data(dp_num*(n-1)+1:dp_num*n,1).*1000,ply_data{n},mand_data,n,E,C_T)./1000000;
	%-----座標、積層データ、マンドレルデータ、翼番号、物性を送って各座標点での曲げ剛性を計算する関数(kgf.m^2)
	%disp("generating S_H");
	S_H(dp_num*(n-1)+1:dp_num*n,1:2)=calc_hardness(M_data(dp_num*(n-1)+1:dp_num*n,1).*1000,ply_data{n},mand_data,n,E,C_T)./1000000;
	%-----各座標での線密度を計算する関数(g/m)
	%disp("generating D_mat");
	if n<div_c
		D_mat(dp_num*(n-1)+1:dp_num*n,1)=calc_condens(M_data(dp_num*(n-1)+1:dp_num*n,1).*1000,ply_data{n},ply_data{n+1},mand_data,n,C_D,C_T).*1000;
	else
		D_mat(dp_num*(n-1)+1:dp_num*n,1)=calc_ldens(M_data(dp_num*(n-1)+1:dp_num*n,1).*1000,ply_data{n},mand_data,n,C_D,C_T).*1000;

	endif
end

%材料計算に必要なデータを別変数に格納
%disp("Data transfer");
S_mat(:,1)=M_data(:,1);
S_mat(:,6)=D_mat(:,1);%座標における線密度
S_mat(:,7)=M_data(:,2);%コード長
S_mat(:,2)=M_data(:,3);%局所揚力係数
S_mat(:,3)=M_data(:,4);%局所抗力係数
S_mat(:,4)=M_data(:,5)./180.*pi;%誘導角度（rad）
S_mat(:,5)=M_data(:,6);%合成速度
ID_swang(:,1)=interp1(nang_line(1:250,1),line_e(1:250,2),S_mat(:,1),'spline','extrap');%撓んだ時のZ方向座標

%disp("genarating d_L");
%-----局所揚力分布(N/m)
S_mat_dL(:,1)=(S_mat(:,2).*cos(S_mat(:,4))-S_mat(:,3).*sin(S_mat(:,4))).*0.5.*S_mat(:,7).*a_density.*(S_mat(:,5).^2)-S_mat(:,6)./1000.*9.8-Snd_rho.*S_mat(:,7).*9.8;

%disp("generating Q");
for i=1:div_c*dp_num-1
	S_mat_dQ(i,1)=(S_mat_dL(i,1)+S_mat_dL(i+1,1))/2*(S_mat(i+1,1)-S_mat(i,1));
end
S_mat_dQ(div_c*dp_num,1)=interp1(S_mat(1:div_c*dp_num-1,1),S_mat_dQ(:,1),S_mat(div_c*dp_num,1),'linear','extrap');

S_mat_Q(div_c*dp_num+1,1)=0;
for i=div_c*dp_num:-1:1
	S_mat_Q(i,1)=S_mat_Q(i+1,1)+S_mat_dQ(i,1);
end
%-----せん断力(kgf)
S_mat_Q=S_mat_Q./9.8;

%disp("generating M");
for i=1:div_c*dp_num-1
	S_mat_dM(i,1)=(S_mat_Q(i,1)+S_mat_Q(i+1,1))/2*(S_mat(i+1,1)-S_mat(i,1))*1000;
end
%-----モーメント(kgf.m)
S_mat_M(div_c*dp_num,1)=0;
for i=div_c*dp_num-1:-1:1
	S_mat_M(i,1)=S_mat_M(i+1,1)+S_mat_dM(i,1);
end

%disp("generating Zz");
for n=1:div_c
	Zz(dp_num*(n-1)+1:dp_num*n,1)=S_I(dp_num*(n-1)+1:dp_num*n,2)./ply_data{n}(rows(ply_data{n}),1).*1000;%断面係数
	S_max_Z(1:dp_num*n,1)=S_mat_M(1:dp_num*n,1)./Zz(1:dp_num*n,1).*9.8./1000;%最大曲げ応力
end
		
for i=1:div_c*dp_num-1
	S_mat_dtheta(i,1)=(S_mat_M(i,1)/S_H(i,2)+S_mat_M(i+1,1)/S_H(i+1,2))/2*(S_mat(i+1,1)-S_mat(i,1));%連続する二つの曲率の平均にその間の距離をかける
	%接線の角度変化
end
S_mat_theta(1,1)=0;%接線と水平軸の角度
for i=2:div_c*dp_num-1
	S_mat_theta(i,1)=S_mat_theta(i-1,1)+S_mat_dtheta(i,1);
end
S_mat_theta(div_c*dp_num,1)=interp1(S_mat(1:div_c*dp_num-1,1),S_mat_theta(:,1),S_mat(div_c*dp_num,1),'linear','extrap');

for i=1:div_c*dp_num-1
	S_mat_dswang(i,1)=(S_mat_theta(i,1)+S_mat_theta(i+1,1))./2.*(S_mat(i+1,1)-S_mat(i,1));%sin x=xと近似できるのでZ方向変化が算出できる
end
S_mat_swang(1,1)=S_mat_dswang(1,1);
for i=2:div_c*dp_num-1
	S_mat_swang(i,1)=S_mat_swang(i-1,1)+S_mat_dswang(i,1);%Z方向位置
end
S_mat_swang(div_c*dp_num,1)=interp1(S_mat(1:div_c*dp_num-1,1),S_mat_swang(:,1),S_mat(div_c*dp_num,1),'linear','extrap');
S_mat_swang(:,1)=S_mat_swang.*10^(-3);

W_Spar=trapz(S_mat(:,1),S_mat(:,6)./1000)*2;
W_wing=W_Spar+trapz(S_mat(:,1),Snd_rho.*S_mat(:,7))*2;

%disp("end bend_disp");

calc_swang=S_mat_swang;