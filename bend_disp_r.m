for n=1:div_c
	%-----座標、積層データ、マンドレルデータ、翼番号、物性を送って各座標点での断面二次モーメントを計算する関数(kgf.m^2)
	S_I(dp_num*(n-1)+1:dp_num*n,1:2)=calc_i(M_data(dp_num*(n-1)+1:dp_num*n,1).*1000,ply_data{n},mand_data,n,E,C_T)./1000000;
	%-----座標、積層データ、マンドレルデータ、翼番号、物性を送って各座標点での曲げ剛性を計算する関数(kgf.m^2)
	S_H(dp_num*(n-1)+1:dp_num*n,1:2)=calc_hardness(M_data(dp_num*(n-1)+1:dp_num*n,1).*1000,ply_data{n},mand_data,n,E,C_T)./1000000;
	%-----各座標での線密度を計算する関数(kg/m)
	D_mat(dp_num*(n-1)+1:dp_num*n,1)=calc_ldens(M_data(dp_num*(n-1)+1:dp_num*n,1).*1000,ply_data{n},mand_data,n,C_D,C_T).*1000;
end

S_mat(:,6)=D_mat(:,1);

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
	S_mat_dswang(i,1)=(S_mat_theta(i,1)+S_mat_theta(i+1,1))/2*(S_mat(i+1,1)-S_mat(i,1));%sin x=xと近似できるのでZ方向変化が算出できる
end
S_mat_swang(1,1)=S_mat_dswang(1,1);
for i=2:div_c*dp_num-1
	S_mat_swang(i,1)=S_mat_swang(i-1,1)+S_mat_dswang(i,1);%Z方向位置
end
S_mat_swang(div_c*dp_num,1)=interp1(S_mat(1:div_c*dp_num-1,1),S_mat_swang(:,1),S_mat(div_c*dp_num,1),'linear','extrap');
S_mat_swang(:,1)=S_mat_swang.*10^(-3);

W_Spar=trapz(S_mat(:,1),S_mat(:,6)./1000)*2;
W_wing=W_Spar+trapz(S_mat(:,1),Snd_rho.*S_mat(:,7))*2;