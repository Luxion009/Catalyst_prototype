%各種パラメータをもとに断面二次モーメントを求める

function S_I=calc_i(sp,ply_data,mandrel_data,cw,E,C_T)
%座標、積層構成、マンドレルデータ、翼番号、ヤング率、
	
	for i=1:rows(sp)
	%座標点ループ
		S_I(i,1:2)=0;
		for j=1:rows(ply_data)
		%積層番号ループ
			if cw==1
			%中央翼ならば
				if sp(i,1)<=ply_data(j,6)
				%座標点に積層が乗っているなら
					if j==1
					%一層目なら
						S_I(i,1:2)=S_I(i,1:2)+calc_m2(ply_data(j,3),C_T(lookup(C_T(:,1),ply_data(j,2)),2)*2,mandrel_data(cw,4));
					else
						S_I(i,1:2)=S_I(i,1:2)+calc_m2(ply_data(j,3),C_T(lookup(C_T(:,1),ply_data(j,2)),2)*2,ply_data(j-1,1));
					end
			
				end
			else
				if sp(i,1)<=ply_data(j,6)+sum(mandrel_data(1:cw-1,1))*1000
				%座標点に積層が乗っているなら
					if j==1
					%一層目なら
						S_I(i,1:2)=S_I(i,1:2)+calc_m2(ply_data(j,3),C_T(lookup(C_T(:,1),ply_data(j,2)),2)*2,mandrel_data(cw,4));
					else
						S_I(i,1:2)=S_I(i,1:2)+calc_m2(ply_data(j,3),C_T(lookup(C_T(:,1),ply_data(j,2)),2)*2,ply_data(j-1,1));
					end
				end
			end
			if S_I(i,1)==0
				S_I(i,1:2)=interp1(sp(1:i-1,1),S_I(1:i-1,1:2),sp(i,1),'nearest','extrap');
			end	
		end
	end
endfunction