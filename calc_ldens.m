%各種パラメータをもとに局所線密度を求める
function S_LD=calc_ldens(Spar_dS,ply_data,mandrel_data,cw,C_D,C_T)
	
	for i=1:rows(Spar_dS)
	%座標点ループ
		S_LD(i,1)=0;
		for j=1:rows(ply_data)
		%積層ループ
			if Spar_dS(i,1)<=ply_data(j,6)+(sum(mandrel_data(1:cw,1))-mandrel_data(cw,1))*1000
			%
				if j==1
					S_LD(i,1)=S_LD(i,1)+C_D(lookup(C_D(:,1),ply_data(j,2)),2).*(ply_data(1,1)^2-mandrel_data(cw,4)^2)/8*ply_data(j,3)*pi/180*4;
				else
					S_LD(i,1)=S_LD(i,1)+C_D(lookup(C_D(:,1),ply_data(j,2)),2).*(ply_data(j,1)^2-ply_data(j-1,1)^2)/8*ply_data(j,3)*pi/180*4;
				end
			end
		end
	end
endfunction