%各種パラメータをもとに局所線密度を求める
function S_LD=calc_condens(Spar_dS,ply_data,next_ply_data,mandrel_data,cw,C_D,C_T)
	
	nextjunc=mandrel_data(cw+1,3);
	junc_c=round(nextjunc./(mandrel_data(cw,1).*1000).*rows(Spar_dS));

	for i=1:rows(Spar_dS)
	%座標点ループ
		S_LD(i,1)=0;
		for j=1:rows(ply_data)
		%積層ループ
			if Spar_dS(i,1)<=ply_data(j,6)+(sum(mandrel_data(1:cw,1))-mandrel_data(cw,1))*1000
			%積層が座標点で存在するかどうか
				if j==1
					S_LD(i,1)=S_LD(i,1)+C_D(lookup(C_D(:,1),ply_data(j,2)),2).*(ply_data(1,1)^2-mandrel_data(cw,4)^2)/8*ply_data(j,3)*pi/180*4;
				else
					S_LD(i,1)=S_LD(i,1)+C_D(lookup(C_D(:,1),ply_data(j,2)),2).*(ply_data(j,1)^2-ply_data(j-1,1)^2)/8*ply_data(j,3)*pi/180*4;
				end
			end
		end
	end

	for i=rows(Spar_dS)-junc_c+1:rows(Spar_dS)
		
		for j=1:rows(next_ply_data)
		%積層ループ	
			if j==1
				S_LD(i,1)=S_LD(i,1)+C_D(lookup(C_D(:,1),next_ply_data(j,2)),2).*(next_ply_data(1,1)^2-mandrel_data(cw+1,4)^2)./4.*pi;
			else
				S_LD(i,1)=S_LD(i,1)+C_D(lookup(C_D(:,1),next_ply_data(j,2)),2).*(next_ply_data(j,1)^2-next_ply_data(j-1,1)^2)./4.*pi;
			endif

		endfor
	endfor

		

endfunction