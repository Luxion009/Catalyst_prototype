%�e��p�����[�^�����Ƃɒf�ʓ񎟃��[�����g�����߂�

function S_I=calc_i(sp,ply_data,mandrel_data,cw,E,C_T)
%���W�A�ϑw�\���A�}���h�����f�[�^�A���ԍ��A�����O���A
	
	for i=1:rows(sp)
	%���W�_���[�v
		S_I(i,1:2)=0;
		for j=1:rows(ply_data)
		%�ϑw�ԍ����[�v
			if cw==1
			%�������Ȃ��
				if sp(i,1)<=ply_data(j,6)
				%���W�_�ɐϑw������Ă���Ȃ�
					if j==1
					%��w�ڂȂ�
						S_I(i,1:2)=S_I(i,1:2)+calc_m2(ply_data(j,3),C_T(lookup(C_T(:,1),ply_data(j,2)),2)*2,mandrel_data(cw,4));
					else
						S_I(i,1:2)=S_I(i,1:2)+calc_m2(ply_data(j,3),C_T(lookup(C_T(:,1),ply_data(j,2)),2)*2,ply_data(j-1,1));
					end
			
				end
			else
				if sp(i,1)<=ply_data(j,6)+sum(mandrel_data(1:cw-1,1))*1000
				%���W�_�ɐϑw������Ă���Ȃ�
					if j==1
					%��w�ڂȂ�
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