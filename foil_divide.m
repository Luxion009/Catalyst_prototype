function divide=foil_divide(x1,x2,alpha)
	divide=(-2+sqrt(4-4*(2*alpha*alpha*(x2-x1)*(-alpha*alpha*(x1*x1+x2*x2)*(x2-x1)-1*(x1+x2)))))/2/(2*alpha*alpha*(x2-x1));
endfunction
%荷重曲線に沿って線素を等分する関数