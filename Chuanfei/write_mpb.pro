@idlrcc

;Rho Ux Uy Uz Bx By Bz Pe P HpRho
;HpUx HpUy HpUz HpP O2pRho O2pUx O2pUy O2pUz O2pP OpRho
;OpUx OpUy OpUz OpP CO2pRho CO2pUx CO2pUy CO2pUz CO2pP jx
;jy jz gradpex gradpey gradpez gradper b1x b1y b1z

filename='../z=0_HallOff.out'
npict=10
.r getpict
xmn =0 & xmx=2 & !x.range=[xmn,xmx]
ymn =-1 & ymx=1 & !y.range=[ymn,ymx]

xx =x(*,0,0) & yy=x(*,0,1)
rho = w(*,0,0) & HpRho = w(*,0,9) & OpRho = w(*,0,19) & O2pRho = w(*,0,14) & CO2pRho = w(*,0,24)
vel_Hx = w(*,0,10) & vel_Hy = w(*,0,11) & vel_Hz = w(*,0,12)
Bx = w(*,0,4) & By = w(*,0,5) & Bz = w(*,0,6)
b1x = w(*,0,36) & b1y = w(*,0,37) & b1z = w(*,0,38)
Pe = w(*,0,7) & Pp = w(*,0,8) & PH = w(*,0,13) & PO2 = w(*,0,18) & PO = w(*,0,23) & PCO2 = w(*,0,28)
jx = w(*,0,29) & jy = w(*,0,30) & jz = w(*,0,31)
epx = w(*,0,32) & epy = w(*,0,33) & epz = w(*,0,34)

p = where(xx gt xmn and xx lt xmx and yy gt ymn and yy lt ymx)

xxx = xx(p) & yyy = yy(p)
rrho = rho(p) & Hrho = HpRho(p) & Orho = OpRho(p) & O2rho = O2pRho(p) & CO2rho = CO2pRho(p)
vx = vel_Hx(p) & vy = vel_Hy(p) & vz = vel_Hz(p)
Bbx = Bx(p) & Bby = By(p) & Bbz = Bz(p)
bb1x = b1x(p) & bb1y = b1y(p) & bb1z = b1z(p)
ppe = Pe(p) & ppp = Pp(p) & pph = PH(p) & ppO = PO(p) & ppO2 = PO2(p) & ppCO2 = PCO2(p)
jjx = jx(p) & jjy = jy(p) & jjz = jz(p)
eepx = epx(p) & eepy = epy(p) & eepz = epz(p)

N=(size(xxx))(1)

openw,1,'outputs/posicion.csv'
for i=0,N-1 do printf,1,xxx(i),yyy(i)
close,1

openw,1,'outputs/rho.csv'
for i=0,N-1 do printf,1,rrho(i),Hrho(i),Orho(i),O2rho(i),CO2rho(i)
close,1

openw,1,'outputs/velocidad_h.csv'
for i=0,N-1 do printf,1,vx(i),vy(i),vz(i)
close,1

openw,1,'outputs/campo.csv'
for i=0,N-1 do printf,1,Bbx(i),Bby(i),Bbz(i),bb1x(i),bb1y(i),bb1z(i)
close,1

openw,1,'outputs/presion.csv'
for i=0,N-1 do printf,1,ppe(i),ppp(i),pph(i),ppO(i),ppO2(i),ppCO2(i)
close,1

openw,1,'outputs/corrientes.csv'
for i=0,N-1 do printf,1,jjx(i),jjy(i),jjz(i)
close,1

openw,1,'outputs/grad_p.csv'
for i=0,N-1 do printf,1,eepx(i),eepy(i),eepz(i)
close,1
