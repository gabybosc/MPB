FUNCTION DIFF,x,Y
COMMON SHARE,par
eps=par(0)
vy0=par(1)
vz0=par(2)
dy0=eps*(vy0*Y(2)-vz0*Y(1))/Y(0)
dy1=vz0/(eps*Y(0))+(eps/Y(0)-1./eps)*Y(2)
dy2=-vy0/(eps*Y(0))-(eps/Y(0)-1./eps)*Y(1)
return,[dy0,dy1,dy2]
end
;
	PRO mpb,siono,vin
;
COMMON SHARE,par
eps=.1 ; eps=vA/u0 (see draft by Gomez 2020, May 8)
vy0=vin & vz0=vin
par=[eps,vy0,vz0]
;
; Integrates dY/dx=F(x,Y,par) using RK4
; Y: N-dim vector of unknowns
; x: position in Mars-Sun direction x0 and xf
; par: k-dim array of parameters
; Y0: N-dim array of initial values
; res=RK4(Y,Dydx,x,dx,Derivs)
;
; Solve mpb equations for 1D, stationary two-fluid eqs
x0=5. & xf=0. & Nx=29000
N=3 & Y=fltarr(N)
; 
; Define parameters and initial values
Y(0)=1.             ; u(0)
Y(1)=0.             ; by(0)
Y(2)=0.             ; bz(0)
dx=(xf-x0)/float(Nx) & x=x0
;
YY=fltarr(N,Nx) & YY(*,0)=Y
xx=fltarr(Nx) & xx(0)=x0
;
for ix=1,Nx-1 do begin
x=x+dx
xx(ix)=x
dydx=diff(x,Y)
Y=RK4(Y,dydx,x,dx,'diff')
if Y(0) lt 0. then goto,out
YY(*,ix)=Y
endfor
out:
;
; Trim array off trailing zeros
Nx=ix-1
xxx=fltarr(Nx) & YYY=fltarr(N,Nx)
xxx=xx(0:Nx-1)
YYY=YY(*,0:Nx-1)
xx=xxx & YY=YYY
;
; Plot and Overplot by vs bz
if siono eq 0 then begin
xmx=sqrt(2.)/eps
window,0,xs=800,ys=800
plot,yy(1,*),yy(2,*),xr=xmx*[-1,1],yr=xmx*[-1,1],xsty=1,ysty=1,chars=.001
endif
oplot,yy(1,*),yy(2,*)
end
