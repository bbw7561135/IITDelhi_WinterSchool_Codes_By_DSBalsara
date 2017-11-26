
pro onedplot

nx = 400
nplanes = 1

loadct, 13
set_plot, 'PS'

x = fltarr(nx)
y = fltarr(nx)

twodarray = fltarr(nx,nplanes)

namearray=['rhoa','prsa','vlxa','vlya','vlza','bfxa','bfya','bfza','sp01']


for index = 0, 7 do begin

datamin = 1.0d30
datamax = -1.0d30


for iplane = 1, nplanes do begin

unit = string(iplane,format='(i4.4)')
ifname = namearray[index]+unit
ofname = ifname+".ps"

close, 2

openr, 2, ifname

readf, 2, y

readf, 2, x

close, 2

print, "read ifname = ", ifname

for ix = 0,nx-1 do begin

twodarray ( ix, iplane - 1) = y ( ix)

endfor

endfor



;this if condition is only if log scaling is desired

;IF (index EQ 0) OR (index EQ 1) THEN BEGIN

;for iplane = 0, nplanes-1 do begin
;for ix = 0,nx-1 do begin

;twodarray ( ix, iplane) = alog10 ( twodarray ( ix, iplane))

;endfor
;endfor

;ENDIF



datamax = MAX(twodarray, ix)
datamin = MIN(twodarray, ix)

tempaa = datamax - datamin

datamax = datamax + 0.08 * tempaa
datamin = datamin - 0.08 * tempaa

xmax = MAX(x, ix)
xmin = MIN(x, ix)

print, "datamax, datamin = ", datamax, datamin
print, "after ifname", ifname

for iplane = 1, nplanes do begin

unit = string(iplane,format='(i4.4)')
ifname = namearray[index]+unit
ofname = ifname+".ps"

for ix = 0,nx-1 do begin

y ( ix) = twodarray ( ix, iplane - 1)

endfor

device, filename=ofname, xsize=18, ysize=18, xoffset=1, yoffset=7, $
 /COLOR, BITS = 8

print, "writing to ofname = ", ofname

;This sets the axis range and style

!X.STYLE = 1
!Y.STYLE = 1
!X.RANGE=[xmin,xmax]
!Y.RANGE=[datamin,datamax]

;This sets the line thickness, style and color.
;To set linestyle, use LINESTYLE=0 for solid; 1 dotted; 2 dashed; 
;3 dash-dot; 4 dash-dot-dot; 5 long-dash

!P.THICK=6
!P.LINESTYLE = 0
!P.COLOR = 1


;For log scaling of x and/or y-axis, put /XLOG and/or /YLOG in "PLOT..."


;PLOT,x,y,XTITLE='x',YTITLE=namearray[index],xthick=6,ythick=6,charthick=5,xcharsize=2.0,ycharsize=2.0,/device


;To overplot, set new line thickness, style and color and set a new "y" & call:

;OPLOT,x,y

;To explicitly plot out data points, uncomment following line and comment out
;the previous line. 
;To set symbols, use PSYM=1 for +; 2 for *; 3 for .; 4 for diamond; 
;5 for triangle; 6 for square; 7 for X;


PLOT,x,y,PSYM=1,XTITLE='x',YTITLE=namearray[index],TITLE='ADER 2nd order, HLLE, CFL.8',/device


;To overplot, set new line thickness, symbol and color and set a new "y" & call:

;OPLOT,x,y,PSYM=4

device, /close

endfor



endfor

end










