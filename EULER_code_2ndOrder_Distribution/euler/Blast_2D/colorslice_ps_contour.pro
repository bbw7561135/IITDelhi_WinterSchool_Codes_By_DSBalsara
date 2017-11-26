
pro colorslice_ps_contour

nx = 96
ny = 96
nplanes = 13

xmin = -0.5
xmax = 0.5

ymin = -0.5
ymax = 0.5

numcontours = 20

set_plot, 'PS'
;loadct, 13

;TVLCT,R,G,B,/GET
;r[255]=0
;g[255]=0
;b[255]=0
;r[0]=255
;g[0]=255
;b[0]=255
;TVLCT,R,G,B

slicearray = fltarr(nx,ny)

intslicearray = intarr(nx,ny)

threedarray = fltarr(nx,ny,nplanes)

namearray = ['rhoa', 'prsa', 'vlxa', 'vlya', 'vlza', 'bfxa', 'bfya', 'bfza', 'divb', 'levl', 'sp01', 'sp02']


for index = 0, 4 do begin

datamin = 1.0d30
datamax = -1.0d30


for iplane = 1, nplanes do begin

unit = string(iplane,format='(i4.4)')
ifname = namearray[index]+unit
ofname = ifname+".ps"

slicearray[*,*] = 0.0d0

close, 2

openr, 2, ifname

readf, 2, slicearray

close, 2

print, "read ifname = ", ifname

for iy = 0,ny-1 do begin
for ix = 0,nx-1 do begin

threedarray ( ix, iy, iplane - 1) = slicearray ( ix, iy)

endfor
endfor

endfor



;this if condition is only if log scaling is desired

;IF (index EQ 0) OR (index EQ 1) THEN BEGIN

;for iplane = 0, nplanes-1 do begin
;for iy = 0,ny-1 do begin
;for ix = 0,nx-1 do begin

;threedarray ( ix, iy, iplane) = alog10 ( threedarray ( ix, iy, iplane))

;endfor
;endfor
;endfor

;ENDIF



datamax = MAX(threedarray, ix)
datamin = MIN(threedarray, ix)


;print, "datamax, datamin = ", datamax, datamin
;print, "after ifname", ifname

for iplane = 1, nplanes do begin

unit = string(iplane,format='(i4.4)')
ifname = namearray[index]+unit
ofname = ifname+".ps"

for iy = 0,ny-1 do begin
for ix = 0,nx-1 do begin

slicearray ( ix, iy) = threedarray ( ix, iy, iplane - 1)

intslicearray ( ix, iy) = 245.0 * ( slicearray ( ix, iy) - datamin) $
                        /(datamax-datamin) + 4.0

endfor
endfor

device, filename=ofname
print, "writing to ofname = ", ofname


; This is where we make the contour plots with "numcontours" contours from
; (datamin, datamax). If the below lines are uncommented then each
; image will be scaled to its own internal min and max. This may be
; desirable in some instances.

; datamin = MIN(slicearray,ix)
; datamax = MAX(slicearray,ix) + 0.00001

print, "datamax, datamin = ", datamax, datamin
print, "after ofname", ofname

; x and y are arrays dimensioned nx and ny that contain position.

x = xmin + (xmax - xmin) * FINDGEN(nx)/FLOAT(nx-1)
y = ymin + (ymax - ymin) * FINDGEN(ny)/FLOAT(ny-1)

levels = datamin $
 + (datamax - datamin) * FINDGEN(numcontours)/FLOAT(numcontours-1)

; This call gives us a global contouring as it applies to all the files.

!P.THICK=4

contour, slicearray, x, y, levels = levels,TITLE ='3d MHD Blast, 129x129 zones',XSTYLE=1,YSTYLE=1,xthick=6,ythick=6,charthick=5,xcharsize=1.4,ycharsize=1.4,/isotropic

; This call gives us a local contouring as it applies to the present file.

;contour, slicearray, x, y, nlevels=numcontours,TITLE ='2d RP 1, 400x400 zones',XSTYLE=1,YSTYLE=1,xthick=6,ythick=6,charthick=5,xcharsize=1.4,ycharsize=1.4,/isotropic

device, /close 

endfor


wait,1.0


endfor

end










