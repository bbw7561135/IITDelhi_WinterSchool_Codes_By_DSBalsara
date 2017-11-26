
pro colorslice_jpg

nx = 96
ny = 96
nplanes = 15

set_plot, 'X'
loadct, 13

TVLCT,R,G,B,/GET
r[255]=0
g[255]=0
b[255]=0
r[0]=255
g[0]=255
b[0]=255
TVLCT,R,G,B

slicearray = fltarr(nx,ny)

intslicearray = intarr(nx,ny)

threedarray = fltarr(nx,ny,nplanes)

namearray = ['rhoa', 'prsa', 'vlxa', 'vlya', 'vlza', 'bfxa', 'bfya', 'bfza', 'divb', 'levl', 'sp01', 'sp02']


for index = 0, 7 do begin

datamin = 1.0d30
datamax = -1.0d30


for iplane = 1, nplanes do begin

unit = string(iplane,format='(i4.4)')
ifname = namearray[index]+unit
ofname = ifname+".jpg"

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


print, "datamax, datamin = ", datamax, datamin
print, "after ifname", ifname

for iplane = 1, nplanes do begin

unit = string(iplane,format='(i4.4)')
ifname = namearray[index]+unit
ofname = ifname+".jpg"

for iy = 0,ny-1 do begin
for ix = 0,nx-1 do begin

slicearray ( ix, iy) = threedarray ( ix, iy, iplane - 1)

endfor
endfor

; If you want each image to be scaled to its own min and max then uncomment
; the next two lines. Otherwise, the images will be scaled to a global minmax.
; datamax = MAX(slicearray, ix) + 0.000001
; datamin = MIN(slicearray, ix)

for iy = 0,ny-1 do begin
for ix = 0,nx-1 do begin

intslicearray ( ix, iy) = 245.0 * ( slicearray ( ix, iy) - datamin) $
                        /(datamax-datamin) + 4.0

endfor
endfor

device, retain=2,decomposed=0
window,0,xsize=nx,ysize=ny*1.25
print, "writing to ofname = ", ofname

tv, intslicearray, /device
colorbar,ncolors=255,position=[0.15,0.90,0.90,0.98], range=[datamin,datamax], format='(F4.1)',DIVISIONS=4
wait,0.1

write_jpeg,ofname,tvrd(true=3),true=3,quality=100


endfor


wait,1.0


endfor

end










