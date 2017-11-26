pro draw_sph

; To run do "idl"
; .compile draw_sph
; draw_sph


npoints = 250000 ;79205 ; 93440 ; 98304 ; read this number of draw.sph.out

nplane_end = 4

nplane_begin = 1


outval = fltarr ( npoints)
xout = fltarr ( npoints)
yout = fltarr ( npoints)

; loadct == 13 gives a color table with a black zero point,
; loadct == 25 gives a color table with a blue zero point which is
;    easier on the eyes.
loadct,  25
set_plot, 'PS'


namearray = ['rhol', 'prsl','vx01','vy01', 'vz01', 'sp01', 'sp02', 'rhoa', 'prsa', 'vlya', 'vlza']
titlearray = ['v_radial', 'log_temperature','temperature','log_density', 'log_pressure', 'species_fraction_1', 'species_fraction_2', 'density',  'pressure',   'v_lattitude', 'v_toroidal']

if ( SIZE(namearray,/N_ELEMENTS) NE SIZE(titlearray,/N_ELEMENTS) ) then stop

;index = 0

; ---------------------------------------------------------------------

; Start of the major loop that loops over flow variables whose name begins
; with the starting tag in "namearray".

for index = 0, 4 do begin

; -------------------------------------------------

globalmin = 1.0e9
globalmax =  -1.0e9

print, "Finding global ymin, ymax for: ", namearray[index]

for iplane = nplane_begin, nplane_end do begin

unit = string(iplane,format='(i4.4)')
ifname = namearray[index]+unit
ofname = ifname+".ps"

print, "reading ifname = ", ifname

close, 2
openr, 2, ifname
readf, 2, outval, xout, yout
close, 2

if globalmin GT MIN(outval,npoints) then globalmin = MIN(outval,npoints)
if globalmax LT MAX(outval,npoints) then globalmax = MAX(outval,npoints)

endfor

print, "Found global ymin, ymax = ",globalmin,globalmax

; -------------------------------------------------

print, "Started making globally autoscaled images for: ", namearray[index]

for iplane = nplane_begin, nplane_end do begin

unit = string(iplane,format='(i4.4)')
ifname = namearray[index]+unit
ofname = ifname+".ps"

print, "reading ifname = ", ifname

close, 2

openr, 2, ifname

readf, 2, outval, xout, yout

close, 2


print, "writing ofname = ", ofname


; Rescale the color range to fall between 2 and 254
outval = (outval - globalmin) / (globalmax-globalmin)
outval = outval * 252.0
outval = outval + 2.0

xmax = MAX(xout,npoints)
xmin = MIN(xout,npoints)

ymax = MAX(yout,npoints)
ymin = MIN(yout,npoints)


; This is to provide a good aspect ratio of the data when displayed on
; a screen
monitoraspect = 3.5/3.0 ; for a normal monitor
;monitoraspect = 16.0/9.0 ; for a widescreen monitor

; Postscript files have arbitrary scaling which is usually set near the top. 
; Here we provide a number which is good on a 15" monitor
; The 1.0 determines how much screen real estate is used.
scaling = 8.0/((xmax-xmin)*(ymax-ymin))

; This is a good dimensioning if one has a filled circle
; device,filename=ofname, xsize=20, ysize=20, xoffset=5, yoffset=7, $
;   /COLOR, BITS=8

; This is a good dimensioning if one has a half circle
; device,filename=ofname, xsize=12.5, ysize=25, xoffset=3.5, yoffset=1.5, $
;     /COLOR, BITS=8

; This is a good dimensioning if one has a 30 degree wedge about the equator
;device,filename=ofname, xsize=10.5, ysize=20.5, xoffset=1.5, yoffset=1.5, $
;    /COLOR, BITS=8

; This prescription replaces all of the above
device,filename=ofname, xsize=scaling*monitoraspect*(xmax-xmin), ysize=scaling*(ymax-ymin), xoffset=1.5, yoffset=1.5, /COLOR, BITS=8



!p.position=[0.035,0.0,1.0,0.965]

; Use this for putting about 30 black contours
; CONTOUR, outval, xout, yout, NLEVELS=30

; Use this for putting about 30 colored contours with colors and no filling.
; CONTOUR, outval, xout, yout, NLEVELS=30, $
;    C_COLOR=INDGEN(30)*8+2


; Use this for putting about 30 colored contours with colors and use filling.
; CONTOUR, outval, xout, yout, NLEVELS=30, /FILL, $
;    C_COLOR=INDGEN(30)*8+2

; Use this for putting about 256 colored contours with color fills; i.e.
; an 8 bit raster plot.
CONTOUR, outval, xout, yout, LEVELS=FINDGEN(252)+2,  $
  /FILL, /IRREGULAR, CHARSIZE=0.1


;colorbar,ncolors=255,Range=[globalmin,globalmax], $
;  format='(f6.2)',vertical=1,left=1,position=[0.0,0.0,0.03,0.965], $
;  charsize=0.8, TITLE=titlearray[index]


; VELOVECT, velx, vely, x_cart, y_cart, LENGTH = 0.6, /OVERPLOT, $
;   DOTS=0, MISSING = 1000.0


device, /close


endfor

print, "Finished making globally autoscaled images for: ", namearray[index]

; --------------------------------------------------

endfor

; End of the major loop that loops over flow variables whose name begins
; with the starting tag in "namearray".

; ----------------------------------------------------------------------

end
