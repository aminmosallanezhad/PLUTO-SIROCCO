;+
;
; NAME: WRITE_VTK
;
; AUTHOR:  Andrea Mignone (andrea.mignone@unito.it)
;
; PURPOSE:  Create a VTK data file out of a PLUTO binary file assuming
;           the grid information previously obtained with pload. 
;           A scalar VTK file is created if only one argument is given;
;           otherwise a VTK vector file is created if 3 arguments are given.
;   
; SYNTAX: WRITE_VTK, QX[, QY, QZ][,FILENAME=string][,VARNAME=string]
;                      [/APPEND][REBIN_ARRAY=REBIN_ARRAY]
;
; ARGUMENTS:
;
;   QX    = a 3D single-precision (float) binary data to be written
;           in the VTK file.
;
;   QY,QZ = the y- and z- components of the vector fields to be 
;           written.
;
;   FILENAME = the name of the output VTK file. Default is "data.vtk". 
;
;   VARNAME  = the name of the variable appearing in the VTK header.
;              Default is "IDL_var".
;
; KEYWORDS:
;
;  /APPEND:       open file and append data
;  REBIN_ARRAY:   an array of integers giving the dimension of the rebinned
;                 array. For instance, setting rebin_array=[nx1/2,nx2/2,nx3/2]
;                 will produce an array with halved size.
;                 The original array and coordinates will be replaced by
;                 new ones.
;
; LAST MODIFIED:
;
;  Dec 19, 2013 by A. Mignone (andrea.mignone@unito.it)
;   
;- 
PRO WRITE_VTK, qx, qy, qz, filename=filename,varname=varname,append=append,$
               rebin_array=rebin_array

  COMMON PLUTO_GRID
  COMMON PLUTO_RUN

  IF (NOT KEYWORD_SET(filename)) THEN filename = "data.vtk"

; --------------------------------
; Check if scalars or vector 
; must be written
; --------------------------------

  write_scalar = 0
  write_vector = 0

  sqx = size(qx)
  sqy = size(qy)
  sqz = size(qz)

  nx1 = sqx[1]
  nx2 = sqx[2]
  nx3 = sqx[3]

  IF ((sqx(1) EQ sqy(1)) AND (sqx(1) EQ sqz(1))) THEN $
  write_vector = 1 ELSE write_scalar = 1

; ---------------------------------------
;   Rebin data files if necessary
; ---------------------------------------

  IF (KEYWORD_SET(rebin_array)) THEN BEGIN
    IF (write_scalar) THEN qx = REBIN(qx, rebin_array)
    IF (write_vector) THEN BEGIN
      qx = REBIN(qx, rebin_array)
      qy = REBIN(qy, rebin_array)
      qz = REBIN(qz, rebin_array)
    ENDIF
   
    nx1b = rebin_array[0] & x1b = FLTARR(nx1b) & dx1b = FLTARR(nx1b)
    nx2b = rebin_array[1] & x2b = FLTARR(nx2b) & dx2b = FLTARR(nx2b)
    nx3b = rebin_array[2] & x3b = FLTARR(nx3b) & dx3b = FLTARR(nx3b)

   ; ----------------------------------------------------------------
   ;  New grid must have interfaces coincident with the old grid:
   ;   
   ;   i  =  | 0 |  1  |  2  |  3  |  4  |   5   |   6   | 7 |  8  | ..
   ;   ib =  |         0           |             1           |  ...  
   ; ----------------------------------------------------------------

    r = nx1/REBIN_ARRAY[0]
print,nx1, nx2, nx3, r
    FOR i = 0, nx1b-1 DO BEGIN     
      xl = x1[r*i]         - dx1[r*i]*0.5
      xr = x1[r*(i+1)-1] + dx1[r*(i+1)-1]*0.5
      x1b[i]  = 0.5*(xl + xr)
      dx1b[i] = xr - xl
    ENDFOR
   
    r = nx2/nx2b
    FOR i = 0, nx2b-1 DO BEGIN
      xl = x2[r*i]       - dx2[r*i]*0.5
      xr = x2[r*(i+1)-1] + dx2[r*(i+1)-1]*0.5
      x2b[i]  = 0.5*(xl + xr)
      dx2b[i] = xr - xl
    ENDFOR

    r = nx3/nx3b
    FOR i = 0, nx3b-1 DO BEGIN
      xl = x3[r*i]       - dx3[r*i]*0.5
      xr = x3[r*(i+1)-1] + dx3[r*(i+1)-1]*0.5
      x3b[i]  = 0.5*(xl + xr)
      dx3b[i] = xr - xl
    ENDFOR
    x1  = x1b  & x2  = x2b  & x3  = x3b
    dx1 = dx1b & dx2 = dx2b & dx3 = dx3b
    nx1 = nx1b & nx2 = nx2b & nx3 = nx3b
  
  ENDIF; Rebin Array

; ---------------------------------------
;      compute node coordinates 
; ---------------------------------------

  ioff = 1
  joff = 1
  koff = 1
  IF (NX3 EQ 1) THEN koff = 0

  x1n = FLTARR(NX1+ioff)
  x2n = FLTARR(NX2+joff)
  x3n = FLTARR(NX3+koff)

  FOR i = 0, NX1-1 DO x1n(i) = x1(i)-0.5*dx1(i)
  x1n(i) = x1n(i-1) + dx1(i-1)

  FOR j = 0, NX2-1 DO x2n(j) = x2(j)-0.5*dx2(j)
  x2n(j) = x2n(j-1) + dx2(j-1)

  IF (koff EQ 1) THEN BEGIN;  we're in 3D
    FOR k = 0, NX3-1 DO x3n(k) = x3(k)-0.5*dx3(k)
    x3n(k) = x3n(k-1) + dx3(k-1)
  ENDIF ELSE BEGIN
    x3n(0) = 0.0;
  ENDELSE

; --------------------------------
;   obtain the variable name
; --------------------------------

  IF (NOT KEYWORD_SET(varname)) THEN varname = "IDL_data"
; varname = "IDL_var"; ** default name **
; n = strpos(filename,".")
; IF (n GT 1) THEN varname = strmid(filename,0,n)

; --------------------------------
;   open file for writing
; --------------------------------

  sn1 = strcompress(string(NX1 + ioff,format='(i4)'))
  sn2 = strcompress(string(NX2 + joff,format='(i4)'))
  sn3 = strcompress(string(NX3 + koff,format='(i4)'))

  IF (KEYWORD_SET (APPEND)) THEN BEGIN
    OPENW, U, filename, /GET_LUN, /APPEND
  ENDIF ELSE BEGIN
    OPENW, U, filename, /GET_LUN

 ; ------------------------------------------------
 ;   write header section and coordinates
 ; ------------------------------------------------

   PRINTF,U,"# vtk DataFile Version 2.0"
   PRINTF,U,"Generated by IDL with WRITE_VTK (A. Mignone)"
   PRINTF,U,"BINARY"

   printf,U,"DATASET RECTILINEAR_GRID"
   printf,U,"DIMENSIONS "+sn1+"  "+sn2+"  "+sn3
;   printf,U,''

   printf,U,"X_COORDINATES "+sn1+" float"
   WRITEU,U,SWAP_ENDIAN(FLOAT(x1n),/swap_if_little_endian)
   printf,U,''
   printf,U,"Y_COORDINATES "+sn2+" float"
   WRITEU,U,SWAP_ENDIAN(FLOAT(x2n),/swap_if_little_endian)
   printf,U,''
   printf,U,"Z_COORDINATES "+sn3+" float"
   WRITEU,U,SWAP_ENDIAN(FLOAT(x3n),/swap_if_little_endian)

   scrh = double(ULONG64(NX1)*ULONG64(NX2)*ULONG64(NX3))
   printf,U,''
   printf,U,"CELL_DATA "+strcompress(string(scrh,format='(I16)'))
   printf,U,''
 ENDELSE

; ----------------------------------------------
;                 VTK Scalar
; ----------------------------------------------

 IF (write_scalar) THEN BEGIN
   print,"> WRITE_VTK: writing scalar VTK data [",filename,"]..."
   printf,U,"SCALARS "+varname+" float"
   printf,U,"LOOKUP_TABLE default"
   WRITEU,U,SWAP_ENDIAN(float(qx),/swap_if_little_endian)
 ENDIF

; ----------------------------------------------
;                 VTK Vector
;
;  we need to create data structure like 
;  [v1,v2,v3] (i,j,k). We reform 2D slices on 
;  1-D column vector of length NX1*NX2 and then
;  join them together with [V1, V2, V3].
;  This is *MUCH* faster than creating 
;  point-by-point triplets like 
;  q1d = [qx(i,j,k), qy(i,j,k), qz(i,j,k)]
; ----------------------------------------------

 IF (write_vector) THEN BEGIN
   print,"> WRITE_VTK: writing vector VTK data [",filename,"]..."
   printf,U,"VECTORS "+varname+" float"

   FOR k = 0, NX3-1 DO BEGIN
     qx2d = REFORM(qx(*,*,k))
     qy2d = REFORM(qy(*,*,k))
     qz2d = REFORM(qz(*,*,k))
     AX = REFORM(qx2d, 1, ULONG64(NX1)*ULONG64(NX2))
     AY = REFORM(qy2d, 1, ULONG64(NX1)*ULONG64(NX2))
     AZ = REFORM(qz2d, 1, ULONG64(NX1)*ULONG64(NX2))
     C = [AX, AY, AZ] 
     WRITEU,U,SWAP_ENDIAN(float(C),/swap_if_little_endian)
   ENDFOR

;
;  OLD, MUCH SLOWER VERSION
;
;   FOR k = 0, NX3-1 DO BEGIN
;   FOR j = 0, NX2-1 DO BEGIN
;   FOR i = 0, NX1-1 DO BEGIN
;     q1d = [qx(i,j,k), qy(i,j,k), qz(i,j,k)]
;     WRITEU,U,SWAP_ENDIAN(float(q1d),/swap_if_little_endian)
;   ENDFOR
;   ENDFOR
;   ENDFOR
 ENDIF

 CLOSE,U
END

