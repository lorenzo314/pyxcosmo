; IDL Version 8.4 (linux x86_64 m64)
; Journal File for lfacciol@sappcr106
; Working directory: /dsm/sappccosmos0/work/lfacciol/icosmo/pyxcosmo/test_data/data_from_andrea/grid_like_res/grid_omega_m_sigma8_theo_3d_scatter_all
; Date: Thu Mar 17 14:46:04 2016
 
x=mrdfits('grid.fits')
;MRDFITS: Null image, NAXIS=0
x=mrdfits('grid.fits',1)
;MRDFITS: Binary table.  3 columns by  1 rows.
help, x, /str
index = indgen(15)*6+8
minmax(x.x)
;      0.11500000      0.34500000
minmax(x.x[index])
;      0.13339999      0.32659999
andreadir = icosmo + '/pyxcosmo/test_data/data_from_andrea/'
; % Variable is undefined: ICOSMO.
andreadir = getenv('ICOSMO') + '/pyxcosmo/test_data/data_from_andrea/'
data=mrdfits(andreadir + 'grid_andrea.fits')
;MRDFITS: Null image, NAXIS=0
omega_m=data['X']
; % Type conversion error: Unable to convert given STRING to Long64.
sigma8=data['Y']
; % Type conversion error: Unable to convert given STRING to Long64.
omega_m=data['X']
; % Type conversion error: Unable to convert given STRING to Long64.
omega_m=data.X
; % Illegal variable attribute: X.
help, data
data=mrdfits(andreadir + 'grid_andrea.fits')
;MRDFITS: Null image, NAXIS=0
data=mrdfits(andreadir + 'grid_andrea.fits',1)
;MRDFITS: Binary table.  3 columns by  1 rows.
help, data
minmax(data.x)
;      0.13000000      0.33000001
minmax(x.x[index])
;      0.13339999      0.32659999
minmax(data.y)
;      0.43000001       1.2300000
minmax(x.y[index])
;      0.48139998       1.1786001
index = indgen(15)*6+8
