drizzle.outnx=167
drizzle.outny=167
drizzle.pixfrac=0.85
drizzle.scale=0.8
drizzle.kernel="gaussian"
drizzle non_drizzled_psf-1.fits dripsf.fits outweig = "dripsfw.fits" xsh=-0.375 ysh=-0.375
drizzle non_drizzled_psf-2.fits dripsf.fits outweig = "dripsfw.fits" xsh=-0.375 ysh=0.125
drizzle non_drizzled_psf-3.fits dripsf.fits outweig = "dripsfw.fits" xsh=0.125 ysh=-0.375
drizzle non_drizzled_psf-4.fits dripsf.fits outweig = "dripsfw.fits" xsh=0.125 ysh=0.125
drizzle non_drizzled_psf-5.fits dripsf.fits outweig = "dripsfw.fits" xsh=0.375 ysh=-0.125
drizzle non_drizzled_psf-6.fits dripsf.fits outweig = "dripsfw.fits" xsh=0.375 ysh=0.375
drizzle non_drizzled_psf-7.fits dripsf.fits outweig = "dripsfw.fits" xsh=-0.125 ysh=-0.125
drizzle non_drizzled_psf-8.fits dripsf.fits outweig = "dripsfw.fits" xsh=-0.125 ysh=0.375