#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 12 15:32:43 2021

@author: Xuheng Ding
"""
import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as pyfits
import glob
import pickle

Bla_type = 'BL_Lac'
Bla_type = 'FSRQ' #k 27

filename = 'model_result/{0}_host_result.txt'.format(Bla_type)
#Read files:
f = open(filename,"r")
string = f.read()
results = string.split('\n')   # Split in to \n
ids_file = glob.glob('../../data/{1}/*-{0}.fits'.format('I', Bla_type))
ids = [ids_file[i].split('/')[-1].split('_')[0] for i in range(len(ids_file))]

for k in range(27,28):
# for k in [23]:    
    bands = ['G', 'R', 'I' ,'Z' , 'Y']
    mags = []
    chisqs = []
    residuals = []
    for band in bands:
        line = [i for i in range(len(results)) if ids[k]+'_HSC'+'-'+band in results[i]]
        if line != []:
            fit_result = results[line[0]].split("'magnitude': ")[1]
            if len(fit_result)>5:
                mags.append(float(fit_result[:-3]))
            else:
                mags.append(30.0)
        else:
            mags.append(30.0)
        pklfile_name = 'model_result/'+'{0}_HSC-{1}.pkl'.format(ids[k], band)
        if glob.glob(pklfile_name) != []:
            fitting_run_class = pickle.load(open(pklfile_name,'rb'))
            chisqs.append(fitting_run_class.reduced_Chisq)
            residuals.append(fitting_run_class.fitting_specify_class.kwargs_data['image_data'] - fitting_run_class.image_ps_list[0])#/fitting_run_class.fitting_specify_class.kwargs_data
        else:
            chisqs.append(99)
            residuals.append(np.zeros((60,60))-10 )
    # mags = [round(mags[i],2) for i in range(len(mags))]
    # chisqs = [round(chisqs[i],2) for i in range(len(chisqs))]
    # if np.max(chisqs) < 10 and np.min(mags)>0:
    if np.max(chisqs) > 0:
        print(ids[k], k) 
        print(mags)
        mag_err = [0.4, 0.4, 0.2, 0.4, 0.4]
        # for i in range(len(bands)):
        #     band = bands[i]
        #     print("Chisq", chisqs[i])
        #     print("min residual", np.sum(residuals[i]<-0.5))
        #     pklfile_name = 'model_result/'+'{0}_HSC-{1}.pkl'.format(ids[k], band) 
        #     fitting_run_class = pickle.load(open(pklfile_name,'rb'))
        #     fitting_run_class.plot_final_qso_fit(target_ID = ids[k]+'-'+band)
        # fnu = [10 ** ((mags[i]-25)/(-2.5)) for i in range(len(mags))]
        # fnu_up = [10 ** ((mags[i]-mag_err[i]-25)/(-2.5)) for i in range(len(mags))]
        # fnu_dw = [10 ** ((mags[i]+mag_err[i]-25)/(-2.5)) for i in range(len(mags))]
        # fnu_err = [(fnu_up[i]-fnu_dw[i])/2 for i in range(len(mags))]
        # lam = np.array([4798.2, 6218.4, 7727.0, 8908.2, 9775.1])
        # mags = np.array(mags)
        # flambda = 10**(-0.4*(mags+2.402+5.0*np.log10(lam))) * 10**19
        # lam = lam/10000.
        # plt.scatter(lam, flambda, c='r', zorder =100)
        # plt.show()
        
        ID_org = ids[k]
        ID = str("".join(filter(str.isdigit, ID_org)))
        
        mag = mags
        fnu = [10 ** ((mag[i]-25)/(-2.5)) for i in range(len(mag))]
        fnu_up = [10 ** ((mag[i]-mag_err[i]-25)/(-2.5)) for i in range(len(mag))]
        fnu_dw = [10 ** ((mag[i]+mag_err[i]-25)/(-2.5)) for i in range(len(mag))]
        fnu_err = [(fnu_up[i]-fnu_dw[i])/2 for i in range(len(mag))]
        for i in range(len(bands)):
            if chisqs[i]>2:
                if np.sum(residuals[i]<-0.5)>25:
                    fnu_err[i] = 150
            if mags[i] > 23:
                fnu_err[i] = 150
                
        print("Start Photo-z fitting:")
        print("mags:", mag)
        print("errs:", fnu_err)
        #Generate the cat file as SED photo-z input
        sed_fnu_name = 'SED_file/' + ID+'_fnu.cat'
        write_file = open(sed_fnu_name,'w') 
        write_file.write("# id F314 E314 F315 E315 F316 E316 F317 E317 F318 E318\n")
        write_file.write(ID) 
        for i in range(len(fnu)):
            write_file.write(" {0} {1}".format(round(fnu[i],4), round(fnu_err[i],4)))
        write_file.write("\n")    
        # write_file.write("113245 0.001 10 28.195 2.601 106.401 9.814 133.879 12.348 188.171 17.356")    
        write_file.close()
        
        #Take takahiro's template
        # print([(round(fnu[i],3), round(fnu_err[i],3)) for i in range(len(mag))])
        # ## Full spectral fitting of a galaxy;
        # 
        import os
        from astropy import __version__ as asver
        asver
        # In[3]:
        # https://github.com/mtakahiro/gsf/tree/version1.4
        import gsf
        print('gsf version',gsf.__version__)
        
        from gsf.function import get_input
        from gsf.gsf import run_gsf_template
        from gsf.plot_sed_logA import plot_sed, plot_corner_physparam_frame, plot_corner_physparam_summary
        from gsf.plot_sfh_logA import plot_sfh
        # ### Setup gsf
        
        # In[6]:
        # Initial setup for gsf.
        # Auto load input dictionary;
        inputs = get_input()
        
        # change Z;
        # Flag;
        fplt = 0
        inputs['DIR_TEMP'] = './SED_templates/'
        # Output directory;
        inputs['DIR_OUT'] = './SED_file/'
        
        # # If templates exit already, then let's save time.
        # # (But if you changed metallicity range or age pixels, fplt needs to be 0.)
        if os.path.exists('%s/spec_all.asdf'%inputs['DIR_TEMP']):
            fplt = 1
        
        inputs['ID'] = ID
        
        # Initial guess of redshift, or the true value if known. 
        # We will later do redshift fit later, though.
        inputs['ZGAL'] = 0.4        #!!!
        
        # Redshift as a free parameter?
        inputs['ZMC'] = 1
        
        # Metallicity range, in logZsun;
        inputs['ZMIN'] = 0.
        inputs['ZMAX'] = 0.
        inputs['DELZ'] = 0.
        # You can fix metallicity;
        inputs['ZFIX'] = 0.0
        # Dust attenuation
        inputs['AVMIN'] = 0.
        inputs['AVMAX'] = 0.
        # You can fix;
        inputs['AVFIX'] = 0.0
        # Templates;
        inputs['BPASS'] = 0
        inputs['AGE'] = '0.01,0.03,0.1,0.3,0.5,0.7,1.0,2.0,3.0'
        # You can fix age;
        #inputs['AGEFIX'] = '0.3' # '0.1,0.3,0.5'
        inputs['NIMF'] = 1 #
        # Data;
        DIR_EXTR = './'
        spec_file = '' #'./l3_nis_G150C_s00003_1d_cont_fnu.txt'
        inputs['DIR_EXTR'] = DIR_EXTR
        #inputs['SPEC_FILE'] = spec_file
        inputs['SPEC_FILE'] =  None # If no spec, then leave this None, or '', so gsf does broadband only SED fitting.
        inputs['DIR_FILT'] = '/Users/Dartoon/Astro/Packages/gsf_v1.3/gsf/example/filter/'
        #inputs['DIR_FILT'] = '/Users/Dartoon/Astro/Packages/gsf/gsf/example/filter/'
        inputs['CAT_BB'] = sed_fnu_name
        
        # Filters;
        # Each number corresponds to EAZY's filter ids. See also filter/filt_Sep20.lst
        # These numbers need to be found in inputs['CAT_BB'] file.
        inputs['FILTER'] = '314,315,316,317,318'
        # Morphology convolution; Necessary for NIRISS spectra;
        #filt = 'f200w'
        #inputs['MORP'] = 'moffat'
        #inputs['MORP_FILE'] = './l3_nis_f200w_G150C_s00003_moffat.txt'
        
        # MCMC part;
        inputs['NCPU'] = 1 # For notebook, somehow multiprocessing causes error. So set to 1.
        inputs['NMC'] = 1000 # NMC for the main SED fit
        inputs['NMCZ'] = 30 # NMC for the redshift fit
        
        
        # Visual inspection;
        # Set to 0 (False), as Notebook cannot show actively iterating plot;
        inputs['ZVIS'] = 0
        
        # Emission line masking;
        #LW = [3727, 4341, 4861, 4960, 5008, 6563, 6717, 6731] # in AA, rest.
        #inputs['LINE'] = LW
        
        # Initial fit:
        inputs['FNELD'] = 0
        
        # In[7]:
        # Then, run template generate function;
        mb = run_gsf_template(inputs, fplt=0)
        fplt = 1
        # In[ ]:
        
        # # You can write down the input file in an ascii file.
        # from gsf.function import write_input
        # write_input(inputs, file_out='gsf.input')
        
        # In[ ]:
        
        
        # Do a quick fit at z=z_guess;
        mb.zprev = mb.zgal
        out, fm_tmp, xm_tmp = mb.quick_fit(mb.zgal, mb.Cz0, mb.Cz1)
        # ### Now, let's improve the fit by finding the true redshift;
        
        # In[ ]:
        # Preparing Fitting Spectral Template from the library generated above;
        # Here, we use 5 templates for find redshift;
        dict = mb.read_data(mb.Cz0, mb.Cz1, mb.zgal)
        ages = [0.01,0.03,0.1,0.3,1.0]
        ntmp = len(ages)
        
        for nn in range(ntmp):
            # For simplicity, no dust attenuation (Av=0), Z fixed to solar (Z=0).
            flux_all, wave_all = mb.fnc.get_template(mb.lib_all, Amp=1.0, T=ages[nn], Av=0.0, Z=0.0, zgal=mb.zgal)
            
            con_tmp = (1000 < wave_all / (1.+mb.zgal)) & (wave_all / (1.+mb.zgal) < 60000)
        
            # Don't forget to blueshift the template.
            xm_tmp = wave_all[con_tmp] / (1.+mb.zgal)
            fm_tmp = flux_all[con_tmp]
        
            if nn == 0:
                fm_tmps = np.zeros((ntmp,len(xm_tmp)),'float')
        
            fm_tmps[nn,:] = fm_tmp[:]
        
        # In[ ]:
        
        
        # Then, run redshift fitting; 
        # dict : dictionary that includes a lot of things, including data.
        # zliml, zlimu : Redshift search range, lower and upper limits.
        
        # This should not be too small, if z-distribution is used as prior.
        delzz = 0.1
        zspace, chi2s = mb.search_redshift(dict, xm_tmp, fm_tmps, zliml=0.1, zlimu=1., delzz=delzz)
        # In[ ]:
        
        
        # Plot;
        plt.plot(zspace,chi2s[:,1])
        plt.ylabel('$\chi^2$',fontsize=18)
        plt.xlabel('$z$',fontsize=18)
        plt.title('Redshift Fitting Result')
        #plt.ylim(0,120)
        # In[ ]:
        # Since the result above looked suspicious:
        if False:
            # Get z at the chi2 minimum.
            izfit = np.argmin(chi2s[:,1])
            zfit = zspace[izfit]
            print('zfit is %.2f'%(zfit))
        
        # In[ ]:
        # Use chi2 as a prior
        # User can provide phot-z prob by EAZY too.
        prior = {}
        prior['z'] = zspace
        prior['chi2'] = chi2s[:,1]
        prior['prob'] = np.exp(-0.5 * prior['chi2'])
        
        #prior
        #chi2s.shape
        # ## Since phot-z is not constraining, let's set a prior..
        
        # In[ ]:
        # Since the result above looked suspicious,
        # use an arbitrary, flat prior
        if True:
            # Or define a new prior:
            zspace_tmp = np.arange(0,13,0.01)
            chi2s_tmp = zspace_tmp * 0 + 99
            con_tmp = (zspace_tmp>0.2) & (zspace_tmp<0.4)
            chi2s_tmp[con_tmp] = 1.0
        
            prior = {}
            prior['z'] = zspace_tmp
            prior['chi2'] = chi2s_tmp
            
            plt.plot(prior['z'],np.exp(-.5 * prior['chi2']), color='cyan', linestyle='--')
            plt.xlim(0.1,0.6)
            
        # In[ ]:
        
        # prior
        
        # In[ ]:
        
        
        # Repeat the quick fit at the proposed redshift;
        #inputs['ZGAL'] = zfit
        inputs['NMCZ'] = 30
        
        # Update with a new z input
        mb = run_gsf_template(inputs, fplt=fplt)
        
        mb.zprev = mb.zgal
        out, fm_tmp, xm_tmp = mb.quick_fit(mb.zgal, mb.Cz0, mb.Cz1)
        # ### Now the result looks good
        
        # In[ ]:
        
        plt.close()
        
        # Plot the result;
        flux_all, wave_all = mb.fnc.tmp04_val(out, mb.zgal, mb.lib_all)
        
        # Template
        plt.errorbar(wave_all, flux_all, ls='-', color='b', zorder=0, label='Fit')
        
        # plot;
        plt.scatter(dict['xbb'], dict['fybb'], marker='o', c='orange', edgecolor='k', s=150, zorder=2, alpha=1, label='Broadband')
        
        if True: # Spec data;
            plt.errorbar(dict['x'], dict['fy'], yerr=dict['ey'], ls='', color='gray', zorder=1, alpha=0.3)
            plt.scatter(dict['x'], dict['fy'], marker='o', color='r',edgecolor='r', s=10, zorder=1, alpha=1, label='Spectrum')
        
        plt.xlim(1000,80000)
        #plt.ylim(0,25)
        plt.xscale('log')
        
        plt.legend(loc=2)
        plt.xlabel('Wavelength (AA)')
        plt.ylabel('Flux (MJy/sr)')
        # ### Now fit redshift in more details;
        
        # In[ ]:
        dict = mb.read_data(mb.Cz0, mb.Cz1, mb.zgal)
        
        # By usinng the bbest fit template above;
        con_tmp = () #(7000 < wave_all) & (wave_all < 25000)
        xm_tmp = wave_all[con_tmp]
        fm_tmp = flux_all[con_tmp]
        
        # Update inputs; 
        inputs['NMCZ'] = 500
        inputs['NWALKZ'] = 30
        mb.update_input(inputs)
        
        # Do you continue Redshift fit when only BB photometry is available?
        f_bb_zfit = True
        
        # This works only when spectrum is provided.
        mb.fit_redshift(dict, xm_tmp, fm_tmp, delzz=0.01, zliml=2., zlimu=2.55, ezmin=0.01, snlim=0,                 f_bb_zfit=f_bb_zfit, priors=prior)
        
        
        # In[ ]:
        # This is normalization;
        # Should be ~1, as we have already normalized the spectra to BB fluxes.
        print('Redshift 16/50/84th percentile range :', mb.z_cz)
        print(mb.Czrec0)
        print(mb.Czrec1)
        # In[ ]:
        # Take a look at z distribution
        # fig_zdist, ax1 = mb.get_zdist(f_interact=True)
        # In[ ]:
        
        # ax1.set_xlim(0, 1)
        # fig_zdist.savefig('%szprob_tmp.png'%inputs['DIR_OUT'])
        # fig_zdist.show()
        
        # In[ ]:
        from IPython.display import Image 
        Image('%szprob_tmp.png'%inputs['DIR_OUT'])
        # ### Now, run the whole SED fitting;
        
        # In[ ]:
        # No interactive fit;
        inputs['ZMC'] = 1
        inputs['ZVIS'] = 0
        inputs['NMC'] = 3000
        inputs['ZGAL'] = mb.z_cz[1]
        
        # Update inputs; 
        mb.update_input(inputs)
        # Since already z-fit done, we can skip z-fit;
        skip_fitz = True
        
        # Main;
        flag_suc = mb.main(cornerplot=True, specplot=1, sigz=1.0, ezmin=0.01, ferr=0, f_move=False, skip_fitz=skip_fitz)
        
        
        # # In[ ]:
        # # Plot SFH;
        
        # # Plot Starforming Main Sequence from Speagle+14?
        f_SFMS = True
        #%%
        try:
            plot_sfh(mb, f_comp=mb.ftaucomp, fil_path=mb.DIR_FILT,
                     inputs=mb.inputs, dust_model=mb.dust_model, DIR_TMP=mb.DIR_TMP, f_SFMS=f_SFMS, f_fill=False)
        except:
            print('\nFile is missing : _param.fits\n')
        #%%
        try:            
            plot_sed(mb, fil_path=mb.DIR_FILT,
                     figpdf=False, save_sed=True, inputs=mb.inputs, mmax=300,
                     f_fill=True, dust_model=mb.dust_model, DIR_TMP=mb.DIR_TMP, f_label=True)
        except:
            print('\nFile is missing : _param.fits\n')
            
        # In[ ]:
        # Physical parameters;
        # plot_corner_physparam_summary(mb)
        # In[ ]:
        #To estimate the AB mag for a given stellar template
        
        import matplotlib as mat
        mat.rcParams['font.family'] = 'STIXGeneral'
        import glob
        # from specutils import Spectrum1D
        # filename_p  = './SED_file/summary_{0}_PA00.fits'.format(ID)
        # filename_p = glob.glob(filename_p)[0]
        # hdul = pyfits.open(filename_p)
        # table = hdul[1].data
        # name = hdul[1].columns
        # z = table['zmc'][1]
        # smass = table['ms'][1]
        filename_p  = './SED_file/SFH_{0}_PA00_param.fits'.format(ID)
        filename_p = glob.glob(filename_p)[0]
        hdul = pyfits.open(filename_p)
        table = hdul[1].data
        name = hdul[1].columns
        age_idx = [i for i in range(len(name)) if 'T_MW' in str(name[i])][0]
        age = str(round(10**table[1][age_idx],3)) # 'AGE:' Gyr
        z_idx = [i for i in range(len(name)) if 'zmc' in str(name[i])][0]
        z = str(round(table[1][z_idx],2))
        mel_idx = [i for i in range(len(name)) if 'Z_MW' in str(name[i])][0]
        mel = str(round(10**table[1][mel_idx],3)) # 'Mel:' Z*/Z_sun
        
        spec1d = pyfits.open("./SED_file/gsf_spec_{0}.fits".format(ID))  
        name_spec = spec1d[1].columns
        table_spec = spec1d[1].data
        
        plt.figure(figsize=(10, 6))
        f = table_spec['f_model_50']            #This is the best-fit spectra (i.e. 50% by MCMC).
        wave = table_spec['wave_model']/10000.
        f_max = f.max()
        plt.plot(wave, f)
        
        f = table_spec['f_model_16']
        wave = table_spec['wave_model']/10000.
        plt.plot(wave, f)
        
        f = table_spec['f_model_84']
        wave = table_spec['wave_model']/10000.
        plt.plot(wave, f)
        
        HSC_g_fil = np.loadtxt('/Users/Dartoon/Astro/Packages/gsf_v1.3/gsf/example/filter/314.fil')
        HSC_g_fil[:,2] = HSC_g_fil[:,2]/HSC_g_fil[:,2].max()* f_max/6
        plt.plot(HSC_g_fil[:,1]/10000., HSC_g_fil[:,2], label='HSC_G response')
        
        HSC_r_fil = np.loadtxt('/Users/Dartoon/Astro/Packages/gsf_v1.3/gsf/example/filter/315.fil')
        HSC_r_fil[:,2] = HSC_r_fil[:,2]/HSC_r_fil[:,2].max()* f_max/6
        plt.plot(HSC_r_fil[:,1]/10000., HSC_r_fil[:,2], label='HSC_R response')
        
        HSC_i_fil = np.loadtxt('/Users/Dartoon/Astro/Packages/gsf_v1.3/gsf/example/filter/316.fil')
        HSC_i_fil[:,2] = HSC_i_fil[:,2]/HSC_i_fil[:,2].max()* f_max/6
        plt.plot(HSC_i_fil[:,1]/10000., HSC_i_fil[:,2], label='HSC_I response')
        
        HSC_z_fil = np.loadtxt('/Users/Dartoon/Astro/Packages/gsf_v1.3/gsf/example/filter/317.fil')
        HSC_z_fil[:,2] = HSC_z_fil[:,2]/HSC_z_fil[:,2].max()* f_max/6
        plt.plot(HSC_z_fil[:,1]/10000., HSC_z_fil[:,2], label='HSC_Z response')
        
        HSC_y_fil = np.loadtxt('/Users/Dartoon/Astro/Packages/gsf_v1.3/gsf/example/filter/318.fil')
        HSC_y_fil[:,2] = HSC_y_fil[:,2]/HSC_y_fil[:,2].max()* f_max/6
        plt.plot(HSC_y_fil[:,1]/10000., HSC_y_fil[:,2], label='HSC_Y response')
        
        yerr_l = np.array(fnu_err)/np.array(fnu)
        
        lam = np.array([4798.2, 6218.4, 7727.0, 8908.2, 9775.1])
        mag = -2.5*np.log10(fnu) + 25
        flambda = 10**(-0.4*(mag+2.402+5.0*np.log10(lam))) * 10**19
        lam = lam/10000.
        
        plt.scatter(lam, flambda, c='r', zorder =100)
        plt.errorbar(lam, flambda, yerr=yerr_l*flambda,fmt='.',color='gray',markersize=1, zorder =90)
        plt.xlim(0.25, 3)
        plt.ylim(-50,700)
        xmin, xmax, ymin, ymax = plt.axis()
        plt.title("photo-z of "+ID_org, fontsize=25)
        plt.text( (xmax-xmin)*0.7, (ymax-ymin)*0.6, 'stellar population with:', fontsize=17)
        plt.text( (xmax-xmin)*0.8, (ymax-ymin)*0.53, 'z={0}'.format(z), fontsize=17)
        plt.text( (xmax-xmin)*0.8, (ymax-ymin)*0.45, 'age={0:.3f} Gyr'.format(float(age)), fontsize=17)
        plt.text( (xmax-xmin)*0.8, (ymax-ymin)*0.37, 'Z$_*$/Z$_\odot$={0}'.format(mel), fontsize=17)
        plt.legend(prop={'size':15}, ncol=3, loc = 1)
        plt.tick_params(labelsize=15)
        plt.xlabel("um",fontsize=27)
        plt.ylabel(r"f$_\lambda$ 10$^{\rm -19}$ erg/s/cm$^2$/$\AA$",fontsize=27)
        plt.savefig('./SED_file/photo-z_' + ID_org + '.pdf')
        plt.show()
