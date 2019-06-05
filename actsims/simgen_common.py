from soapack import interfaces as sints

def get_save_paths(model,version,coadd,season=None,patch=None,array=None,mkdir=False,overwrite=False,mask_patch=None):
    paths = sints.dconfig['actsims']

    try: assert paths['plot_path'] is not None
    except: paths['plot_path'] = "./"
    assert paths['covsqrt_path'] is not None
    assert paths['trial_sim_path'] is not None

    # Prepare output dirs
    pdir = "%s/%s/" % (paths['plot_path'] ,version)
    cdir = "%s/%s/" % (paths['covsqrt_path'] ,version)
    sdir = "%s/%s/" % (paths['trial_sim_path'] ,version)

    if mkdir:
        exists1 = util.mkdir(pdir)
        exists2 = util.mkdir(cdir)
        exists3 = util.mkdir(sdir)
        if any([exists1,exists2,exists3]):
            if not(overwrite): raise IOError
            warnings.warn("Version directory already exists. Overwriting.")

    if model=='planck_hybrid':
        assert season is None
        suff = '_'.join([model,patch,array,"coadd_est_"+str(coadd)])
    else:
        suff = '_'.join([model,season,patch,array,"coadd_est_"+str(coadd)])


    pout = pdir + suff
    cout = cdir + suff
    sout = sdir

    if mask_patch is not None:
        if mask_patch != patch:
            pout = pout+"_"+mask_patch
            cout = cout+"_"+mask_patch
            sout = sout+mask_patch+"_"

    return pout,cout,sout

