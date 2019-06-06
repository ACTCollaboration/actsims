from actsims import signal, simgen
from orphics import io, maps
from pixell import enmap

def make_plots(prefix, emaps):
    for idx, cmb_type in enumerate(['I', 'Q', 'U']):
        file_name = "{}_{}.png".format(prefix, cmb_type)
        io.high_res_plot_img(emaps[idx], file_name, down=3)

# our test data set

version = 'v4.0_mask_version_mr3c_20190215_pickupsub_190301'
season, array, patch, freq = ('s13', 'pa1', 'deep1', 'f150')
shape, wcs = simgen.get_default_geometry(version, season, patch, array, freq)


SG      = signal.SignalGen()
emaps   = SG.get_signal_sim(season, patch, array, freq, 0, 0, oshape=shape, owcs=wcs)
cmbmaps = SG.get_cmb_sim(season, patch, array, freq, 0, 0, oshape=shape, owcs=wcs)
fgmaps  = SG.get_fg_sim(season, patch, array, freq, 0, 0, oshape=shape, owcs=wcs) 

make_plots('cmb', cmbmaps)
make_plots('cmb_fg', emaps)
make_plots('fg', fgmaps)
make_plots('diff_fg', emaps-cmbmaps)

template   = cmbmaps[0]
template   = enmap.pad(template, 100)


shape, wcs = template.shape, template.wcs
SG      = signal.SignalGen(extract_region_shape=shape, extract_region_wcs=wcs)
cmbmaps = SG.get_cmb_sim('s15', 'pa3', 'deep56', 'f090', 0, 0)
make_plots('cmb_ext1', cmbmaps)

extract_region = enmap.zeros(shape, wcs)
SG      = signal.SignalGen(extract_region=extract_region)
cmbmaps = SG.get_cmb_sim('s15', 'pa3', 'deep56', 'f090', 0, 0)
make_plots('cmb_ext2', cmbmaps)
