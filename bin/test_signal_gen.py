from actsims import signal
from orphics import io, maps
from pixell import enmap

def make_plots(prefix, emaps):
    for idx, cmb_type in enumerate(['I', 'Q', 'U']):
        file_name = "{}_{}.png".format(prefix, cmb_type)
        io.high_res_plot_img(emaps[idx], file_name, down=3)

SG      = signal.SignalGen()
emaps   = SG.get_signal_sim('s15', 'pa3', 'deep56', 'f090', 0, 0)
cmbmaps = SG.get_cmb_sim('s15', 'pa3', 'deep56', 'f090', 0, 0)
fgmaps  = SG.get_fg_sim('s15', 'pa3', 'deep56', 'f090', 0, 0) 

make_plots('cmb', cmbmaps)
make_plots('cmb_fg', emaps)
make_plots('fg', fgmaps)
make_plots('diff_fg', emaps-cmbmaps)

template   = emaps[0,:100, :100]

shape, wcs = template.shape, template.wcs
SG      = signal.SignalGen(extract_region_shape=shape, extract_region_wcs=wcs)
cmbmaps = SG.get_cmb_sim('s15', 'pa3', 'deep56', 'f090', 0, 0)
make_plots('cmb_ext1', cmbmaps)

extract_region = enmap.zeros(shape, wcs)
SG      = signal.SignalGen(extract_region=extract_region)
cmbmaps = SG.get_cmb_sim('s15', 'pa3', 'deep56', 'f090', 0, 0)
make_plots('cmb_ext2', cmbmaps)
