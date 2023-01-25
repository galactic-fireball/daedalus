from astropy.io import fits
import astropy.units as u
import numpy as np
import pathlib
import plotly.subplots
import sys


LINES = [
    ('FeVI_1', 1.01089*u.um), # 75 eV
    ('FeXIII', 1.07462*u.um), # 331 eV
    ('SiX_1', 1.43008*u.um), # 351 eV
    ('SiXI_1', 1.93446*u.um), # 523 eV
    ('SiVI_1', 1.96247*u.um), # 167 eV
    ('AlIX_2', 2.04444*u.um), # 285 eV
    ('CaVIII_2', 2.32117*u.um), # 127 eV
    ('SiVII_2', 2.48071*u.um), # 205 eV
    ('SiIX_2', 2.580*u.um),
    ('AlV_2', 2.9045*u.um), # 120 eV
    ('MgVIII_3', 3.027950*u.um),
    ('CaIV_3', 3.20610*u.um), # 51 eV
    ('AlVI_3', 3.65932*u.um), # 154 eV
    ('SiIX_3', 3.92820*u.um), # 304 eV
    ('MgIV_4', 4.48712*u.um), # 80 eV
    ('ArVI_4', 4.52800*u.um), # 75 eV
    ('MgVII_5', 5.503*u.um), # 187 eV
    ('MgV_5', 5.60700*u.um), # 109 eV
    ('ArIII_8', 8.9910*u.um),
    ('MgVII_9', 9.0090*u.um),
]



def dered(wave, z):
    return wave / (1 + z)


def plot_plotly(wave, spec, lines=None, out=None):
    fig = plotly.subplots.make_subplots(rows=1, cols=1)
    fig.add_trace(plotly.graph_objects.Scatter(x=wave, y=spec))

    fig.update_layout(
        yaxis_title='Flux ({})'.format(spec.unit),
        xaxis_title='Wavelength ({})'.format(wave.unit),
        title=out.stem if out else '',
        hovermode='x',
        template='plotly_white'
    )

    fig.update_xaxes(
        range=(np.nanmin(wave.value), np.nanmax(wave.value)),
        constrain='domain'
    )
    fig.update_yaxes(
        range=(np.nanmin(spec.value)-.3, np.nanmax(spec.value)+.3),
        constrain='domain'
    )

    if lines:
        for name, line in lines:
            fig.add_vline(x=line.to(wave.unit).value, line_dash='dash', line_color='red', annotation_text=name)

    if out:
        fig.write_html(out, include_mathjax="cdn")
    else:
        fig.show()


def plot_fits(fits_file, z=0.0):
    hdu = fits.open(fits_file)

    wave = hdu[1].data['WAVE'] * u.Unit(hdu[1].header['WAVEUNIT'])
    wave = dered(wave, z)
    spec = hdu[1].data['SPEC'] * u.Unit(hdu[1].header['SPECUNIT'])
    plot_plotly(wave, spec, lines=LINES, out=fits_file.parent.joinpath('%s_aperture.html'%fits_file.stem))


def main():
    if len(sys.argv) < 2:
        raise Exception('expected argument')

    fits_file = pathlib.Path(sys.argv[1]).resolve()
    z = float(sys.argv[2]) if len(sys.argv) > 2 else 0.0

    plot_fits(fits_file, z=z)


if __name__ == '__main__':
    main()
