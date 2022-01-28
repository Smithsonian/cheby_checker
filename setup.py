import setuptools

setuptools.setup(
    name="cheby",
    version="0.0.1",
    install_requires=['numpy', 'pytest', 'astropy', 'astropy_healpix', 'astroquery', 'scipy', 'novas', 'novas_de405', 'jplephem'],
    packages=setuptools.find_packages()
)
