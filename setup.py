import setuptools

setuptools.setup(
    name="cheby",
    version="0.0.1",
    packages=setuptools.find_packages(),
    package_data={'cheby_checker': ['../../dev_data/*']},
    include_package_data=True
)
