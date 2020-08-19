import setuptools

def readme():
    with open("README.rst", 'r') as f:
        return f.read()


setuptools.setup(
    name='geonomics',
    # version num.: MAJOR.MINOR.PATCH
    version='0.0.25',
    author='Drew Ellison Hart',
    author_email='drew.ellison.hart@gmail.com',
    description='A package for landscape genomic simulation',
    long_description=readme(),
    long_description_content_type='text/x-rst',
    url='https://github.com/drewhart/geonomics',
    # include the download URL, from the latest release on Github
    download_url='https://github.com/drewhart/geonomics/archive/0.0.25.tar.gz',
    include_package_data=True,
    # packages=setuptools.find_packages(),
    packages=['geonomics',
              'geonomics.sim',
              'geonomics.utils',
              'geonomics.structs',
              'geonomics.ops',
              'geonomics.help',
              'geonomics.demos'],
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
        'Intended Audience :: Science/Research',
        'Intended Audience :: Education',
        'Natural Language :: English',
        'Topic :: Scientific/Engineering :: Bio-Informatics'
    ],
    keywords=('landscape genomics genetics ecology evolution simulation model '
              'environmental model agent-based'),
    project_urls={
        'Documentation': ('https://htmlpreview.github.io/?https://github.com/'
                          'drewhart/geonomics/blob/master/doc/built/doc.html'),
        # 'Methods Paper': 'PAPER URL HERE!',
        'Source': 'https://github.com/drewhart/geonomics',
        # 'Tracker': 'BUGTRACKERSITHERE',
    },
    install_requires=['numpy', 'matplotlib>=3.0.0', 'pandas>=0.23.4', 'geopandas',
                      'scipy>=1.3.1', 'scikit-learn', 'statsmodels>=0.9.0',
                      'shapely', 'bitarray', 'pyvcf', 'rasterio',
                      'msprime>=0.7.4', 'tskit>=0.2.3'],
    extras_require={'simulation on neutral landscape models': ['nlmpy'],
                   '3d plots for Yosemite demo': ['pykrige']},
    python_requires='>=3.6',
    # package_data={
    # 'demos': ['geonomics/demos/yosemite/yosemite_30yr_normals_90x90.tif']
    # },
)
