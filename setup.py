import setuptools


def readme():
    with open("README.rst", 'r') as f:
        return f.read()


setuptools.setup(
    name='geonomics',
    version='0.0',
    author='Drew Ellison Hart',
    author_email='drew.ellison.hart@gmail.com',
    description='A package for landscape genomic simulation',
    long_description=readme(),
    long_description_content_type='text/x-rst',
    url='https://github.com/drewhart/geonomics',
    include_package_data=True,
    packages=setuptools.find_packages(),
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
        'Documentation': 'PUTDOCURLHERE',
        'Methods Paper': 'URLTOMETHODSPAPERHERE',
        'Source': 'https://github.com/drewhart/geonomics',
        'Tracker': 'BUGTRACKERSITHERE',
    },
    install_requires=['numpy', 'matplotlib', 'pandas', 'scipy', 'scikit-learn',
                      'statsmodels', 'shapely', 'bitarray', 'pyvcf'],
    extras_require={
        'simulation on neutral landscape models': ['nlmpy'],
        'reading and writing of common raster data formats': ['GDAL']
    },
    python_requires='>=3.5',
    packages_data={
        'example': ['yosemite_30yr_normals_90x90.tif']
    },
)
