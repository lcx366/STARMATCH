from setuptools import setup,find_packages 

setup(
    name='starmatch',
    version='0.2.8',
    description='A package to handle the Star chart matching and astrometric positioning',
    author='Chunxiao Li',
    author_email='lcx366@126.com',
    url='https://github.com/lcx366/STARMATCH',
    license='MIT',
    long_description_content_type='text/markdown',
    long_description=open('README.md', 'rb').read().decode('utf-8'),
    keywords = ['Star Map Matching','Astrometric positioning','Distortion Calibration','Astronomical Correction'],
    python_requires = '>=3.10',
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Education',
        'Intended Audience :: Science/Research',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
        'Programming Language :: Python :: 3.12',
        'License :: OSI Approved :: MIT License',
        ],
    packages = find_packages(),
    include_package_data=True,
    install_requires=[
        'numpy',
        'scipy',
        'matplotlib',
        'scikit-image',
        'astropy',
        'GPy',
        'statsmodels',
        'loess',
        'starcatalogquery'
        ],
)
