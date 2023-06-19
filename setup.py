from setuptools import setup,find_packages 

setup(
    name='starmatch',
    version='0.1.1',
    description='A package to handle the Star map matching and astronomical positioning',
    author='Chunxiao Li',
    author_email='lcx366@126.com',
    url='https://github.com/lcx366/STARMATCH',
    license='MIT',
    long_description_content_type='text/markdown',
    long_description=open('README.md', 'rb').read().decode('utf-8'),
    keywords = ['Star Map Matching','Astronomical positioning','Distortion Correction','pointing calibration'],
    python_requires = '>=3.8',
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Education',
        'Intended Audience :: Science/Research',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
        'License :: OSI Approved :: MIT License',
        ],
    packages = find_packages(),
    include_package_data=True,
    install_requires=[
        'numpy',
        'scipy',
        'matplotlib',
        'scikit-image',
        'photutils',
        'astropy>=4.3.1',
        'GPy',
        'Pillow',
        'colorama',
        'starcatalogquery'
        ],
)
