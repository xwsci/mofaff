from setuptools import setup, find_packages

setup(
    name='mofaff',
    version='0.1.23',
    packages=find_packages(),
    description='A package to calculate MOF adsorption with framework flexibility',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    author='Xijun Wang, Randall Q. Snurr*',
    author_email='xijun@northwestern.edu / wangxijun1016@gmail.com',
    url='https://github.com/XXX',
    install_requires=[
        'numpy>=1.24.4',
        'ase>=3.22.1',
    ],
    python_requires='>=3.8',
    entry_points={
        'console_scripts': [
            'mofaff=mofaff.mofaff_main:main',
        ],
    },
)

