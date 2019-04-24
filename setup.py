from setuptools import setup, find_packages

with open('README') as fh:
    long_description = fh.read()

requirements = [
    'biopython==1.72',
    'numpy==1.15.4',
    'scipy==1.2.0rc1',
    'numba==0.40.1',
    'matplotlib==3.0.3',
]

setup(
    name="pydca",
    version="0.1",
    python_requires='>=3.5',
    packages=find_packages(
        exclude=["*.tests","*.tests.*","tests.*", "tests",
            "*.extras", "*.extras.*", "extras.*", "extras",
        ],
    ),
    install_requires= requirements,
    tests_require = requirements,
    entry_points={
        'console_scripts':[
            'mfdca = pydca.main:run_meanfield_dca',
        ],
    },
    test_suite="tests",
    author='Mehari B. Zerihun',
    author_email='mbzerihun@gmail.com',
    description='mean-field implementation of DCA for protein, and RNA sequences',
    long_description=long_description,

)
