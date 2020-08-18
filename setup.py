from setuptools import setup, find_packages, Extension
from setuptools.command.build_ext import build_ext 

with open("README.md") as fh:
    long_description = fh.read()

requirements = [
    "scipy==1.3.1",
    "biopython==1.74",
    "numpy>=1.13.3, <=1.15.4",
    'llvmlite==0.30.0',
    "numba==0.46.0",
    "matplotlib==3.0.0",
    "requests>=2.22.0",
]

plmdca_compile_args = ["-fopenmp", "-std=c++11", "-O3"]  
plmdca_link_args = ["-fopenmp", "-O3"] 


plmdca_ext = Extension(
    'pydca.plmdca._plmdcaBackend',
    [   
        'pydca/plmdca/lbfgs/lib/lbfgs.cpp',
        'pydca/plmdca/plmdca_numerics.cpp',
        'pydca/plmdca/plmdcaBackend.cpp', 
    ],
    include_dirs=[
        'pydca/plmdca/include/',
        'pydca/plmdca/lbfgs/include/',
    ],
    extra_compile_args = plmdca_compile_args,
    extra_link_args = plmdca_link_args,
    language = "c++",  
)

setup(
    name="pydca",
    version="1.23",
    author="Mehari B. Zerihun",
    author_email="mbzerihun@gmail.com",
    python_requires=">=3.5",
    description="Direct couplings analysis (DCA) for protein and RNA sequences",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/KIT-MBS/pydca",
    download_url="https://pypi.org/project/pydca/",
    packages=find_packages(
        exclude=["*.tests","*.tests.*","tests.*", "tests",
            "*.extras", "*.extras.*", "extras.*", "extras",
            "examples", "*.examples", "examples.*", "*.examples.*",
            "install.sh",
        ],
    ),
    ext_modules = [plmdca_ext],
    classifiers=[
        "Programming Language :: Python :: 3",
        "Programming Language :: C++",
        "Programming Language :: C",
        "License :: OSI Approved :: MIT License",
        "Development Status :: 4 - Beta",
        "Operating System :: POSIX :: Linux",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    install_requires= requirements,
    tests_require = requirements,
    entry_points={
        "console_scripts":[
            "mfdca = pydca.mfdca_main:run_meanfield_dca",
            "plmdca = pydca.plmdca_main:run_plm_dca",
            "pydca = pydca.main:run_pydca",
        ],
    },
    test_suite="tests",
)