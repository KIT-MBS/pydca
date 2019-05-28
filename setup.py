from setuptools import setup, find_packages

with open("README.md") as fh:
    long_description = fh.read()

requirements = [
    "biopython==1.72",
    "numpy==1.15.4",
    "scipy==1.2.0rc1",
    "numba==0.40.1",
    "matplotlib==3.0.3",
]

setup(
    name="pydca",
    version="0.2",
    author="Mehari B. Zerihun",
    author_email="mbzerihun@gmail.com",
    python_requires=">=3.5",
    description="Direct couplings analysis (DCA) for protein and RNA sequences",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/KIT-MBS/pydca",
    packages=find_packages(
        exclude=["*.tests","*.tests.*","tests.*", "tests",
            "*.extras", "*.extras.*", "extras.*", "extras",
        ],
    ),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Development Status :: 4 - Beta",
    ],
    install_requires= requirements,
    tests_require = requirements,
    entry_points={
        "console_scripts":[
            "mfdca = pydca.main:run_meanfield_dca",
        ],
    },
    test_suite="tests",
)
