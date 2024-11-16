from setuptools import setup, find_packages

setup(
    name="csdl_integrator",
    version="0.1.0",
    description="Numerical integrators for ODEs using CSDL.",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    author="Sankalp Kaushik",
    author_email="sskaushi@ucsd.edu",
    url="https://github.com/butterboy265/csdl_integrator.git",
    packages=find_packages(),
    install_requires=[
        "numpy>=1.20.0",
        "csdl-alpha>=0.0.1"
    ],
    dependency_links=[
        "git+https://github.com/LSDOlab/CSDL_alpha.git"
    ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.7',
)
