from setuptools import setup
from gumpy import __version__

setup(
    name='Seq&Treat data processing',
    version=__version__,
    author='Philip W Fowler',
    author_email="philip.fowler@ndm.ox.ac.uk",
    description="Simple scripts and functions to validate and process genetic and AST data donated to the Seq&Treat project.",
    url="https://github.com/oxfordmmm/seqtreat",
    packages=['seqtreat'],
    package_data={''},
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent"],
    python_requires='>=3.5',
    install_requires=[
        "numpy >= 1.18",
        "pandas >= 1.0"
    ],
    license="MIT",
    scripts=['bin/seqtreat-spreadsheet-validate.py'],\
    zip_safe=False
)
