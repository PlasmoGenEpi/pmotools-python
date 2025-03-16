from setuptools import setup, find_packages

setup(
    name='pmotools',
    version='1.0.0',
    packages=find_packages(),
    install_requires=[
        'pandas==2.2.2',
        'openpyxl==3.1.4',
        'biopython==1.83',
        'plotly==5.22.0',
        'jsonschema==4.23.0'
    ]
)
