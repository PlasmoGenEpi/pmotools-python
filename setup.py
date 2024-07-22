from setuptools import setup, find_packages

setup(
    name='pmotools',
    version='0.1',
    packages=find_packages(),
    install_requires=[
        'pandas==2.2.2',
        'openpyxl==3.1.4',
        'biopython==1.83',
        'plotly==5.22.0'
    ],
    entry_points={
        'console_scripts': [
            'microhaplotype_table_to_json=scripts.microhaplotype_table_to_json:main',
        ],
    },
)
