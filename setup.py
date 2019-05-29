from setuptools import setup, find_packages

setup(
    name="tapestry",
    version="0.9.1",
    author="John Davey",
    author_email="johnomics@gmail.com",
    description="Validate and edit small eukaryotic genome assemblies",
    url="https://github.com/johnomics/tapestry",
    packages=['tapestry'],
    package_data={
        'tapestry': ['report/template.html', 'report/static/*.js'],
    },
    python_requires='>=3.6',
)
