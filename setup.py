from setuptools import setup, find_packages

setup(
    name="tapestry",
    version="0.9.3",
    author="John Davey",
    author_email="johnomics@gmail.com",
    description="Validate and edit small eukaryotic genome assemblies",
    url="https://github.com/johnomics/tapestry",
    packages=['tapestry'],
    package_data={
        'tapestry': ['report/template.html', 'report/static/*.js', 'report/static/*.css'],
    },
    test_suite = 'test',
    scripts=['weave', 'clean'],
    python_requires='>=3.6',
)
