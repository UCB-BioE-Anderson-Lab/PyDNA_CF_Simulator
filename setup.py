from setuptools import setup, find_packages

setup(
    name='PyDNA_CF_Simulator',
    version='0.1',
    packages=find_packages(),
    install_requires=[
        'flask',
        'flask-cors',
        'pydna'
    ],
    author='J. Christopher Anderson',
    author_email='jcanderson@berkeley.edu',
    description='A ChatGPT plugin wrapper of pydna to simulate construction file.',
    license='MIT',
    keywords='synthetic biology, DNA assembly, simulation',
)