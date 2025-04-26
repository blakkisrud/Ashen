from setuptools import setup, find_packages

setup(
    name="ashen",
    version="0.1",
    packages=find_packages(),
    install_requires=[
        "click",
    ],
    entry_points={
        "console_scripts": [
            "ashen=ashen.cli:cli",
        ],
    },
)
