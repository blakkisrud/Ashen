from setuptools import setup, find_packages

setup(
    name="ashen",
    version="1.0.8",
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
