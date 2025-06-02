from setuptools import setup, find_packages

setup(
    name="ashen",
    version="0.1.6",
    packages=find_packages(),
    install_requires=[
        "click",
    ],
    entry_points={
        "console_scripts": [
            "ashen=ashen.cli:cli",
        ],
    },
    package_data={
        "ashen": ["resources/*.TXT", "resources/*.RAD"],
    }
)
