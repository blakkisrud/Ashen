from setuptools import setup, find_packages

setup(
    name="ashen",
    version="0.1.2",
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
        "ashen_utils": ["resources/*.TXT", "resources/*.RAD"],
    }
)
