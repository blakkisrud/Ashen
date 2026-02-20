from setuptools import setup, find_packages

setup(
    name="ashen",
    version="1.0.8",
    packages=find_packages(),
    install_requires=[
        "click",
    ],
    #include_package_data=True,
    #package_data={"ashen": ["resources/FULL_RAD_LIST.RAD"]},
    package_data={
        "ashen": ["resources/*.TXT", "resources/*.RAD", "resources/*.BET"],
    },
    entry_points={
        "console_scripts": [
            "ashen=ashen.cli:cli",
        ],
    },
)
