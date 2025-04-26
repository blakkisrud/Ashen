# Ashen

 <p align="center"> <img src="resources/logo.png" alt="Ashen Logo" width="200"/> </p> 

Ashen is a Python resource designed to extract and work with information from ICRP Publication 107, a comprehensive collection of physical properties of radionuclides. It provides simple access to radionuclide data for use in research, clinical applications, or educational purposes.

This was mostly made as a pet project in accociation with absorbed dose calculation for radionuclide therapy. The goal is to provide a simple interface to access the data in ICRP 107, and to make it easy to use in Python.

This is still under development and contributions are welcome, please make a pull request or open an issue if you have any suggestions or improvements.

## Features

- Access radionuclide physical properties from ICRP 107

- Easy-to-use Python functions

- Demonstration script to get started

- Simple command-line tool for quick access to data

- Planned support for:

    - Simple graphical user interface (GUI)

## Installation

Clone the repository and install the required dependencies:
```bash
git clone https://github.com/yourusername/ashen.git
cd ashen
pip install -r requirements.txt
```

To install the ashen command line tool, run the following command in the root directory of the project:
```bash
pip install -e .
```
## Usage

- A simple example script is provided in the demo.py file to demonstrate how to use the functions in Ashen.

### Command line tool
You can also use the command line tool to access the data. The command line tool is called `ashen` and can be run from the terminal. The basic usage is as follows:

```bash
ashen [options] <radionuclide_name>
```

Where `<radionuclide_name>` is the name of the radionuclide you want to access. For example, to get the half life of I-131, you can run:

```bash
ashen get-half-life I-131
```

This will return the half life of I-131.

To get the summed energy emitted by I-131, you can run:

```bash
ashen energy I-131
```

Where you can specify the type of energy you want to get. The available options are:
- `--(b)eta`: Get the beta energy emitted by the radionuclide.
- `--(a)lpha`: Get the alpha energy emitted by the radionuclide.
- `--(n)on-penetrative`: Non penetrative energy emitted by the radionuclide.

like this:

```bash
ashen energy I-131 --b
```

All energy values are in MeV.

You can also use the `--help` option to see a list of available options and commands:

```bash
ashen --help
```

Future versions will include:

-     A lightweight graphical interface for browsing radionuclide properties.

## License

This work make use of the ICRP Publication 107 data. Please refer to the ICRP website for more information on the license and usage of the data.

Ashen is licensed under the MIT License. See the LICENSE file for more information.

