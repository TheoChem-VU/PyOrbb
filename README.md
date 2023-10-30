[![Documentation](https://github.com/TheoChem-VU/PyOrb/actions/workflows/build_docs.yml/badge.svg)](https://github.com/TheoChem-VU/PyOrb/actions/workflows/build_docs.yml) 
[![Testing](https://github.com/TheoChem-VU/PyOrb/actions/workflows/testing.yml/badge.svg)](https://github.com/TheoChem-VU/PyOrb/actions/workflows/testing.yml)
[![Publishing to PyPI](https://github.com/TheoChem-VU/PyOrb/actions/workflows/pypi_publish.yml/badge.svg?branch=main)](https://github.com/TheoChem-VU/PyOrb/actions/workflows/pypi_publish.yml)

# PyOrb
Our group aims to empower chemists with state-of-the-art quantum-chemical analysis tools. The latest project is PyOrb, used for the automated analysis of complex molecular orbital interactions. PyOrb makes expert orbital analysis available to all chemists, regardless of their theoretical background, ensuring that knowledge knows no bounds!

## Welcome 
Welcome to the PyOrb project repository. This project is currently led by members of the TheoCheM group at the Vrije Universiteit Amsterdam (Tori Gijzen, Yuman Hordijk, and Trevor Hamlin). We aim to empower computational and non-computational chemists to perform advanced orbital analyses. By presentation of computational data in 

This README file serves as a central hub providing project-related information. Feel free to navigate directly to specific sections below or simply scroll down to learn more.

- [Motivation](#Motivation)
- [Installation](#Installation)
- [Workflow](#workflow)
- [Examples](#Examples)
- [Get involved](#involved)
- [Who are we?](#Who)
- [Contact us](#Contact)
- [Find out more](#more)

We are striving to optimise our application to facilitate anyone’s journey into the realm of advanced orbital analysis. So if you see any points where we can improve, don’t hesitate to contact us or open an issue in this repository! 

## Motivation 
The PyOrb project aims to assist in advanced orbital analysis for computational chemistry calculations. PyOrb will generate important orbital interaction schemes and related orbital figures automatically, as well as summarise important information such as SFO coefficients, populations, orbital levels and more. These properties are commonly used in computational chemical research analyses.

It is often said 
- Hard to get a clear overview of relevant data.
- Often done manually, which can be a chaotic and laborious task.
- The information overload in the output of quantum chemical programs can be intimidating.

The aim of the PyOrb project:
- To provide an efficient way to isolate all the relevant data
- Clear projection  of the data that is relevant (color and size)
- Let you choose what data you want to include in your analysis

Previous work by our group has focused on the automation of fragment approaches and their quantitative analyses along a reaction coordinate. This allowed users to perform the Activation Strain Model (AMS) and the canonical energy decomposition analysis (EDA) for a whole reaction profile with one single input (link).

## Installation 
*Currently, the installation can be done by cloning this repository and having all necessary dependencies.*

The following is for people who would like to install the repository themselves. For example, to edit and/or contribute code to the project.

First clone this repository:
```
git clone https://github.com/TheoChem-VU/Pyorb.git
```

Then move into the new directory and install the package:

```
cd PyOrb
python -m pip install --upgrade build 
python -m build 
python -m pip install -e .
```

To get new updates, simply run:
```
git pull
```
## Workflow <a name=workflow></a> 
[AMS rkf file]
[should be included once this is known]
[Figure horizontal si better than vertical]

## Examples
[insert example outputs]

## Get involved
At the moment, contribution to this project is limited to our research group. In the foreseeable future, we want to open up to a greater community. This group will include users who a not necessarily familiar with computational analysis of molecular systems or fragment approaches. At the point of the first public rollout, we want to have such groups involved in the process of enhancing this repository. 

## Who are we?
A preliminary version of PyOrb was roled out by Xiaobo Sun, who was a PostDoc in the TheoCheM group at the Vrije Universiteit Amsterdam.

Currently, the project is adopted by a team of active members of the TheoCheM group consisting of Yuman Hordijk (@YumanHordijk), Tori Gijzen (@ToriGijzen), and dr. Trevor A. Hamlin (@TrevorAHamlin). 

Additionally, we are being consulted by Nadine Spychala (@NadineSpychala) from the University of Sussex as a part of cohort 8 in the Open Life Sciences Scheme (https://openlifesci.org/openseeds/ols-8/projects-participants.html, @openlifesci).

## Contact us
If you want to report a problem or suggest an enhancement we'd encourage for you to open an issue in the GitHub repository.  Otherwise you can reach us by email (????) or twitter (????)
[Add contact page from theochem website]

## Find out more
- The link to the theochem group - https://www.theochem.nl/ 
- Link to pry-frag - https://www.theochem.nl/pyfrag2019
- TC-utility - https://github.com/TheoChem-VU/TCutility
- Reference paper for theory 
    - https://doi.org/10.1038/s41596-019-0265-0 
    - Should be more
