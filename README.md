[![Documentation](https://github.com/TheoChem-VU/PyOrb/actions/workflows/build_docs.yml/badge.svg)](https://github.com/TheoChem-VU/PyOrb/actions/workflows/build_docs.yml) 
[![Testing](https://github.com/TheoChem-VU/PyOrb/actions/workflows/testing.yml/badge.svg)](https://github.com/TheoChem-VU/PyOrb/actions/workflows/testing.yml)
[![Publishing to PyPI](https://github.com/TheoChem-VU/PyOrb/actions/workflows/pypi_publish.yml/badge.svg?branch=main)](https://github.com/TheoChem-VU/PyOrb/actions/workflows/pypi_publish.yml)

# PyOrb
Our group aims to empower chemists with state-of-the-art quantum-chemical analysis tools. The latest project is PyOrb, used for the automated analysis of complex molecular orbital interactions. PyOrb makes expert orbital analysis available to all chemists, regardless of their theoretical background, ensuring that knowledge knows no bounds!

## Welcome 
Welcome to the PyOrb project repository. This project is currently led by members of the TheoCheM group at the Vrije Universiteit Amsterdam (Tori Gijzen, Yuman Hordijk, and Trevor Hamlin). We aim to empower computational and non-computational chemists to perform advanced orbital analyses. By presentation of computational data in 

This README file serves as a central hub providing project-related information. Feel free to navigate directly to specific sections below or simply scroll down to learn more.

- [Motivation](#motivation)
- [Installation](#installation)
- [Workflow](#workflow)
- [Examples](#examples)
- [Get involved](#involved)
- [Who are we?](#members)
- [Contact us](#contact)
- [Find out more](#more_info)

We are striving to optimise our application to facilitate anyone’s journey into the realm of advanced orbital analysis. So if you see any points where we can improve, don’t hesitate to contact us or open an issue in this repository! 

## Motivation <a name=motivation></a>
The PyOrb project aims to assist in advanced orbital analysis for computational chemistry calculations. PyOrb will generate important orbital interaction schemes and related orbital figures automatically, as well as summarise important information such as SFO coefficients, populations, orbital levels and more. These properties are commonly used in computational chemical research analyses.

People often run into the following challenges:
- Hard to get a clear overview of relevant data.
- Data acquisition is often done manually, which can be a chaotic and laborious task.
- The information overload in the output of quantum chemical programs can be intimidating.

The aim of the PyOrb project is:
- Enhance accesability to advanced orbital analyses
- To provide an efficient way to isolate all the relevant data
- Create a clear graphical representation  of the relevant data that is relevant 

Previously, our group has published the [PyFrag 2019](https://onlinelibrary.wiley.com/doi/10.1002/jcc.25871) program for the automation of fragment approaches and their quantitative analyses along a reaction coordinate. This allowed users to perform the Activation Strain Model (AMS) and the canonical energy decomposition analysis (EDA) for a whole reaction profile with one single input.

## Installation <a name=installation></a>
*Currently, the installation can be done by cloning this repository and having all necessary dependencies.*

The following is for people who would like to install the repository themselves. For example, to edit and/or contribute code to the project.

First clone this repository:
```
git clone https://github.com/TheoChem-VU/PyOrb.git
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

## Examples <a name=examples></a>
[insert example outputs]

## Get involved <a name=involved></a>
At the moment, contribution to this project is limited to our research group. In the foreseeable future, we want to open up to a greater community. This group will include users who a not necessarily familiar with computational analysis of molecular systems or fragment approaches. At the point of the first public rollout, we want to have such groups involved in the process of enhancing this repository. 

## Who are we? <a name=members></a>
A preliminary version of PyOrb was roled out by Xiaobo Sun assited by Laurens Groot, who was a PostDoc in the TheoCheM group at the Vrije Universiteit Amsterdam.

Currently, the project has been adopted by a team of active TheoCheM group members consisting of Yuman Hordijk ([@YumanHordijk](https://twitter.com/YumanHordijk)), Tori Gijzen ([@ToriGijzen](https://twitter.com/ToriGijzen)), and dr. Trevor A. Hamlin ([@TrevorAHamlin](https://twitter.com/TrevorAHamlin)). 

Additionally, we are being consulted by Nadine Spychala ([@NadineSpychala](https://twitter.com/NadineSpychala)) from the University of Sussex as a part of cohort 8 in the Open Life Sciences Scheme (<https://openlifesci.org/openseeds/ols-8/projects-participants.html>, [@openlifesci](https://twitter.com/openlifesci)).

## Contact us <a name=contact></a>
If you want to report a problem or suggest an enhancement we'd encourage for you to open an issue in the GitHub repository.  Otherwise you can reach us by email (<https://www.theochem.nl/contact>) or twitter ([@VU_TheoCheM](https://twitter.com/VU_TheoCheM))


## Find out more <a name=more_info></a>
You Might be interested in"
- The [TheoCheM group](https://www.theochem.nl/)
- The [Vrije Universiteit Amsterdam](https://vu.nl/en)
- The [PyFrag 2019](https://github.com/TheoChem-VU/PyFrag) program 
- Our python library [TC-utility](https://github.com/TheoChem-VU/TCutility)
- For a thorough overview of the theoretical background we recommend the following articles:
    - P. Vermeeren, S.C.C. van der Lubbe, C. Fonseca Guerra, F.M. Bickelhaupt, T.A.
Hamlin, _Nature Protoc._ **2020**, _15_, 649-667. (<https://doi.org/10.1038/s41596-019-0265-0>)
    - T. A. Hamlin, P. Vermeeren, C. Fonseca Guerra, F. M. Bickelhaupt.
In: _Complementary Bonding Analyses_; S. Grabowski, Ed.; De Gruyter: Berlin,
2021, pp 199-212. (<https://doi.org/10.1515/9783110660074-008>)
